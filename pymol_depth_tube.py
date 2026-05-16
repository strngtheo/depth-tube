"""Camera-distance-driven backbone tube for PyMOL.

Uses PyMOL's built-in ``cartoon putty`` mode, which already varies the tube
radius per residue from the B-factor. We hijack the B-factor channel and write
projected camera-depth into it on a 10 Hz timer.

Two effects compound, matching the ChimeraX implementation:

    proximity taper   per-residue B = depth-rank along the view axis; putty
                      maps that onto [putty_scale_min, putty_scale_max]
    zoom scale        chain mean depth maps onto [zoom_radius_min,
                      zoom_radius_max] between [zoom_d_min Å, zoom_d_max Å]
                      and is fed into cartoon_putty_radius

Load into a PyMOL session::

    run /path/to/pymol_depth_tube.py

then::

    depthtube chain A and polymer
    depthtube chain B, r_proximal=0.55, r_distal=0.06
    depthtube_off chain A      # disable one selection
    depthtube_off             # disable everything
"""
from __future__ import annotations

from collections import OrderedDict

import numpy as np
from pymol import cmd, stored

try:
    from pymol.Qt import QtCore
    _HAS_QT = True
except Exception:
    _HAS_QT = False
    import threading


_ENABLED: "OrderedDict[str, dict]" = OrderedDict()
_CAM_DPOS = 0.25
_CAM_DDIR = 0.003


def _camera():
    v = cmd.get_view()
    R = np.array(v[0:9], dtype=float).reshape(3, 3)
    origin = np.array(v[12:15], dtype=float)
    cam_offset = np.array(v[9:12], dtype=float)
    view_dir = R[2]
    n = float(np.linalg.norm(view_dir)) or 1.0
    view_dir = view_dir / n
    cam_pos = origin - cam_offset[2] * view_dir
    return cam_pos, view_dir


def _depths_for(selection, cam_pos, view_dir):
    coords, ids = [], []
    stored._coords = coords
    stored._ids = ids
    cmd.iterate_state(
        1, f"({selection}) and name CA and polymer",
        "stored._ids.append(ID); stored._coords.append((x,y,z))",
    )
    if not coords:
        return None, None
    coords = np.asarray(coords, dtype=float)
    depths = (coords - cam_pos) @ view_dir
    return np.asarray(ids, dtype=int), depths


def enable(selection,
           r_proximal=0.50, r_distal=0.05,
           zoom_radius_min=0.15, zoom_radius_max=0.55,
           zoom_d_min=25.0, zoom_d_max=160.0,
           tick_hz=10):
    """Enable depth-tube for one selection."""
    sel = selection.strip() or "polymer"
    if cmd.count_atoms(f"({sel}) and name CA") < 2:
        print(f"[depth-tube] {sel!r}: <2 Cα; skipped")
        return

    cmd.cartoon("putty", sel)
    cmd.show("cartoon", sel)
    cmd.set("cartoon_putty_scale_min", float(r_distal), quiet=1)
    cmd.set("cartoon_putty_scale_max", float(r_proximal), quiet=1)
    cmd.set("cartoon_putty_transform", 0, quiet=1)
    cmd.set("cartoon_putty_quality", 11, quiet=1)

    key = sel
    _disable_one(key)

    state = {
        "sel": sel,
        "r_proximal": float(r_proximal),
        "r_distal": float(r_distal),
        "zoom_radius_min": float(zoom_radius_min),
        "zoom_radius_max": float(zoom_radius_max),
        "zoom_d_min": float(zoom_d_min),
        "zoom_d_max": float(zoom_d_max),
        "last_pos": None,
        "last_view": None,
        "tick_hz": int(tick_hz),
    }
    _ENABLED[key] = state
    _start_timer_if_needed()
    _update_one(state, force=True)
    print(f"[depth-tube] enabled on {sel!r}")


def _update_one(state, force=False):
    sel = state["sel"]
    if cmd.count_atoms(f"({sel}) and name CA") < 2:
        _disable_one(sel)
        return
    cam_pos, view_dir = _camera()

    if not force and state["last_pos"] is not None:
        dp = float(np.linalg.norm(cam_pos - state["last_pos"]))
        dd = 1.0 - float(np.dot(view_dir, state["last_view"]))
        if dp < _CAM_DPOS and dd < _CAM_DDIR:
            return
    state["last_pos"] = cam_pos.copy()
    state["last_view"] = view_dir.copy()

    ids, depths = _depths_for(sel, cam_pos, view_dir)
    if ids is None:
        return

    d_min, d_max = float(depths.min()), float(depths.max())
    d_range = max(1e-6, d_max - d_min)
    mean_d = 0.5 * (d_min + d_max)

    zfrac = max(0.0, min(1.0,
                         (mean_d - state["zoom_d_min"])
                         / max(1e-6, state["zoom_d_max"] - state["zoom_d_min"])))
    putty_radius = (state["zoom_radius_min"]
                    + zfrac * (state["zoom_radius_max"] - state["zoom_radius_min"]))

    # Smaller depth = closer to camera in this projection → higher b → fatter.
    fracs = (depths - d_min) / d_range  # 0 = nearest, 1 = farthest
    b_values = 1.0 - fracs              # 0..1, near=1, far=0
    stored._b_lookup = {int(i): float(b) for i, b in zip(ids, b_values)}
    cmd.alter(
        f"({sel}) and polymer",
        "b = stored._b_lookup.get(ID, b)",
    )
    cmd.set("cartoon_putty_radius", float(putty_radius), quiet=1)
    cmd.rebuild(sel, "cartoon")


def _disable_one(key):
    state = _ENABLED.pop(key, None)
    if state is None:
        return
    try:
        cmd.cartoon("automatic", state["sel"])
    except Exception:
        pass


def disable(selection=None):
    if selection is None or not selection.strip():
        for k in list(_ENABLED.keys()):
            _disable_one(k)
    else:
        _disable_one(selection.strip())
    if not _ENABLED:
        _stop_timer()


# ─── timer machinery ─────────────────────────────────────────────────────────
_TIMER = None
_THREAD_TIMER = None


def _tick():
    for state in list(_ENABLED.values()):
        try:
            _update_one(state)
        except Exception as e:
            print(f"[depth-tube] update error on {state['sel']!r}: "
                  f"{type(e).__name__}: {e}")
            _disable_one(state["sel"])


def _start_timer_if_needed():
    global _TIMER, _THREAD_TIMER
    if _ENABLED and _TIMER is None and _THREAD_TIMER is None:
        if _HAS_QT:
            _TIMER = QtCore.QTimer()
            _TIMER.setInterval(100)
            _TIMER.timeout.connect(_tick)
            _TIMER.start()
        else:
            def _loop():
                global _THREAD_TIMER
                _tick()
                if _ENABLED:
                    _THREAD_TIMER = threading.Timer(0.1, _loop)
                    _THREAD_TIMER.daemon = True
                    _THREAD_TIMER.start()
            _loop()


def _stop_timer():
    global _TIMER, _THREAD_TIMER
    if _TIMER is not None:
        try:
            _TIMER.stop(); _TIMER.deleteLater()
        except Exception:
            pass
        _TIMER = None
    if _THREAD_TIMER is not None:
        try:
            _THREAD_TIMER.cancel()
        except Exception:
            pass
        _THREAD_TIMER = None


# ─── PyMOL command registration ──────────────────────────────────────────────
def _cmd_enable(selection, r_proximal=0.50, r_distal=0.05,
                zoom_radius_min=0.15, zoom_radius_max=0.55,
                zoom_d_min=25.0, zoom_d_max=160.0):
    enable(selection,
           r_proximal=float(r_proximal), r_distal=float(r_distal),
           zoom_radius_min=float(zoom_radius_min),
           zoom_radius_max=float(zoom_radius_max),
           zoom_d_min=float(zoom_d_min), zoom_d_max=float(zoom_d_max))


def _cmd_disable(selection=""):
    disable(selection if selection else None)


cmd.extend("depthtube", _cmd_enable)
cmd.extend("depthtube_off", _cmd_disable)

print("[depth-tube] PyMOL commands registered: depthtube, depthtube_off")
