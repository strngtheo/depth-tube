"""Camera-distance-driven backbone tube for ChimeraX.

The selected chain's cartoon is replaced by a synthetic pseudobond tube whose
per-segment radius is recomputed every frame from the camera:

    proximity taper   nearest segment gets `r_proximal`, farthest `r_distal`,
                      linear ramp on depth along the view axis
    zoom scale        chain mean depth maps onto [zoom_scale_min, zoom_scale_max]
                      between [zoom_d_min Å, zoom_d_max Å]

Cα-Cβ stubs and intra-residue side-chain bonds piggy-back at branch_scale (0.50)
× the tube radius so the whole chain reads with consistent depth cueing.

Load into a ChimeraX session::

    runscript "/path/to/chimerax_depth_tube.py"

then enable / disable per chain::

    depthtube #1/A
    depthtube #1/B r_proximal 0.7 r_distal 0.04 zoom_d_max 200
    ~depthtube #1/A           # disable one chain
    ~depthtube                # disable all
"""
from __future__ import annotations

import numpy as np

from chimerax.atomic import AtomicStructure, AtomsArg
from chimerax.core.commands import CmdDesc, FloatArg, register, run


_DYN_BB: dict = {}


def _disable(session, key=None):
    keys = [key] if key is not None else list(_DYN_BB.keys())
    for k in keys:
        entry = _DYN_BB.pop(k, None)
        if entry is None:
            continue
        t = entry.get("timer")
        if t is not None:
            try:
                t.stop(); t.deleteLater()
            except Exception:
                pass
        h = entry.get("handler")
        if h is not None:
            try:
                session.triggers.remove_handler(h)
            except Exception:
                pass
        tube = entry.get("tube")
        if tube is not None and not tube.deleted:
            try:
                session.models.close([tube])
            except Exception:
                pass
        DEFAULT_RADIUS = 0.2
        for triplet in entry.get("branch_bonds", []):
            b = triplet[0]
            if not b.deleted:
                b.radius = DEFAULT_RADIUS


def _ordered_cas(struct, chain_id):
    chain_obj = None
    for c in struct.chains:
        if c.chain_id == chain_id:
            chain_obj = c
            break
    if chain_obj is not None and hasattr(chain_obj, "existing_residues"):
        residues = list(chain_obj.existing_residues)
    else:
        residues = sorted(
            (r for r in struct.residues if r.chain_id == chain_id),
            key=lambda r: (r.number, getattr(r, "insertion_code", "") or ""),
        )
    cas = []
    for r in residues:
        a = r.find_atom("CA") if hasattr(r, "find_atom") else None
        if a is None:
            for aa in r.atoms:
                if aa.name == "CA":
                    a = aa
                    break
        if a is not None:
            cas.append(a)
    return residues, cas


def enable(session, struct, chain_id,
           r_proximal=0.50, r_distal=0.05,
           zoom_scale_min=0.7, zoom_scale_max=2.5,
           zoom_d_min=25.0, zoom_d_max=160.0,
           max_ca_ca=5.0, branch_scale=0.50,
           tick_hz=10):
    """Enable depth-tube for one chain. See module docstring for params."""
    key = (struct.id_string, chain_id)
    _disable(session, key=key)

    residues, cas = _ordered_cas(struct, chain_id)
    if len(cas) < 2:
        session.logger.info(f"[depth-tube] {chain_id}: <2 Cα found; skipped")
        return

    run(session, f"~cartoon #{struct.id_string}/{chain_id}")

    ca_coords = np.array([a.coord for a in cas], dtype=float)
    runs, cur = [], [0]
    for i in range(1, len(cas)):
        if float(np.linalg.norm(ca_coords[i] - ca_coords[i - 1])) > max_ca_ca:
            if len(cur) >= 2:
                runs.append(cur)
            cur = [i]
        else:
            cur.append(i)
    if len(cur) >= 2:
        runs.append(cur)
    if not runs:
        return

    tube = AtomicStructure(session, name=f"depth_tube_{struct.id_string}_{chain_id}")
    session.models.add([tube])

    all_atoms, all_bonds = [], []
    for run_ix, ridxs in enumerate(runs):
        res = tube.new_residue("TUB", "X", run_ix + 1)
        run_atoms = []
        for j in ridxs:
            atom = tube.new_atom(f"C{j}", "C")
            res.add_atom(atom)
            atom.coord = ca_coords[j]
            atom.draw_mode = atom.STICK_STYLE
            atom.display = True
            atom.radius = 0.01
            run_atoms.append(atom)
            all_atoms.append(atom)
        for i in range(len(run_atoms) - 1):
            b = tube.new_bond(run_atoms[i], run_atoms[i + 1])
            b.radius = float(r_proximal)
            b.halfbond = False
            all_bonds.append(b)

    if not all_bonds:
        _disable(session, key=key)
        return

    backbone = {"N", "CA", "C", "O", "OXT", "H", "HA"}
    branch_bonds = []
    seen = set()
    for r in residues:
        ca = cb = None
        for aa in r.atoms:
            if aa.name == "CA":
                ca = aa
            elif aa.name == "CB":
                cb = aa
        if ca is not None and cb is not None:
            for b in cb.bonds:
                if b.other_atom(cb) is ca:
                    branch_bonds.append((b, ca, cb))
                    break
        for a in r.atoms:
            if a.name in backbone:
                continue
            for b in a.bonds:
                other = b.other_atom(a)
                if other is None or other.residue is not r:
                    continue
                if other.name in backbone:
                    continue
                if id(b) in seen:
                    continue
                seen.add(id(b))
                branch_bonds.append((b, a, other))

    mids_bb = np.stack(
        [0.5 * (np.asarray(b.atoms[0].coord, dtype=float)
                + np.asarray(b.atoms[1].coord, dtype=float))
         for b in all_bonds], axis=0,
    )
    mids_branch = (
        np.stack([0.5 * (np.asarray(a1.coord, dtype=float)
                         + np.asarray(a2.coord, dtype=float))
                  for (_, a1, a2) in branch_bonds], axis=0)
        if branch_bonds else np.zeros((0, 3), dtype=float)
    )

    last = {"pos": None, "view": None}
    CAM_DPOS = 0.25
    CAM_DDIR = 0.003

    def _update(*_, **__):
        try:
            if struct.deleted or tube.deleted:
                _disable(session, key=key)
                return
            cam = session.main_view.camera
            cam_pos = np.asarray(cam.position.origin(), dtype=float)
            view_dir = np.asarray(cam.view_direction(), dtype=float)
            view_dir = view_dir / (float(np.linalg.norm(view_dir)) or 1.0)
            if last["pos"] is not None:
                if (float(np.linalg.norm(cam_pos - last["pos"])) < CAM_DPOS
                        and 1.0 - float(np.dot(view_dir, last["view"])) < CAM_DDIR):
                    return
            last["pos"] = cam_pos.copy(); last["view"] = view_dir.copy()

            d_bb = (mids_bb - cam_pos) @ view_dir
            d_min, d_max = float(d_bb.min()), float(d_bb.max())
            d_range = max(1e-6, d_max - d_min)
            mean_d = 0.5 * (d_min + d_max)
            zfrac = max(0.0, min(1.0, (mean_d - zoom_d_min)
                                 / max(1e-6, zoom_d_max - zoom_d_min)))
            zoom = zoom_scale_min + zfrac * (zoom_scale_max - zoom_scale_min)

            fracs = (d_bb - d_min) / d_range
            radii_bb = ((r_proximal + fracs * (r_distal - r_proximal)) * zoom
                        ).astype(np.float32)
            try:
                if len(tube.bonds) == len(all_bonds):
                    tube.bonds.radii = radii_bb
                else:
                    raise RuntimeError("bond count mismatch")
            except Exception:
                for b, r_ in zip(all_bonds, radii_bb):
                    b.radius = float(r_)
            try:
                if len(tube.atoms) == len(all_atoms):
                    a_r = np.empty(len(all_atoms), dtype=np.float32)
                    a_r.fill(float(radii_bb.mean()))
                    tube.atoms.radii = a_r
            except Exception:
                pass

            if branch_bonds:
                d_br = (mids_branch - cam_pos) @ view_dir
                f_br = np.clip((d_br - d_min) / d_range, 0.0, 1.0)
                radii_br = ((r_proximal + f_br * (r_distal - r_proximal))
                            * zoom * branch_scale).astype(np.float32)
                for (b, _, _), r_ in zip(branch_bonds, radii_br):
                    if not b.deleted:
                        b.radius = float(r_)
        except Exception as e:
            session.logger.warning(
                f"[depth-tube] {chain_id} update error: "
                f"{type(e).__name__}: {e}")
            _disable(session, key=key)

    h = None
    for trig in ("new frame", "frame drawn"):
        try:
            h = session.triggers.add_handler(trig, _update)
            break
        except Exception:
            pass

    from Qt.QtCore import QTimer
    timer = QTimer()
    timer.setInterval(max(20, int(1000 / max(1, tick_hz))))
    timer.timeout.connect(_update)
    timer.start()

    _DYN_BB[key] = {
        "tube": tube, "bonds": all_bonds, "atoms": all_atoms,
        "branch_bonds": branch_bonds, "handler": h, "timer": timer,
    }
    _update()
    session.logger.info(
        f"[depth-tube] {chain_id} enabled — {len(all_bonds)} bb segments + "
        f"{len(branch_bonds)} side-chain branches, {tick_hz} Hz")


def _chain_keys_from_atoms(atoms):
    seen = set()
    out = []
    for a in atoms:
        k = (a.structure, a.residue.chain_id)
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


def _cmd_enable(session, atoms,
                r_proximal=0.50, r_distal=0.05,
                zoom_scale_min=0.7, zoom_scale_max=2.5,
                zoom_d_min=25.0, zoom_d_max=160.0):
    for s, cid in _chain_keys_from_atoms(atoms):
        enable(session, s, cid,
               r_proximal=r_proximal, r_distal=r_distal,
               zoom_scale_min=zoom_scale_min, zoom_scale_max=zoom_scale_max,
               zoom_d_min=zoom_d_min, zoom_d_max=zoom_d_max)


def _cmd_disable(session, atoms=None):
    if atoms is None or len(atoms) == 0:
        _disable(session)
        return
    for s, cid in _chain_keys_from_atoms(atoms):
        _disable(session, key=(s.id_string, cid))


register(
    "depthtube",
    CmdDesc(
        required=[("atoms", AtomsArg)],
        keyword=[
            ("r_proximal", FloatArg),
            ("r_distal", FloatArg),
            ("zoom_scale_min", FloatArg),
            ("zoom_scale_max", FloatArg),
            ("zoom_d_min", FloatArg),
            ("zoom_d_max", FloatArg),
        ],
        synopsis="Enable camera-depth-driven backbone tube on selection",
    ),
    _cmd_enable,
)

register(
    "~depthtube",
    CmdDesc(
        optional=[("atoms", AtomsArg)],
        synopsis="Disable depth-tube on selection (or all if omitted)",
    ),
    _cmd_disable,
)

print("[depth-tube] commands registered: depthtube, ~depthtube")
