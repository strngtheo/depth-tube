# depth-tube

Camera-distance-driven backbone tube for **ChimeraX** and **PyMOL**.

The selected chain's backbone is rendered as a tube whose per-segment
radius is recomputed on every frame from two effects:

| Effect          | What it does                                                                                                  |
|-----------------|---------------------------------------------------------------------------------------------------------------|
| Proximity taper | The Cα closest to the camera gets `r_proximal`; the farthest gets `r_distal`; linear ramp on projected depth. |
| Zoom scale      | The chain's mean depth maps onto `[zoom_*_min, zoom_*_max]` between `[zoom_d_min, zoom_d_max]` Å.             |

This makes near loops read as fat / saturated and far loops as thin —
strong spatial cueing for crowded interfaces, flexible loops, and
docking poses.

## ChimeraX

```text
runscript /path/to/chimerax_depth_tube.py

depthtube #1/A
depthtube #1/B r_proximal 0.7 r_distal 0.04 zoom_d_max 200
~depthtube #1/A       # disable one chain
~depthtube            # disable all
```

How it works (file: `chimerax_depth_tube.py`):

* On enable, a synthetic `AtomicStructure` is created with one atom per
  Cα and bonds between consecutive Cαs. Chain breaks at >5 Å Cα-Cα
  open a new run.
* The original chain's cartoon is hidden (`~cartoon #id/chain`).
* A `QTimer` ticks at 10 Hz; each tick reads `session.main_view.camera`,
  computes `(midpoint - cam_pos) @ view_dir` per bond, and writes new
  radii via the batched `Structure.bonds.radii = …` assignment.
* Cα-Cβ stubs and intra-residue sidechain bonds on the *real* structure
  piggy-back at 0.50× the tube radius.
* Camera-delta gating skips redundant ticks (`< 0.25 Å` translation AND
  `< 0.003` cos-delta in view direction).

## PyMOL

```text
run /path/to/pymol_depth_tube.py

depthtube chain A and polymer
depthtube chain B and polymer, r_proximal=0.55, r_distal=0.06
depthtube_off chain A
depthtube_off
```

How it works (file: `pymol_depth_tube.py`):

* Uses PyMOL's built-in `cartoon putty` mode, which already varies the
  cartoon tube radius per residue from the B-factor (smooth, GPU-accelerated).
* On a 10 Hz QTimer (falls back to `threading.Timer` if Qt is unavailable):
  * reads `cmd.get_view()` and reconstructs camera position + view direction;
  * projects each Cα onto the view axis;
  * normalises to `[0, 1]` (near=1, far=0);
  * writes those values into the B-factor channel via `cmd.alter` + a
    `stored._b_lookup` dict (single `alter` call, no Python loop);
  * adjusts `cartoon_putty_radius` from the chain's mean depth;
  * `cmd.rebuild(sel, "cartoon")` to refresh.
* The proximity taper is realised by `cartoon_putty_scale_min` /
  `cartoon_putty_scale_max` (set to `r_distal` and `r_proximal`).

The PyMOL version actually looks *smoother* than the ChimeraX one
because putty interpolates the radius along the spline; the ChimeraX
version draws discrete per-segment cylinders.

## Tuning

| Param                              | Default | Effect                                            |
|------------------------------------|---------|---------------------------------------------------|
| `r_proximal`                       | 0.50    | Tube radius (Å) at the chain's nearest Cα         |
| `r_distal`                         | 0.05    | Tube radius (Å) at the chain's farthest Cα        |
| `zoom_scale_min` / `_max` (CXC)    | 0.7 / 2.5 | Multiplier on `r_proximal`/`r_distal` from zoom |
| `zoom_radius_min` / `_max` (PML)   | 0.15 / 0.55 | Final `cartoon_putty_radius` band             |
| `zoom_d_min` / `_max`              | 25 / 160 | Mean-depth Å range that drives the zoom factor   |
| `max_ca_ca` (CXC only)             | 5.0     | Cα-Cα threshold for chain-break detection         |
| `tick_hz`                          | 10      | Update rate                                       |

## License

MIT — see `LICENSE`.
