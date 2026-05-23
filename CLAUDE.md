# CLAUDE.md

Project context for Claude (and any other AI coding assistant) working in
this repository.

## What this is

`depth-tube` is a minimal two-file rendering plugin for **ChimeraX** and
**PyMOL** that rebinds backbone tube radius to camera distance, every
frame. The goal is to make 3D depth localisation graspable from a *still*
image without relying on rotation, fog, or stereo cues. See `README.md`
for the value proposition and tuning table.

The repo intentionally stays small and dependency-light. Don't add
package managers, build systems, or framework wrappers unless there's a
concrete reason — the whole point is "drop one .py into your viewer and
it works."

## Repo layout

```
depth-tube/
├── chimerax_depth_tube.py   ChimeraX implementation
├── pymol_depth_tube.py      PyMOL implementation
├── examples/
│   ├── load_chimerax.cxc    one-liner demo loader for ChimeraX
│   └── load_pymol.pml       one-liner demo loader for PyMOL
├── README.md                user-facing docs + motivation
├── CLAUDE.md                this file
├── LICENSE                  MIT
└── .gitignore
```

No tests yet — both viewers require their own runtime to load these
scripts, so the only realistic smoke-test is loading a PDB and running
the demo `.cxc` / `.pml` files (see "Testing" below).

## Implementation overview

### ChimeraX (`chimerax_depth_tube.py`)

Activation: `runscript chimerax_depth_tube.py` registers two commands
in the ChimeraX shell — `depthtube <atomspec>` and `~depthtube
[<atomspec>]`.

Core algorithm in `enable()`:

1. Collect ordered Cα atoms for the chain (chain order via
   `chain.existing_residues`, fall back to residue-number sort).
2. Split into "runs" wherever consecutive Cα-Cα distance exceeds
   `max_ca_ca=5.0 Å` (chain breaks / icode jumps).
3. Build a synthetic `AtomicStructure` named
   `depth_tube_<source_id>_<chain>`, one atom per Cα, bonds between
   consecutive atoms in each run.
4. Hide the source chain's cartoon with `~cartoon #id/chain`.
5. Cache static bond midpoints for the tube, the Cα-Cβ stubs, and
   intra-residue side-chain bonds.
6. Start a `Qt.QtCore.QTimer` at 10 Hz running `_update()`.

Per tick (`_update`):

- Read `session.main_view.camera` → `cam_pos`, normalised `view_dir`.
- Skip if `|Δcam_pos| < 0.25 Å` AND `1 - cos(Δview) < 0.003` (idle gating).
- `depths = (mids - cam_pos) @ view_dir` — projection along view axis.
- Linear ramp: `r_proximal` at `depths.min()`, `r_distal` at `depths.max()`.
- Multiply by `zoom_scale`, derived from the chain's mean depth via
  the `[zoom_d_min, zoom_d_max]` → `[zoom_scale_min, zoom_scale_max]` map.
- Batched write: `tube.bonds.radii = new_radii` (fast path) with a per-bond
  Python loop fallback. Same scheme for Cα-Cβ and side-chain bonds at
  `branch_scale=0.50` × tube radius.

State is held in module-level dict `_DYN_BB`, keyed by
`(source_struct.id_string, chain_id)`. Multiple chains can run
independently. `_disable()` closes the synthetic structure, stops the
timer, removes the trigger handler, and resets touched side-chain bond
radii to `DEFAULT_RADIUS=0.2`.

### PyMOL (`pymol_depth_tube.py`)

Activation: `run pymol_depth_tube.py` registers `depthtube
<selection>` and `depthtube_off [<selection>]`.

Key insight: **PyMOL already does width-from-channel rendering** via
`cartoon putty` mode, which reads B-factor per residue. We don't draw
any custom geometry — we just write camera-projected depth into the
B-factor channel.

Core algorithm:

1. `cmd.cartoon("putty", sel)`, set `cartoon_putty_scale_min/max` to
   `r_distal` / `r_proximal`, `cartoon_putty_transform=0` (B-factors
   used directly, no normalisation).
2. Start a `pymol.Qt.QtCore.QTimer` at 10 Hz (falls back to
   `threading.Timer` if `pymol.Qt` is unavailable, e.g. headless).

Per tick (`_update_one`):

- Reconstruct camera position + view direction from `cmd.get_view()`
  (see `_camera()`).
- Project each Cα onto the view axis.
- Normalise depths to `[0, 1]` (near=1, far=0), drop into
  `stored._b_lookup = {ID: b}`.
- Single `cmd.alter("(sel) and polymer", "b = stored._b_lookup.get(ID, b)")`
  — no per-atom Python loop.
- Update `cartoon_putty_radius` from chain mean depth (zoom effect).
- `cmd.rebuild(sel, "cartoon")` to refresh the spline.

The PyMOL output looks smoother than the ChimeraX one because putty
interpolates the radius continuously along the spline; ChimeraX draws
discrete per-bond cylinders.

## Conventions

- **Two-platform parity.** Behaviour-affecting changes (e.g. tweaking
  the `zoom_d_min/max` defaults, adjusting the camera-delta gating)
  should land in *both* `chimerax_depth_tube.py` and
  `pymol_depth_tube.py` in the same PR so the two viewers don't drift.
- **Defaults are tuned together.** The tuning table in `README.md`
  must match the function signatures of `enable()` in both files.
- **No external dependencies** beyond what the host viewer already
  ships with (NumPy on both, Qt on both). Don't introduce a `requirements.txt`.
- **Self-contained modules.** Each `.py` is loadable on its own; no
  shared utility module. Duplication between the two files is fine
  and expected — they target different APIs.
- **Comments**: only when the *why* is non-obvious. Skip docstring
  paragraphs that just describe what the function does.
- **Logging** goes through `session.logger` (ChimeraX) or `print()`
  (PyMOL) — both viewers route those to their consoles.

## Testing

There is no automated test suite. To smoke-test changes:

```bash
# ChimeraX (from the ChimeraX command line, not bash):
runscript /path/to/chimerax_depth_tube.py
open 1bvk
cartoon
depthtube #1/A
# Rotate the model — radius should track depth, idle holds steady.
~depthtube

# PyMOL (from the PyMOL command line):
run /path/to/pymol_depth_tube.py
fetch 1bvk, async=0
show cartoon
depthtube chain A and polymer
# Rotate — putty width should track depth.
depthtube_off
```

For syntax-only checks without launching a viewer:

```bash
python -m py_compile chimerax_depth_tube.py pymol_depth_tube.py
```

## Common change patterns

| Change                                  | Touch                                                                 |
|-----------------------------------------|-----------------------------------------------------------------------|
| Tune default radius / zoom              | function signature of `enable()` in *both* `.py` + README tuning table |
| Change tick rate                        | `QTimer.setInterval` in both files                                    |
| Add a new render param                  | add to `enable()` signature, expose via the registered command in both |
| Improve PyMOL putty smoothness          | `cartoon_putty_quality` / `cartoon_putty_radius` in `pymol_depth_tube.py` |
| Add a new viewer port (e.g. Mol*)       | new top-level file `<viewer>_depth_tube.<ext>`; keep parity            |

## License

MIT, see `LICENSE`. PRs implicitly under the same license unless the
PR description says otherwise.
