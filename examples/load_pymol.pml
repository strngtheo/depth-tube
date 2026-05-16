# Example: PyMOL depth-tube on chain A.
# Edit the paths before running.

fetch 1bvk, async=0
hide everything
show cartoon
color spectrum, chain A
bg_color white
set ray_shadows, 0

# Load the script (registers `depthtube` and `depthtube_off`).
run ~/repos/depth-tube/pymol_depth_tube.py

# Enable on chain A.
depthtube chain A and polymer

# Custom parameters (named-arg style):
# depthtube chain B and polymer, r_proximal=0.55, r_distal=0.06, zoom_d_max=200

# Disable everything:
# depthtube_off
