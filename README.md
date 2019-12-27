GBB | Gradient-Based Boundary
===

- registration between anatomy and functional image is difficult
- especially for high-resolution applications
- BBR registration shows robust registration between different contrasts [[1]](#1).
- similar approach
- assume surface and image in same space (e.g. after first rigid alignment)
- vertices are randomly chosen
- position of vertex and neighborhood is updated with (shifted towards the GM/WM) boundary (highest contrast)
- contrast is only checked along one axis
- this should coincide to phase-encoding direction
- especially, in t2* images, signal dropouts can lead to falses boundary detections
- therefore, a vein mask should be given as well
- so, we avoid shifting vertices into veins
- a cost function is computed equivalent to the BBR paper
- refinement stops if maximum of iterations are reached
- or if cost function reaches steady-state
- exemplary data can is given

## Prerequisites
- Freesurfer should be included in the search path
- used python packages: os, sys, shutil, subprocess, numpy, nibabel, matplotlib, scipy, cv2
- one function from my scripts repository

## Explanation of input parameters
```python
# input files
input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
path_output = "/home/daniel/Schreibtisch/parameters1"
name_output = "lh.layer10_refined"

# parameters
t2s = True # underlying image contrast
line_dir = 2 # line axis in ras convention
line_length = 3 # line length in one direction in mm
r_size = [5, 2.5, 1] # neighborhood radius in mm
l_rate = [0.1, 0.1, 0.1] # learning rate
max_iterations = [100000, 250000, 500000] # maximum iterations
cost_threshold = [1e-3,5e-4,1e-4] # cost function threshold
cleanup = True

# gradient preparation
sigma = 1 # gaussian filter
kernel_size = 3 # kernel size used by gradient calculation

# output
show_cost = True # show temporary cost function
write_gradient = True # write gradient image
write_step = 1000 # step size to write intermediate surfaces (if set > 0)
```

- number of steps should be larger than number of vertices -> hard to achieve
- for approximately the first 1000 steps, learning rate should start with a value that is close to unity -> we do not do that because we have already a good starting point
- normally, learning rate is decreased -> we keep the same small learning rate over all iterations
- we decrease neighborhood size over time (three steps; couplint strength between neighbors)
- at the end: neighborhood size and learning rate is small


## References
<a id="1">[1]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63&ndash;72 (2009).