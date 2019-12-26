GBB | Gradient-Based Boundary
===

- general description

## Main idea
- choose randomly one vertex
- get gradient along one direction
- update neighborhood
- update cost function [[1]](#1).
- stop after maximum of iterations or cost function saturates

### something which should be explained
- number of steps should be larger than number of vertices -> hard to achieve
- for approximately the first 1000 steps, learning rate should start with a value that is close to unity -> we do not do that because we have already a good starting point
- normally, learning rate is decreased -> we keep the same small learning rate over all iterations
- we decrease neighborhood size over time (three steps; couplint strength between neighbors)
- at the end: neighborhood size and learning rate is small

## Prerequisites
- Freesurfer should be included in the search path
- python packages: os, sys, subprocess, random, numpy, nibabel, matplotlib, scipy, cv2
- one function from my scripts repository

## Important parameters
```python
input_surf = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
path_output = "/home/daniel/Schreibtisch/parameters13"
name_output = "lh.layer10_refined"

# parameters
t2s = True # underlying image contrast
line_dir = 2 # line axis in ras convention
line_length = 3 # line length in one direction in mm
r_size = [5, 2.5, 1] # neighborhood radius in mm
l_rate = [0.1, 0.1, 0.1] # learning rate
max_iterations = [50000, 50000, 50000] # maximum iterations
cost_threshold = [1e-8,1e-8,1e-8] # cost function threshold

# gradient preparation
sigma = 1
kernel_size = 3

# output
show_line = False
show_cost = True
write_gradient = True
write_intermediate = True
write_step = 1000
```

## References
<a id="1">[1]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63&ndash;72 (2009).