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
x += 2
```

## References
<a id="1">[1]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63--72 (2009).