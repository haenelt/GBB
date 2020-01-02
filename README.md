GBB | Gradient-Based Boundary
===

The gradient-based boundary surface mesh refinement is based on the idea used in boundary-based registration (BBR) which is a method for robust alignment of different image contrasts [[1]](#1). In BBR, the clear GM/WM border of the target image is used to find the transformation of the source image which maximizes its contrast along the tissue boundary of the target image. This method only finds a global rigid-body transformation between both images. Different distortions between source and target image can therefore not adequately taken into account. The proposed method turns the idea of BBR around and locally deforms the surface boundary of the target image to maximize the tissue contrast in the source image.

## Method
We assume that we have a (distorted) mean epi image and a surface mesh from a (non-distorted) anatomy in the same space after applying a rigid-body transformation to the GM/WM surface boundary based on another registration method (see example data). We now randomly select single vertices in an iterative procedure. In each step, we calculate the gradient along the phase-encoding direction of the epi acquisition around the selected vertex point. The vertex is then shifted a small amount (learning rate) towards the point of maximum GM/WM contrast. Neighboring vertices are also shifted depending on their distance to the considered vertex point. During the iterative procedure, we shrink the neighborhood size and the learning rate. A vein mask should also be applied to prevent the surface to fall into large intenstiy holes. A cost function is computed equivalent to the function used in the BBR method. The procedure stops if the maximum of iterations is reached or if the cost function saturates.

## Prerequisites
- Freesurfer should be included in the search path
- used python packages: os, sys, shutil, subprocess, numpy, nibabel, matplotlib, scipy, cv2
- one function from my scripts repository is used

## Example data
A time series mean from a GE-EPI acquisition is used as target image. A binary vein mask and surface meshs are given as well. Surface meshs were computed from a separete MP2RAGE acquisition using FreeSurfer and were rigidly aligned to the epi image. To see the effect of the proposed method, a functional contrast image from an visual experiment can be found as well.

## References
<a id="1">[1]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63&ndash;72 (2009).