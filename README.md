GBB | Gradient-Based Boundary | not finished
===

[![Python](https://img.shields.io/badge/Python-3-blue)](https://github.com/haenelt/GBB)
[![License](https://img.shields.io/badge/license-GPL--3.0-orange)](https://github.com/haenelt/GBB)

<p align="center">
  <img src="https://github.com/haenelt/GBB/blob/master/gbb_logo.gif?raw=true" width=75% height=75% alt="Illustration of GBB"/>
</p>

The gradient-based boundary surface mesh refinement is based on the idea used in boundary-based registration (BBR), which is a method for robust alignment of images with different image contrasts [[1]](#1). In BBR, the clear GM/WM border of the target image is used to find the transformation of the source image which maximizes its contrast along the tissue boundary of the target image. This method only finds a global rigid-body transformation between both images. Different distortions between source and target image can therefore not adequately be taken into account. The proposed method turns the idea of BBR around and deforms the surface boundary of the target image locally to maximize the tissue contrast around the surface boundary in the source image.

## Methods
We assume that we have a (distorted) mean epi image and a surface mesh from a (non-distorted) anatomy in the same space after applying a rigid-body transformation to the GM/WM surface boundary based on another registration method (see example data). We now randomly select single vertices in an iterative procedure. In each step, we calculate the gradient along defined direction (along one axis, e.g. phase-encoding direction or along the surface normal) around the selected vertex point. The vertex is then shifted a small amount (learning rate) towards the point of maximum GM/WM contrast. Neighboring vertices are also shifted depending on their distance to the considered vertex point (neighborhood size). During the iterative procedure, we shrink the neighborhood size and the learning rate. A vein mask should also be applied to prevent the surface to fall into large intenstiy holes. A cost function is computed equivalent to the function used in the BBR method. The procedure stops when the maximum of iterations is reached or the cost function saturates.

To improve the starting mesh for the GBB method, two processes can optionally be applied beforehand, called *deveining* and *anchoring*. In deveining, all surface mesh vertex points are pulled out of masked veins. In anchoring, anchor points are manually set and the surface mesh is locally shifted to these points. Vertices now located within these anchor points are locked and not moved during the GBB method.

## Installation
- freesurfer dependency

### Stable release (PyPI)

GBB can be installed via the `pip` command. I recommend to use `Anaconda to create a new python environment with Python >= 3 and then simply run the following line from a terminal with the environment being activated:
```
pip install gbb
```

### Development version (github)

Alternatively, it is possible to clone this repository and run one of the following lines from the directory in which the repository was cloned:
```
python setup.py install
```
or
```
pip install .[<options>]
```

## Example data
A time series mean from a GE-EPI acquisition is used as target image. A binary vein mask and surface meshs are given as well. Surface meshs were computed from a separete MP2RAGE acquisition using FreeSurfer and were rigidly aligned to the epi image.

## References
<a id="1">[1]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63&ndash;72 (2009).

## Contact
If you have questions, problems or suggestions regarding the GBB package, please feel free to contact [me](<daniel.haenelt@gmail.com>).

# TODO
- general information and purpose
- URL of the main source of the software
- freesurfer dependency (and other?)
- basic credit information
- howto install
- Example
- contact
