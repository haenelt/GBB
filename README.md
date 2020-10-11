GBB | Gradient-Based Boundary
===

[![Python](https://img.shields.io/badge/Python-3.6%7C3.7%7C3.8-blue)](https://github.com/haenelt/GBB)
[![License](https://img.shields.io/badge/license-GPL--3.0-orange)](https://github.com/haenelt/GBB)

<p align="center">
  <img src="https://github.com/haenelt/GBB/blob/master/gbb_logo.gif?raw=true" width=75% height=75% alt="Illustration of GBB"/>
</p>

High-resolution fMRI offers the great potential to differentiate signal processing within the cerebral cortex. Since the cortex is organized into separate layers [[1]](#1) and these layers differ in their respective connectivity pattern [[2]](#2), this promises the gain of new insights into the functioning of the living human brain. However, precise sampling of fMRI data within the cortex is still challenging because of limitations in segmentation and registration. This package proposes a method to improve the alignment between the position of a gray matter/white matter (GM/WM) boundary surface mesh from a separate segmentation and the same border found in the functional image.

The gradient-based boundary surface mesh refinement is based on the idea used in boundary-based registration (BBR), which is a method for robust alignment of images with different image contrasts [[3]](#3). In BBR, the clear GM/WM border of the target image is used to find the transformation of the source image which maximizes the source image contrast along the tissue boundary of the target image. This method only finds a global rigid-body transformation between both images. Different distortions between source and target image can therefore not adequately be taken into account. The proposed method turns the idea of BBR around and deforms the surface boundary of the target image locally to maximize the tissue contrast around the surface boundary in the source image.

## Methods
We assume that we have a distorted image from an fMRI scan and a cortical boundary surface mesh from a non-distorted anatomical MR scan (see [example\_data](https://github.com/haenelt/GBB/tree/master/example_data)). We further assume that both are in the same space without being perfectly aligned due to inaccuracies in segmentation and registration. The method starts by randomly selecting one vertex. For the selected vertex, the gradient along a defined direction (along one axis, e.g. phase-encoding direction, or along the surface normal) is calculated. The vertex is then shifted a small amount (learning rate) towards the point of maximum GM/WM contrast. Neighboring vertices are also shifted depending on their distance to the considered vertex (neighborhood size). A new iteration starts by randomly selecting a new vertex. During an iterative procedure, neighborhood size and learning rate can be changed. A vein mask should be used as additional input to prevent the surface to fall into large intenstiy holes in T2*-weighted functional images. A cost function is computed equivalent to the function defined in the BBR method. The procedure stops when the maximum number of iterations is reached or the cost function saturates.

An alternative explanation can be found in this [blog post](https://haenelt.github.io/blog/blog_1/).

To improve the initial mesh for the GBB method, two processes can optionally be applied beforehand, called *deveining* and *anchoring*. In deveining, all vertices are pulled out of masked veins. In anchoring, control points are manually set and the surface mesh is locally shifted to those points. Anchored vertices are locked and not moved during the GBB method.

## Installation

### Stable release (PyPI)
GBB can be installed via the `pip` command. I recommend to use `Anaconda` to create a new python environment with `Python >= 3.6`. Then, run the following line from a terminal with the environment being activated:

```
pip install gbb
```

### Development version (github)
Alternatively, it is possible to clone this repository and run the following line from the directory in which the repository was cloned:

```
python setup.py install
```

### Other dependencies (FreeSurfer)
The package needs a FreeSurfer installation on the used machine [[4]](#4). The package was tested with version 6.0.0p1.

## Run the package
After installation, you should be able to import the package with `import gbb`. The main processing steps can be run by calling the function `main.py` in the package. Alternatively, they can be run from the command line with `python -m gbb <args>`. The main function uses a list of parameters which are stored in a configuration file which can be found [here](https://github.com/haenelt/GBB/tree/master/gbb/config). E.g., you can define in this file if you want to run deveining or not. To change parameters, I recommend that you save a copy of the configuration file and apply all changes there. You can provide the filename of this customized configuration file if you run the package. If no file is provided, the default configuration file will be read. If you run the main module, a JSON file will be created which stores all the used parameters.

## Example data
Data is provided for testing purposes and to illustrate the data format which is needed to run the package. You can find the data [here](https://github.com/haenelt/GBB/tree/master/example_data).

## References
<a id="1">[1]</a> Brodmann K, Vergleichende Lokalisationslehre der Grosshirnrinde, Barth-Verlag (1909).<br/>
<a id="2">[2]</a> Felleman DJ, Van Essen DC, Distributed hierarchical processing in the primate cerebral cortex, Cerebral Cortex 1(1), 1&ndash;47 (1991).<br/>
<a id="2">[3]</a> Greve DN, Fischl B, Accurate and robust brain image alignment using boundary-based registration, Neuroimage 48(1), 63&ndash;72 (2009).<br/>
<a id="3">[4]</a> https://surfer.nmr.mgh.harvard.edu/ (accessed 05-10-2020).

## Contact
If you have questions, problems or suggestions regarding the GBB package, please feel free to contact [me](mailto:daniel.haenelt@gmail.com).
