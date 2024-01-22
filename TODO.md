## general
- [ ] new logo
- [ ] include docs (sphinx)
- [ ] clean example_data
- [ ] example_data -> osf?
- [ ] example_data -> test_data?
- [ ] update .gitignore
- [ ] update CHANGELOG.md
- [ ] update README.md
- [ ] installation with poetry
- [ ] update requirements.txt
- [ ] update VERSION
- [ ] new module structure
- [ ] add example scripts
- [ ] add tests
- [ ] add docker
- [ ] docs on readthedocs
- [ ] github: about, release, package, workflow

## modules
- [ ] use artificial data set (sphere and cube; add noise)
- [ ] include spatial filters: gaussian filter, median filter (edge-preserving), bilateral filter (edge-preserving)
- [ ] VOX and FREESURFER in config
- [ ] print info in io which coordinate system is used
- [ ] identification of “bad vertex regions” by intensity histogram of BBR cost percentiles
- [ ] learning rate like in DNNs: step size dependent on slope
- [ ] include regularization: less change for high change rate or something like that (to prevent that the mesh gets to messy)
- [ ] include regularization to be robust against veins (mean in t2* time series of std phase)
- [ ] make own deformation field (natural neighbor interpolation)
- [ ] check preservation of smoothness
- [ ] check preservation of topology (prevent mesh to intersect)
- [ ] include own vox2ras to be independent to freesurfer
- [ ] after alignment, rearrange vertices to homogenesouly cover the surface
- [ ] induce deformation field to artificially deform mesh (for testing purposes)
- [ ] include api for BBR, rBBR and fieldmap undistortion
- [ ] include api for initial registration
- [ ] include api for phase unwrapping
- [ ] include code for preprocessing
- [ ] weight vertices with probablity function based on distance to zero crossing
- [ ] weight shift direction with normal direction (less in normal direction at later stages)
- [ ] include preprocessing using old_preprocessing.py

## Kriegeskorte2019
- gradient descent: we adjust each weight in the direction that reduces the cost (the error) and by an amount proportional to the derivative of the cost with respect to the weight.
- batch size: in practice, the best solution is to use small batches of training examples to estimate the gradient before taking a step. Compared to the single-example approach, this gives us a more stable sense of direction.

```
.
└── GBB/
    ├── scripts/
    ├── gbb/
    │   ├── __init__
    │   ├── __main__
    │   ├── config
    │   ├── data
    │   ├── prepare
    │   ├── io
    │   ├── meshutils
    │   ├── meshdeform
    │   ├── testvolume
    │   ├── utils
    │   ├── plot
    │   ├── interpolation(?)
    │   └── neighbor(?)
    ├── docs/
    ├── tests/
    ├── .gitignore
    ├── CHANGELOG.md
    ├── LICENSE
    ├── README.md
    ├── poetry.lock
    ├── requirements.txt
    └── pyproject.toml
```