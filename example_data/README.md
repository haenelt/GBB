GBB/example_data
===

The folder contains data files which can be used as input arguments for testing purposes:

#### mean\_epi.nii and mean\_epi\_enhanced.nii
mean\_epi.nii or mean\_epi\_enhanced.nii can be set as reference volume (argument `-r`). The image is the arithmetic mean of an fMRI time series. A coronal slab covering the occipital lobe of a human volunteer was imaged with an isotropic voxel size of 0.8 mm (single-shot 2D GE-EPI @ 7 T). It is recommended to use mean\_epi\_enhanced.nii as reference volume which shows an improved GM/WM contrast. Contrast enhancement was achieved by weighting the magnitude image by its phase. The explicit procedure will be explained elsewhere.

#### lh.white
The GM/WM boundary surface mesh was computed using FreeSurfer [[1]](#1) based on an acquired anatomical whole-brain MR scan (MP2RAGE @ 7 T). The resulting GM/WM boundary surface mesh was linearly transformed to mean\_epi.nii. This surface should be used as white surface (argument `-w`).

#### lh.pial
THE GM/CSF boundary surface can be used as additional surface mesh (argument `-p`). The surface resulted from the same segmentation as explained for `lh.white`. If set, the final deformation field will be applied to that surface as well.

#### vein.nii
Optionally, a vein mask can be set (argument `-v`). This mask is mandatory if `devein_params["run"] = True` in the configuration file `config.py`. It will also be used if reg_params["run"] = True` but is not mandatory here. The binary mask was computed by thresholding the mean_epi.nii and a corresponding tSNR image.

#### control_points.dat
For anchoring, control points have to be set (argument `-a`). The file was generated using the FreeSurfer [[1]](#1) viewer FreeView [[2]](#2).

#### ignore.nii
Optionally, a binary mask can be used (argument `-i`). If set, all masked voxels will be excluded from further processing.

## References
<a id="1">[1]</a> https://surfer.nmr.mgh.harvard.edu/ (version 6.0.0p1)<br/>
<a id="2">[2]</a> https://surfer.nmr.mgh.harvard.edu/fswiki/FreeviewGuide/FreeviewIntroduction
