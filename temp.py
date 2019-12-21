import random
import nibabel as nb
from lib_gbb.utils import get_gradient
from nibabel.freesurfer.io import read_geometry
from lib_gbb.normal import get_normal_direction

# get current vertex point
n_coords = np.arange(0,len(vtx),1)
n_vertex = random.choice(n_coords)

#%%

input_geometry = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_ref = "/home/daniel/projects/GBB/test_data/mean_data.nii"
input_vein = "/home/daniel/projects/GBB/test_data/vein.nii"
ind = 27000
mm_line = 3 # in mm
n_line = 1000
show_plot = True
delta_dir = 2
vol_xmax = 100
vol_ymax = 100
vol_zmax = 100

# get transformation
vox2ras, ras2vox = vox2ras(input_ref)

# get vol and gradient
vol_array = nb.load(input_ref).get_fdata()
grad_array = get_gradient(input_ref, input_vein, sigma_gaussian=1, grad_dir=1, kernel_size=3)

# get surf
vtx, fac = read_geometry(input_geometry)

# get normals
n = get_normal_direction(vtx, fac, diff_dir