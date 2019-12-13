def get_normal_fsurf(input_surf):
    """
    This function computes the vertex-wise surface normals of an input surface using a freesurfer
    function. An ascii surface file is written as output.
    Inputs:
        *input_surf: filename of input surface.
        
    created by Daniel Haenelt
    Date created: 13-12-2019
    Last modified: 13-12-2019
    """
    import nipype.interfaces.freesurfer as fs
    
    # get surface normas
    mris = fs.MRIsConvert()
    mris.inputs.in_file = input_surf
    mris.inputs.out_file = input_surf+"_normal.asc"
    mris.inputs.normal = True
    mris.run()