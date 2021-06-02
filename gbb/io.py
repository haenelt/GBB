import nibabel as nb
import numpy as np
import pandas as pd
import csv
import os



def load_volume(volume):
    """
    Load volumetric data into a
    `Nibabel SpatialImage <http://nipy.org/nibabel/reference/nibabel.spatialimages.html#nibabel.spatialimages.SpatialImage>`_
    Parameters
    ----------
    volume: niimg
        Volumetric data to be loaded, can be a path to a file that nibabel can
        load, or a Nibabel SpatialImage
    Returns
    ----------
    image: Nibabel SpatialImage
    Notes
    ----------
    Originally created as part of Laminar Python [1]_ .
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    """  # noqa

    # if input is a filename, try to load it
    # python 2 version if isinstance(volume, basestring):
    if isinstance(volume, str):
        # importing nifti files
        image = nb.load(volume)
    # if volume is already a nibabel object
    elif isinstance(volume, nb.spatialimages.SpatialImage):
        image = volume
    else:
        raise ValueError('Input volume must be a either a path to a file in a '
                         'format that Nibabel can load, or a nibabel'
                         'SpatialImage.')
    return image


def save_volume(filename, volume, dtype='float32', overwrite_file=True):
    """
    Save volumetric data that is a
    `Nibabel SpatialImage <http://nipy.org/nibabel/reference/nibabel.spatialimages.html#nibabel.spatialimages.SpatialImage>`_
    to a file
    Parameters
    ----------
    filename: str
        Full path and filename under which volume should be saved. The
        extension determines the file format (must be supported by Nibabel)
    volume: Nibabel SpatialImage
        Volumetric data to be saved
    dtype: str, optional
        Datatype in which volumetric data should be stored (default is float32)
    overwrite_file: bool, optional
        Overwrite existing files (default is True)
    Notes
    ----------
    Originally created as part of Laminar Python [1]_ .
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    """  # noqa
    if dtype is not None:
        volume.set_data_dtype(dtype)
    if os.path.isfile(filename) and overwrite_file is False:
        print("\nThis file exists and overwrite_file was set to False, "
              "file not saved.")
    else:
        try:
            volume.to_filename(filename)
            print("\nSaving {0}".format(filename))
        except AttributeError:
            print('\nInput volume must be a Nibabel SpatialImage.')



















def load_mesh(surf_mesh):
    '''
    Load a mesh into a dictionary with entries
    "points", "faces" and "data"
    Parameters
    ----------
    surf_mesh:
        Mesh to be loaded, can be a path to a file
        (currently supported formats are freesurfer geometry formats,
        gii and ASCII-coded vtk, ply or obj) or a dictionary with the
        keys "points", "faces" and (optionally) "data"
    Returns
    ----------
    dict
        Dictionary with a numpy array with key "points" for a Numpy array of
        the x-y-z coordinates of the mesh vertices and key "faces" for a
        Numpy array of the the indices (into points) of the mesh faces.
        Optional "data" key is a Numpy array of values sampled on the "points".
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''

    if surf_mesh.endswith('vtk'):
        points, faces, data = _read_vtk(surf_mesh)
        return {'points': points, 'faces': faces, 'data': data}
    else:
        geom = load_mesh_geometry(surf_mesh)
        return geom


def save_mesh(filename, surf_dict):
    '''
    Saves surface mesh to file
    Parameters
    ----------
    filename: str
        Full path and filename under which surfaces data should be saved. The
        extension determines the file format. Currently supported are
        freesurfer geometry formats, gii and ASCII-coded vtk, obj, ply. Note
        that only ASCII-coded vtk currently saves data, the others only save
        the geometry.
    surf_dict: dict
        Surface mesh geometry to be saved. Dictionary with a numpy array with
        key "points" for a Numpy array of the x-y-z coordinates of the mesh
        vertices and key "faces" for a Numpy array of the the indices
        (into points) of the mesh faces. Optional "data" key is a Numpy array
        of values sampled on the "points"
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''
    if filename.endswith('vtk'):
        _write_vtk(filename,
                   surf_dict['points'],
                   surf_dict['faces'],
                   surf_dict['data'])
    else:
        save_mesh_geometry(filename,
                           surf_dict)


def load_mesh_geometry(surf_mesh):
    '''
    Load a mesh geometry into a dictionary with entries
    "points" and "faces"
    Parameters
    ----------
    surf_mesh:
        Mesh geometry to be loaded, can be a path to a file
        (currently supported formats are freesurfer geometry formats,
        gii and ASCII-coded vtk, ply or obj) or a dictionary with the
        keys "points" and "faces"
    Returns
    ----------
    dict
        Dictionary with a numpy array with key "points" for a Numpy array of
        the x-y-z coordinates of the mesh vertices and key "faces" for a
        Numpy array of the the indices (into points) of the mesh faces
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''
    # if input is a filename, try to load it with nibabel
    if isinstance(surf_mesh, str):
        if (surf_mesh.endswith('orig') or surf_mesh.endswith('pial') or
                surf_mesh.endswith('white') or surf_mesh.endswith('sphere') or
                surf_mesh.endswith('inflated')):
            points, faces = nb.freesurfer.io.read_geometry(surf_mesh)
        elif surf_mesh.endswith('vtk'):
            points, faces, _ = _read_vtk(surf_mesh)
        else:
            raise ValueError('Currently supported file formats are freesurfer '
                             'geometry formats and gii, vtk, ply, obj')
    elif isinstance(surf_mesh, dict):
        if ('faces' in surf_mesh and 'points' in surf_mesh):
            points, faces = surf_mesh['points'], surf_mesh['faces']
        else:
            raise ValueError('If surf_mesh is given as a dictionary it '
                             'must contain items with keys "points" and '
                             '"faces"')
    else:
        raise ValueError('Input surf_mesh must be a either filename or a '
                         'dictionary containing items with keys "points" '
                         'and "faces"')
    return {'points': points, 'faces': faces}


def load_mesh_data(surf_data):
    '''
    Loads mesh data into a Numpy array
    Parameters
    ----------
    surf_data:
        Mesh data to be loaded, can be a Numpy array or a path to a file.
        Currently supported formats are freesurfer data formats (mgz, curv,
        sulc, thickness, annot, label), nii, gii, ASCII-coded vtk and txt
    Returns
    ----------
    np.ndarray
        Numpy array containing the data
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''
    # if the input is a filename, load it
    if isinstance(surf_data, str):
        if (surf_data.endswith('nii') or surf_data.endswith('nii.gz') or
                surf_data.endswith('mgz')):
            data = np.squeeze(nb.load(surf_data).get_data())
        elif (surf_data.endswith('curv') or surf_data.endswith('sulc') or
                surf_data.endswith('thickness')):
            data = nb.freesurfer.io.read_morph_data(surf_data)
        elif surf_data.endswith('annot'):
            data = nb.freesurfer.io.read_annot(surf_data)[0]
        elif surf_data.endswith('label'):
            data = nb.freesurfer.io.read_label(surf_data)
        # check if this works with multiple indices (if dim(data)>1)
        elif surf_data.endswith('vtk'):
            _, _, data = _read_vtk(surf_data)
        elif surf_data.endswith('txt'):
            data = np.loadtxt(surf_data)
        else:
            raise ValueError('Format of data file not recognized. Currently '
                             'supported formats are freesurfer data formats '
                             '(mgz, sulc, curv, thickness, annot, label)'
                             'nii', 'gii, ASCII-coded vtk and txt')
    elif isinstance(surf_data, np.ndarray):
        data = np.squeeze(surf_data)
    return data


def save_mesh_data(filename, surf_data):
    '''
    Saves surface data that is a Numpy array to file
    Parameters
    ----------
    filename: str
        Full path and filename under which surfaces data should be saved. The
        extension determines the file format. Currently supported are
        freesurfer formats curv, thickness, sulc and ASCII-coded txt'
    surf_data: np.ndarray
        Surface data to be saved
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''
    if isinstance(filename, str) and isinstance(surf_data, np.ndarray):
        if (filename.endswith('curv') or filename.endswith('thickness') or
                filename.endswith('sulc')):
            nb.freesurfer.io.write_morph_data(filename, surf_data)
            print("\nSaving {0}".format(filename))
        elif filename.endswith('txt'):
            np.savetxt(filename, surf_data)
            print("\nSaving {0}".format(filename))
        else:
            raise ValueError('File format not recognized. Currently supported '
                             'are freesurfer formats curv, sulc, thickness '
                             'and ASCII coded vtk and txt')
    else:
        raise ValueError('Filename must be a string')


def save_mesh_geometry(filename, surf_dict):
    '''
    Saves surface mesh geometry to file
    Parameters
    ----------
    filename: str
        Full path and filename under which surfaces data should be saved. The
        extension determines the file format. Currently supported are
        freesurfer geometry formats, gii and ASCII-coded vtk, obj, ply'
    surf_dict: dict
        Surface mesh geometry to be saved. Dictionary with a numpy array with
        key "points" for a Numpy array of the x-y-z coordinates of the mesh
        vertices and key "faces2 for a Numpy array of the the indices
        (into points) of the mesh faces
    Notes
    ----------
    Originally created as part of Laminar Python [1]_
    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    '''
    if isinstance(filename, str) and isinstance(surf_dict, dict):
        if (filename.endswith('orig') or filename.endswith('pial') or
                filename.endswith('white') or filename.endswith('sphere') or
                filename.endswith('inflated')):
            nb.freesurfer.io.write_geometry(filename, surf_dict['points'],
                                            surf_dict['faces'])
            print("\nSaving {0}".format(filename))
        elif filename.endswith('vtk'):
            if 'data' in surf_dict.keys():
                _write_vtk(filename, surf_dict['points'], surf_dict['faces'],
                           surf_dict['data'])
                print("\nSaving {0}".format(filename))
            else:
                _write_vtk(filename, surf_dict['points'], surf_dict['faces'])
                print("\nSaving {0}".format(filename))
    else:
        raise ValueError('Filename must be a string and surf_dict must be a '
                         'dictionary with keys "points" and "faces"')


# function to read vtk files
# ideally use pyvtk, but it didn't work for our data, look into why
def _read_vtk(file):
    '''
    Reads ASCII coded vtk files using pandas,
    returning vertices, faces and data as three numpy arrays.
    '''
    # read full file while dropping empty lines
    try:
        vtk_df = pd.read_csv(file, header=None, engine='python')
    except csv.Error:
        raise ValueError(
            'This vtk file appears to be binary coded currently only ASCII '
            'coded vtk files can be read')
    vtk_df = vtk_df.dropna()
    # extract number of vertices and faces
    number_vertices = int(vtk_df[vtk_df[0].str.contains(
                                            'POINTS')][0].iloc[0].split()[1])
    number_faces = int(vtk_df[vtk_df[0].str.contains(
                                            'POLYGONS')][0].iloc[0].split()[1])
    # read vertices into df and array
    start_vertices = (vtk_df[vtk_df[0].str.contains(
                                            'POINTS')].index.tolist()[0]) + 1
    vertex_df = pd.read_csv(file, skiprows=range(start_vertices),
                            nrows=number_vertices, delim_whitespace=True,
                            header=None, engine='python')
    if np.array(vertex_df).shape[1] == 3:
        vertex_array = np.array(vertex_df)
    # sometimes the vtk format is weird with 9 indices per line,
    # then it has to be reshaped
    elif np.array(vertex_df).shape[1] == 9:
        vertex_df = pd.read_csv(file, skiprows=range(start_vertices),
                                nrows=int(number_vertices / 3) + 1,
                                delim_whitespace=True, header=None,
                                engine='python')
        vertex_array = np.array(vertex_df.iloc[0:1, 0:3])
        vertex_array = np.append(vertex_array, vertex_df.iloc[0:1, 3:6],
                                 axis=0)
        vertex_array = np.append(vertex_array, vertex_df.iloc[0:1, 6:9],
                                 axis=0)
        for row in range(1, (int(number_vertices / 3) + 1)):
            for col in [0, 3, 6]:
                vertex_array = np.append(vertex_array, np.array(
                    vertex_df.iloc[row:(row + 1), col:(col + 3)]), axis=0)
        # strip rows containing nans
        vertex_array = vertex_array[~np.isnan(vertex_array)].reshape(
                                                            number_vertices, 3)
    else:
        print("vertex indices out of shape")
    # read faces into df and array
    start_faces = (vtk_df[vtk_df[0].str.contains(
                                            'POLYGONS')].index.tolist()[0]) + 1
    face_df = pd.read_csv(file, skiprows=range(start_faces),
                          nrows=number_faces, delim_whitespace=True,
                          header=None, engine='python')
    face_array = np.array(face_df.iloc[:, 1:4])
    # read data into df and array if exists
    if vtk_df[vtk_df[0].str.contains('POINT_DATA')].index.tolist() != []:
        start_data = (vtk_df[vtk_df[0].str.contains(
                                        'POINT_DATA')].index.tolist()[0]) + 3
        number_data = number_vertices
        data_df = pd.read_csv(file, skiprows=range(start_data),
                              nrows=number_data, delim_whitespace=True,
                              header=None, engine='python')
        data_array = np.array(data_df)
    else:
        data_array = None

    return vertex_array, face_array, data_array


def _write_vtk(filename, vertices, faces, data=None, comment=None):
    '''
    Creates ASCII coded vtk file from numpy arrays using pandas.
    Inputs:
    -------
    (mandatory)
    * filename: str, path to location where vtk file should be stored
    * vertices: numpy array with vertex coordinates,  shape (n_vertices, 3)
    * faces: numpy array with face specifications, shape (n_faces, 3)
    (optional)
    * data: numpy array with data points, shape (n_vertices, n_datapoints)
        NOTE: n_datapoints can be =1 but cannot be skipped (n_vertices,)
    * comment: str, is written into the comment section of the vtk file
    Usage:
    ---------------------
    _write_vtk('/path/to/vtk/file.vtk', v_array, f_array)
    '''
    # infer number of vertices and faces
    number_vertices = vertices.shape[0]
    number_faces = faces.shape[0]
    if data is not None:
        number_data = data.shape[0]
    # make header and subheader dataframe
    header = ['# vtk DataFile Version 3.0',
              '%s' % comment,
              'ASCII',
              'DATASET POLYDATA',
              'POINTS %i float' % number_vertices
              ]
    header_df = pd.DataFrame(header)
    sub_header = ['POLYGONS %i %i' % (number_faces, 4 * number_faces)]
    sub_header_df = pd.DataFrame(sub_header)
    # make dataframe from vertices
    vertex_df = pd.DataFrame(vertices)
    # make dataframe from faces, appending first row of 3's
    # (indicating the polygons are triangles)
    triangles = np.reshape(3 * (np.ones(number_faces)), (number_faces, 1))
    triangles = triangles.astype(int)
    faces = faces.astype(int)
    faces_df = pd.DataFrame(np.concatenate((triangles, faces), axis=1))
    # write dfs to csv
    header_df.to_csv(filename, header=None, index=False)
    with open(filename, 'a') as f:
        vertex_df.to_csv(f, header=False, index=False, float_format='%.3f',
                         sep=' ')
    with open(filename, 'a') as f:
        sub_header_df.to_csv(f, header=False, index=False)
    with open(filename, 'a') as f:
        faces_df.to_csv(f, header=False, index=False, float_format='%.0f',
                        sep=' ')
    # if there is data append second subheader and data
    if data is not None:
        if len(data.shape)>1:
            datapoints = data.shape[1]
            sub_header2 = ['POINT_DATA %i' % (number_data),
                           'SCALARS Scalars float %i' % (datapoints),
                           'LOOKUP_TABLE default']
        else:
            datapoints = 1
            sub_header2 = ['POINT_DATA %i' % (number_data),
                           'SCALARS Scalars float',
                           'LOOKUP_TABLE default']
        sub_header_df2 = pd.DataFrame(sub_header2)
        data_df = pd.DataFrame(data)
        with open(filename, 'a') as f:
            sub_header_df2.to_csv(f, header=False, index=False)
        with open(filename, 'a') as f:
            data_df.to_csv(f, header=False, index=False, float_format='%.16f',
                           sep=' ')
