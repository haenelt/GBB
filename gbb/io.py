# -*- coding: utf-8 -*-
"""I/O functions."""

# python standard library inputs
import os
import csv

# external inputs
import numpy as np
import nibabel as nb
import pandas as pd

__all__ = ['get_filename', 'load_volume', 'save_volume']


def get_filename(file_in):
    """Get path, basename and file extension for an input file. The loop checks
    for some given extension names (nii, mgh, mgz). Otherwise it stops after the
    last found dot in the string.

    Parameters
    ----------
    file_in : str
        Filename.

    Returns
    -------
    path_file : str
        Path of filename.
    name_file : str
        Basename of filename.
    ext_file : str
        File extension of filename.

    """

    # get path and basename
    path_file = os.path.dirname(file_in)
    name_file = os.path.basename(file_in)

    # split basename and file extension
    ext_file = ""
    exit_loop = 0
    ext_key = [".nii", ".mgh", ".mgz"]
    while exit_loop == 0:
        name_file, ext_temp = os.path.splitext(name_file)
        ext_file = ext_temp + ext_file

        if not ext_temp:
            exit_loop = 1

        if ext_file in ext_key:
            exit_loop = 1

    return path_file, name_file, ext_file


def load_volume(file_in):
    """Load nifti volume.

    Parameters
    ----------
    file_in : str
        File name of nifti volume.

    Returns
    -------
    arr : (N, M, O, ...) np.ndarray
        Image array.
    affine : (4, 4) np.ndarray
        Affine transformation matrix.
    header : nb.nifti1.Nifti1Header
        Header information.

    """

    if isinstance(file_in, str):
        image = nb.load(file_in)
        arr = image.get_fdata()
        affine = image.affine
        header = image.header
    else:
        raise ValueError('No valid file name!')

    return arr, affine, header


def save_volume(file_out, arr, affine=None, header=None):
    """Save nifti volume.

    Parameters
    ----------
    file_out : str
        File name of saved nifti volume.
    arr : (N, M, O, ...) np.ndarray
        Image array.
    affine : (4, 4) np.ndarray, optional
        Affine transformation matrix.
    header : nb.nifti1.Nifti1Header, optional
        Header information

    Returns
    -------
    None

    """

    if isinstance(file_out, str):

        # make output folder
        dir_out = os.path.dirname(file_out)
        if not os.path.exists(dir_out):
            os.makedirs(dir_out)

        # use identity matrix if no affine transformation matrix is set
        if affine is None:
            affine = np.eye(4)

        # use empty header if no header information is set
        if header is None:
            header = nb.Nifti1Header()

        # write to disk
        image = nb.Nifti1Image(arr, affine, header)
        image.to_filename(file_out)
    else:
        raise ValueError('No valid file name!')



















def load_mesh(surf_mesh):

    if surf_mesh.endswith('vtk'):
        points, faces, data = _read_vtk(surf_mesh)
        return {'points': points, 'faces': faces, 'data': data}
    else:
        geom = load_mesh_geometry(surf_mesh)
        return geom



def save_mesh(filename, surf_dict):

    if filename.endswith('vtk'):
        _write_vtk(filename,
                   surf_dict['points'],
                   surf_dict['faces'],
                   surf_dict['data'])
    else:
        save_mesh_geometry(filename,
                           surf_dict)


def load_mesh_geometry(surf_mesh):

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


def save_mesh_geometry(filename, surf_dict):

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















def load_mesh_data():
    """mgh etc."""
    pass

def save_mesh_data():
    pass









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






