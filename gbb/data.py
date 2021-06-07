# -*- coding: utf-8 -*-
"""Download example data."""

# python standard library inputs
import os

# external inputs
import gdown

__all__ = ['download_data']


def download_data(dir_out):
    """Download data.

    This function downloads the example data from my personal google drive [1]_
    which enables you to directly execute the code of this repository. Data is
    taken from a human 7 T fMRI study. Further information about the dataset
    can be found in the provided documentation.

    Parameters
    ----------
    dir_out : str
        Directory in which downloaded files are stored.

    Raises
    ------
    FileExistsError
        If `dir_out` already exists.

    Returns
    -------
    dict
        Dictionary collecting the output under the following keys

        * control_points (str) : Path to control points.
        * ignore (str) : Path to ignore mask.
        * mean_epi_enhanced (str) : Path to mean epi with enhanced contrast.
        * mean_epi (str) : Path to mean epi.
        * vein (str) : Path to vein mask.
        * pial (str) : path to pial surface.
        * white (str) : path to white surface.

    References
    ----------
    .. [1] https://drive.google.com/drive/folders/1dwra3FDRIQ-W3iO81vfxvYRHKqGgqSWn

    """

    # make output folder
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    # base url to google drive folder
    url = 'https://drive.google.com/uc?export=download&id='

    # file IDs of source files
    file_id = ['1Ljotc8rgbJ33Uial708OSCPA5jzosGiR',
               '1w-iZlmKk6133uLxpl38LoVOkJel2kpWl',
               '1x1KO_ZQsTiht0IERLck4uNDRweTFiQYy',
               '1r6yUQpxmnALp6MWfB-pCvytB11off5gN',
               '1j02nHLrLX-sO5NnU5QoNt39sytdhPprh',
               '1ZY2Ld8VXocSfbxJIxIqm9hX1rBDB_99X',
               '1us0lNqNZe-YuO8SodsLbSe5gOmUE52cU',
               '1_FR17jOaEZpCtYL36nG0TMe9SYhK52Dg',
               ]

    # file names of output files
    filename = ['control_points.dat',
                'ignore.nii',
                'mean_epi_enhanced.nii',
                'mean_epi.nii',
                'vein.nii',
                'lh.pial',
                'lh.white',
                'Readme.md',
                ]

    file_sources = [url + id_ for id_ in file_id]
    file_targets = [os.path.join(dir_out, file) for file in filename]

    for source, target in zip(file_sources, file_targets):

        # check if file already exists
        if os.path.exists(target):
            raise FileExistsError("The file " + target + " already exists!")
        else:
            gdown.download(source, target, quiet=False)

    return dict(control_points=file_targets[0],
                ignore=file_targets[1],
                mean_epi_enhanced=file_targets[2],
                mean_epi=file_targets[3],
                vein=file_targets[4],
                pial=file_targets[5],
                white=file_targets[6])
