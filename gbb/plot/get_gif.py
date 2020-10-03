# -*- coding: utf-8 -*-
import os
import imageio
import numpy as np


def get_gif(img_file, path_output, name_output, nsteps, duration):
    """
    Create gif movie from n input images with fading from one image to the next 
    image.
    Inputs:
        *img_file (list): list of input images.
        *path_output (str): path where output is saved.
        *name_output (str): base name of output file.
        *nsteps (int): number of generated transition images.
        *duration (float): duration for each frame in seconds.

    created by Daniel Haenelt
    Date created: 17-11-2018
    Last modified: 03-10-2020
    """

    # append first list item to the end of the list
    img_file.append(img_file[0])

    images = []
    for i in range(len(img_file)-1):
    
        # input
        img1 = imageio.imread(img_file[i])
        img2 = imageio.imread(img_file[i+1])
    
        # generates a list with transitions from img1 to img2
        for j in range(nsteps):
            img = np.multiply(img1, (nsteps-j)/nsteps) + np.multiply(img2, j/nsteps)
            img = img.astype("uint8")
            images.append(img)

    # save output gif
    imageio.mimwrite(os.path.join(path_output,name_output+".gif"), images, duration=duration)