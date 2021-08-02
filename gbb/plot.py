# -*- coding: utf-8 -*-
"""I/O functions."""

# external inputs
import numpy as np
import matplotlib.pyplot as plt

__all__ = ['plot_moco']


def plot_moco(params_in, file_out, mtype="trans"):
    """Plot translational or rotational motion parameters from afni's motion
    parameter file (1Dfile). The file is expected to contain 6 ASCII formatted
    columns (roll, pitch, yaw, dS, dL, dP). Parameters are plotted in RAS
    convention.

    Parameters
    ----------
    params_in : str
        File name of afni's motion parameter file (1Dfile).
    file_out : str
        File name of saved figure (svg format).
    mtype : str, optional
        Plot translation or rotation parameters (trans or rot).

    Returns
    -------
    None

    """

    # load afni 1Dparams file
    mp = np.loadtxt(params_in)

    # make plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_axes([0.12, 0.15, 0.8, 0.77])

    # afni convention
    # ---
    # roll: rotation about the I - S axis (degrees ccw)
    # pitch: rotation about the R - L axis (degrees ccw)
    # yaw: rotation about the A - P axis (degrees ccw)
    # dS: displacement in the superior direction (mm)
    # dL: displacement in the left direction (mm)
    # dP: displacement in the posterior direction (mm)

    # plotting convention
    # ---
    # roll: rotation about the P - A axis(degrees ccw)
    # pitch: rotation about the L - R axis (degrees ccw)
    # yaw: rotation about the I - S axis (degrees ccw)
    # dR: displacement in the right direction (mm)
    # dA: displacement in the anterior direction (mm)
    # dS: displacement in the superior direction (mm)

    if mtype == "rot":
        ax.plot(-1*mp[:, 2], lw=0.75, label="roll")
        ax.plot(-1*mp[:, 1], lw=0.75, label="pitch")
        ax.plot(mp[:, 0], lw=0.75, label="yaw")
        ax.set_ylabel("Rotation in deg (ccw)")
    elif mtype == "trans":
        ax.plot(-1*mp[:, 4], lw=0.75, label="x")
        ax.plot(-1*mp[:, 5], lw=0.75, label="y")
        ax.plot(mp[:, 3], lw=0.75, label="z")
        ax.set_ylabel("Translation in mm")
    else:
        raise ValueError("some statement")

    # make plot nice
    ax.legend(loc=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Volume")

    if file_out:
        plt.savefig(file_out,
                    dpi=300,
                    bbox_inches="tight",
                    transparent=True,
                    format="svg")
    else:
        plt.show()


def plot_cost():
    pass
