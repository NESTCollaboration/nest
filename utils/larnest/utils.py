"""
Various utilities for plotting, etc.
"""
import numpy as np
from matplotlib import pyplot as plt


def generate_plot_grid(
    num_plots,
    **kwargs,
):
    """
    Generate a grid of N plots, excluding extra
    ones which are not used.
    """
    nrows = int(np.floor(np.sqrt(num_plots)))
    ncols = int(np.ceil(num_plots/nrows))
    fig, axs = plt.subplots(
        nrows=nrows, ncols=ncols,
        **kwargs
    )
    nplots = nrows * ncols
    nextra = nplots - num_plots
    for ii in range(nextra):
        axs.flat[-(ii+1)].set_visible(False)
    return fig, axs