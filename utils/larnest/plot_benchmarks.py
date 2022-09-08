"""
Benchmark functions for LArNEST
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

from parameters import *
from lar_dataset import LArDataset

if __name__ == "__main__":

    dataset_file = "data/lar_data.npz"
    mean_yields_benchmarks_file = "data/mean_yields_benchmarks.csv"
    fluctuations_benchmarks_file = "data/fluctuation_benchmarks.csv"
    lar_datasets = []

    # create the dataset object
    lar_datasets.append(LArDataset(
        dataset_file=dataset_file,
        benchmarks_file=mean_yields_benchmarks_file,
        plot_config=default_plot_config,
        fluctuations=False,
    ))
    lar_datasets.append(LArDataset(
        dataset_file=dataset_file,
        benchmarks_file=fluctuations_benchmarks_file,
        plot_config=default_plot_config,
        fluctuations=True,
    ))

    # make benchmark plots
    lar_datasets[0].plot_all_yields()
    lar_datasets[1].plot_all_fluctuations()

    for ii, lar_dataset in enumerate(lar_datasets):
        # plot with data
        # nr yields
        lar_dataset.plot_data('nr_total')
        lar_dataset.plot_data_grid('nr_total')
        lar_dataset.plot_data('nr_charge')
        lar_dataset.plot_data_grid('nr_charge')
        lar_dataset.plot_data('nr_light')
        lar_dataset.plot_data_grid('nr_light')

        # er yields
        lar_dataset.plot_data('er_charge')
        lar_dataset.plot_data_grid('er_charge')
        lar_dataset.plot_data('er_light')
        lar_dataset.plot_data_grid('er_light')