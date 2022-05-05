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
    benchmarks_file = "data/benchmarks.csv"
    
    # create the dataset object
    lar_dataset = LArDataset(
        dataset_file=dataset_file,
        benchmarks_file=benchmarks_file,
        plot_config=default_plot_config,
    )

    # make benchmark plots
    lar_dataset.plot_all_yields()