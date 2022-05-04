"""
Functions for constructing training datasets from 
LAr data files.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import csv

from parameters import *

class LArDataset:
    """
    """
    def __init__(self,
        dataset_dir:    str='',
    ):
        self.dataset_dir = dataset_dir
        self.nr_total = pd.read_csv(
            self.dataset_dir + "nr_total_alt.csv"
        )
        self.nr_charge = pd.read_csv(
            self.dataset_dir + "nr_charge_alt.csv"
        )
        self.nr_light = pd.read_csv(
            self.dataset_dir + "nr_light_alt.csv"
        )
        

        # for ii, dataset in enumerate(datasets):
        #     energy, energy_sl, energy_sh = [], [], []
        #     field, field_sl, field_sh = [], [], []
        #     yields, yields_sl, yields_sh = [], [], []
        #     for data in datasets[dataset][0]:
        #         df = pd.read_csv(
        #             self.dataset_dir + data,
        #             usecols=datasets[dataset][1].values(),
        #         )
        #         name.append(df[datasets[dataset][1]['n']])
        #         energy.append(df[datasets[dataset][1]['energy']])
        #         energy_sl.append(df[datasets[dataset][1]['energy_sl']])
        #         energy_sh.append(df[datasets[dataset][1]['energy_sh']])

        #         field.append(df[datasets[dataset][1]['field']])
        #         field_sl.append(df[datasets[dataset][1]['field_sl']])
        #         field_sh.append(df[datasets[dataset][1]['field_sh']])

        #         yields.append(df[datasets[dataset][1]['yield']])
        #         yields_sl.append(df[datasets[dataset][1]['yield_sl']])
        #         yields_sh.append(df[datasets[dataset][1]['yield_sh']])

        #     energy = np.concatenate(energy)
        #     energy_sl = np.concatenate(energy_sl)
        #     energy_sh = np.concatenate(energy_sh)
            
        #     field = np.concatenate(field)
        #     field_sl = np.concatenate(field_sl)
        #     field_sh = np.concatenate(field_sh)
            
        #     yields = np.concatenate(yields)
        #     yields_sl = np.concatenate(yields_sl)
        #     yields_sh = np.concatenate(yields_sh)
            
        #     unique_field = np.unique(field)
        #     energy_s = np.array(list(zip(energy_sl,energy_sh))).T
        #     field_s = np.array(list(zip(field_sl,field_sh))).T
        #     yields_s = np.array(list(zip(yields_sl,yields_sh))).T

        #     self.datasets[dataset] = [
        #         energy, energy_sl, energy_sh,
        #         field, field_sl, field_sh,
        #         yields, yields_sl, yields_sh,
        #         energy_s, field_s, yields_s,
        #         unique_field,
        #     ]
    
    def plot_dataset(self,
        dataset
    ):
        fig = plt.figure(figsize=(10,6))
        axs = fig.add_subplot(projection='3d')
        axs.scatter(
            self.datasets[dataset][0], self.datasets[dataset][3], self.datasets[dataset][6],
            s=self.datasets[dataset][0], c=self.datasets[dataset][3]
        )
        for name in datasets[dataset][2]:
            axs.scatter([],[],[],label=name)
        axs.set_xlabel(datasets[dataset][1]['energy'])
        axs.set_ylabel(datasets[dataset][1]['field'])
        axs.set_zlabel(datasets[dataset][1]['yield'])
        axs.set_title(f"{dataset}"+r"$(E,\vec{E})$")
        plt.legend(bbox_to_anchor=(1.10, 1.0), loc='upper left')
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":

    lar_dataset = LArDataset()
    lar_dataset.plot_dataset("NR Ly")