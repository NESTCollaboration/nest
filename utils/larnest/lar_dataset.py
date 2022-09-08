"""
Functions for constructing training datasets from 
LAr data files.
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import csv
import os

from parameters import *
from utils import generate_plot_grid

class LArDataset:
    """
    """
    def __init__(self,
        dataset_file: str='data/lar_data.npz',
        benchmarks_file: str='',
        plot_config: dict=default_plot_config,
        fluctuations: bool=False,
    ):
        self.dataset_file = dataset_file
        self.benchmarks_file = benchmarks_file
        self.plot_config = plot_config
        self.fluctuations = fluctuations
        # set up datasets
        self.data = dict(np.load(self.dataset_file, allow_pickle=True))
        self.datasets = {key: self.data[key].item() for key in self.data}

        self.benchmarks = pd.read_csv(
            self.benchmarks_file,
            names=[
                'type','energy', 'efield' , 
                'TotalYields', 'QuantaYields', 'LightYields', 
                'Nph' , 'Ne', 'Nex', 'Nion',
                'TotalYields_std', 'QuantaYields_std', 'LightYields_std',
                'Nph_std', 'Ne_std', 'Nex_std', 'Nion_std',
            ],
            header=0,
            #dtype=float
        )
        self.benchmarks.replace([np.inf, -np.inf, np.nan], 0.0)
        self.ylabels = {
            'nr_total':  r"$Y_q$ [quanta/keV]",
            'nr_charge': r"$Y_{e^-}$ [electron/keV]",
            'nr_light':  r"$Y_{\gamma}$ [photon/keV]",
            'nr_nph':    r"$\langle N_{ph}\rangle$ [photon]",
            'nr_ne':     r"$\langle N_{e^{-}}\rangle$ [electron]",
            'nr_nex':    r"$\langle N_{ex}\rangle$ [exciton]",
            'nr_nion':   r"$\langle N_{ion}\rangle$ [ion]",
            'er_total':  r"$Y_q$ [quanta/keV]",
            'er_charge': r"$Y_{e^-}$ [electron/keV]",
            'er_light':  r"$Y_{\gamma}$ [photon/keV]",
            'er_nph':    r"$\langle N_{ph}\rangle$ [photon]",
            'er_ne':     r"$\langle N_{e^{-}}\rangle$ [electron]",
            'er_nex':    r"$\langle N_{ex}\rangle$ [exciton]",
            'er_nion':   r"$\langle N_{ion}\rangle$ [ion]",
            'alpha_total':  r"$Y_q$ [quanta/keV]",
            'alpha_charge': r"$Y_{e^-}$ [electron/keV]",
            'alpha_light':  r"$Y_{\gamma}$ [photon/keV]",
            'alpha_nph':    r"$\langle N_{ph}\rangle$ [photon]",
            'alpha_ne':     r"$\langle N_{e^{-}}\rangle$ [electron]",
            'alpha_nex':    r"$\langle N_{ex}\rangle$ [exciton]",
            'alpha_nion':   r"$\langle N_{ion}\rangle$ [ion]",
        }
        self.titles = {
            'nr_total':  "LArNEST Nuclear Recoil Total Yields",
            'nr_charge': "LArNEST Nuclear Recoil Charge Yields",
            'nr_light':  "LArNEST Nuclear Recoil Light Yields",
            'nr_nph':    "LArNEST Nuclear Recoil Photon Yields",
            'nr_ne':     "LArNEST Nuclear Recoil Electron Yields",
            'nr_nex':    "LArNEST Nuclear Recoil Exciton Yields",
            'nr_nion':   "LArNEST Nuclear Recoil Ion Yields",
            'er_total':  "LArNEST Electron Recoil Total Yields",
            'er_charge': "LArNEST Electron Recoil Charge Yields",
            'er_light':  "LArNEST Electron Recoil Light Yields",
            'er_nph':    "LArNEST Electron Recoil Photon Yields",
            'er_ne':     "LArNEST Electron Recoil Electron Yields",
            'er_nex':    "LArNEST Electron Recoil Exciton Yields",
            'er_nion':   "LArNEST Electron Recoil Ion Yields",
            'alpha_total':  "LArNEST Alpha Total Yields",
            'alpha_charge': "LArNEST Alpha Charge Yields",
            'alpha_light':  "LArNEST Alpha Light Yields",
            'alpha_nph':    "LArNEST Alpha Photon Yields",
            'alpha_ne':     "LArNEST Alpha Electron Yields",
            'alpha_nex':    "LArNEST Alpha Exciton Yields",
            'alpha_nion':   "LArNEST Alpha Ion Yields",
        }
        self.benchmark_types = {
            'nr_total':  ['NR','TotalYields'],
            'nr_charge': ['NR','QuantaYields'],
            'nr_light':  ['NR','LightYields'],
            'nr_nph':    ['NR','Nph'],
            'nr_ne':     ['NR','Ne'],
            'nr_nex':    ['NR','Nex'],
            'nr_nion':   ['NR','Nion'],
            'er_total':  ['ER','TotalYields'],
            'er_charge': ['ER','QuantaYields'],
            'er_light':  ['ER','LightYields'],
            'er_nph':    ['ER','Nph'],
            'er_ne':     ['ER','Ne'],
            'er_nex':    ['ER','Nex'],
            'er_nion':   ['ER','Nion'],
            'alpha_total':  ['Alpha','TotalYields'],
            'alpha_charge': ['Alpha','QuantaYields'],
            'alpha_light':  ['Alpha','LightYields'],
            'alpha_nph':    ['Alpha','Nph'],
            'alpha_ne':     ['Alpha','Ne'],
            'alpha_nex':    ['Alpha','Nex'],
            'alpha_nion':   ['Alpha','Nion'],
        }

        self.fluctuation_types = [
            'Nph', 'Ne', 'Nex', 'Nion'
        ]

        # create plotting directory
        if not os.path.isdir("plots/mean_yields/"):
            os.makedirs("plots/mean_yields")
        if not os.path.isdir("plots/fluctuations/"):
            os.makedirs("plots/fluctuations")
        if not os.path.isdir("plots/data/"):
            os.makedirs("plots/data")
        if not os.path.isdir("plots/data_fluctuations/"):
            os.makedirs("plots/data_fluctuations")

    def plot_data_grid(self,
        dataset_type:   str='nr_total'
    ):
        if (len(self.plot_config[dataset_type].keys()) == 1):
            return self.plot_data(dataset_type)
        fig, axs = generate_plot_grid(
            len(self.plot_config[dataset_type].keys()),
            figsize=(20,12)
        )
        for ii, efield in enumerate(self.plot_config[dataset_type]):
            # plot the benchmark values for the given efield
            benchmark_mask = (self.benchmarks['type'] == self.benchmark_types[dataset_type][0]) & (self.benchmarks['efield'] == efield)
            energy = self.benchmarks['energy'][benchmark_mask]
            yields = self.benchmarks[self.benchmark_types[dataset_type][1]][benchmark_mask]
            yields[(yields<0)] = 0
            axs.flat[ii].plot(energy, yields, label=f"NEST - {efield} V/cm")
            if self.fluctuations:
                errors = self.benchmarks[self.benchmark_types[dataset_type][1]+"_std"][benchmark_mask]
                low_errors = yields - errors
                high_errors = yields + errors
                low_errors[(low_errors)<0] = 0
                axs.flat[ii].fill_between(energy, low_errors, high_errors)
            # plot the corresponding data points for this efield
            for jj, dataset in enumerate(self.plot_config[dataset_type][efield]):
                if (len(self.plot_config[dataset_type][efield][dataset]) < 4):
                    for kk, field in enumerate(self.plot_config[dataset_type][efield][dataset]):
                        dataset_mask = (self.datasets[dataset_type]['Dataset'] == dataset) & (self.datasets[dataset_type]['field'] == field)
                        dataset_energy = self.datasets[dataset_type]['energy'][dataset_mask]
                        dataset_energy_sl = self.datasets[dataset_type]['energy_sl'][dataset_mask]
                        dataset_energy_sh = self.datasets[dataset_type]['energy_sh'][dataset_mask]
                        dataset_yields = self.datasets[dataset_type]['yield'][dataset_mask]
                        dataset_yields_sl = self.datasets[dataset_type]['yield_sl'][dataset_mask]
                        dataset_yields_sh = self.datasets[dataset_type]['yield_sh'][dataset_mask]
                        axs.flat[ii].errorbar(
                            dataset_energy, dataset_yields,
                            xerr=[dataset_energy_sl, dataset_energy_sh],
                            yerr=[dataset_yields_sl, dataset_yields_sh],
                            linestyle='', marker='.',
                            label=f"{dataset} - {field} V/cm"
                        )
                else:
                    dataset_mask = (self.datasets[dataset_type]['Dataset'] == str(dataset)) & np.isin(np.array(self.datasets[dataset_type]['field']),self.plot_config[dataset_type][efield][dataset])
                    dataset_energy = self.datasets[dataset_type]['energy'][dataset_mask]
                    dataset_energy_sl = self.datasets[dataset_type]['energy_sl'][dataset_mask]
                    dataset_energy_sh = self.datasets[dataset_type]['energy_sh'][dataset_mask]
                    dataset_yields = self.datasets[dataset_type]['yield'][dataset_mask]
                    dataset_yields_sl = self.datasets[dataset_type]['yield_sl'][dataset_mask]
                    dataset_yields_sh = self.datasets[dataset_type]['yield_sh'][dataset_mask]
                    axs.flat[ii].errorbar(
                        dataset_energy, dataset_yields,
                        xerr=[dataset_energy_sl, dataset_energy_sh],
                        yerr=[dataset_yields_sl, dataset_yields_sh],
                        linestyle='', marker='.',
                        label=f"{dataset} - {min(self.plot_config[dataset_type][efield][dataset])}-{max(self.plot_config[dataset_type][efield][dataset])} V/cm"
                    )
            axs.flat[ii].set_xscale("log")
            if ("_n" in dataset_type):
                axs.flat[ii].set_yscale("log")
            axs.flat[ii].legend()
        fig.supxlabel("Energy [keV]")
        fig.supylabel(self.ylabels[dataset_type])
        plt.suptitle(self.titles[dataset_type])
        plt.tight_layout()
        if self.fluctuations:
            plt.savefig(f"plots/data_fluctuations/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}_Data.png")
        else:
            plt.savefig(f"plots/data/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}_Data.png")
        plt.close()

    def plot_data(self,
        dataset_type:   str='nr_total'
    ):
        for ii, efield in enumerate(self.plot_config[dataset_type]):
            fig, axs = plt.subplots()
            # plot the benchmark values for the given efield
            benchmark_mask = (self.benchmarks['type'] == self.benchmark_types[dataset_type][0]) & (self.benchmarks['efield'] == efield)
            energy = self.benchmarks['energy'][benchmark_mask]
            yields = self.benchmarks[self.benchmark_types[dataset_type][1]][benchmark_mask]
            yields[(yields<0)] = 0
            axs.plot(energy, yields, label=f"NEST - {efield} V/cm")
            if self.fluctuations:
                errors = self.benchmarks[self.benchmark_types[dataset_type][1]+"_std"][benchmark_mask]
                low_errors = yields - errors
                high_errors = yields + errors
                low_errors[(low_errors)<0] = 0
                axs.fill_between(energy, low_errors, high_errors)
            # plot the corresponding data points for this efield
            for jj, dataset in enumerate(self.plot_config[dataset_type][efield]):
                if (len(self.plot_config[dataset_type][efield][dataset]) < 4):
                    for kk, field in enumerate(self.plot_config[dataset_type][efield][dataset]):
                        dataset_mask = (self.datasets[dataset_type]['Dataset'] == str(dataset)) & (self.datasets[dataset_type]['field'] == field)
                        dataset_energy = self.datasets[dataset_type]['energy'][dataset_mask]
                        dataset_energy_sl = self.datasets[dataset_type]['energy_sl'][dataset_mask]
                        dataset_energy_sh = self.datasets[dataset_type]['energy_sh'][dataset_mask]
                        dataset_yields = self.datasets[dataset_type]['yield'][dataset_mask]
                        dataset_yields_sl = self.datasets[dataset_type]['yield_sl'][dataset_mask]
                        dataset_yields_sh = self.datasets[dataset_type]['yield_sh'][dataset_mask]
                        axs.errorbar(
                            dataset_energy, dataset_yields,
                            xerr=[dataset_energy_sl, dataset_energy_sh],
                            yerr=[dataset_yields_sl, dataset_yields_sh],
                            linestyle='', marker='.',
                            label=f"{dataset} - {field} V/cm"
                        )
                else:
                    dataset_mask = (self.datasets[dataset_type]['Dataset'] == str(dataset)) & np.isin(np.array(self.datasets[dataset_type]['field']),self.plot_config[dataset_type][efield][dataset])
                    dataset_energy = self.datasets[dataset_type]['energy'][dataset_mask]
                    dataset_energy_sl = self.datasets[dataset_type]['energy_sl'][dataset_mask]
                    dataset_energy_sh = self.datasets[dataset_type]['energy_sh'][dataset_mask]
                    dataset_yields = self.datasets[dataset_type]['yield'][dataset_mask]
                    dataset_yields_sl = self.datasets[dataset_type]['yield_sl'][dataset_mask]
                    dataset_yields_sh = self.datasets[dataset_type]['yield_sh'][dataset_mask]
                    axs.errorbar(
                        dataset_energy, dataset_yields,
                        xerr=[dataset_energy_sl, dataset_energy_sh],
                        yerr=[dataset_yields_sl, dataset_yields_sh],
                        linestyle='', marker='.',
                        label=f"{dataset} - {min(self.plot_config[dataset_type][efield][dataset])}-{max(self.plot_config[dataset_type][efield][dataset])} V/cm"
                    )
            axs.set_xscale("log")
            if ("_n" in dataset_type):
                axs.set_yscale("log")
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            axs.set_xlabel("Energy [keV]")
            axs.set_ylabel(self.ylabels[dataset_type])
            plt.title(self.titles[dataset_type])
            plt.tight_layout()
            if self.fluctuations:
                plt.savefig(f"plots/data_fluctuations/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}_{efield}_Data.png")
            else:
                plt.savefig(f"plots/data/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}_{efield}_Data.png")
            plt.close()
    
    def plot_mean_yields(self,
        dataset_type:   str='nr_total'
    ):
        fig, axs = plt.subplots(figsize=(10,6))
        # plot the benchmark values for the given efield
        benchmark_mask = (self.benchmarks['type'] == self.benchmark_types[dataset_type][0])
        unique_efields = np.unique(self.benchmarks['efield'][benchmark_mask])
        for efield in unique_efields:
            efield_mask = (self.benchmarks['efield'][benchmark_mask] == efield)
            energy = self.benchmarks['energy'][benchmark_mask][efield_mask]
            yields = self.benchmarks[self.benchmark_types[dataset_type][1]][benchmark_mask][efield_mask]
            yields[(yields<0)] = 0
            axs.plot(energy, yields, label=f"NEST - {efield} V/cm")
        axs.set_xscale("log")
        if ("_n" in dataset_type):
            axs.set_yscale("log")
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        axs.set_xlabel("Energy [keV]")
        axs.set_ylabel(self.ylabels[dataset_type])
        plt.title(self.titles[dataset_type])
        plt.tight_layout()
        plt.savefig(f"plots/mean_yields/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}.png")
        plt.close()
    
    def plot_all_yields(self,
    ):
        for benchmark_type in self.benchmark_types:
            self.plot_mean_yields(benchmark_type)

    def plot_fluctuations(self,
        dataset_type: str='nr_total'
    ):
        # plot the benchmark values for the given efield
        benchmark_mask = (self.benchmarks['type'] == self.benchmark_types[dataset_type][0])
        unique_efields = np.unique(self.benchmarks['efield'][benchmark_mask])
        for efield in unique_efields:
            fig, axs = plt.subplots(figsize=(10,6))
            efield_mask = (self.benchmarks['efield'][benchmark_mask] == efield)
            energy = self.benchmarks['energy'][benchmark_mask][efield_mask]
            yields = self.benchmarks[self.benchmark_types[dataset_type][1]][benchmark_mask][efield_mask]
            yields[(yields<0)] = 0
            errors = self.benchmarks[self.benchmark_types[dataset_type][1]+"_std"][benchmark_mask][efield_mask]
            low_errors = yields - errors
            high_errors = yields + errors
            low_errors[(low_errors)<0] = 0
            axs.fill_between(energy, low_errors, high_errors, color='y', alpha=0.7)
            axs.plot(energy, yields, label=f"NEST - {efield} V/cm", color='k', linestyle='-')
            axs.set_xscale("log")
            if ("_n" in dataset_type):
                axs.set_yscale("log")
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            axs.set_xlabel("Energy [keV]")
            axs.set_ylabel(self.ylabels[dataset_type])
            plt.title(self.titles[dataset_type])
            plt.tight_layout()
            plt.savefig(f"plots/fluctuations/{self.benchmark_types[dataset_type][0]}_{self.benchmark_types[dataset_type][1]}_{efield}Vcm.png")
            plt.close()
    
    def plot_all_fluctuations(self,
    ):
        for benchmark_type in self.benchmark_types:
            self.plot_fluctuations(benchmark_type)