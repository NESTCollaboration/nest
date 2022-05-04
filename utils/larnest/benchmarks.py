"""
Benchmark functions for LArNEST
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

from parameters import *
from lar_dataset import LArDataset

def generate_plot_grid(
    num_plots,
    **kwargs,
):
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

class LArNESTBenchmarks:
    """
    """
    def __init__(self,
        input_file: str,
        dataset_dir: str='',
    ):
        self.input_file = input_file
        self.dataset_dir = dataset_dir
        self.data = pd.read_csv(
            self.input_file,
            names=['energy', 'efield' , 'TotalYields', 'QuantaYields', 'LightYields', 'Nph' , 'Ne', 'Nex', 'Nion'],
            header=0,
            dtype=float
        )
        self.data.replace([np.inf, -np.inf, np.nan], 0.0)
        self.unique_efield = np.unique(self.data['efield'])
        self.unique_energy = np.unique(self.data['energy'])

        if dataset_dir != '':
            self.dataset = LArDataset(self.dataset_dir)

        if not os.path.isdir("plots/"):
            os.makedirs("plots/")

    def get_dataset(self,
        yields: str='TotalYields',
        recoil: str='NR'
    ):
        if recoil == "NR": 
            if yields == "TotalYields":
                return self.dataset.nr_total
            elif yields == "QuantaYields":
                return self.dataset.nr_charge
            elif yields == "LightYields":
                return self.dataset.nr_light
            else:
                return []
        elif recoil == "ER":
            if yields == "TotalYields":
                return []
            elif yields == "QuantaYields":
                return self.dataset.er_charge
            elif yields == "LightYields":
                return self.dataset.er_light   
            else:
                return []

    def plot_yields(self,
        yields:  str='Ne',
        recoil:  str='NR',
        include_data:   bool=False,
        label_by_efield:bool=True,
        title:  str='',
        save:   str='',
        show:   bool=False,
    ):
        fig, axs = plt.subplots(figsize=(10,6))
        for efield in self.unique_efield:
            mask = (self.data['efield'] == efield)            
            axs.plot(
                self.unique_energy, 
                np.array(self.data[yields][mask],dtype=float),
                label=f"NEST - {efield} V/cm"
            )
        if include_data:
            data = self.get_dataset(yields,recoil)
            if len(data) != 0:
                unique_names = np.unique(data['Dataset'])
                for name in unique_names:
                    name_mask = (data['Dataset'] == name)
                    if label_by_efield:
                        unique_efield = np.unique(data['field'][name_mask])
                        for efield in unique_efield:
                            efield_mask = (data['field'][name_mask] == efield)
                            energy = data['energy'][name_mask][efield_mask]
                            energy_sl = data['energy_sl'][name_mask][efield_mask]
                            energy_sh = data['energy_sh'][name_mask][efield_mask]
                            results = data['yield'][name_mask][efield_mask]
                            yields_sl = data['yield_sl'][name_mask][efield_mask]
                            yields_sh = data['yield_sh'][name_mask][efield_mask]
                            axs.errorbar(
                                energy, results, 
                                xerr=[energy_sl, energy_sh], 
                                yerr=[yields_sl, yields_sh], 
                                linestyle='', marker='.',
                                label=f"{name} - {efield:.2f} V/cm"
                            )
                    else:
                        energy = data['energy'][name_mask]
                        energy_sl = data['energy_sl'][name_mask]
                        energy_sh = data['energy_sh'][name_mask]
                        results = data['yield'][name_mask]
                        yields_sl = data['yield_sl'][name_mask]
                        yields_sh = data['yield_sh'][name_mask]
                        axs.errorbar(
                            energy, results, 
                            xerr=[energy_sl, energy_sh], 
                            yerr=[yields_sl, yields_sh], 
                            linestyle='', marker='.',
                            label=f"{name}"
                        )
        axs.set_xlabel("Energy [keV]")
        if yields == 'TotalYields':
            axs.set_ylabel("Total Yield [quanta/keV]")
            axs.set_yscale("log")
        elif yields == 'QuantaYields':
            axs.set_ylabel("Quanta Yield [quanta/keV]")
            axs.set_yscale("linear")
        elif yields == 'LightYields':
            axs.set_ylabel("Light Yield [quanta/keV]")
            axs.set_yscale("linear")
        elif yields == 'Nph':
            axs.set_ylabel("Number of photons [photons]")
            axs.set_yscale("log")
        elif yields == 'Ne':
            axs.set_ylabel("Number of electrons [electrons]")
            axs.set_yscale("log")
        elif yields == 'Nex':
            axs.set_ylabel("Number of excitons [quanta]")
            axs.set_yscale("log")
        else:
            axs.set_ylabel("Number of ions [quanta]")
            axs.set_yscale("log")
        axs.set_xscale("log")
        if title != '':
            plt.title(title)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        if save != '':
            plt.savefig(f"plots/{save}.png")
        if show:
            plt.show()
        plt.close()
    
    def plot_yields_separated(self,
        yields:  str='Ne',
        recoil:  str='NR',
        include_data:   bool=False,
        title:  str='',
        save:   str='',
        show:   bool=False,
    ):
        fig, axs = generate_plot_grid(
            len(self.unique_efield), 
            figsize=(10,6)
        )
        for ii, efield in enumerate(self.unique_efield):
            mask = (self.data['efield'] == efield)            
            axs.flat[ii].plot(
                self.unique_energy, 
                np.array(self.data[yields][mask],dtype=float),
                label=f"NEST - {efield} V/cm"
            )
            axs.flat[ii].legend()
            axs.flat[ii].set_xlabel("Energy [keV]")
            if yields == 'TotalYields':
                axs.flat[ii].set_ylabel("Total Yield [quanta/keV]")
                axs.flat[ii].set_yscale("log")
            elif yields == 'QuantaYields':
                axs.flat[ii].set_ylabel("Quanta Yield [quanta/keV]")
                axs.flat[ii].set_yscale("linear")
            elif yields == 'LightYields':
                axs.flat[ii].set_ylabel("Light Yield [quanta/keV]")
                axs.flat[ii].set_yscale("linear")
            elif yields == 'Nph':
                axs.flat[ii].set_ylabel("Number of photons [photons]")
                axs.flat[ii].set_yscale("log")
            elif yields == 'Ne':
                axs.flat[ii].set_ylabel("Number of electrons [electrons]")
                axs.flat[ii].set_yscale("log")
            elif yields == 'Nex':
                axs.flat[ii].set_ylabel("Number of excitons [quanta]")
                axs.flat[ii].set_yscale("log")
            else:
                axs.flat[ii].set_ylabel("Number of ions [quanta]")
                axs.flat[ii].set_yscale("log")
            axs.flat[ii].set_xscale("log")
        if include_data:
            data = self.get_dataset(yields,recoil)
            if len(data) != 0:
                unique_names = np.unique(data['Dataset'])
                for name in unique_names:
                    name_mask = (data['Dataset'] == name)
                    unique_efield = np.unique(data['field'][name_mask])
                    for efield in unique_efield:
                        efield_mask = (data['field'][name_mask] == efield)
                        energy = data['energy'][name_mask][efield_mask]
                        energy_sl = data['energy_sl'][name_mask][efield_mask]
                        energy_sh = data['energy_sh'][name_mask][efield_mask]
                        results = data['yield'][name_mask][efield_mask]
                        yields_sl = data['yield_sl'][name_mask][efield_mask]
                        yields_sh = data['yield_sh'][name_mask][efield_mask]
                        axs.errorbar(
                            energy, results, 
                            xerr=[energy_sl, energy_sh], 
                            yerr=[yields_sl, yields_sh], 
                            linestyle='', marker='.',
                            label=f"{name} - {efield:.2f} V/cm"
                        )
                else:
                    energy = data['energy'][name_mask]
                    energy_sl = data['energy_sl'][name_mask]
                    energy_sh = data['energy_sh'][name_mask]
                    results = data['yield'][name_mask]
                    yields_sl = data['yield_sl'][name_mask]
                    yields_sh = data['yield_sh'][name_mask]
                    axs.errorbar(
                        energy, results, 
                        xerr=[energy_sl, energy_sh], 
                        yerr=[yields_sl, yields_sh], 
                        linestyle='', marker='.',
                        label=f"{name}"
                    )
        if title != '':
            plt.title(title)
        #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        if save != '':
            plt.savefig(f"plots/{save}.png")
        if show:
            plt.show()
        plt.close()


    def plot_all_yields(self,
        recoil: str="NR",
        include_data:   bool=False,
        label_by_efield:bool=True,
    ):
        for y in ['TotalYields', 'QuantaYields', 'LightYields', 'Nph' , 'Ne', 'Nex', 'Nion']:
            self.plot_yields(
                yields=y,
                recoil=recoil,
                include_data=include_data,
                label_by_efield=label_by_efield,
                title=f'{recoil} {y}',
                save=f'{recoil}_{y}'
            )

if __name__ == "__main__":


    # first make the NR plots
    nr_file = "../../build/larnest_benchmarks_nr.csv"
    larnest_benchmarks = LArNESTBenchmarks(
        nr_file,
        "data/",
    )
    larnest_benchmarks.plot_yields_separated(
        "TotalYields", "NR", show=True
    )
    # larnest_benchmarks.plot_all_yields(
    #     "NR",
    #     include_data=True,
    #     label_by_efield=False
    # )

    # # now the ER plots
    # er_file = "../../build/larnest_benchmarks_er.csv"
    # larnest_benchmarks = LArNESTBenchmarks(
    #     er_file
    # )
    # larnest_benchmarks.plot_all_yields("ER")

    # # Alpha plots
    # er_file = "../../build/larnest_benchmarks_alpha.csv"
    # larnest_benchmarks = LArNESTBenchmarks(
    #     er_file
    # )
    # larnest_benchmarks.plot_all_yields("Alpha")

    # # legacy plots
    # plot_legacy = True
    # if plot_legacy:
    #     # first make the NR plots
    #     nr_file = "../../build/legacy_larnest_benchmarks_2112.csv"
    #     larnest_benchmarks = LArNESTBenchmarks(
    #         nr_file
    #     )
    #     larnest_benchmarks.plot_all_yields("LegacyNR")

    #     # now the ER plots
    #     er_file = "../../build/legacy_larnest_benchmarks_13.csv"
    #     larnest_benchmarks = LArNESTBenchmarks(
    #         er_file
    #     )
    #     larnest_benchmarks.plot_all_yields("LegacyER")