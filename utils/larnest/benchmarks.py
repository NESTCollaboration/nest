"""
Benchmark functions for LArNEST
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

class LArNESTBenchmarks:
    """
    """
    def __init__(self,
        input_file: str,
    ):
        self.input_file = input_file
        self.data = pd.read_csv(
            self.input_file,
            names=['energy', 'efield' , 'TotalYields', 'QuantaYields', 'LightYields', 'Nph' , 'Ne', 'Nex', 'Nion'],
            header=0,
            dtype=float
        )
        self.data.replace([np.inf, -np.inf, np.nan], 0.0)
        self.unique_efield = np.unique(self.data['efield'])
        self.unique_energy = np.unique(self.data['energy'])
        if not os.path.isdir("plots/"):
            os.makedirs("plots/")

    def plot_yields(self,
        yields:  str='Ne',
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
                label=f"{efield} V/cm"
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
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig(f"plots/{save}.png")
        if show:
            plt.show()

    def plot_all_yields(self,
        recoil: str="NR",
    ):
        for y in ['TotalYields', 'QuantaYields', 'LightYields', 'Nph' , 'Ne', 'Nex', 'Nion']:
            self.plot_yields(
                yields=y,
                title=f'{recoil} {y}',
                save=f'{recoil}_{y}'
            )

if __name__ == "__main__":

    # first make the NR plots
    nr_file = "../../build/larnest_benchmarks_nr.csv"
    larnest_benchmarks = LArNESTBenchmarks(
        nr_file
    )
    larnest_benchmarks.plot_all_yields("NR")

    # now the ER plots
    er_file = "../../build/larnest_benchmarks_er.csv"
    larnest_benchmarks = LArNESTBenchmarks(
        er_file
    )
    larnest_benchmarks.plot_all_yields("ER")

    # legacy plots
    plot_legacy = True
    if plot_legacy:
        # first make the NR plots
        nr_file = "../../build/legacy_larnest_benchmarks_2112.csv"
        larnest_benchmarks = LArNESTBenchmarks(
            nr_file
        )
        larnest_benchmarks.plot_all_yields("LegacyNR")

        # now the ER plots
        er_file = "../../build/legacy_larnest_benchmarks_13.csv"
        larnest_benchmarks = LArNESTBenchmarks(
            er_file
        )
        larnest_benchmarks.plot_all_yields("LegacyER")