"""
Benchmark functions for LArNEST
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv

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
        self.unique_efield = np.unique(self.data['efield'])
        self.unique_energy = np.unique(self.data['energy'])

    def plot_yields(self,
        yields:  str='Ne',
        title:  str='',
        save:   str='',
        show:   bool=False,
    ):
        fig, axs = plt.subplots(figsize=(10,6))
        for efield in self.unique_efield:
            mask = (self.data['efield'] == efield)
            print(
                np.array(self.data[yields][mask],dtype=float),self.unique_energy)
            axs.plot(
                self.unique_energy, 
                np.array(self.data[yields][mask],dtype=float) / self.unique_energy,
                label=f"{efield} V/cm"
            )
        axs.set_xlabel("Energy (keV)")
        if yields == 'TotalYields':
            axs.set_ylabel("Total Yield (quanta)")
        elif yields == 'QuantYields':
            axs.set_ylabel("Quanta Yield (quanta)")
        elif yields == 'LightYields':
            axs.set_ylabel("Light Yield (quanta)")
        elif yields == 'Nph':
            axs.set_ylabel("Number of photons (phe)")
        elif yields == 'Ne':
            axs.set_ylabel("Number of electrons (e-)")
        elif yields == 'Nex':
            axs.set_ylabel("Number of excitons (quanta)")
        else:
            axs.set_ylabel("Number of ions (quanta)")
        axs.set_xscale("log")
        if title != '':
            plt.title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig(f"{save}.png")
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

    # larnest_benchmarks = LArNESTBenchmarks(
    #     "../../build/LArNESTBenchmarks_output_NR.csv"
    # )
    # larnest_benchmarks.plot_all_yields("NR")
    larnest_benchmarks = LArNESTBenchmarks(
        "../../build/LArNESTBenchmarks_output_ER.csv"
    )
    larnest_benchmarks.plot_all_yields("ER")