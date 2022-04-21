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
            names=['energy', 'efield' , 'photons' , 'electrons'],
            header=0,
        )
        self.unique_efield = np.unique(self.data['efield'])
        self.unique_energy = np.unique(self.data['energy'])

    def plot_yields(self,
        yields:  str='electrons',
        title:  str='',
        save:   str='',
        show:   bool=False,
    ):
        fig, axs = plt.subplots(figsize=(10,6))
        for efield in self.unique_efield:
            means = []
            for energy in self.unique_energy:
                mask = (self.data['efield'] == efield) & (self.data['energy'] == energy)
                means.append(np.mean(self.data[yields][mask]))
            axs.plot(
                self.unique_energy, 
                np.array(means) / self.unique_energy,
                label=f"{efield} V/cm"
            )
        axs.set_xlabel("Energy (keV)")
        if yields == 'electrons':
            axs.set_ylabel("Number of electrons (e-)")
        else:
            axs.set_ylabel("Number of photons (phe)")
        axs.set_xscale("log")
        if title != '':
            plt.title(title)
        plt.legend()
        plt.tight_layout()
        if save != '':
            plt.savefig(f"{save}.png")
        if show:
            plt.show()

if __name__ == "__main__":

    larnest_benchmarks = LArNESTBenchmarks(
        "../../build/larnest_output_2112.csv"
    )
    larnest_benchmarks.plot_yields(
        yields='photons',
        title='NR Light Yield',
        save='nr_photons'
    )
    larnest_benchmarks.plot_yields(
        yields='electrons',
        title='NR Qy Yield',
        save='nr_electrons'
    )
