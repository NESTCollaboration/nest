"""
Benchmark functions for LArNEST
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv

class LegacyBenchmarks:
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
    
    def plot_photons(self,
    ):
        fig, axs = plt.subplots(figsize=(10,6))
        for efield in self.unique_efield:
            means = []
            for energy in self.unique_energy:
                mask = (self.data['efield'] == efield) & (self.data['energy'] == energy)
                means.append(np.mean(self.data['photons'][mask]))
            axs.plot(
                self.unique_energy * 1000, 
                np.array(means) / 1000,
                label=f"{efield} V/cm"
            )
        axs.set_xlabel("Energy (keV)")
        axs.set_ylabel("Number of photons (phe)")
        axs.set_xscale("log")
        plt.legend()
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":

    legacy_benchmarks = LegacyBenchmarks("../../build/larnest_output_2112.csv")
    legacy_benchmarks.plot_photons()
