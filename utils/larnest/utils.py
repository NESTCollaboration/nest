"""
Various utilities for plotting, etc.
"""
import numpy as np
from matplotlib import pyplot as plt
import csv


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

if __name__ == "__main__":

    nr_total = {
        "Dataset":  [],
        "energy":   [],
        "energy_sl":[],
        "energy_sh":[],
        "field":    [],
        "field_sl": [],
        "field_sh": [],
        "yield":   [],
        "yield_sl": [],
        "yield_sh": []
    }
    with open("data/nr_total_alt.csv","r") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader)
        for row in reader:
            nr_total["Dataset"].append(row[0])
            nr_total["energy"].append(float(row[1]))
            nr_total["energy_sl"].append(float(row[2]))
            nr_total["energy_sh"].append(float(row[3]))
            nr_total["field"].append(float(row[4]))
            nr_total["field_sl"].append(float(row[5]))
            nr_total["field_sh"].append(float(row[6]))
            nr_total["yield"].append(float(row[7]))
            nr_total["yield_sl"].append(float(row[8]))
            nr_total["yield_sh"].append(float(row[9]))
    nr_total = {key: np.array(nr_total[key]) for key in nr_total}
    nr_charge = {
        "Dataset":  [],
        "energy":   [],
        "energy_sl":[],
        "energy_sh":[],
        "field":    [],
        "field_sl": [],
        "field_sh": [],
        "yield":   [],
        "yield_sl": [],
        "yield_sh": []
    }
    with open("data/nr_charge_alt.csv","r") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader)
        for row in reader:
            nr_charge["Dataset"].append(row[0])
            nr_charge["energy"].append(float(row[1]))
            nr_charge["energy_sl"].append(float(row[2]))
            nr_charge["energy_sh"].append(float(row[3]))
            nr_charge["field"].append(float(row[4]))
            nr_charge["field_sl"].append(float(row[5]))
            nr_charge["field_sh"].append(float(row[6]))
            nr_charge["yield"].append(float(row[7]))
            nr_charge["yield_sl"].append(float(row[8]))
            nr_charge["yield_sh"].append(float(row[9]))
    nr_charge = {key: np.array(nr_charge[key]) for key in nr_charge}
    nr_light = {
        "Dataset":  [],
        "energy":   [],
        "energy_sl":[],
        "energy_sh":[],
        "field":    [],
        "field_sl": [],
        "field_sh": [],
        "yield":   [],
        "yield_sl": [],
        "yield_sh": []
    }
    with open("data/nr_light_alt.csv","r") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader)
        for row in reader:
            nr_light["Dataset"].append(row[0])
            nr_light["energy"].append(float(row[1]))
            nr_light["energy_sl"].append(float(row[2]))
            nr_light["energy_sh"].append(float(row[3]))
            nr_light["field"].append(float(row[4]))
            nr_light["field_sl"].append(float(row[5]))
            nr_light["field_sh"].append(float(row[6]))
            nr_light["yield"].append(float(row[7]))
            nr_light["yield_sl"].append(float(row[8]))
            nr_light["yield_sh"].append(float(row[9]))
    nr_light = {key: np.array(nr_light[key]) for key in nr_light}
    er_charge = {
        "Dataset":  [],
        "energy":   [],
        "energy_sl":[],
        "energy_sh":[],
        "field":    [],
        "field_sl": [],
        "field_sh": [],
        "yield":   [],
        "yield_sl": [],
        "yield_sh": []
    }
    with open("data/er_charge.csv","r") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader)
        for row in reader:
            er_charge["Dataset"].append(row[0])
            er_charge["energy"].append(float(row[1]))
            er_charge["energy_sl"].append(float(row[2]))
            er_charge["energy_sh"].append(float(row[3]))
            er_charge["field"].append(float(row[4]))
            er_charge["field_sl"].append(float(row[5]))
            er_charge["field_sh"].append(float(row[6]))
            er_charge["yield"].append(float(row[7]))
            er_charge["yield_sl"].append(float(row[8]))
            er_charge["yield_sh"].append(float(row[9]))
    er_charge = {key: np.array(er_charge[key]) for key in er_charge}
    er_light = {
        "Dataset":  [],
        "energy":   [],
        "energy_sl":[],
        "energy_sh":[],
        "field":    [],
        "field_sl": [],
        "field_sh": [],
        "yield":   [],
        "yield_sl": [],
        "yield_sh": []
    }
    with open("data/er_light.csv","r") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader)
        for row in reader:
            er_light["Dataset"].append(row[0])
            er_light["energy"].append(float(row[1]))
            er_light["energy_sl"].append(float(row[2]))
            er_light["energy_sh"].append(float(row[3]))
            er_light["field"].append(float(row[4]))
            er_light["field_sl"].append(float(row[5]))
            er_light["field_sh"].append(float(row[6]))
            er_light["yield"].append(float(row[7]))
            er_light["yield_sl"].append(float(row[8]))
            er_light["yield_sh"].append(float(row[9]))
    er_light = {key: np.array(er_light[key]) for key in er_light}
    
    np.savez(
        "data/lar_data.npz",
        nr_total=nr_total,
        nr_charge=nr_charge,
        nr_light=nr_light,
        er_charge=er_charge,
        er_light=er_light,
    )