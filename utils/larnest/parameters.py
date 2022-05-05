"""
Collections of datasets and information.
"""
nr_charge_yield_default = {
    100: {
        "SCENE 2015": [96.4],
    },
    250: {
        "SCENE 2015": [193, 293],
        "ARIS 2018":  [200],
        "Joshi 2014": [240],
    },
    500: {
        "SCENE 2015": [486],
        "Joshi 2014": [640],
    },
    1750: {
        "Joshi 2014": [1600, 2130],
    },
}
nr_light_yield_default = {
    50: {
        "SCENE 2013": [50],
        "ARIS 2018":  [50],
    },
    100: {
        "SCENE 2015": [100],
        "ARIS 2018":  [100],
    },
    250: {
        "SCENE 2015": [193, 293],
        "ARIS 2018":  [200],
        "SCENE 2013": [300],
    },
    500: {
        "ARIS 2018":  [500],
    },
}

default_plot_config = {
    #"nr_total":     nr_total_yield_default,
    "nr_charge":    nr_charge_yield_default,
    "nr_light":     nr_light_yield_default,
    #"er_charge":    er_charge_yield_default,
    #"er_light":     er_light_yield_default,
}