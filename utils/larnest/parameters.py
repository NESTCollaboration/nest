"""
Collections of datasets and information.
"""
nr_total_yield_default = {
    100: {
        "SCENE 2015": [96.4, 193, 293],
        "ARIS 2018":  [200],
    },
}
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
    1: {
        "ARIS 2018":  [0],
        "SCENE 2015": [0],
        "Regenfus 2012":[0],
        "WARP 2005":    [0],
        "CREUS 2015":   [0],
        "MicroClean 2012": [0],
    },
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
er_charge_yield_default = {
    1: {
        "ARIS 2018": [0],
        "DarkSide10 2013":  [0],
        "Lippincott 2010":  [0],
        "WARP 2005":    [0],
        "Kimura 2019":  [0],
        "Doke 2002":    [0]
    },
    100: {
        "ARIS 2018":    [50,100],
        "Scalettar 1982": [84,89,94,101,110,119,128,129,130,139,145,148],
    },
    200: {
        "ARIS 2018":    [200],
        "Joshi 2014":   [200],
        "Scalettar 1982":[152,160,176,180,185,188,192,200,205,210,240,243,248,268,270,276,285,309,313,351,371,375,388,391]
    },
    600: {
        "ARIS 2018":    [500],
        "Joshi 2014":   [550],
        "Scalettar 1982":   [410,411,427,441,478,481,490,510,531,536,601,620,626,661],
        "Bondar 2016":  [600],
        #"Aprile ": [572]
    },
    1000: {
        "Joshi 2014":   [1200],
        "Scalettar 1982":    [801,820,832,899,904,943],
        #"Aprile ": []
    },
    1500: {
        "Joshi 2014":   [1600,1750],
        "Scalettar 1982":   [1004,1012,1064,1244,1403,1455],
        "Bondar 2016":  [1750],
    },
    2500: {
        "Joshi 2014":   [2150,2400,3000],
        "Scalettar 1982":[2009,2064,2806,2913],
        "Bondar 2016":  [2400],
        "Sangiorgio":   [2400],
        "Doke 2002":    [2010],
    },
    6000: {
        "Scalettar 1982":[4600,5704,6586,6693],
        "Doke 2002":    [4020,5000,6010],
        #"Aprile",
    },
    9500: {
        "Scalettar 1982":[8490,8673,9481,9681],
        "Doke 2002":    [8000,9000,9990],
        #"Aprile",
    },
}
er_light_yield_default = {
    1: {
        "ARIS 2018": [0],
        "DarkSide10 2013":  [0],
        "Lippincott 2010":  [0],
        "WARP 2005":    [0],
        "Kimura 2019":  [0],
    },
    100: {
        "ARIS 2018":    [50,100],
        "Scalettar 1982": [84,89,94,101,110,119,128,129,130,139,145,148],
    },
    200: {
        "ARIS 2018":    [200],
        "Joshi 2014":   [200],
        "Scalettar 1982":[152,160,176,180,185,188,192,200,205,210,240,243,248,268,270,276,285,309,313,351,371,375,388,391]
    },
    600: {
        "ARIS 2018":    [500],
        "Joshi 2014":   [550],
        "Scalettar 1982":   [410,411,427,441,478,481,490,510,531,536,601,620,626,661],
        "Bondar 2016":  [600],
        #"Aprile ": [572]
    },
    1000: {
        "Joshi 2014":   [1200],
        "Scalettar 1982":    [801,820,832,899,904,943],
        #"Aprile ": []
    },
    1500: {
        "Joshi 2014":   [1600,1750],
        "Scalettar 1982":   [1004,1012,1064,1244,1403,1455],
        "Bondar 2016":  [1750],
    },
    2500: {
        "Joshi 2014":   [2150,2400,3000],
        "Scalettar 1982":[2009,2064,2806,2913],
        "Bondar 2016":  [2400],
        "Sangiorgio":   [2400],
        "Doke 2002":    [2010],
    },
    6000: {
        "Scalettar 1982":[4600,5704,6586,6693],
        "Doke 2002":    [4020,5000,6010],
        #"Aprile",
    },
    9500: {
        "Scalettar 1982":[8490,8673,9481,9681],
        "Doke 2002":    [8000,9000,9990],
        #"Aprile",
    },
}

default_plot_config = {
    "nr_total":     nr_total_yield_default,
    "nr_charge":    nr_charge_yield_default,
    "nr_light":     nr_light_yield_default,
    "er_charge":    er_charge_yield_default,
    "er_light":     er_light_yield_default,
}