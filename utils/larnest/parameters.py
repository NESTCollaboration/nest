"""
Collections of datasets and information.
"""
"""----------------Alpha datasets----------------"""
alpha_charge_datasets = [
    "alpha/alpha_charge_hitachi_1987.csv",   
    "alpha/alpha_charge_scalettar_1982.csv"
]
alpha_light_datasets = [
    "alpha/alpha_light_agnes_2016.csv",
    "alpha/alpha_light_hitachi_1987.csv"
]

"""----------------ER Charge datasets----------------"""
er_charge_datasets = [
    "er/er_charge_aprile_1987.csv",
    "er/er_charge_bondar_2016.csv",
    "er/er_charge_doke_2002.csv",
    "er/er_charge_joshi_2014.csv",
    "er/er_charge_sangiorgio_2013.csv",
    "er/er_charge_scalettar_1982.csv",
    "er/er_charge_thomasimel_1987.csv",
]
er_charge_cols = {
    'energy':       'Energy [keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field [V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield [e/keV]',
    'yield_sl':     'Yield SL[e/keV]',
    'yield_sh':     'Yield SH[e/keV]',
}
er_charge_names = [
    "Aprile 1987",
    "Bondar 2016",
    "Doke 2002",
    "Joshi 2014",
    "Sangiorgio 2013",
    "Scalettar 1982",
    "Thomas-Imel 1987",
]
er_charge_converted_datasets = [
    "er/er_charge_aris_2018.csv",
    "er/er_charge_darkside10_2013.csv",
    "er/er_charge_doke_converted_2002.csv",
    "er/er_charge_kimura_2019.csv",
    "er/er_charge_lippincott_2010.csv",
    "er/er_charge_warp_2014.csv",
]

"""----------------ER Charge fits datasets----------------"""
er_charge_fits_datasets = [
    "er/er_charge_aris_2018.csv",
    "er/er_charge_aprile_1987.csv",
    "er/er_charge_bondar_2016.csv",
    "er/er_charge_darkside10_2013.csv",
    "er/er_charge_doke_2002.csv",
    "er/er_charge_joshi_2014.csv",
    "er/er_charge_kimura_2019.csv",
    "er/er_charge_lippincott_2010.csv",
    "er/er_charge_sangiorgio_2013.csv",
    "er/er_charge_scalettar_1982.csv",
    "er/er_charge_thomasimel_1987.csv",
    "er/er_charge_warp_2014.csv",
]
er_charge_fits_cols = {
    'energy':       'Energy [keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field [V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield [e/keV]',
    'yield_sl':     'Yield SL[e/keV]',
    'yield_sh':     'Yield SH[e/keV]',
}
er_charge_fits_names = [
    "ARIS 2018",
    "Aprile 1987",
    "Bondar 2016",
    "DarkSide10 2013",
    "Doke 2002",
    "Joshi 2014",
    "Kimura 2019",
    "Lippincott 2010",
    "Sangiorgio 2013",
    "Scalettar 1982",
    "Thomas-Imel 1987",
    "WARP 2014",
]

"""----------------ER Light datasets----------------"""
er_light_datasets = [
    "er/er_light_aris_2018.csv",
    "er/er_light_darkside10_2013.csv",
    "er/er_light_doke_2002.csv",
    "er/er_light_kimura_2019.csv",
    "er/er_light_kimuralidine_2019.csv",
    "er/er_light_lippincott_2010.csv",
    "er/er_light_warp_2014.csv",
]
er_light_cols = {
    'energy':       'Energy [keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field [V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield [photons/keV]',
    'yield_sl':     'Yield SL[photons/keV]',
    'yield_sh':     'Yield SH[photons/keV]',
}
er_light_names = [
    "ARIS 2018",
    "DarkSide10 2013",
    "Doke 2002",
    "Kimura 2019",
    "Kimura LIDINE 2019",
    "Lippincott 2010",
    "WARP 2014",
]
er_light_converted_datasets = [
    "er/er_light_bondar_2016.csv",
    "er/er_light_doke_converted_2002.csv",
    "er/er_light_joshi_2014.csv",
    "er/er_light_sangiorgio_2013.csv",
    "er/er_light_scalettar_1982.csv",
]
"""----------------ER Light fits datasets----------------"""
er_light_fits_datasets = [
    "er/er_light_aris_2018.csv",
    "er/er_light_bondar_2016.csv",
    "er/er_light_darkside10_2013.csv",
    "er/er_light_doke_2002.csv",
    "er/er_light_joshi_2014.csv",
    "er/er_light_kimura_2019.csv",
    "er/er_light_kimuralidine_2019.csv",
    "er/er_light_lippincott_2010.csv",
    "er/er_light_sangiorgio_2013.csv",
    "er/er_light_scalettar_1982.csv",
    "er/er_light_warp_2014.csv",
]
er_light_fits_cols = {
    'energy':       'Energy [keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field [V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield [photons/keV]',
    'yield_sl':     'Yield SL[photons/keV]',
    'yield_sh':     'Yield SH[photons/keV]',
}
er_light_fits_names = [
    "ARIS 2018",
    "Bondar 2016",
    "DarkSide10 2013",
    "Doke 2002",
    "Joshi 2014",
    "Kimura 2019",
    "Kimura LIDINE 2019",
    "Lippincott 2010",
    "Sangiorgio 2013",
    "Scalettar 1982",
    "WARP 2014",
]

"""----------------NR Charge datasets----------------"""
nr_charge_datasets = [
    "nr/nr_charge_bondar_2017.csv",
    "nr/nr_charge_joshi_2015.csv",
    "nr/nr_charge_scene_2015.csv",
]
nr_charge_names = [
    "Bondar 2017",
    "Joshi 2015",
    "SCENE 2015",
]
nr_charge_converted_datasets = [
    "nr/nr_charge_aris_2018.csv",
    "nr/nr_charge_bondar_2017.csv",
    "nr/nr_charge_creus_2015.csv",
    "nr/nr_charge_joshi_2015.csv",
    "nr/nr_charge_microclean_2012.csv",
    "nr/nr_charge_regenfus_2012.csv",
    "nr/nr_charge_scene_2013.csv",
    "nr/nr_charge_scene_converted_2015.csv",
    "nr/nr_charge_warp_2005.csv",
]
nr_charge_cols = {
    'energy':       'Energy[keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field[V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield[e/keV]',
    'yield_sl':     'Yield SL[e/keV]',
    'yield_sh':     'Yield SH[e/keV]',
}
"""----------------NR Charge fits datasets----------------"""
nr_charge_fits_datasets = [
    "nr/nr_charge_aris_2018.csv",
    "nr/nr_charge_bondar_2017.csv",
    "nr/nr_charge_joshi_2015.csv",
    "nr/nr_charge_scene_2015.csv",
]
nr_charge_fits_names = [
    "ARIS 2018",
    "Bondar 2017",
    "Joshi 2015",
    "SCENE 2015",
]
nr_charge_fits_cols = {
    'energy':       'Energy[keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field[V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield[e/keV]',
    'yield_sl':     'Yield SL[e/keV]',
    'yield_sh':     'Yield SH[e/keV]',
}

"""----------------NR Light datasets----------------"""
nr_light_datasets = [
    "nr/nr_light_aris_2018.csv",
    "nr/nr_light_creus_2015.csv",
    "nr/nr_light_kimura_2019.csv",
    "nr/nr_light_microclean_2012.csv",
    "nr/nr_light_regenfus_2012.csv",
    "nr/nr_light_scene_2013.csv",
    "nr/nr_light_scene_2015.csv",
    "nr/nr_light_warp_2005.csv",
]
nr_light_converted_datasets = [
]
nr_light_cols = {
    'energy':       'Energy[keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field[V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield[photons/keV]',
    'yield_sl':     'Yield SL[photons/keV]',
    'yield_sh':     'Yield SH[photons/keV]',
}
nr_light_names = [
    "ARIS 2018",
    "CREUS 2015",
    "Kimura 2019",
    "MicroCLEAN 2012",
    "Regenfus 2012",
    "SCENE 2013",
    "SCENE 2015",
    "WArP 2005",
]

"""----------------NR Light fits datasets----------------"""
nr_light_fits_datasets = [
    "nr/nr_light_aris_2018.csv",
    "nr/nr_light_creus_2015.csv",
    "nr/nr_light_kimura_2019.csv",
    "nr/nr_light_microclean_2012.csv",
    "nr/nr_light_regenfus_2012.csv",
    "nr/nr_light_scene_2013.csv",
    "nr/nr_light_scene_2015.csv",
    "nr/nr_light_warp_2005.csv",
]
nr_light_fits_cols = {
    'energy':       'Energy[keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field[V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield[photons/keV]',
    'yield_sl':     'Yield SL[photons/keV]',
    'yield_sh':     'Yield SH[photons/keV]',
}
nr_light_fits_names = [
    "ARIS 2018",
    "CREUS 2015",
    "Kimura 2019",
    "MicroCLEAN 2012",
    "Regenfus 2012",
    "SCENE 2013",
    "SCENE 2015",
    "WArP 2005",
]


"""----------------NR Total datasets----------------"""
nr_total_datasets = [
    "nr_alt/nr_total_aris_2018_alt.csv",
    "nr_alt/nr_total_scene_2015_alt.csv",
]
nr_total_cols = {
    'energy':       'Energy[keV]',
    'energy_sl':    'Energy SL[keV]',
    'energy_sh':    'Energy SH[keV]',
    'field':        'Field[V/cm]',
    'field_sl':     'Field SL[V/cm]',
    'field_sh':     'Field SH[V/cm]',
    'yield':        'Yield[quanta/keV]',
    'yield_sl':     'Yield SL[quanta/keV]',
    'yield_sh':     'Yield SH[quanta/keV]',
}
nr_total_names = [
    "ARIS 2018",
    "SCENE 2015",
]

"""----------------NR Charge alternate datasets----------------"""
nr_charge_alt_datasets = [
    "nr/nr_charge_bondar_2015_alt.csv",
    "nr/nr_charge_joshi_2015_alt.csv",
    "nr/nr_charge_scene_2015_alt.csv",
]
nr_charge_alt_converted_datasets = [
    "nr/nr_charge_aris_2018_alt.csv",
    "nr/nr_charge_creus_2015_alt.csv",
    "nr/nr_charge_microclean_2012.csv",
    "nr/nr_charge_regenfus_2012_alt.csv",
    "nr/nr_charge_scene_2013_alt.csv",
    "nr/nr_charge_scene_converted_2015_alt.csv",
    "nr/nr_charge_warp_2005_alt.csv",
]

"""----------------NR Light alternate datasets----------------"""
nr_light_alt_datasets = [
    "nr/nr_light_aris_2018_alt.csv",
    "nr/nr_light_creus_2015_alt.csv",
    "nr/nr_light_microclean_2012_alt.csv",
    "nr/nr_light_regenfus_2012_alt.csv",
    "nr/nr_light_scene_2013_alt.csv",
    "nr/nr_light_scene_2015_alt.csv",
    "nr/nr_light_warp_2005_alt.csv",
    "nr/nr_total_aris_2018_alt.csv",
    "nr/nr_total_scene_2015_alt.csv",
]
nr_light_alt_converted_datasets = [
]

"""----------------Drift datasets----------------"""
drift_datasets = [
    "drift/gushchin_1982.csv",
    "drift/halpern_1967.csv",
    "drift/icarus_2004.csv",
    "drift/li_2016.csv",
    "drift/microboone_2009.csv",
    "drift/miller_1968.csv",
    "drift/walkowiak_2000.csv",
    "drift/yoshino_1976.csv",
]


datasets = {
    'NR Qy': [
        nr_charge_datasets, 
        nr_charge_cols,
        nr_charge_names,
    ],
    'NR Ly': [
        nr_light_datasets,
        nr_light_cols,
        nr_light_names,
    ],
    'ER Qy': [
        er_charge_datasets, 
        er_charge_cols,
        er_charge_names,
    ],
    'ER Ly': [
        er_light_datasets,
        er_light_cols,
        er_light_names,
    ],
    'NR TotalYields fits': [
        nr_total_datasets,
        nr_total_cols,
        nr_total_names
    ],
    'NR Qy fits': [
        nr_charge_fits_datasets,
        nr_charge_fits_cols,
        nr_charge_fits_names
    ],
    'NR Ly fits': [
        nr_light_fits_datasets,
        nr_light_fits_cols,
        nr_light_fits_names
    ],
    # 'ER Ty fits': [
    #     er_total_datasets,
    #     er_total_cols,
    #     er_total_names
    # ],
    # 'ER Qy fits': [
    #     er_charge_fits_datasets,
    #     er_charge_fits_cols,
    #     er_charge_fits_names
    # ],
    # 'ER Ly fits': [
    #     er_light_fits_datasets,
    #     er_light_fits_cols,
    #     er_light_fits_names
    # ],
    # 'Alpha Ty fits': [
    #     alpha_total_datasets,
    #     alpha_total_cols,
    #     alpha_total_names
    # ],
    # 'Alpha Qy fits': [
    #     alpha_charge_fits_datasets,
    #     alpha_charge_fits_cols,
    #     alpha_charge_fits_names
    # ],
    # 'Alpha Ly fits': [
    #     alpha_light_fits_datasets,
    #     alpha_light_fits_cols,
    #     alpha_light_fits_names
    # ],
}