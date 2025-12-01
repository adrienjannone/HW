"""
LELME2150 - Thermal cycles
Project - CO2 battery

Test code for your CO2 battery model

@authors : Igor Grégoire, Adrien Jannone
@date: November 30th, 2025
"""


from CO2_battery import CO2_battery

# From the slides
p_storage_co2_liquid = 70 * 1e5     # Pa
T_storage_water = 22 + 273.15       # K
T_storage_TES = 450 + 273.15        # K
Pe = 20e6                           # W

# From the Paper
T_amb = 15 + 273.15                 # K
p_dome = 1e5                        # Pa
k_dome = 0.5 * 1e-2                 # Pressure drop coefficient dome (inlet, outlet)
k_TES = 4e-2                        # Pressure drop coefficient TES (inlet, outlet)
eta_comp = 0.85                     # Compressor adiabatic efficiency (isentropic)
eta_turb = 0.9                      # Turbine adiabatic efficiency (isentropic)
eta_pump = 0.85                     # Pump adiabatic efficiency (isentropic)
eta_mec = 0.99                      # Mechanical efficiency (shaft)
eta_elec = 0.985                    # Electrical efficiency (motor/generator)
pinch_TES0 = 2                      # K
pinch_TES = 7.5                     # K
p_amb = 1e5                         # Pa
inputs = Pe
params = {
    'p_storage_co2_liquid': p_storage_co2_liquid,
    'T_storage_water': T_storage_water,
    'T_storage_TES': T_storage_TES,
    'T_amb': T_amb,
    'p_amb': p_amb,
    'p_dome': p_dome,
    'k_dome': k_dome,
    'eta_comp': eta_comp,
    'eta_turb': eta_turb,
    'eta_pump': eta_pump,
    'eta_mec': eta_mec,
    'eta_elec': eta_elec,
    'pinch_TES0': pinch_TES0,
    'pinch_TES': pinch_TES,
    'fluid': 'CO2'
}


dome = CO2_battery(inputs, params, True)
dome.evaluate()
