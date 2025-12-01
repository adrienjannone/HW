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
eta_comp = 0.85                     # Compressor adiabatic efficiency (isentropic)
eta_turb = 0.9                      # Turbine adiabatic efficiency (isentropic)
eta_pump = 0.85                     # Pump adiabatic efficiency (isentropic)
pinch = 2 + 273.15                  # K (jsp où tu l'as trouvé)
inputs = Pe
params = {
    'p_storage_co2_liquid': p_storage_co2_liquid,
    'T_storage_water': T_storage_water,
    'T_storage_TES': T_storage_TES,
    'T_amb': T_amb,
    'p_dome': p_dome,
    'k_dome': k_dome,
    'eta_comp': eta_comp,
    'eta_turb': eta_turb,
    'eta_pump': eta_pump
}


dome = CO2_battery(inputs, params, True)
dome.evaluate()
