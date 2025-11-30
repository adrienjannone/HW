"""
LELME2150 - Thermal cycles
Project - CO2 battery

Test code for your CO2 battery model

@authors : Igor Grégoire, Adrien Jannone
@date: November 30th, 2025
"""

#
# ===IMPORT PACKAGES============================================================

from CO2_battery import CO2_battery


p_0 = 1e+5  # Pa
T_0 = 15 + 273.15  # K
p_2 = 70e+5  # Pa
T_4 = T_0
eta_is_comp = 0.85
eta_is_turb = 0.93
eta_mec = 0.99
eta_elec = 0.985
pinch = 2  # K #pinch of single heat exchanger

params = {'p_0': p_0,   
          'T_0': T_0,
          'p_2': p_2,
          'eta_is_comp': eta_is_comp,
          'eta_is_turb': eta_is_turb,
          'eta_mec': eta_mec,
          'eta_elec': eta_elec,
          'pinch': pinch
          }

my_battery = CO2_battery(params, True)
my_battery.evaluate()
