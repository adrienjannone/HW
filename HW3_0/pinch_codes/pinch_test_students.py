"""
LELME2150 - Thermal cycles
Homework 3.0 - Pinch model heat exchanger

Test code for your function

@author: Mattéo Hauglustaine
@date: October 10th, 2025
"""

#
# ===IMPORT PACKAGES============================================================
#

from pinch_group_zoe import heat_exchanger

p_hs, T_hs_su, T_hs_ex = 2e5, 220 + 273.15, 150 + 273.15  # [Pa], [K]
T_cs_su, T_cs_ex = 125 + 273.15, 200 + 273.15  # [K]

pinch_target = 5  # [K]
wf_hot_side = "INCOMP::T72"  # Therminol 72 thermal oil
wf_cold_side = "Heptane"  # Refrigerant
p_cs_guess = 5e5  # Pa (5 bar)

inputs = p_hs, T_hs_su, T_hs_ex, T_cs_su, T_cs_ex
params = {
    "pinch_target": pinch_target,
    "wf_hot_side": wf_hot_side,
    "wf_cold_side": wf_cold_side,
    "p_cs_guess": p_cs_guess,
}
my_HX = heat_exchanger(inputs, params, True)
my_HX.evaluate()

eta_ex = my_HX.eta_transex
