"""
LELME2150 - Thermal cycles
Homework 2 - Gas turbines

Test code for your function

@author: Antoine Laterre
@date: September 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#

from gas_turbine_group_09 import gas_turbine

p_1, T_1 = 1e+5, 288.15 # [Pa], [K]
r_c = 18
T_3 = 1673.15 # [K]
eta_pi = 0.9 # [-]
k_mec = 0.015 # [-]
k_cc = 0.95
P_e = 230e+6

inputs = p_1,T_1,P_e
params =  {'T_3':       T_3,
           'r_c':       r_c,
           'eta_pi_c':  eta_pi,
           'eta_pi_t':  eta_pi,
           'k_cc':      k_cc,
           'k_mec':     k_mec,
           'air':       ['N2','O2'],
           'air_prop':  [0.79,0.21],
           'alkane':    [1,4]
           }
my_GT = gas_turbine(inputs,params,True)
my_GT.evaluate()

eta_en = my_GT.eta_toten