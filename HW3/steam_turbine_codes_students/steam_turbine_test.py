"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Test code for your function

@author: Antoine Laterre
@date: October 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#

import os
import sys
HW1 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HW1,'corrections'))

from steam_turbine_group_09 import steam_turbine

p_1, p_3, p_4   = 100e+5, 310e+5, 70e+5
p_ref, T_ref    = 1e+5,   288.15
T_max           = 838.15
T_cd_out        = 305.15
T_cd_subcool    = 1
T_pinch_sc      = 0
T_pinch_ex      = 0
T_pinch_cd      = 0
T_drum          = 421.85
x_6min          = 0.88
eta_is_HP       = 0.92
eta_is_LP       = 0.90
eta_mec         = 0.98
eta_pump        = 0.85
k_heat          = 31/35
P_e             = 288e+6

inputs = P_e
params =  {'p_1':          p_1,
           'p_3':          p_3,
           'p_4':          p_4,
           'p_ref':        p_ref,
           'T_ref':        T_ref,
           'T_max':        T_max,
           'T_cd_out':     T_cd_out,
           'T_cd_subcool': T_cd_subcool,
           'T_pinch_sc':   T_pinch_sc,
           'T_pinch_ex':   T_pinch_ex,
           'T_pinch_cd':   T_pinch_cd,
           'T_drum':       T_drum,
           'x_6min':       x_6min,
           'eta_is_HP':    eta_is_HP,
           'eta_is_LP':    eta_is_LP,
           'eta_mec':      eta_mec,
           'eta_pump':     eta_pump,
           'k_heat':       k_heat
           }

my_ST = steam_turbine(inputs,params,True)
my_ST.evaluate()

eta_en = my_ST.eta_cyclen