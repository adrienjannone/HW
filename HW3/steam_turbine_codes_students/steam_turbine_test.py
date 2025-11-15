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
import numpy as np
HW1 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HW1,'corrections'))
import matplotlib.pyplot as plt

from steam_turbine_group_09 import steam_turbine

p_1, p_3, p_4   = 100e+5, 310e+5, 70e+5
p_ref, T_ref    = 1e+5,   288.15
T_max           = 838.15
T_cd_out        = 305.15 # Condenseur
T_cd_subcool    = 1 # Condeuseur sous-refroidissement
T_pinch_sc      = 0 # Pinch parfait
T_pinch_ex      = 0 # Pinch parfait
T_pinch_cd      = 0 # Pinch parfait
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

# my_ST = steam_turbine(inputs,params,False)
# my_ST.evaluate()

# eta_en = my_ST.eta_cyclen


# Impact de p4
# p4s = np.linspace(50e5, 90e5)
# x6s = []
# eta_cyclen = []
# eta_cyclex = []
# for p4 in p4s:
#     params['p_4'] = p4
#     my_ST = steam_turbine(inputs,params,False)
#     my_ST.evaluate()
#     eta_cyclen.append(my_ST.eta_cyclen)
#     eta_cyclex.append(my_ST.eta_cyclex)
#     x6s.append(my_ST.x_6)

# plt.figure()
# plt.plot(p4s, eta_cyclen, label=r'$\eta_{cyclen}$')
# plt.plot(p4s, x6s, label=r'$x_6$')
# plt.plot(p4s, eta_cyclex, label=r'$\eta_{cyclex}$')
# plt.legend()
# plt.axhline(y=0.88, color='r', linestyle='--')
# plt.xlabel('p4 (Pa)')
# plt.ylabel(r'\eta, x6')
# plt.title(r'Variation of $\eta_{cyclen}$, $x_6$ and $\eta_{cyclex}$ with $p_4$')
# plt.grid()
# plt.savefig('impact_p4.svg')

Tmaxs = np.linspace(700.15, 2000.15)
cyclen = []
cyclex = []
flowrates = []
x6s = []
for Tmax in Tmaxs:
    params['T_max'] = Tmax
    my_ST = steam_turbine(inputs,params,False)
    my_ST.evaluate()
    cyclen.append(my_ST.eta_cyclen)
    cyclex.append(my_ST.eta_cyclex)
    flowrates.append(my_ST.dotm_tot)
    x6s.append(my_ST.x_6)
plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()
ax1.plot(Tmaxs, cyclen, 'r', label=r'$\eta_{cyclen}$')
ax1.plot(Tmaxs, cyclex, 'b', label=r'$\eta_{cyclex}$')
ax1.plot(Tmaxs, x6s, 'orange', label=r'$x_6$')
ax2.plot(Tmaxs, flowrates, 'g', label=r'$\dot{m}_{tot}$')
ax1.set_xlabel(r'$T_{max}$ (K)')
ax1.set_ylabel(r'$\eta$')
ax2.set_ylabel(r'$\dot{m}_{tot}$ (kg/s)')
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc=0)
plt.title(r'Variation of $\eta_{cyclen}$, $\eta_{cyclex}$ and $\dot{m}_{tot}$ with $T_{max}$')
plt.grid()
plt.savefig('impact_Tmax.svg')
