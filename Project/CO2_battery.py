"""
LELME2150 - Thermal cycles
Project - CO2 battery

Test code for your CO2 battery model

@authors : Igor Grégoire, Adrien Jannone
@date: November 30th, 2025
"""


import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

"""
D1: Outlet of the liquid CO2 storage tank and inlet of the evaporator
D2: Outlet of the evaporator and inlet of the HXD
D3: Outlet of the HXD and inlet of the compressor
D7: Outlet of the compressor and inlet of the turbine
D8: Outlet of the turbine and inlet of the PCHX
D9: Outlet of the PCHX and inlet of the CO2 gas storage dome
C1: Storage dome 
"""

class CO2_battery(object):
    def __init__(self, inputs, params, plot):
        self.p_storage_co2_liquid = params['p_storage_co2_liquid']
        self.T_storage_water = params['T_storage_water']
        self.T_storage_TES = params['T_storage_TES']
        self.T_amb = params['T_amb']
        self.p_amb = params['p_amb']
        self.p_dome = params['p_dome']
        self.k_dome = params['k_dome']
        self.eta_comp = params['eta_comp']
        self.eta_turb = params['eta_turb']
        self.eta_pump = params['eta_pump']
        self.fluid = params['fluid']
        self.Pe = inputs
        self.plot = plot

        # Outside states
        self.p_0 = self.p_amb
        self.T_0 = self.T_amb
        self.h_0 = CP.PropsSI('H', 'P', self.p_0, 'T', self.T_0, self.fluid)
        self.s_0 = CP.PropsSI('S', 'P', self.p_0, 'T', self.T_0, self.fluid)

        # Dome states
        self.p_C1 = self.p_dome
        self.T_C1 = self.T_amb
        self.h_C1 = CP.PropsSI('H', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.s_C1 = CP.PropsSI('S', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.e_C1 = self.h_C1 - self.h_0 - self.T_0 * (self.s_C1 - self.s_0)

    def discharge_phase(self):
        # Liquid CO2 storage
        self.p_D1 = self.p_storage_co2_liquid
        self.T_D1 = self.T_amb
        self.h_D1 = CP.PropsSI('H', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.s_D1 = CP.PropsSI('S', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.e_D1 = self.h_D1 - self.h_0 - self.T_0 * (self.s_D1 - self.s_0)

        # PCHX outlet - Dome inlet
        self.p_D9 = self.p_dome / self.k_dome
        


        
      
    def evaluate(self):
        pass

