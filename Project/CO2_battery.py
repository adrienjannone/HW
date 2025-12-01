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
D2: Outlet of the evaporator and inlet of the HX
D6: Outlet of the HX and inlet of the compressor
D7: Outlet of the compressor and inlet of the turbine
D8: Outlet of the turbine and inlet of the PCHX
D9: Outlet of the PCHX and inlet of the CO2 gas storage dome (/!\ no PCHX in the study)
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
        self.k_TES = params['k_TES']
        self.eta_comp = params['eta_comp']
        self.eta_turb = params['eta_turb']
        self.eta_pump = params['eta_pump']
        self.eta_mec = params['eta_mec']
        self.eta_elec = params['eta_elec']
        self.pinch_TES0 = params['pinch_TES0']
        self.pinch_TES = params['pinch_TES']
        self.fluid = params['fluid']
        self.Pe = inputs
        self.plot = plot

        # Outside states
        self.p_ref = self.p_amb
        self.T_ref = self.T_amb
        self.h_ref = CP.PropsSI('H', 'P', self.p_ref, 'T', self.T_ref, self.fluid)
        self.s_ref = CP.PropsSI('S', 'P', self.p_ref, 'T', self.T_ref, self.fluid)

        # Dome states
        self.p_C1 = self.p_dome
        self.T_C1 = self.T_amb
        self.h_C1 = CP.PropsSI('H', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.s_C1 = CP.PropsSI('S', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.e_C1 = self.h_C1 - self.h_ref - self.T_ref * (self.s_C1 - self.s_ref)
        self.x_C1 = CP.PropsSI('Q', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
    
    def evaporator(self):
        # Counter-flow heat exchanger between CO2 and water
        # J'ai besoin de la temperature de l'eau en sortie du condensateur  
        pass

    def TES_charge(self):
        # TES pinch = 7.5 K
        self.T_D6 = self.T_storage_TES - self.pinch_TES
        self.p_D6 = self.p_D2 * self.k_TES
        self.h_D6 = CP.PropsSI('H', 'P', self.p_D6, 'T', self.T_D6, self.fluid)
        self.s_D6 = CP.PropsSI('S', 'P', self.p_D6, 'T', self.T_D6, self.fluid)
        self.e_D6 = self.h_D6 - self.h_ref - self.T_ref * (self.s_D6 - self.s_ref)
        self.x_D6 = CP.PropsSI('Q', 'P', self.p_D6, 'T', self.T_D6, self.fluid)


    def discharge_phase(self):
        # State 1 - Liquid CO2 storage
        self.p_D1 = self.p_storage_co2_liquid
        self.T_D1 = self.T_amb
        self.h_D1 = CP.PropsSI('H', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.s_D1 = CP.PropsSI('S', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.e_D1 = self.h_D1 - self.h_ref - self.T_ref * (self.s_D1 - self.s_ref)
        self.x_D1 = CP.PropsSI('Q', 'P', self.p_D1, 'T', self.T_D1, self.fluid)

        # State 2 - Evaporator outlet - HX inlet
        self.evaporator()

        # State 6 - HX outlet - Compressor inlet
        self.TES_charge()

        # State 9 - PCHX outlet - Dome inlet
        self.p_D9 = self.p_dome / self.k_dome

        # State 8 - Turbine outlet
        self.p_D8 = self.p_D9
        self.h_8DS = CP.PropsSI('H', 'S', self.s_D7, 'P', self.p_D8, self.fluid)
        self.h_D8 = self.h_D7 - self.eta_turb * (self.h_D7 - self.h_8DS)
        self.T_D8 = CP.PropsSI('T', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.s_D8 = CP.PropsSI('S', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.e_D8 = self.h_D8 - self.h_ref - self.T_ref * (self.s_D8 - self.s_ref)
        self.x_D8 = CP.PropsSI('Q', 'P', self.p_D8, 'T', self.T_D8, self.fluid)





        
      
    def evaluate(self):
        pass

