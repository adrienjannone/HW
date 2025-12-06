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
#from scipy.optimize import fsolve

"""
D1: Outlet of the liquid CO2 storage tank and inlet of the evaporator
D2: Outlet of the evaporator and inlet of the HX
D6: Outlet of the HX and inlet of the turbine
D8: Outlet of the turbine and inlet of the PCHX
D9: Outlet of the PCHX and inlet of the CO2 gas storage dome (/!\ no PCHX in the study)
C1: Storage dome 
"""


class CO2_battery(object):
    def __init__(self, inputs, params, plot):
        self.p_storage_co2_liquid = params['p_storage_co2_liquid']# = T_D0
        self.T_storage_water = params['T_storage_water']# = T_w_cool
        self.T_storage_TES = params['T_storage_TES'] # = T_TES_hot
        self.T_amb = params['T_amb'] # = T_D0 =  T_C1
        self.p_amb = params['p_amb'] # = p_C1
        self.k_dome = params['k_dome'] #= p_D9/p_C1 = 1-k_dome
        self.k_TES = params['k_TES'] # = ?
        self.eta_comp = params['eta_comp']
        self.eta_turb = params['eta_turb']
        self.eta_pump = params['eta_pump']
        self.eta_mec = params['eta_mec']
        self.eta_elec = params['eta_elec']
        self.pinch_TES0 = params['pinch_TES0']
        self.pinch_TES = params['pinch_TES']
        self.fluid = params['fluid']
        self.m_dot_CO2 = params['m_dot_cycle']
        self.Pe = inputs
        self.plot = plot

        # Outside states
        self.p_ref = self.p_amb
        self.T_ref = self.T_amb
        self.h_ref = CP.PropsSI('H', 'P', self.p_ref, 'T', self.T_ref, self.fluid)
        self.s_ref = CP.PropsSI('S', 'P', self.p_ref, 'T', self.T_ref, self.fluid)

        # Dome states
        self.p_C1 = self.p_amb
        self.T_C1 = self.T_amb
        self.h_C1 = CP.PropsSI('H', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.s_C1 = CP.PropsSI('S', 'P', self.p_C1, 'T', self.T_C1, self.fluid)
        self.e_C1 = self.h_C1 - self.h_ref - self.T_ref * (self.s_C1 - self.s_ref)
        self.x_C1 = CP.PropsSI('Q', 'P', self.p_C1, 'T', self.T_C1, self.fluid)

        #evaporator
        self.T_D2 = 20 + 273.15
        self.T_w_in = 25+273.15
        self.T_w_out = self.T_storage_water
        self.p_w = 1.1e5
        self.m_dot_r = 0
        self.measured_pinch = 0

    def mass_ratio(self, p_evap):
        h_hs_su = CP.PropsSI("H", "T", self.T_w_in, "P", self.p_w, "water")  
        h_hs_ex = CP.PropsSI("H", "T", self.T_w_out, "P", self.p_w, "water") 

        h_cs_su = CP.PropsSI("H", "T", self.T_D1, "P", p_evap, self.fluid)  
        h_cs_ex = CP.PropsSI("H", "T", self.T_D2, "P", p_evap, self.fluid)  

        self.h_hs = (h_hs_su - h_hs_ex)  
        self.h_cs = (h_cs_ex - h_cs_su)  
        self.m_dot_r = self.h_hs / self.h_cs
    
    def get_pinch_SAT(self, p_evap):
        self.mass_ratio(p_evap)
        T_c = CP.PropsSI("T", "P", p_evap, "Q", 0, self.fluid)  
        h_c = CP.PropsSI("H", "P", p_evap, "Q", 0, self.fluid)  
        h_cs_su = CP.PropsSI("H", "T", self.T_D1, "P", p_evap, self.fluid) 
        dh = h_c - h_cs_su  
        h_hs_ex = CP.PropsSI("H", "T", self.T_w_out, "P", self.p_w, "water")  
        
        h_h = h_hs_ex + dh*self.m_dot_r  
        T_h = CP.PropsSI("T", "H", h_h, "P", self.p_w, "water")  
        pinch = T_h - T_c 
        self.measured_pinch = pinch
        return pinch
    
    def get_pinch_exit(self, p_evap): 
        return self.T_w_in - self.T_D2

    def pinch_objective(self, p_evap):
        pinch = min(self.get_pinch_SAT(p_evap), self.get_pinch_exit(p_evap))
        return pinch - self.pinch_TES0
    
    def evaporator(self):
        self.T_D1 = self.T_D0
        p_cs_guess = 57e5 
        self.mass_ratio(p_cs_guess)
        p_evap_solution = opt.fsolve(self.pinch_objective, p_cs_guess)[0]
        self.p_D1 = p_evap_solution
        print(p_evap_solution)
        return 
        

    def TES_discharge(self):
        self.T_D7 = self.T_storage_TES - self.pinch_TES
        self.p_D7 = self.p_D2 * (1-self.k_TES)
        self.h_D7 = CP.PropsSI('H', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.s_D7 = CP.PropsSI('S', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.e_D7 = self.h_D7 - self.h_ref - self.T_ref * (self.s_D7 - self.s_ref)
        self.x_D7 = CP.PropsSI('Q', 'P', self.p_D7, 'T', self.T_D7, self.fluid)


    def discharge_phase(self):
        # State 0 - Liquid CO2 storage
        self.p_D0 = self.p_storage_co2_liquid
        self.T_D0 = self.T_amb
        self.h_D0 = CP.PropsSI('H', 'P', self.p_D0, 'T', self.T_D0, self.fluid)
        self.s_D0= CP.PropsSI('S', 'P', self.p_D0, 'T', self.T_D0, self.fluid)
        self.e_D0 = self.h_D0 - self.h_ref - self.T_ref * (self.s_D0 - self.s_ref)
        self.x_D0 = CP.PropsSI('Q', 'P', self.p_D0, 'T', self.T_D0, self.fluid)

        # State 1 - Evaporator inlet - TS0 inlet
        self.evaporator()
        self.h_D1 = CP.PropsSI('H', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.s_D1= CP.PropsSI('S', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.e_D1 = self.h_D1 - self.h_ref - self.T_ref * (self.s_D1 - self.s_ref)
        self.x_D1 = CP.PropsSI('Q', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        

        # State 2 - Evaporator outlet - TS0 outlet/TES inlet
        self.p_D2 = self.p_D1 # TSO isobaric
        self.h_D2 = CP.PropsSI('H', 'P', self.p_D2, 'T', self.T_D2, self.fluid)
        self.s_D2= CP.PropsSI('S', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.e_D2 = self.h_D2 - self.h_ref - self.T_ref * (self.s_D2 - self.s_ref)
        self.x_D2 = CP.PropsSI('Q', 'P', self.p_D2, 'T', self.T_D2, self.fluid)

         # State 7 - TES outlet  
        self.TES_discharge()

        # State 9 - PCHX outlet - Dome inlet
        self.p_D9 = self.p_amb *(1+self.k_dome)
        self.h_D9 = self.h_C1
        self.T_D9 = CP.PropsSI('T', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        self.s_D9 = CP.PropsSI('S', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        self.x_D9 = CP.PropsSI('Q', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        self.e_D9 = self.h_D9 - self.h_ref - self.T_ref * (self.s_D9 - self.s_ref)

        # State 8 - Turbine outlet
        self.p_D8 = self.p_D9
        self.h_8DS = CP.PropsSI('H', 'S', self.s_D7, 'P', self.p_D8, self.fluid)
        self.h_D8 = self.h_D7 - self.eta_turb * (self.h_D7 - self.h_8DS)
        self.T_D8 = CP.PropsSI('T', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.s_D8 = CP.PropsSI('S', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.e_D8 = self.h_D8 - self.h_ref - self.T_ref * (self.s_D8 - self.s_ref)
        self.x_D8 = CP.PropsSI('Q', 'P', self.p_D8, 'T', self.T_D8, self.fluid)

    
        # Mass flow rate
        self.m_dot_CO2 = self.Pe / ((self.h_D7 - self.h_D8) * self.eta_mec * self.eta_elec)
        print(self.m_dot_CO2) # !!!! = 54 kg/s dans le paper

        self.m_dot_TS0 = self.m_dot_CO2/self.m_dot_r # !!!! = 700 kg/s dans le paper
        print(self.m_dot_TS0)

        
      
    def evaluate(self):
        self.discharge_phase()
        pass

