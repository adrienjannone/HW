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
        self.pinch_PCHX = params['pinch_PCHX']
        self.fluid = params['fluid']
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

        #evaporator TS0
        self.T_D2 = 20 + 273.15

        self.T_w_hot = 25+273.15
        self.T_w_cold = self.T_storage_water

        self.p_w = 1.1e5
        self.m_dot_r = 0
        self.measured_pinch = 0

        #TES
        self.T_D7 = 430 + 273.15
        self.T_TES_out = 25 + 273.15

        

        #TS0_in
        self.TS0_in = self.T_w_hot 
        self.p_TS0 = self.p_w
        self.h_TS0_in = CP.PropsSI('H', 'T', self.TS0_in, 'P', self.p_TS0, 'water')
        self.s_TS0_in = CP.PropsSI('S', 'T', self.TS0_in, 'P', self.p_TS0, 'water')
        self.e_TS0_in = self.h_TS0_in - self.h_ref - self.T_ref * (self.s_TS0_in - self.s_ref)
        
        #TS0_out
        self.TS0_out = self.T_w_cold
        self.h_TS0_out = CP.PropsSI('H', 'T', self.TS0_out, 'P', self.p_TS0, 'water')
        self.s_TS0_out = CP.PropsSI('S', 'T', self.TS0_out, 'P', self.p_TS0, 'water')
        self.e_TS0_out = self.h_TS0_out - self.h_ref - self.T_ref * (self.s_TS0_out - self.s_ref)

        self.set_ref()

    def set_ref(self):
        # CO2 H2O
        #s=0, h=0 @ 1 bar, 15°C = 288.15 K
        DCO2 = CP.PropsSI('DMOLAR','P',1e+5,'T',288.15,'CO2')
        CP.set_reference_state('CO2', 288.15, DCO2, 0, 0)

        Dw = CP.PropsSI('DMOLAR','Q',1,'T',273.15,'H2O')
        CP.set_reference_state('H2O', 273.15, Dw, 0, 0)

    def vaporator_efficiency(self):
        self.eta_transex_TES0 = (self.m_dot_CO2*(self.e_D2-self.e_D1))/(self.m_dot_TS0*(self.e_TS0_in - self.e_TS0_out))

    def TES_efficiency(self):
        # Si TES est de l'eau (improbable)
        #self.m_dot_TES_r = (self.h_D7 - self.h_D2) / (4180*(self.T_storage_TES - self.T_TES_out)) # J'ai inversé ?
        #T1_m = (self.T_storage_TES-self.T_TES_out)/(np.log(self.T_storage_TES/self.T_TES_out))
        #T2_m = (self.T_D7-self.T_D2)/(np.log(self.T_D7/self.T_D2)) 
        #self.eta_transex_TES = self.m_dot_TES_r* ((T2_m - self.T_ref)/T2_m)*(T1_m/(T1_m - self.T_ref))

        self.eta_transex_TES = (self.m_dot_CO2*(self.e_D7 - self.e_D2))/(self.m_dot_TES*(self.e_TES_in - self.e_TES_out))

    def PCHX_efficiency(self):
        #T1_m = (self.T_D8-self.T_D9)/(np.log(self.T_D8/self.T_D9))
        #T2_m = (self.T_w_hot-self.T_storage_water)/(np.log(self.T_w_hot/self.T_storage_water))
        #self.eta_transex = self.m_dot_r* ((T2_m - self.T_ref)/T2_m)*(T1_m/(T1_m - self.T_ref))
        self.eta_transex_PCHX = (self.m_dot_TS0*(self.e_TS0_in - self.e_TS0_out))/(self.m_dot_CO2*(self.e_D8 - self.e_D9))

    def mass_ratio_TS0(self, p_evap):
        h_hs_su = CP.PropsSI("H", "T", self.T_w_hot, "P", self.p_w, "water")  
        h_hs_ex = CP.PropsSI("H", "T", self.T_w_cold, "P", self.p_w, "water") 

        h_cs_su = CP.PropsSI("H", "T", self.T_D1, "P", p_evap, self.fluid)  
        h_cs_ex = CP.PropsSI("H", "T", self.T_D2, "P", p_evap, self.fluid)  

        self.h_hs = (h_hs_su - h_hs_ex)  
        self.h_cs = (h_cs_ex - h_cs_su)  
        self.m_dot_r = self.h_hs / self.h_cs
    
    def get_pinch_SAT_TS0(self, p_evap):
        self.mass_ratio_TS0(p_evap)
        T_c = CP.PropsSI("T", "P", p_evap, "Q", 0, self.fluid)  
        h_c = CP.PropsSI("H", "P", p_evap, "Q", 0, self.fluid)  
        h_cs_su = CP.PropsSI("H", "T", self.T_D1, "P", p_evap, self.fluid) 
        dh = h_c - h_cs_su  
        h_hs_ex = CP.PropsSI("H", "T", self.T_w_cold, "P", self.p_w, "water")  
        
        h_h = h_hs_ex + dh*self.m_dot_r  
        T_h = CP.PropsSI("T", "H", h_h, "P", self.p_w, "water")  
        pinch = T_h - T_c 
        self.measured_pinch = pinch
        return pinch
    
    def get_pinch_exit_TS0(self, p_evap): 
        return self.T_w_hot - self.T_D2

    def pinch_objective_TS0(self, p_evap):
        self.mass_ratio_TS0(p_evap)
        pinch = min(self.get_pinch_SAT_TS0(p_evap), self.get_pinch_exit_TS0(p_evap))
        return pinch - self.pinch_TES0
    
    def evaporator(self):
        p_cs_guess = 57e5 
        p_evap_solution = opt.fsolve(self.pinch_objective_TS0, p_cs_guess)[0]
        self.mass_ratio_TS0(p_evap_solution)
        self.p_D1 = p_evap_solution
        return
    
    def mass_ratio_TES(self, p_hs):
        h_hs_su = CP.PropsSI('H', 'T', self.T_storage_TES, 'P', p_hs,"water")
        h_hs_ex = CP.PropsSI('H', 'T', self.T_TES_out, 'P', p_hs, 'water')

        h_cs_su = CP.PropsSI('H', 'P', self.p_D2, 'T', self.T_D2, self.fluid)
        h_cs_ex = CP.PropsSI('H','P', self.p_D2, 'T', self.T_D7, self.fluid)
        self.h_hs = h_hs_su - h_hs_ex
        self.h_cs = h_cs_ex - h_cs_su
        self.m_dot_TES_r = self.h_hs / self.h_cs
    
    def get_pinch_SAT_TES(self, p_hs):
        self.mass_ratio_TES(p_hs)
        T_h = CP.PropsSI('T', 'P', p_hs, 'Q', 1, "water")
        h_h = CP.PropsSI('H', 'P', p_hs, 'Q', 1, "water")
        h_cs_su = CP.PropsSI('H', 'P', self.p_D2, 'T', self.T_D2, self.fluid)
        h_hs_ex = CP.PropsSI('H', 'T', self.T_TES_out, 'P', p_hs, 'water')
        dh = h_h - h_hs_ex
        h_c = h_cs_su + dh/self.m_dot_TES_r
        T_c = CP.PropsSI('T', 'H', h_c, 'P', self.p_D2, self.fluid)
        pinch = T_h - T_c

        return pinch
    
    def get_pinch_inlet_TES(self):
        return self.T_TES_out - self.T_D2
    
    def get_pinch_exit_TES(self):
        return self.T_storage_TES - self.T_D7
    
    def pinch_objective_TES(self, p_hs):
        self.mass_ratio_TES(p_hs)
        pinch = min(self.get_pinch_SAT_TES(p_hs), self.get_pinch_exit_TES(), self.get_pinch_inlet_TES())
        return pinch - self.pinch_TES
    

    def plotTQTES(self):
        """
        Hot source supply: TES storage
        Hot source exit: TES outlet
        Cold source supply: D2
        Cold source exit: D7
        """
        T_hs = np.linspace(self.T_TES_out, self.T_storage_TES, 100)
        T_cs = np.linspace(self.T_D2, self.T_D7, 100)
        h_hs = CP.PropsSI("H", "T", T_hs, "P", self.p_TES, "water") * 1e-3
        h_cs = CP.PropsSI("H", "T", T_cs, "P", self.p_D1, "CO2") *1e-3

        h_cs = h_cs * self.m_dot_TES_r # Scale
        h_hs = h_hs - np.ones_like(h_hs)*h_hs[0] # put h_hs[0] = 0
        h_cs = h_cs - np.ones_like(h_cs)*h_cs[0] # put h_cs[0] = 0
        h_hs /= h_hs[-1] # put h_hs[-1] = 1
        h_cs /= h_cs[-1] # put h_cs[-1] = 1

        plt.figure()
        plt.plot(h_hs, T_hs-273.15, label="Water side at p = {:.2f} bar".format(self.p_TES/1e5))
        plt.plot(h_cs, T_cs-273.15, label="CO2 side at p = {:.2f} bar".format(self.p_D1/1e5))
        plt.xlabel("Normalized cumulative heat transfer [-]")
        plt.ylabel("Temperature [°C]")
        plt.title("TES TQ diagram")
        plt.legend()
        plt.grid()
        plt.show()

    def plotTQTES0(self):
        """
        Hot source supply: water hot
        Hot source exit: water cold
        Cold source supply: D1
        Cold source exit: D2
        """
        T_hs = np.linspace(self.T_w_cold, self.T_w_hot, 100)
        T_cs = np.linspace(self.T_D1, self.T_D2, 100)
        h_hs = CP.PropsSI("H", "T", T_hs, "P", self.p_w, "water") * 1e-3
        h_cs = CP.PropsSI("H", "T", T_cs, "P", self.p_D1,"CO2") *1e-3

        h_cs = h_cs * self.m_dot_r # Scale
        h_hs = h_hs - np.ones_like(h_hs)*h_hs[0] # put h_hs[0] = 0
        h_cs = h_cs - np.ones_like(h_cs)*h_cs[0] # put h_cs[0] = 0
        h_hs /= h_hs[-1] # put h_hs[-1] = 1
        h_cs /= h_cs[-1] # put h_cs[-1] = 1

        plt.figure()
        plt.plot(h_hs, T_hs-273.15, label="Water side at p = {:.2f} bar".format(self.p_w/1e5))
        plt.plot(h_cs, T_cs-273.15, label="CO2 side at p = {:.2f} bar".format(self.p_D1/1e5))
        plt.xlabel("Normalized cumulative heat transfer [-]")
        plt.ylabel("Temperature [°C]")
        plt.title("Evaporator (TS0) TQ diagram")
        plt.legend()
        plt.grid()
        plt.show()

    def TES_discharge(self):

        #self.T_D7 = self.T_storage_TES - self.pinch_TES
        self.p_D7 = self.p_D2 * (1-self.k_TES)
        self.h_D7 = CP.PropsSI('H', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.s_D7 = CP.PropsSI('S', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.e_D7 = self.h_D7 - self.h_ref - self.T_ref * (self.s_D7 - self.s_ref)
        self.x_D7 = CP.PropsSI('Q', 'P', self.p_D7, 'T', self.T_D7, self.fluid)

        p_hs_guess = 166e5
        p_TES_solution = opt.fsolve(self.pinch_objective_TES, p_hs_guess)[0]
        self.mass_ratio_TES(p_TES_solution)
        self.p_TES = p_TES_solution
        return


    def discharge_phase(self):
        # State 0 - Liquid CO2 storage
        self.p_D0 = self.p_storage_co2_liquid
        self.T_D0 = self.T_amb
        self.h_D0 = CP.PropsSI('H', 'P', self.p_D0, 'T', self.T_D0, self.fluid)
        self.s_D0= CP.PropsSI('S', 'P', self.p_D0, 'T', self.T_D0, self.fluid)
        self.e_D0 = self.h_D0 - self.h_ref - self.T_ref * (self.s_D0 - self.s_ref)
        self.x_D0 = CP.PropsSI('Q', 'P', self.p_D0, 'T', self.T_D0, self.fluid)

        # State 1 - Evaporator inlet - TS0 inlet
        self.T_D1 = self.T_D0
        self.evaporator()
        # self.mass_ratio(self.p_D1)

        self.h_D1 = CP.PropsSI('H', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.s_D1= CP.PropsSI('S', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        self.e_D1 = self.h_D1 - self.h_ref - self.T_ref * (self.s_D1 - self.s_ref)
        self.x_D1 = CP.PropsSI('Q', 'P', self.p_D1, 'T', self.T_D1, self.fluid)
        

        # State 2 - Evaporator outlet - TS0 outlet/TES inlet
        self.p_D2 = self.p_D1 # TSO isobaric
        self.h_D2 = CP.PropsSI('H', 'P', self.p_D2, 'T', self.T_D2, self.fluid)
        self.s_D2= CP.PropsSI('S', 'P', self.p_D2, 'T', self.T_D2, self.fluid)
        self.e_D2 = self.h_D2 - self.h_ref - self.T_ref * (self.s_D2 - self.s_ref)
        self.x_D2 = CP.PropsSI('Q', 'P', self.p_D2, 'T', self.T_D2, self.fluid)

         # State 7 - TES outlet  
        self.TES_discharge()


        print("p_TES (bar):", self.p_TES*1e-5) 

        #TES_in
        self.T_TES_in = self.T_storage_TES
        self.h_TES_in = CP.PropsSI('H', 'T', self.T_TES_in, 'P', self.p_TES, 'water')
        self.s_TES_in = CP.PropsSI('S', 'T', self.T_TES_in, 'P', self.p_TES, 'water')
        self.e_TES_in = self.h_TES_in - self.h_ref - self.T_ref * (self.s_TES_in - self.s_ref)

        #TES_out
        self.h_TES_out = CP.PropsSI('H', 'T', self.T_TES_out, 'P', self.p_TES, 'water')
        self.s_TES_out = CP.PropsSI('S', 'T', self.T_TES_out, 'P', self.p_TES, 'water')
        self.e_TES_out = self.h_TES_out - self.h_ref - self.T_ref * (self.s_TES_out - self.s_ref)

        # State 9 - PCHX outlet - Dome inlet
        self.p_D9 = self.p_amb *(1+self.k_dome)
        # self.h_D9 = self.h_C1
        # self.T_D9 = CP.PropsSI('T', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        # self.s_D9 = CP.PropsSI('S', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        # self.x_D9 = CP.PropsSI('Q', 'P', self.p_D9, 'H', self.h_D9, self.fluid)
        # self.e_D9 = self.h_D9 - self.h_ref - self.T_ref * (self.s_D9 - self.s_ref)

        # State 8 - Turbine outlet
        self.p_D8 = self.p_D9
        self.h_8DS = CP.PropsSI('H', 'S', self.s_D7, 'P', self.p_D8, self.fluid)
        self.h_D8 = self.h_D7 - self.eta_turb * (self.h_D7 - self.h_8DS)
        self.T_D8 = CP.PropsSI('T', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.s_D8 = CP.PropsSI('S', 'P', self.p_D8, 'H', self.h_D8, self.fluid)
        self.e_D8 = self.h_D8 - self.h_ref - self.T_ref * (self.s_D8 - self.s_ref)
        self.x_D8 = CP.PropsSI('Q', 'P', self.p_D8, 'T', self.T_D8, self.fluid)

        # State 9 - PCHX outlet - Dome inlet
        self.T_D9 = self.T_storage_water + self.pinch_PCHX
        self.h_D9 = CP.PropsSI('H', 'P', self.p_D9, 'T', self.T_D9, self.fluid)
        self.s_D9 = CP.PropsSI('S', 'P', self.p_D9, 'T', self.T_D9, self.fluid)
        self.x_D9 = CP.PropsSI('Q', 'P', self.p_D9, 'T', self.T_D9, self.fluid)
        self.e_D9 = self.h_D9 - self.h_ref - self.T_ref * (self.s_D9 - self.s_ref)

    
        # Mass flow rate
        self.m_dot_CO2 = self.Pe / ((self.h_D7 - self.h_D8) * self.eta_mec * self.eta_elec) 
        self.m_dot_TS0 = self.m_dot_CO2/self.m_dot_r 
        #self.m_dot_TES_r = (self.h_TES_in-self.h_TES_out)/(self.h_D7 - self.h_D2)
        self.m_dot_TES = self.m_dot_CO2/self.m_dot_TES_r
        print(self.m_dot_TES)

        # Efficiencies
        self.eta_rotex = (self.h_D7 - self.h_D8) / (self.e_D7 - self.e_D8)
        self.vaporator_efficiency()
        self.TES_efficiency()
        self.PCHX_efficiency()

        # Losses
        # Wrong losses
        self.loss_rotex = self.m_dot_CO2 * ((self.e_D7 - self.e_D8) - (self.h_D7 - self.h_D8))
        self.loss_evaporator = self.m_dot_TS0 * (self.e_TS0_in- self.e_TS0_out) - self.m_dot_CO2 * (self.e_D2 - self.e_D1)
        self.loss_TES = self.m_dot_TES * (self.e_TES_in - self.e_TES_out) - self.m_dot_CO2 * (self.e_D7 - self.e_D2)
        self.loss_PCHX = self.m_dot_CO2 * (self.e_D8 - self.e_D9) - self.m_dot_TS0 * (self.e_TS0_in - self.e_TS0_out)
        W_shaft = self.Pe / (self.eta_mec * self.eta_elec)
        self.loss_mec = W_shaft * (1 - self.eta_mec)
        self.loss_elec = W_shaft * self.eta_mec * (1 - self.eta_elec)

      
    def evaluate(self):
        self.discharge_phase()
        if self.plot:
            self.fig_pie_ex()
            self.plotTQTES()
            self.plotTQTES0()
        pass

    def fig_pie_ex(self):
        label = ['Rotex losses', 'Evaporator losses', 'TES losses', 'PCHX losses', 'Mechanical losses', 'Electrical losses', "Effective Power"]
        sizes = [self.loss_rotex, self.loss_evaporator, self.loss_TES, self.loss_PCHX, self.loss_mec, self.loss_elec, self.Pe]
        plt.figure(figsize=(8, 8))
        plt.pie(sizes, labels=label, autopct='%1.1f%%', startangle=140)
        plt.title('Exergy Losses Distribution in CO2 Battery Discharge Phase')
        plt.show()

    def print_results(self):
        print("Results of the CO2 battery discharge phase:")
        print("Energy produced: {:.2f} kW".format(self.Pe*1e-3))
        print("--- States ---")
        print("State D0: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D0*1e-5, self.T_D0-273.15, self.h_D0*1e-3, self.s_D0*1e-3, self.x_D0))
        print("State D1: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D1*1e-5, self.T_D1-273.15, self.h_D1*1e-3, self.s_D1*1e-3, self.x_D1))
        print("State D2: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D2*1e-5, self.T_D2-273.15, self.h_D2*1e-3, self.s_D2*1e-3, self.x_D2))
        print("State D7: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D7*1e-5, self.T_D7-273.15, self.h_D7*1e-3, self.s_D7*1e-3, self.x_D7))
        print("State D8: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D8*1e-5, self.T_D8-273.15, self.h_D8*1e-3, self.s_D8*1e-3, self.x_D8))
        print("State D9: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D9*1e-5, self.T_D9-273.15, self.h_D9*1e-3, self.s_D9*1e-3, self.x_D9))
        print("State C1: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_C1*1e-5, self.T_C1-273.15, self.h_C1*1e-3, self.s_C1*1e-3, self.x_C1))
        print("Water storage temperature: {:.2f} °C".format(self.T_storage_water-273.15))
        print("Water heated temperature: {:.2f} °C".format(self.T_w_hot-273.15))
        print("TES hot outlet temperature: {:.2f} °C".format(self.T_storage_TES-273.15))
        print("TES cold outlet temperature: {:.2f} °C".format(self.T_TES_out-273.15))
        print("--- Flow rates ---")
        print("Mass flow rate of CO2: {:.2f} kg/s".format(self.m_dot_CO2))
        print("Mass flow rate of TS0: {:.2f} kg/s".format(self.m_dot_TS0))
        print("--- Efficiencies ---")
        print("Rotex efficiency: {:.2f} %".format(self.eta_rotex*100))
        print("Heat exchanger efficiency TES0: {:.2f} %".format(self.eta_transex_TES0*100))
        print("Heat exchanger efficiency TES: {:.2f} %".format(self.eta_transex_TES*100))
        print("Heat exchanger efficiency PCHX: {:.2f} %".format(self.eta_transex_PCHX*100))
        print("--- Losses ---")
        print("Rotex losses: {:.2f} kW".format(self.loss_rotex*1e-3))
        print("Evaporator losses: {:.2f} kW".format(self.loss_evaporator*1e-3))
        print("TES losses: {:.2f} kW".format(self.loss_TES*1e-3))
        print("PCHX losses: {:.2f} kW".format(self.loss_PCHX*1e-3))
        print("Mechanical losses: {:.2f} kW".format(self.loss_mec*1e-3))
        print("Electrical losses: {:.2f} kW".format(self.loss_elec*1e-3))
