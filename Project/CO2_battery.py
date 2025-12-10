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

        #evaporator TS0
        self.T_D2 = 20 + 273.15

        self.T_w_hot = 25+273.15
        self.T_w_cold = self.T_storage_water

        self.p_w = 1.1e5
        self.m_dot_r = 0
        self.measured_pinch = 0
        
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

        #TES
        self.T_TES1_in = 42+273.15
        self.T_TES1_out = 25+273.15
        self.p_TES1    = 1.1e5
        self.T_TES2_in = 97+273.15
        self.T_TES2_out = 42+273.15
        self.p_TES2    = 1.1e5
        self.T_TES3_in = 300+273.15
        self.T_TES3_out = 97+273.15
        self.p_TES3    = 86e5
        self.T_TES4_in = self.T_storage_TES
        self.T_TES4_out = 300+273.15
        self.p_TES4    = 1.1e5

    def vaporator_efficiency(self):
        self.eta_transex_TES0 = (self.m_dot_CO2*(self.e_D2-self.e_D1))/(self.m_dot_TS0*(self.e_TS0_in - self.e_TS0_out))

    def TES_efficiency(self):
        self.eta_transex_TES1 = (self.m_dot_CO2*(self.e_D3 - self.e_D2))/(self.m_dot_TSE1*(self.e_TES1_in - self.e_TES1_out))
        self.eta_transex_TES2 = (self.m_dot_CO2*(self.e_D4 - self.e_D3))/(self.m_dot_TSE2*(self.e_TES2_in - self.e_TES2_out))
        self.eta_transex_TES3 = (self.m_dot_CO2*(self.e_D5 - self.e_D4))/(self.m_dot_TSE3*(self.e_TES3_in - self.e_TES3_out))
        self.eta_transex_TES4 = (self.m_dot_CO2*(self.e_D7 - self.e_D5))/(self.m_dot_TSE4*(self.e_TES4_in - self.e_TES4_out))


    def PCHX_efficiency(self):
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
    
    
    def TES1(self, TES1_in, TES1_out, p_TES1):
        self.h_TES1_in = CP.PropsSI('H', 'T', TES1_in, 'P', p_TES1, 'water')
        print(self.h_TES1_in)
        self.s_TES1_in = CP.PropsSI('S', 'T', TES1_in, 'P', p_TES1, 'water')
        print(self.s_TES1_in)
        self.e_TES1_in = self.h_TES1_in - self.h_ref - self.T_ref * (self.s_TES1_in - self.s_ref)
        self.h_TES1_out = CP.PropsSI('H', 'T', TES1_out, 'P', p_TES1, 'water')
        self.s_TES1_out = CP.PropsSI('S', 'T', TES1_out, 'P', p_TES1, 'water')
        self.x_TES1_out = CP.PropsSI('Q', 'T', TES1_out, 'P', p_TES1, 'water')
        print(self.x_TES1_out)
        self.e_TES1_out = self.h_TES1_out - self.h_ref - self.T_ref * (self.s_TES1_out - self.s_ref)

        self.T_D3 = TES1_in - self.pinch_TES
        self.p_D3 = self.p_D2
        self.h_D3 = CP.PropsSI('H', 'P', self.p_D3, 'T', self.T_D3, self.fluid)
        self.s_D3= CP.PropsSI('S', 'P', self.p_D3, 'T', self.T_D3, self.fluid)
        self.e_D3 = self.h_D3 - self.h_ref - self.T_ref * (self.s_D3 - self.s_ref)
        self.x_D3 = CP.PropsSI('Q', 'P', self.p_D3, 'T', self.T_D3, self.fluid)

        return

    def TES2(self, TES2_in, TES2_out, p_TES2):
        self.h_TES2_in = CP.PropsSI('H', 'T', TES2_in, 'P', p_TES2, 'water')
        self.s_TES2_in = CP.PropsSI('S', 'T', TES2_in, 'P', p_TES2, 'water')
        self.e_TES2_in = self.h_TES2_in - self.h_ref - self.T_ref * (self.s_TES2_in - self.s_ref)
        self.h_TES2_out = CP.PropsSI('H', 'T', TES2_out, 'P', p_TES2, 'water')
        self.s_TES2_out = CP.PropsSI('S', 'T', TES2_out, 'P', p_TES2, 'water')
        self.e_TES2_out = self.h_TES2_out - self.h_ref - self.T_ref * (self.s_TES2_out - self.s_ref)

        self.T_D4 = TES2_in - self.pinch_TES
        self.p_D4 = self.p_D3
        self.h_D4 = CP.PropsSI('H', 'P', self.p_D4, 'T', self.T_D4, self.fluid)
        self.s_D4= CP.PropsSI('S', 'P', self.p_D4, 'T', self.T_D4, self.fluid)
        self.e_D4 = self.h_D4 - self.h_ref - self.T_ref * (self.s_D4 - self.s_ref)
        self.x_D4 = CP.PropsSI('Q', 'P', self.p_D4, 'T', self.T_D4, self.fluid)

        return
    
    def TES3(self, TES3_in, TES3_out, p_TES3):
        self.h_TES3_in = CP.PropsSI('H', 'T', TES3_in, 'P', p_TES3, 'water')
        self.s_TES3_in = CP.PropsSI('S', 'T', TES3_in, 'P', p_TES3, 'water')
        self.e_TES3_in = self.h_TES3_in - self.h_ref - self.T_ref * (self.s_TES3_in - self.s_ref)
        self.h_TES3_out = CP.PropsSI('H', 'T', TES3_out, 'P', p_TES3, 'water')
        self.s_TES3_out = CP.PropsSI('S', 'T', TES3_out, 'P', p_TES3, 'water')
        self.e_TES3_out = self.h_TES3_out - self.h_ref - self.T_ref * (self.s_TES3_out - self.s_ref)

        self.T_D5 = TES3_in - self.pinch_TES
        self.p_D5 = self.p_D4
        self.h_D5 = CP.PropsSI('H', 'P', self.p_D5, 'T', self.T_D5, self.fluid)
        self.s_D5= CP.PropsSI('S', 'P', self.p_D5, 'T', self.T_D5, self.fluid)
        self.e_D5 = self.h_D5 - self.h_ref - self.T_ref * (self.s_D5 - self.s_ref)
        self.x_D5 = CP.PropsSI('Q', 'P', self.p_D5, 'T', self.T_D5, self.fluid)

        return
    
    def TES4(self, TES4_in, TES4_out, p_TES4):
        self.h_TES4_in = CP.PropsSI('H', 'T', TES4_in, 'P', p_TES4, 'INCOMP::NaK')
        self.s_TES4_in = CP.PropsSI('S', 'T', TES4_in, 'P', p_TES4, 'INCOMP::NaK')
        self.e_TES4_in = self.h_TES4_in - self.h_ref - self.T_ref * (self.s_TES4_in - self.s_ref)
        self.h_TES4_out = CP.PropsSI('H', 'T', TES4_out, 'P', p_TES4, 'INCOMP::NaK')
        self.s_TES4_out = CP.PropsSI('S', 'T', TES4_out, 'P', p_TES4, 'INCOMP::NaK')
        self.e_TES4_out = self.h_TES4_out - self.h_ref - self.T_ref * (self.s_TES4_out - self.s_ref)

        self.T_D7 = TES4_in - self.pinch_TES
        self.p_D7 = self.p_D5
        self.h_D7 = CP.PropsSI('H', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.s_D7= CP.PropsSI('S', 'P', self.p_D7, 'T', self.T_D7, self.fluid)
        self.e_D7 = self.h_D7 - self.h_ref - self.T_ref * (self.s_D7 - self.s_ref)
        self.x_D7 = CP.PropsSI('Q', 'P', self.p_D7, 'T', self.T_D7, self.fluid)

        return
    

    def TES_discharge(self):
        self.TES1(self.T_TES1_in, self.T_TES1_out, self.p_TES1)
        self.TES2(self.T_TES2_in, self.T_TES2_out, self.p_TES2)
        self.TES3(self.T_TES3_in, self.T_TES3_out, self.p_TES3)
        self.TES4(self.T_TES4_in, self.T_TES4_out, self.p_TES4)
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


        
        # State 9 - PCHX outlet - Dome inlet
        self.p_D9 = self.p_amb *(1+self.k_dome)

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
        self.m_dot_TES1_r = (self.h_TES1_in-self.h_TES1_out)/(self.h_D3 - self.h_D2)
        self.m_dot_TES2_r = (self.h_TES2_in-self.h_TES2_out)/(self.h_D4 - self.h_D3)
        self.m_dot_TES3_r = (self.h_TES3_in-self.h_TES3_out)/(self.h_D5 - self.h_D4)
        self.m_dot_TES4_r = (self.h_TES4_in-self.h_TES4_out)/(self.h_D7 - self.h_D5)
        self.m_dot_TSE1 = self.m_dot_CO2/self.m_dot_TES1_r
        self.m_dot_TSE2 = self.m_dot_CO2/self.m_dot_TES2_r
        self.m_dot_TSE3 = self.m_dot_CO2/self.m_dot_TES3_r
        self.m_dot_TSE4 = self.m_dot_CO2/self.m_dot_TES4_r
        

        # Efficiencies
        self.eta_rotex = (self.h_D7 - self.h_D8) / (self.e_D7 - self.e_D8)
        self.vaporator_efficiency()
        self.TES_efficiency()
        self.PCHX_efficiency()

        # Losses
        self.loss_rotex = self.m_dot_CO2 * ((self.e_D7 - self.e_D8) - (self.h_D7 - self.h_D8))
        self.loss_evaporator = self.m_dot_TS0 * (self.e_TS0_in- self.e_TS0_out) - self.m_dot_CO2 * (self.e_D2 - self.e_D1)
        self.loss_TES1 = self.m_dot_TSE1 * (self.e_TES1_in - self.e_TES1_out) - self.m_dot_CO2 * (self.e_D3 - self.e_D2)
        self.loss_TES2 = self.m_dot_TSE2 * (self.e_TES2_in - self.e_TES2_out) - self.m_dot_CO2 * (self.e_D4 - self.e_D3)
        self.loss_TES3 = self.m_dot_TSE3 * (self.e_TES3_in - self.e_TES3_out) - self.m_dot_CO2 * (self.e_D5 - self.e_D4)
        self.loss_TES4 = self.m_dot_TSE4 * (self.e_TES4_in - self.e_TES4_out) - self.m_dot_CO2 * (self.e_D7 - self.e_D5)
        self.loss_PCHX = self.m_dot_CO2 * (self.e_D8 - self.e_D9) - self.m_dot_TS0 * (self.e_TS0_in - self.e_TS0_out)
        W_shaft = self.Pe / (self.eta_mec * self.eta_elec)
        self.loss_mec = W_shaft * (1 - self.eta_mec)
        self.loss_elec = W_shaft * self.eta_mec * (1 - self.eta_elec)


        print(CP.PropsSI('Tcrit', 'Water')-273.15)
        print("Water evaporation pressure at T_TES3_in",CP.PropsSI('P', 'T', 300+273.15, 'Q', 0, "water"))
        #print("Water evaporation pressure at T_TES4_in",CP.PropsSI('P', 'T', self.T_storage_TES, 'Q', 0, "water"))

      
    def evaluate(self):
        self.discharge_phase()
        pass

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
        plt.plot(h_hs, T_hs-273.15, label="Hot side at p = {:.2f} bar".format(self.p_w/1e5))
        plt.plot(h_cs, T_cs-273.15, label="Cold side at p = {:.2f} bar".format(self.p_D1/1e5))
        plt.xlabel("Normalized cumulative heat transfer [-]")
        plt.ylabel("Temperature [°C]")
        plt.title("Heat exchanger TQ diagram")
        plt.legend()
        plt.grid()
        plt.show()


    def fig_pie_ex(self):
        label = ['Rotex losses', 'Evaporator losses (TS0)', 'TS1 losses', 'TS2 losses', 'TS3 losses', 'TS4 losses', 'PCHX losses', 'Mechanical losses', 'Electrical losses', "Effective Power"]
        sizes = [self.loss_rotex, self.loss_evaporator, self.loss_TES1, self.loss_TES2, self.loss_TES3, self.loss_TES4, self.loss_PCHX, self.loss_mec, self.loss_elec, self.Pe]
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
        print("State D3: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D3*1e-5, self.T_D3-273.15, self.h_D3*1e-3, self.s_D3*1e-3, self.x_D3))
        print("State D4: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D4*1e-5, self.T_D4-273.15, self.h_D4*1e-3, self.s_D4*1e-3, self.x_D4))
        print("State D5: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D5*1e-5, self.T_D5-273.15, self.h_D5*1e-3, self.s_D5*1e-3, self.x_D5))
        print("State D7: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D7*1e-5, self.T_D7-273.15, self.h_D7*1e-3, self.s_D7*1e-3, self.x_D7))
        print("State D8: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D8*1e-5, self.T_D8-273.15, self.h_D8*1e-3, self.s_D8*1e-3, self.x_D8))
        print("State D9: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_D9*1e-5, self.T_D9-273.15, self.h_D9*1e-3, self.s_D9*1e-3, self.x_D9))
        print("State C1: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K, x = {:.2f}".format(self.p_C1*1e-5, self.T_C1-273.15, self.h_C1*1e-3, self.s_C1*1e-3, self.x_C1))
        print("State TS0 inlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TS0*1e-5, self.TS0_in-273.15, self.h_TS0_in*1e-3, self.s_TS0_in*1e-3))
        print("State TS0 outlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TS0*1e-5, self.TS0_out-273.15, self.h_TS0_out*1e-3, self.s_TS0_out*1e-3))
        print("State TES1 inlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES1*1e-5, self.T_TES1_in-273.15, self.h_TES1_in*1e-3, self.s_TES1_in*1e-3))
        print("State TES1 outlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES1*1e-5, self.T_TES1_out-273.15, self.h_TES1_out*1e-3, self.s_TES1_out*1e-3))
        print("State TES2 inlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES2*1e-5, self.T_TES2_in-273.15, self.h_TES2_in*1e-3, self.s_TES2_in*1e-3))
        print("State TES2 outlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES2*1e-5, self.T_TES2_out-273.15, self.h_TES2_out*1e-3, self.s_TES2_out*1e-3))
        print("State TES3 inlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES3*1e-5, self.T_TES3_in-273.15, self.h_TES3_in*1e-3, self.s_TES3_in*1e-3))
        print("State TES3 outlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES3*1e-5, self.T_TES3_out-273.15, self.h_TES3_out*1e-3, self.s_TES3_out*1e-3))
        print("State TES4 inlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES4*1e-5, self.T_TES4_in-273.15, self.h_TES4_in*1e-3, self.s_TES4_in*1e-3))
        print("State TES4 outlet: p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg.K".format(self.p_TES4*1e-5, self.T_TES4_out-273.15, self.h_TES4_out*1e-3, self.s_TES4_out*1e-3))
        print("Water storage temperature: {:.2f} °C".format(self.T_storage_water-273.15))
        print("Water heated temperature: {:.2f} °C".format(self.T_w_hot-273.15))
        print("TES hot outlet temperature: {:.2f} °C".format(self.T_storage_TES-273.15))
        print("--- Flow rates ---")
        print("Mass flow rate of CO2: {:.2f} kg/s".format(self.m_dot_CO2))
        print("Mass flow rate of TS0: {:.2f} kg/s".format(self.m_dot_TS0))
        print("Mass flow rate of TES1: {:.2f} kg/s".format(self.m_dot_TSE1))
        print("Mass flow rate of TES2: {:.2f} kg/s".format(self.m_dot_TSE2))
        print("Mass flow rate of TES3: {:.2f} kg/s".format(self.m_dot_TSE3))
        print("Mass flow rate of TES4: {:.2f} kg/s".format(self.m_dot_TSE4))
        print("--- Efficiencies ---")
        print("Rotex efficiency: {:.2f} %".format(self.eta_rotex*100))
        print("Heat exchanger efficiency TES0: {:.2f} %".format(self.eta_transex_TES0*100))
        print("Heat exchanger efficiency TES1: {:.2f} %".format(self.eta_transex_TES1*100))
        print("Heat exchanger efficiency TES2: {:.2f} %".format(self.eta_transex_TES2*100))
        print("Heat exchanger efficiency TES3: {:.2f} %".format(self.eta_transex_TES3*100))
        print("Heat exchanger efficiency TES4: {:.2f} %".format(self.eta_transex_TES4*100))
        print("Heat exchanger efficiency PCHX: {:.2f} %".format(self.eta_transex_PCHX*100))
        print("--- Losses ---")
        print("Rotex losses: {:.2f} kW".format(self.loss_rotex*1e-3))
        print("Evaporator losses: {:.2f} kW".format(self.loss_evaporator*1e-3))
        print("TES1 losses: {:.2f} kW".format(self.loss_TES1*1e-3))
        print("TES2 losses: {:.2f} kW".format(self.loss_TES2*1e-3))
        print("TES3 losses: {:.2f} kW".format(self.loss_TES3*1e-3))
        print("TES4 losses: {:.2f} kW".format(self.loss_TES4*1e-3))
        print("PCHX losses: {:.2f} kW".format(self.loss_PCHX*1e-3))
        print("Mechanical losses: {:.2f} kW".format(self.loss_mec*1e-3))
        print("Electrical losses: {:.2f} kW".format(self.loss_elec*1e-3))
        self.fig_pie_ex()
        self.plotTQTES0()


        #QUESTIONS : 
        # Tref, p_ref -> différent pour water, FN2, NaK ?
        # la pression (TES 1 et 2) ne semble pas avoir d'incidence sur le rendement -> puisque l'échange se fait à l'état liquide, une 
        # variation de pression àà la même température ne change pratiquement rien à l'enthalpie/entropie ?  
