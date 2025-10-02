"""
LELME2150 - Thermal cycles
Homework 2 - Gas turbine

Signature of the function for the gas turbine

@author: Antoine Laterre
@date: September 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#
import CoolProp.CoolProp as CP
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

class gas_turbine(object):
    """
    Class for the simulation of gas turbines

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (SI units = K, Pa).
    'parameters':
        -> Dictionary with the model parameters.
    'display':
        -> Boolean to plot the graphs.
        
    Methods
    -------
    'evaluate':
        -> Compute the cycle
    """
    
    def __init__(self,inputs,parameters,display):
        """
        Create a gas turbine object.
        """
        p_1,T_1,P_e             = inputs
        self.p_1                = p_1                       # [Pa] compressor supply pressure
        self.T_1                = T_1                       # [K]  compressor supply temperature
        self.P_e                = P_e                       # [W]  gas turbine power output
        self.T_3                = parameters['T_3']         # [K]  turbine inlet temperature 
        self.r_c                = parameters['r_c']         # [-]  compression ratio
        self.eta_pi_c           = parameters['eta_pi_c']    # [-]  polytropic efficiency of compressor
        self.eta_pi_t           = parameters['eta_pi_t']    # [-]  polytropic efficiency of turbine
        self.k_cc               = parameters['k_cc']        # [-]  pressure   losses coefficient in combustion chamber
        self.k_mec              = parameters['k_mec']       # [-]  mechanical losses coefficient
        self.air                = parameters['air']         # air  species
        self.air_prop           = parameters['air_prop']    # molar fractions of air species
        self.alkane             = parameters['alkane']      # alkane composition ['c','h']
        self.display            = display                   # if True, make the plots

        self.p_2 = 0.0
        self.p_3 = 0.0
        self.p_4 = 0.0
        self.T_2 = 0.0
        self.T_4 = 0.0
        self.h_1 = 0.0
        self.h_2 = 0.0
        self.h_3 = 0.0
        self.h_4 = 0.0
        self.s_1 = 0.0
        self.s_2 = 0.0
        self.s_3 = 0.0
        self.s_4 = 0.0
        self.e_1 = 0.0
        self.e_2 = 0.0
        self.e_3 = 0.0
        self.e_4 = 0.0
        self.LHV = 0.0
        self.e_f = 0.0
        self.excess_air = 0.0
        self.gas = []
        self.gas_prop = []
        self.dotm_a = 0.0
        self.dotm_f = 0.0
        self.dotm_g = 0.0
        self.eta_cyclen = 0.0
        self.eta_toten = 0.0
        self.eta_mec = 0.0
        self.eta_cyclex = 0.0
        self.eta_totex = 0.0
        self.eta_rotex = 0.0
        self.eta_combex = 0.0
        self.loss_mec = 0.0
        self.loss_echen = 0.0
        self.loss_rotex = 0.0
        self.loss_combex = 0.0
        self.loss_echex = 0.0

    def cp_func(self, T, p=1e+5):
        cp = 0
        for i in range(len(self.air)):
            cp += self.air_prop[i]*CP.PropsSI('CPMASS','P',p,'T',T,self.air[i])
        return cp
    
    def cp_avg(self, T1, T2, p):
        return sc.integrate.quad(self.cp_func, T1, T2, args=(p,))[0]/(T2-T1)
    
    def T2_func(self, T2):
        cp = self.cp_avg(self.T_1, T2, (self.p_2+self.p_1)/2)
        return T2 - self.T_1 * (self.p_2 / self.p_1) ** (self.R / (cp * self.eta_pi_c))
    
    def get_R(self):
        R_u = CP.PropsSI("GAS_CONSTANT", "air")
        M_mix = 0
        for i in range(len(self.air)):
            M_mix += self.air_prop[i] * CP.PropsSI('MOLAR_MASS', self.air[i])
        self.R = R_u / M_mix

    def T4_func(self, T4):
        cp = self.cp_avg(self.T_3, T4, (self.p_4+self.p_3)/2)
        return T4 - self.T_3 * (self.p_4 / self.p_3) ** ((self.R*self.eta_pi_t) / cp)

    def evaluate(self):
        """
        This is the main method of the gas turbine class.
        It evaluates the different thermodynamic quantities at each state of 
        the cycle, as well as some KPI's.
        """
    
        # >>>>>             <<<<< #    
        # Replace with your model #
        # >>>>>             <<<<< # 
        self.get_R()
        T0 = 273.15
        p0 = 1e+5
        h_ref = [CP.PropsSI('H','P',p0,'T',T0,'N2'),CP.PropsSI('H','P',p0,'T',T0,'O2')]
        s_ref = [CP.PropsSI('S','P',p0,'T',T0,'N2'),CP.PropsSI('S','P',p0,'T',T0,'O2')]
        Cx_fuel = self.alkane[0]
        Hy_fuel = self.alkane[1]
        m_a1 = ((Cx_fuel + (Hy_fuel/4))*(32+28*3.76))/(12*Cx_fuel + Hy_fuel) # [kg_air/kg_fuel]
        print('Stoichiometric air-to-fuel ratio = ', m_a1, 'kg_air/kg_fuel')
        H_form_CH4 = -74.8e+3 # [J/mol]
        H_form_O2 = 0.0
        H_form_CO2 = -393.5e+3 # [J/mol]
        H_form_H2O = -241.8e+3 # [J/mol]
        LHV_molar = H_form_CH4 + 2*H_form_O2 - H_form_CO2 - 2*H_form_H2O # [J/mol]
        self.LHV = LHV_molar/(1000* CP.PropsSI('M', 'T', T0, 'P', p0, 'CH4')) # [kJ/kg]
        print('LHV = ', self.LHV, 'kJ/kg')
        cp_CH4 = 35.3 #kJ/ kmol.K 
        hf = cp_CH4 * self.T_3
        print('hf = ', hf, 'kJ/kmol')


        

        for i in range(len(self.air)):
            self.h_1 += self.air_prop[i]*(CP.PropsSI('H','P',self.p_1,'T',self.T_1,self.air[i])-h_ref[i])
            self.s_1 += self.air_prop[i]*(CP.PropsSI('S','P',self.p_1,'T',self.T_1,self.air[i])-s_ref[i])
        
        self.p_2 = self.r_c*self.p_1
        self.T_2 = sc.optimize.fsolve(self.T2_func, (self.T_1+self.T_3)/2)[0]
        self.h_2 = self.h_1 + self.cp_avg(self.T_1, self.T_2, (self.p_2+self.p_1)/2)*(self.T_2 - self.T_1)
        self.s_2 = self.s_1 + (1-self.eta_pi_c)* self.cp_avg(self.T_1, self.T_2, (self.p_2+self.p_1)/2)*np.log(self.T_2/self.T_1)
        
        self.p_3 = self.p_2*self.k_cc
        #h_3 = function de lambda
        #s3 = 
        
        self.p_4 = self.p_1
        self.T_4 = sc.optimize.fsolve(self.T4_func, (self.T_1+self.T_3)/2)[0]
        self.h_4 = self.h_3 - self.cp_avg(self.T_3, self.T_4, (self.p_3+self.p_4)/2)*(self.T_3 - self.T_4)
        self.s_4 = self.s_3 + ((self.eta_pi_t-1)/self.eta_pi_t)* self.cp_avg(self.T_3, self.T_4, (self.p_3+self.p_4)/2)*np.log(self.T_4/self.T_3)
        

        # States --------------------------------------------------------------
        self.p           = self.p_1, self.p_2, self.p_3, self.p_4
        self.T           = self.T_1, self.T_2, self.T_3, self.T_4
        self.s           = self.s_1, self.s_2, self.s_3, self.s_4
        self.h           = self.h_1, self.h_2, self.h_3, self.h_4
        self.e           = self.e_1, self.e_2, self.e_3, self.e_4
        self.DAT         = self.p,self.T,self.s,self.h,self.e
        # Combustion paramters ------------------------------------------------
        self.COMBUSTION  = self.LHV,self.e_f,self.excess_air,self.gas,self.gas_prop
        # Mass flow rates -----------------------------------------------------
        self.MASSFLOW    = self.dotm_a,self.dotm_f,self.dotm_g
        # Efficiencies --------------------------------------------------------
        self.ETA         = self.eta_cyclen,self.eta_toten,self.eta_mec,self.eta_cyclex,self.eta_totex,self.eta_rotex,self.eta_combex
        # Energy losses -------------------------------------------------------
        self.DATEN       = self.loss_mec,self.loss_echen
        # Exergy losses -------------------------------------------------------
        self.DATEX       = self.loss_mec,self.loss_rotex,self.loss_combex,self.loss_echex
        # Energy and Exergy pie charts
        #if self.display: self.FIG = self.fig_pie_en,self.fig_pie_ex, self.fig_Ts, self.fig_ph