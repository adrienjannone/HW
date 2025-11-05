"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signature of the function for the Brayton cycle

@author: Antoine Laterre
@date: August 26, 2022
"""

#
#===IMPORT PACKAGES============================================================
#

from CoolProp.CoolProp import PropsSI
import numpy as np
import scipy as sc
# Please do not add any other packages

#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

class Brayton_cycle(object):
    """
    Class for the simulation of Brayton cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (SI units = K, Pa).
    'parameters':
        -> Dictionary with the model parameters.
        
    Methods
    -------
    'evaluate':
        -> Compute the cycle and evaluate some KPI's.
    """
    
    def __init__(self,inputs,parameters):
        """
        Create a Brayton cycle object.
        """
        p_1,T_1,p_3,T_3         = inputs
        self.p_1,self.T_1       = p_1, T_1                  # pressure [Pa] and temperature [K] at the compressor inlet
        self.p_3,self.T_3       = p_3, T_3                  # pressure [Pa] and temperature [K] at the turbine    inlet
        
        self.eta_pi_c           = parameters['eta_pi_c']    # compressor polytropic efficiency     [-]
        self.eta_pi_t           = parameters['eta_pi_t']    # turbine    polytropic efficiency     [-]
        self.eta_mec_c          = parameters['eta_mec_c']   # compressor mechanical efficiency     [-]
        self.eta_mec_t          = parameters['eta_mec_t']   # turbine    mechanical efficiency     [-]
        self.k_cc               = parameters['k_cc']        # pressure loss coeff. in comb. chamb. [-]
        self.fluid              = parameters['fluid']       # fluid species (e.g. ['N2','O2'])
        self.fluid_prop         = parameters['fluid_prop']  # species proportions (e.g. [0.79, 0.21])

        self.h_1 = 0
        self.s_1 = 0
        self.h_2 = 0
        self.s_2 = 0
        self.p_2 = 0
        self.T_2 = 0
        self.h_3 = 0
        self.s_3 = 0
        self.T_4 = 0
        self.p_4 = 0
        self.h_4 = 0
        self.s_4 = 0

        self.R = 0


    def cp_func(self, T, p=1e+5):
        cp = 0
        for i in range(len(self.fluid)):
            cp += self.fluid_prop[i]*PropsSI('CPMASS','P',p,'T',T,self.fluid[i])
        return cp

    def cp_avg(self, T1, T2, p):
        return sc.integrate.quad(self.cp_func, T1, T2, args=(p,))[0]/(T2-T1)

    def get_R(self):
        for i in range(len(self.fluid)):
            self.R += self.fluid_prop[i]*PropsSI('GAS_CONSTANT',self.fluid[i])/PropsSI('MOLAR_MASS',self.fluid[i])
    
    def T2_func(self, T2):
        cp = self.cp_avg(self.T_1, T2, (self.p_2+self.p_1)/2)
        return T2 - self.T_1 * (self.p_2 / self.p_1) ** (self.R / (cp * self.eta_pi_c))
    
    def T4_func(self, T4):
        cp = self.cp_avg(self.T_3, T4, (self.p_4+self.p_3)/2)
        return T4 - self.T_3 * (self.p_4 / self.p_3) ** ((self.R*self.eta_pi_t) / cp)


    def evaluate(self):
        """
        This is the main method of the Brayton cycle class.
        It evaluates the different thermodynamic quantities at each state of 
        the Brayton cycle, as well as some KPI's.
        """
        self.get_R()

        # State 1
        for i in range(len(self.fluid)):
            self.h_1 += self.fluid_prop[i]*PropsSI('H','P',self.p_1,'T',self.T_1,self.fluid[i])
            self.s_1 += self.fluid_prop[i]*PropsSI('S','P',self.p_1,'T',self.T_1,self.fluid[i])

        # State 2
        self.p_2 = self.p_3 / self.k_cc
        self.T_2 = sc.optimize.fsolve(self.T2_func, (self.T_1+self.T_3)/2)[0]
        cp_2 = self.cp_avg(self.T_1, self.T_2, (self.p_1+self.p_2)/2)
        self.h_2 = self.h_1 + cp_2*(self.T_2 - self.T_1)
        self.s_2 = self.s_1 + cp_2*np.log(self.T_2/self.T_1) - self.R*np.log(self.p_2/self.p_1)


        # State 3
        cp_3 = self.cp_avg(self.T_2, self.T_3, (self.p_2+self.p_3)/2)
        self.h_3 = self.h_2 + cp_3*(self.T_3 - self.T_2)
        self.s_3 = self.s_2 + cp_3*np.log(self.T_3/self.T_2) - self.R*np.log(self.p_3/self.p_2)

        # State 4
        self.p_4 = self.p_1
        self.T_4 = sc.optimize.fsolve(self.T4_func, (self.T_2+self.T_3)/2)[0]
        cp_4 = self.cp_avg(self.T_3, self.T_4, (self.p_3+self.p_4)/2)
        self.h_4 = self.h_3 + cp_4*(self.T_4 - self.T_3)
        self.s_4 = self.s_3 + cp_4*np.log(self.T_4/self.T_3) - self.R*np.log(self.p_4/self.p_3)

        
        # Cycle efficiency - do not modify
        self.w_comp = (self.h_2-self.h_1)/self.eta_mec_c
        self.q_comb = (self.h_3-self.h_2)
        self.w_turb = (self.h_3-self.h_4)*self.eta_mec_t
        self.eta_en = (self.w_turb-self.w_comp)/self.q_comb        
        # Final outputs - do not modify
        self.p = (self.p_1, self.p_2, self.p_3, self.p_4)
        self.T = (self.T_1, self.T_2, self.T_3, self.T_4)
        self.s = (self.s_1, self.s_2, self.s_3, self.s_4)
        self.h = (self.h_1, self.h_2, self.h_3, self.h_4)
        self.states = (self.p,self.T,self.s,self.h)
        