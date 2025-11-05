"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signature of the function for the Rankine cycle

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
#===RANKINE CYCLE - TO BE IMPLEMENTED==========================================
#

class Rankine_cycle(object):
    """
    Class for the simulation of the Rankine cycle

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
        Create a Rankine cycle object.
        """
        T_1,p_3,T_3,P_el,LHV = inputs
        self.T_1                = T_1                       # temperature [K] at the pump inlet
        self.p_3,self.T_3       = p_3, T_3                  # pressure [Pa] and temperature [K] at the turbine inlet
        self.P_el, self.LHV     = P_el, LHV                 # demanded power [W] and fuel LHV [J/kg]
        
        self.eta_pump           = parameters['eta_pump']    # pump internal efficiency          [-]
        self.eta_is_t           = parameters['eta_is_t']    # turbine isentropic efficiency     [-]
        self.eta_mec_t          = parameters['eta_mec_t']   # turbine mechanical efficiency     [-]
        self.eta_gen            = parameters['eta_gen']     # steam generator efficiency        [-]
        self.k_gen              = parameters['k_gen']       # pressure losses coeff. in s. gen. [-]
        self.fluid              = parameters['fluid']       # fluid (e.g. 'H2O')       
        
    def evaluate(self):
        """
        This is the main method of the Rankine cycle class.
        It evaluates the different thermodynamic quantities at each state of 
        the Rankine cycle, as well as some KPI's.
        """    
        self.h_3 = PropsSI('H','P',self.p_3,'T',self.T_3,self.fluid) # [J/kg]
        self.s_3 = PropsSI('S','P',self.p_3,'T',self.T_3,self.fluid) # [J/kg/K]
        self.x_3 = PropsSI('Q','P',self.p_3,'T',self.T_3,self.fluid) # [-]

        self.T_4 = self.T_1 + 3
        self.p_4 = PropsSI('P','T',self.T_4,'Q',0,self.fluid) # [Pa]
        H4S = PropsSI('H','P',self.p_4,'S',self.s_3,self.fluid) # [J/kg]
        self.h_4 = self.h_3 + (H4S - self.h_3)*self.eta_is_t
        self.s_4 = PropsSI('S','P',self.p_4,'H',self.h_4,self.fluid) # [J/kg/K]
        self.x_4 = PropsSI('Q','P',self.p_4,'H',self.h_4,self.fluid) # [-]


        self.p_1 = self.p_4 # [Pa]
        self.h_1 = PropsSI('H','P',self.p_1,'T',self.T_1,self.fluid) # [J/kg]
        self.s_1 = PropsSI('S','P',self.p_1,'T',self.T_1,self.fluid) # [J/kg/K]
        self.x_1 = PropsSI('Q','P',self.p_1,'T',self.T_1,self.fluid) # [-]

        v = 1/PropsSI('D','P',self.p_1,'T',self.T_1,self.fluid) # [m3/kg]
        self.p_2 = self.p_3/self.k_gen # [Pa]
        self.h_2 = self.h_1 + v*(self.p_2 - self.p_1)/self.eta_pump
        self.s_2 = PropsSI('S','P',self.p_2,'H',self.h_2,self.fluid) # [J/kg/K]
        self.T_2 = PropsSI('T','P',self.p_2,'H',self.h_2,self.fluid) # [K]
        self.x_2 = PropsSI('Q','P',self.p_2,'H',self.h_2,self.fluid) # [-]


       
        
        # Cycle efficiency - do not modify
        self.W_pump = (self.h_2-self.h_1)
        self.Q_comb = (self.h_3-self.h_2)/self.eta_gen
        self.W_turb = (self.h_3-self.h_4)*self.eta_mec_t
        self.eta_en = (self.W_turb-self.W_pump)/self.Q_comb
        self.dot_m_f= self.P_el/self.LHV/self.eta_en
        # Final outputs - do not modify
        self.p = (self.p_1, self.p_2, self.p_3, self.p_4)
        self.T = (self.T_1, self.T_2, self.T_3, self.T_4)
        self.s = (self.s_1, self.s_2, self.s_3, self.s_4)
        self.h = (self.h_1, self.h_2, self.h_3, self.h_4)
        self.x = (self.x_1, self.x_2, self.x_3, self.x_4)
        self.states = (self.p,self.T,self.s,self.h,self.x)