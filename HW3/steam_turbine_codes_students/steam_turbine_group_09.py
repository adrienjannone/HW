"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Signature of the function for the steam turbine

@author: Antoine Laterre
@date: October 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#
import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

class steam_turbine(object):
    """
    Class for the simulation of steam turbines

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
        Create a steam turbine object.
        """        
        P_e                     = inputs
        self.P_e                = P_e                        # [W]  steam turbine power output
        self.p_1                = parameters['p_1']          # [Pa] pressure at state 1
        self.p_3                = parameters['p_3']          # [Pa] maximum steam pressure
        self.p_4                = parameters['p_4']          # [Pa] reheating pressure
        self.p_ref              = parameters['p_ref']        # [Pa] reference pressure for exergy computation
        self.T_ref              = parameters['T_ref']        # [K] reference temperature for exergy computation
        self.T_max              = parameters['T_max']        # [K] maximum steam temperature
        self.T_cd_out           = parameters['T_cd_out']     # [K] condenser cold outlet temperature
        self.T_cd_subcool       = parameters['T_cd_subcool'] # [K] subcooling at the condenser
        self.T_pinch_sc         = parameters['T_pinch_sc']   # [K] pinch temperature at the subcooler
        self.T_pinch_ex         = parameters['T_pinch_ex']   # [K] pinch temperature at heat exchangers
        self.T_pinch_cd         = parameters['T_pinch_cd']   # [K] pinch temperature at the condenser 
        self.T_drum             = parameters['T_drum']       # [K] minimum drum temperature
        self.x_6min             = parameters['x_6min']       # [-] minimum possible vapor quality after final expansion
        self.eta_is_HP          = parameters['eta_is_HP']    # [-] HP turbine isentropic efficiency (see text book pp. 55)
        self.eta_is_LP          = parameters['eta_is_LP']    # [-] LP turbine isentropic efficiency (see text book pp. 55)
        self.eta_mec            = parameters['eta_mec']      # [-] shafts bearings mechanical efficiency
        self.eta_pump           = parameters['eta_pump']     # [-] internal efficiency of the pump (see text book pp. 53)
        self.k_heat             = parameters['k_heat']       # [-] pressure loss coefficient at the generator 
        self.display            = display                    # if True, make the plots
        
    def evaluate(self):
        """
        This is the main method of the steam turbine class.
        It evaluates the different thermodynamic quantities at each state of 
        the cycle (pp. 89, fig 2.33), as well as some KPI's based on several 
        inputs and for a given electrical production. It can also plot the T-s 
        and h-s as well as the energy and exergy pies.
        """
 
        # >>>>>             <<<<< #    
        # Replace with your model #
        # >>>>>             <<<<< # 
        
        # ...
        # ...
        # ...

        # States --------------------------------------------------------------
        self.p           = self.p_1, self.p_2,self.p_3,self.p_4,self.p_5,self.p_6I,self.p_6II,self.p_6III,self.p_6IV,self.p_6V, self.p_6VI,self.p_6VII,self.p_6VIII,self.p_6,self.p_7,self.p_7I,self.p_7II,self.p_7III,self.p_7IV,self.p_7V,self.p_7VI,self.p_7VII,self.p_7VIII,self.p_8,self.p_90,self.p_9I,self.p_9II,self.p_9III,self.p_9IV,self.p_9V,self.p_9VI,self.p_9VII, self.p_9VIII # [Pa]     tuple containing the pressure at each state
        self.T           = self.T_1, self.T_2,self.T_3,self.T_4,self.T_5,self.T_6I,self.T_6II,self.T_6III,self.T_6IV,self.T_6V, self.T_6VI,self.T_6VII,self.T_6VIII,self.T_6,self.T_7,self.T_7I,self.T_7II,self.T_7III,self.T_7IV,self.T_7V,self.T_7VI,self.T_7VII,self.T_7VIII,self.T_8,self.T_90,self.T_9I,self.T_9II,self.T_9III,self.T_9IV,self.T_9V,self.T_9VI,self.T_9VII, self.T_9VIII # [K]      temperature at each state
        self.s           = self.s_1, self.s_2,self.s_3,self.s_4,self.s_5,self.s_6I,self.s_6II,self.s_6III,self.s_6IV,self.s_6V, self.s_6VI,self.s_6VII,self.s_6VIII,self.s_6,self.s_7,self.s_7I,self.s_7II,self.s_7III,self.s_7IV,self.s_7V,self.s_7VI,self.s_7VII,self.s_7VIII,self.s_8,self.s_90,self.s_9I,self.s_9II,self.s_9III,self.s_9IV,self.s_9V,self.s_9VI,self.s_9VII, self.s_9VIII # [J/kg/K] entropy at each state
        self.h           = self.h_1, self.h_2,self.h_3,self.h_4,self.h_5,self.h_6I,self.h_6II,self.h_6III,self.h_6IV,self.h_6V, self.h_6VI,self.h_6VII,self.h_6VIII,self.h_6,self.h_7,self.h_7I,self.h_7II,self.h_7III,self.h_7IV,self.h_7V,self.h_7VI,self.h_7VII,self.h_7VIII,self.h_8,self.h_90,self.h_9I,self.h_9II,self.h_9III,self.h_9IV,self.h_9V,self.h_9VI,self.h_9VII, self.h_9VIII # [J/kg]   enthalpy at each state
        self.x           = self.x_1, self.x_2,self.x_3,self.x_4,self.x_5,self.x_6I,self.x_6II,self.x_6III,self.x_6IV,self.x_6V, self.x_6VI,self.x_6VII,self.x_6VIII,self.x_6,self.x_7,self.x_7I,self.x_7II,self.x_7III,self.x_7IV,self.x_7V,self.x_7VI,self.x_7VII,self.x_7VIII,self.x_8,self.x_90,self.x_9I,self.x_9II,self.x_9III,self.x_9IV,self.x_9V,self.x_9VI,self.x_9VII, self.x_9VIII # [-]      vapor quality at each state
        self.e           = self.e_1, self.e_2,self.e_3,self.e_4,self.e_5,self.e_6I,self.e_6II,self.e_6III,self.e_6IV,self.e_6V, self.e_6VI,self.e_6VII,self.e_6VIII,self.e_6,self.e_7,self.e_7I,self.e_7II,self.e_7III,self.e_7IV,self.e_7V,self.e_7VI,self.e_7VII,self.e_7VIII,self.e_8,self.e_90,self.e_9I,self.e_9II,self.e_9III,self.e_9IV,self.e_9V,self.e_9VI,self.e_9VII, self.e_9VIII # [J/kg]   exergy at each state (use ref conditions)
        self.DAT         = self.p,self.T,self.s,self.h,self.x,self.e
        # Mass flow rates -----------------------------------------------------
        self.MASSFLOW    = self.dotm_v,self.dotm_tot
        #      o dotm_v         [kg/s]  mass flow rate of water at the condenser
        #      o dotm_tot       [kg/s]  total mass flow rate
        # Efficiencies --------------------------------------------------------
        self.ETA         = self.eta_cyclen,self.eta_cyclex,self.eta_condex,self.eta_rotex
        # -> see text book pp. 53-94
        #      o eta_cyclen    [-]      cycle energy efficiency
        #      o eta_cyclex    [-]      cycle exergy efficiency
        #      o eta_condex    [-]      condenser exergy efficiency
        #      o eta_rotex     [-]      pumps and turbines exergy efficiency
        # Energy losses -------------------------------------------------------
        self.DATEN       = self.loss_mec,self.loss_conden
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_conden   [W]      condenser energy losses
        # Exergy losses -------------------------------------------------------
        self.DATEX       = self.loss_mec,self.loss_rotex,self.loss_condex
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_rotex    [W]      pumps and turbines exergy losses
        #      o loss_condex   [W]      condenser exergy losses
        self.XMASSFLOW  = self.dotm_6I,self.dotm_6II,self.dotm_6III,self.dotm_6IV,self.dotm_6V,self.dotm_6VI,self.dotm_6VII,self.dotm_6VIII        
        # -> massflow rate in each feedheater (w.r.t. fig. 2.33, pp. 91)
        #      o dotm_6I          [kg/s]
        #      o dotm_6II         [kg/s]
        #      o dotm_6...        [kg/s]    
        #      o dotm_6VIII       [kg/s]
        # Energy and Exergy pie charts ----------------------------------------
        if self.display: self.FIG = self.fig_pie_en,self.fig_pie_ex, self.fig_Ts, self.fig_hs
        #      o fig_pie_en: pie chart of energy losses
        #      o fig_pie_ex: pie chart of exergy losses
        #      o fig_Ts_diagram: T-s diagram of the ST cycle
        #      o fig_hs_diagram: h-s diagram of the ST cycle 