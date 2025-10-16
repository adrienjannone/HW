"""
LELME2150 - Thermal cycles
Homework 3.0 - Pinch model heat exchanger

Signature of the function for the heat exchanger pinch model.

@author: Mattéo Hauglustaine
@date: October 10th, 2025
"""

#
# ===IMPORT PACKAGES============================================================

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

#
# ===HEAT EXCHANGER MODEL - TO BE IMPLEMENTED==========================================
#


class heat_exchanger(object):
    """
    Class for the simulation of heat exchangers

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
        -> Compute the heat exchange
    """

    def __init__(self, inputs, parameters, display):
        """
        Create a heat exchanger object.
        """
        p_hs, T_hs_su, T_hs_ex, T_cs_su, T_cs_ex = inputs
        self.p_hs = p_hs                                # [Pa] hot side fluid pressure
        self.T_hs_su = T_hs_su                          # [K]  hot side fluid supply temperature
        self.T_hs_ex = T_hs_ex                          # [K]  hot side fluid exit temperature
        self.T_cs_su = T_cs_su                          # [K]  cold side fluid supply temperature
        self.T_cs_ex = T_cs_ex                          # [K]  cold side fluid exit temperature
        self.pinch_target = parameters["pinch_target"]  # [K]  target for pinch point temperature difference
        self.fluid_hs = parameters["wf_hot_side"]       # hot side working fluid
        self.fluid_cs = parameters["wf_cold_side"]      # cold side working fluid
        self.p_cs_guess = parameters["p_cs_guess"]      # cold side working fluid pressure guess
        self.display = display                          # if True, make the plots

        self.h_hs = 0
        self.h_cs = 0
        self.m_dot_r = 0

        self.measured_pinch = 0
        self.eta_transex = 0


    def mass_ratio(self, p_evap):
        h_hs_su = CP.PropsSI("H", "T", self.T_hs_su, "P", self.p_hs, self.fluid_hs)  # [J/kg]
        h_hs_ex = CP.PropsSI("H", "T", self.T_hs_ex, "P", self.p_hs, self.fluid_hs)  # [J/kg]

        h_cs_su = CP.PropsSI("H", "T", self.T_cs_su, "P", p_evap, self.fluid_cs)  # [J/kg]
        h_cs_ex = CP.PropsSI("H", "T", self.T_cs_ex, "P", p_evap, self.fluid_cs)  # [J/kg]

        self.h_hs = (h_hs_su - h_hs_ex)  # [J/kg] 
        self.h_cs = (h_cs_ex - h_cs_su)  # [J/kg]
        self.m_dot_r = self.h_cs / self.h_hs



    def get_pinch(self, p_evap):
        self.mass_ratio(p_evap)
        T_c = CP.PropsSI("T", "P", p_evap, "Q", 0, self.fluid_cs)  # [K] saturation temperature
        h_c = CP.PropsSI("H", "P", p_evap, "Q", 0, self.fluid_cs)  # [J/kg] saturated liquid enthalpy
        h_cs_su = CP.PropsSI("H", "T", self.T_cs_su, "P", p_evap, self.fluid_cs)  # [J/kg]
        dh = h_c - h_cs_su  # [J/kg]
        h_hs_ex = CP.PropsSI("H", "T", self.T_hs_ex, "P", self.p_hs, self.fluid_hs)  # [J/kg]

        
        h_h = h_hs_ex + dh/self.m_dot_r  # [J/kg]
        T_h = CP.PropsSI("T", "H", h_h, "P", self.p_hs, self.fluid_hs)  # [K]
        pinch = T_h - T_c  # [K]
        return pinch
    
    def pinch_objective(self, p_evap):
        return self.get_pinch(p_evap) - self.pinch_target


    def evaluate(self):

        self.T_0 = 15 + 273.15 #[K] reference temperature
        self.p_0 = 1e5 #[Pa] reference pressure
        """
        This is the main method of the heat exchanger class.
        It finds the working fluid pressure to reach the pinch point and
        computes the exergy efficiency of the heat exchanger.
        """
        self.mass_ratio(self.p_cs_guess)
        self.p_evap_solution = fsolve(self.pinch_objective, self.p_cs_guess)[0]
        T1_m = (self.T_hs_su-self.T_hs_ex)/(np.log(self.T_hs_su/self.T_hs_ex))
        T2_m = (self.T_cs_ex-self.T_cs_su)/(np.log(self.T_cs_ex/self.T_cs_su)) 
        self.eta_transex = ((T2_m - self.T_0)/T2_m)*(T1_m/(T1_m - self.T_0))

        print(f"Optimal evaporation pressure: {self.p_evap_solution/1e5:.2f} bar")
        print(f"Pinch point value measured: {self.measured_pinch} K")
        print(f"Exergy efficiency of HEX: {self.eta_transex:.2%}")

        # Figures -------------------------------------------------------------
        if self.display:
            # --- Tq diagram ------------------------------------------------------
            pass
