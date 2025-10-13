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
        self.p_hs = p_hs  # [Pa] hot side fluid pressure
        self.T_hs_su = T_hs_su  # [K]  hot side fluid supply temperature
        self.T_hs_ex = T_hs_ex  # [K]  hot side fluid exit temperature
        self.T_cs_su = T_cs_su  # [K]  cold side fluid supply temperature
        self.T_cs_ex = T_cs_ex  # [K]  cold side fluid exit temperature
        self.pinch_target = parameters[
            "pinch_target"
        ]  # [K]  target for pinch point temperature difference
        self.fluid_hs = parameters["wf_hot_side"]  # hot side working fluid
        self.fluid_cs = parameters["wf_cold_side"]  # cold side working fluid
        self.p_cs_guess = parameters[
            "p_cs_guess"
        ]  # cold side working fluid pressure guess
        self.display = display  # if True, make the plots

    def evaluate(self):

        self.T_0 = 15 + 273.15 #[K] reference temperature
        self.p_0 = 1e5 #[Pa] reference pressure
        """
        This is the main method of the heat exchanger class.
        It finds the working fluid pressure to reach the pinch point and
        computes the exergy efficiency of the heat exchanger.
        """

        # >>>>>>>>               <<<<<<<<<<
        # >>>> Replace with your model <<<<
        # >>>>>>>>               <<<<<<<<<<

        # ...
        # ...
        # ...

        print(f"Optimal evaporation pressure: {self.p_evap_solution/1e5:.2f} bar")
        print(f"Pinch point value measured: {self.measured_pinch} K")
        print(f"Exergy efficiency of HEX: {self.eta_transex:.2%}")

        # Figures -------------------------------------------------------------
        if self.display:
            # --- Tq diagram ------------------------------------------------------
            pass
