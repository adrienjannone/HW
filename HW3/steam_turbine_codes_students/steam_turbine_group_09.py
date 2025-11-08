"""
LELME2150 - Thermal cycles
Homework 3 - Steam turbine

Signature of the function for the steam turbine

@author: Antoine Laterre
@date: October 30, 2022
"""

"""
Current working time : 2:00
05/11/25 : 14:00 - 16:00


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

        self.R = 8.314/CP.PropsSI('M','Water')               # [J/kg/K] specific gas constant for water vapor
        

    def pressure_bleedings(self,h_bleed, h_prev_bleed, s_prev_bleed):
        h_bleed_s = h_prev_bleed - (h_prev_bleed - h_bleed)/self.eta_is_LP
        p_bleed = CP.PropsSI('P','H',h_bleed_s,'S',s_prev_bleed,'Water')
        T_bleed = CP.PropsSI('T','H',h_bleed,'P',p_bleed,'Water')
        x_bleed = CP.PropsSI('Q','H',h_bleed,'P',p_bleed,'Water')
        s_bleed = CP.PropsSI('S','H',h_bleed,'P',p_bleed,'Water')
        e_bleed = (h_bleed - self.h_ref) - self.T_ref*(s_bleed - self.s_ref)
        return T_bleed, p_bleed, x_bleed, s_bleed, e_bleed  
    
    def saturated_liquid_7(self, pressure):
        x_7 = 0
        T_7 = CP.PropsSI('T','P',pressure,'Q',x_7,'Water')
        h_7 = CP.PropsSI('H','P',pressure,'Q',x_7,'Water')
        s_7 = CP.PropsSI('S','P',pressure,'Q',x_7,'Water')
        e_7 = (h_7 - self.h_ref) - self.T_ref*(s_7 - self.s_ref)
        return T_7, x_7, h_7, s_7, e_7
    
    def heat_exchangers_9(self, pressure, T_out): 
        T_9 = T_out - self.T_pinch_ex
        x_9 = CP.PropsSI('Q','P',pressure,'T',T_9,'Water')
        h_9 = CP.PropsSI('H','P',pressure,'T',T_9,'Water')
        s_9 = CP.PropsSI('S','P',pressure,'T',T_9,'Water')
        e_9 = (h_9 - self.h_ref) - self.T_ref*(s_9 - self.s_ref)
        return T_9, x_9, h_9, s_9, e_9


    def evaluate(self):
        """
        This is the main method of the steam turbine class.
        It evaluates the different thermodynamic quantities at each state of 
        the cycle (pp. 89, fig 2.33), as well as some KPI's based on several 
        inputs and for a given electrical production. It can also plot the T-s 
        and h-s as well as the energy and exergy pies.
        """

        #h_ref and s_ref 
        self.h_ref = CP.PropsSI('H','P',self.p_ref,'T',self.T_ref,'Water')    # [J/kg] enthalpy at reference state
        self.s_ref = CP.PropsSI('S','P',self.p_ref,'T',self.T_ref,'Water')    # [J/kg/K] entropy at reference state

        # Generator 
        self.p_2 = self.p_3/self.k_heat # [Pa] pressure after generator

        self.T_3 = self.T_max                                           # [K] temperature at state 3
        self.h_3 = CP.PropsSI('H','P',self.p_3,'T',self.T_3,'Water')    # [J/kg] enthalpy at state 3
        self.s_3 = CP.PropsSI('S','P',self.p_3,'T',self.T_3,'Water')    # [J/kg/K] entropy at state 3
        self.x_3 = CP.PropsSI('Q','P',self.p_3,'T',self.T_3,'Water')    # [-] vapor quality at state 3
        self.e_3 = (self.h_3 - self.h_ref) - self.T_ref*(self.s_3 - self.s_ref) # [J/kg] exergy at state 3
        print("State 3 : %f %f %f %f %f %f" % (self.T_3,self.p_3,self.x_3,self.h_3,self.s_3,self.e_3))

        self.T_5 = self.T_max                                           # [K] temperature at state 5
        self.p_5 = self.p_4*self.k_heat                                 # [Pa] pressure at state 5
        self.h_5 = CP.PropsSI('H','P',self.p_5,'T',self.T_5,'Water')    # [J/kg] enthalpy at state 5
        self.s_5 = CP.PropsSI('S','P',self.p_5,'T',self.T_5,'Water')    # [J/kg/K] entropy at state 5
        self.x_5 = CP.PropsSI('Q','P',self.p_5,'T',self.T_5,'Water')    # [-] vapor quality at state 5 
        self.e_5 = (self.h_5 - self.h_ref) - self.T_ref*(self.s_5 - self.s_ref) # [J/kg] exergy at state 5  
        print("State 5 : %f %f %f %f %f %f" % (self.T_5,self.p_5,self.x_5,self.h_5,self.s_5,self.e_5))

        h4s = CP.PropsSI('H','P',self.p_4,'S',self.s_3,'Water')         # [J/kg] isentropic enthalpy at state 4
        self.h_4 = self.h_3 - self.eta_is_HP*(self.h_3 - h4s)           # [J/kg] enthalpy at state 4
        self.T_4 = CP.PropsSI('T','P',self.p_4,'H',self.h_4,'Water')    # [K] temperature at state 4
        self.s_4 = CP.PropsSI('S','P',self.p_4,'H',self.h_4,'Water')    # [J/kg/K] entropy at state 4
        self.x_4 = CP.PropsSI('Q','P',self.p_4,'H',self.h_4,'Water')    # [-] vapor quality at state 4
        self.e_4 = (self.h_4 - self.h_ref) - self.T_ref*(self.s_4 - self.s_ref) # [J/kg] exergy at state 4
        print("State 4 : %f %f %f %f %f %f" % (self.T_4,self.p_4,self.x_4,self.h_4,self.s_4,self.e_4))  

        self.T_7  = self.T_cd_out + self.T_pinch_cd
        self.x_7  = 0
        self.p_7  = CP.PropsSI('P','T',self.T_7,'Q',self.x_7,'Water')    # [Pa] pressure at state 7
        self.h_7  = CP.PropsSI('H','P',self.p_7,'Q',self.x_7,'Water')    # [J/kg] enthalpy at state 7
        self.s_7  = CP.PropsSI('S','P',self.p_7,'Q',self.x_7,'Water')    # [J/kg/K] entropy at state 7
        self.e_7  = (self.h_7 - self.h_ref) - self.T_ref*(self.s_7 - self.s_ref) # [J/kg] exergy at state 7
        print("State 7 : %f %f %f %f %f %f" % (self.T_7,self.p_7,self.x_7,self.h_7,self.s_7,self.e_7))

        self.T_6 = self.T_7 + self.T_cd_subcool                  # [K] temperature at state 6
        h6s = CP.PropsSI('H','T',self.T_6,'S',self.s_5,'Water')         # [J/kg] isentropic enthalpy at state 6
        self.h_6 = self.h_5 - self.eta_is_LP*(self.h_5 - h6s)           # [J/kg] enthalpy at state 6
        self.p_6 = CP.PropsSI('P','T',self.T_6,'S',self.s_5,'Water')            # [Pa] pressure at state 6
        self.s_6 = CP.PropsSI('S','H',self.h_6,'P',self.p_6,'Water')            # [J/kg/K] entropy at state 6
        self.x_6 = CP.PropsSI('Q','H',self.h_6,'P',self.p_6,'Water')            # [-] vapor quality at state 6
        self.e_6 = (self.h_6 - self.h_ref) - self.T_ref*(self.s_6 - self.s_ref) # [J/kg] exergy at state 6
        print("State 6 : %f %f %f %f %f %f" % (self.T_6,self.p_6,self.x_6,self.h_6,self.s_6,self.e_6))



        #6_VIII = same than 4
        self.T_6VIII = self.T_4
        self.p_6VIII = self.p_4
        self.h_6VIII = self.h_4
        self.s_6VIII = self.s_4
        self.x_6VIII = self.x_4
        self.e_6VIII = self.e_4


        # Bleedings enthalpies - bleedings occur mid-expansion (p72)
        dh_IP_LP = self.h_5 - self.h_6
        dh_bleed = dh_IP_LP / 8 #8 bleedings
        self.h_6VII = self.h_5 - dh_bleed
        self.h_6VI = self.h_6VII - dh_bleed
        self.h_6V = self.h_6VI - dh_bleed
        self.h_6IV = self.h_6V - dh_bleed   
        self.h_6III = self.h_6IV - dh_bleed
        self.h_6II = self.h_6III - dh_bleed
        self.h_6I = self.h_6II - dh_bleed

        
        # Bleedings states
        self.T_6VII, self.p_6VII, self.x_6VII, self.s_6VII, self.e_6VII = self.pressure_bleedings(self.h_6VII, self.h_5, self.s_5)
        self.T_6VI, self.p_6VI, self.x_6VI, self.s_6VI, self.e_6VI = self.pressure_bleedings(self.h_6VI, self.h_6VII, self.s_6VII)
        self.T_6V, self.p_6V, self.x_6V, self.s_6V, self.e_6V = self.pressure_bleedings(self.h_6V, self.h_6VI, self.s_6VI)
        self.T_6IV, self.p_6IV, self.x_6IV, self.s_6IV, self.e_6IV = self.pressure_bleedings(self.h_6IV, self.h_6V, self.s_6V)
        self.T_6III, self.p_6III, self.x_6III, self.s_6III, self.e_6III = self.pressure_bleedings(self.h_6III, self.h_6IV, self.s_6IV)
        self.T_6II, self.p_6II, self.x_6II, self.s_6II, self.e_6II = self.pressure_bleedings(self.h_6II, self.h_6III, self.s_6III)
        self.T_6I, self.p_6I, self.x_6I, self.s_6I, self.e_6I = self.pressure_bleedings(self.h_6I, self.h_6II, self.s_6II)
        print("State 6_I : %f %f %f %f %f %f" % (self.T_6I,self.p_6I,self.x_6I,self.h_6I,self.s_6I,self.e_6I))
        print("State 6_II : %f %f %f %f %f %f" % (self.T_6II,self.p_6II,self.x_6II,self.h_6II,self.s_6II,self.e_6II))
        print("State 6_III : %f %f %f %f %f %f" % (self.T_6III,self.p_6III,self.x_6III,self.h_6III,self.s_6III,self.e_6III))
        print("State 6_IV : %f %f %f %f %f %f" % (self.T_6IV,self.p_6IV,self.x_6IV,self.h_6IV,self.s_6IV,self.e_6IV))
        print("State 6_V : %f %f %f %f %f %f" % (self.T_6V,self.p_6V,self.x_6V,self.h_6V,self.s_6V,self.e_6V))
        print("State 6_VI : %f %f %f %f %f %f" % (self.T_6VI,self.p_6VI,self.x_6VI,self.h_6VI,self.s_6VI,self.e_6VI))
        print("State 6_VII : %f %f %f %f %f %f" % (self.T_6VII,self.p_6VII,self.x_6VII,self.h_6VII,self.s_6VII,self.e_6VII))
        print("State 6_VIII : %f %f %f %f %f %f" % (self.T_6VIII,self.p_6VIII,self.x_6VIII,self.h_6VIII,self.s_6VIII,self.e_6VIII))


        # Drum
        pDrum = CP.PropsSI('P','T',self.T_drum,'Q',0,'Water')            # [Pa] pressure at the drum
        

        # Right Side exhangers pressures
        self.p_8 = pDrum
        self.p_90 = pDrum
        self.p_9I = pDrum
        self.p_9II = pDrum
        self.p_9III = pDrum
        self.p_7IV = pDrum
        
        # Left Side exchangers pressures
        self.p_9IV = self.p_1
        self.p_9V = self.p_1
        self.p_9VI = self.p_1
        self.p_9VII = self.p_1
        self.p_9VIII = self.p_1

        # Saturated liquid 7 pressures (heat exchangers are isobaric)
        self.p_7I = self.p_6I
        self.p_7II = self.p_6II
        self.p_7III = self.p_6III
        self.p_7V = self.p_6V
        self.p_7VI = self.p_6VI
        self.p_7VII = self.p_6VII
        self.p_7VIII = self.p_6VIII

        self.T_7IV = self.T_drum
        self.x_7IV = 0 # saturated liquid at the drum
        self.h_7IV = CP.PropsSI('H','P',self.p_7IV,'Q',self.x_7IV,'Water')    # [J/kg] enthalpy at state 7IV
        self.s_7IV = CP.PropsSI('S','P',self.p_7IV,'Q',self.x_7IV,'Water')    # [J/kg/K] entropy at state 7IV   
        self.e_7IV = (self.h_7IV - self.h_ref) - self.T_ref*(self.s_7IV - self.s_ref) # [J/kg] exergy at state 7IV  
        print("State 7_IV : %f %f %f %f %f %f" % (self.T_7IV,self.p_7IV,self.x_7IV,self.h_7IV,self.s_7IV,self.e_7IV))

        # Pump Pe
        self.h_8 = self.h_7 + CP.PropsSI('V','P',self.p_7,'T',self.T_7,'Water')*(self.p_8 - self.p_7)/self.eta_pump # [J/kg] enthalpy at state 8
        self.T_8 = CP.PropsSI('T','P',self.p_8,'H',self.h_8,'Water')    # [K] temperature at state 8
        self.s_8 = CP.PropsSI('S','P',self.p_8,'H',self.h_8,'Water')    # [J/kg/K] entropy at state 8
        self.x_8 = CP.PropsSI('Q','P',self.p_8,'H',self.h_8,'Water')    # [-] vapor quality at state 8
        self.e_8 = (self.h_8 - self.h_ref) - self.T_ref*(self.s_8 - self.s_ref) # [J/kg] exergy at state 8
        print("State 8 : %f %f %f %f %f %f" % (self.T_8,self.p_8,self.x_8,self.h_8,self.s_8,self.e_8))

        self.h_9IV = self.h_7IV + CP.PropsSI('V','P',self.p_9IV,'T',self.T_7IV,'Water')*(self.p_9IV - self.p_7IV)/self.eta_pump # [J/kg] enthalpy at state 9IV
        self.T_9IV = CP.PropsSI('T','P',self.p_9IV,'H',self.h_9IV,'Water')    # [K] temperature at state 9IV
        self.s_9IV = CP.PropsSI('S','P',self.p_9IV,'H',self.h_9IV,'Water')    # [J/kg/K] entropy at state 9IV
        self.x_9IV = CP.PropsSI('Q','P',self.p_9IV,'H',self.h_9IV,'Water')    # [-] vapor quality at state 9IV
        self.e_9IV = (self.h_9IV - self.h_ref) - self.T_ref*(self.s_9IV - self.s_ref) # [J/kg] exergy at state 9IV  
        print("State 9_IV : %f %f %f %f %f %f" % (self.T_9IV,self.p_9IV,self.x_9IV,self.h_9IV,self.s_9IV,self.e_9IV))

        # states 7_I to 7_VIII 
        self.T_7I, self.x_7I, self.h_7I, self.s_7I, self.e_7I = self.saturated_liquid_7(self.p_7I)
        self.T_7II, self.x_7II, self.h_7II, self.s_7II, self.e_7II = self.saturated_liquid_7(self.p_7II)
        self.T_7III, self.x_7III, self.h_7III, self.s_7III, self.e_7III = self.saturated_liquid_7(self.p_7III)  
        self.T_7V, self.x_7V, self.h_7V, self.s_7V, self.e_7V = self.saturated_liquid_7(self.p_7V)
        self.T_7VI, self.x_7VI, self.h_7VI, self.s_7VI, self.e_7VI = self.saturated_liquid_7(self.p_7VI)    
        self.T_7VII, self.x_7VII, self.h_7VII, self.s_7VII, self.e_7VII = self.saturated_liquid_7(self.p_7VII)
        self.T_7VIII, self.x_7VIII, self.h_7VIII, self.s_7VIII, self.e_7VIII = self.saturated_liquid_7(self.p_7VIII)
        print("State 7_I : %f %f %f %f %f %f" % (self.T_7I,self.p_7I,self.x_7I,self.h_7I,self.s_7I,self.e_7I))  
        print("State 7_II : %f %f %f %f %f %f" % (self.T_7II,self.p_7II,self.x_7II,self.h_7II,self.s_7II,self.e_7II))
        print("State 7_III : %f %f %f %f %f %f" % (self.T_7III,self.p_7III,self.x_7III,self.h_7III,self.s_7III,self.e_7III))
        print("State 7_V : %f %f %f %f %f %f" % (self.T_7V,self.p_7V,self.x_7V,self.h_7V,self.s_7V,self.e_7V))
        print("State 7_VI : %f %f %f %f %f %f" % (self.T_7VI,self.p_7VI,self.x_7VI,self.h_7VI,self.s_7VI,self.e_7VI))
        print("State 7_VII : %f %f %f %f %f %f" % (self.T_7VII,self.p_7VII,self.x_7VII,self.h_7VII,self.s_7VII,self.e_7VII))    
        print("State 7_VIII : %f %f %f %f %f %f" % (self.T_7VIII,self.p_7VIII,self.x_7VIII,self.h_7VIII,self.s_7VIII,self.e_7VIII))

        # state 9_0 to 9_VIII
        self.T_9I, self.x_9I, self.h_9I, self.s_9I, self.e_9I = self.heat_exchangers_9(self.p_9I, self.T_7I)
        self.T_9II, self.x_9II, self.h_9II, self.s_9II, self.e_9II = self.heat_exchangers_9(self.p_9II, self.T_7II)
        self.T_9III, self.x_9III, self.h_9III, self.s_9III, self.e_9III = self.heat_exchangers_9(self.p_9III, self.T_7III)
        self.T_9V, self.x_9V, self.h_9V, self.s_9V, self.e_9V = self.heat_exchangers_9(self.p_9V, self.T_7V)
        self.T_9VI, self.x_9VI, self.h_9VI, self.s_9VI, self.e_9VI = self.heat_exchangers_9(self.p_9VI, self.T_7VI)
        self.T_9VII, self.x_9VII, self.h_9VII, self.s_9VII, self.e_9VII = self.heat_exchangers_9(self.p_9VII, self.T_7VII)
        self.T_9VIII, self.x_9VIII, self.h_9VIII, self.s_9VIII, self.e_9VIII = self.heat_exchangers_9(self.p_9VIII, self.T_7VIII)
        self.T_1, self.x_1, self.h_1, self.s_1, self.e_1 = self.heat_exchangers_9(self.p_1, self.T_7VIII)
        print("State 9_I : %f %f %f %f %f %f" % (self.T_9I,self.p_9I,self.x_9I,self.h_9I,self.s_9I,self.e_9I))
        print("State 9_II : %f %f %f %f %f %f" % (self.T_9II,self.p_9II,self.x_9II,self.h_9II,self.s_9II,self.e_9II))
        print("State 9_III : %f %f %f %f %f %f" % (self.T_9III,self.p_9III,self.x_9III,self.h_9III,self.s_9III,self.e_9III))
        print("State 9_V : %f %f %f %f %f %f" % (self.T_9V,self.p_9V,self.x_9V,self.h_9V,self.s_9V,self.e_9V))
        print("State 9_VI : %f %f %f %f %f %f" % (self.T_9VI,self.p_9VI,self.x_9VI,self.h_9VI,self.s_9VI,self.e_9VI))   
        print("State 9_VII : %f %f %f %f %f %f" % (self.T_9VII,self.p_9VII,self.x_9VII,self.h_9VII,self.s_9VII,self.e_9VII))
        print("State 9_VIII : %f %f %f %f %f %f" % (self.T_9VIII,self.p_9VIII,self.x_9VIII,self.h_9VIII,self.s_9VIII,self.e_9VIII))
        print("State 1 : %f %f %f %f %f %f" % (self.T_1,self.p_1,self.x_1,self.h_1,self.s_1,self.e_1))
        # 

        # state 9_0
        self.T_90 = self.T_7I - self.T_pinch_sc
        self.x_90 = CP.PropsSI('Q','P',self.p_90,'T',self.T_90,'Water')
        self.h_90 = CP.PropsSI('H','P',self.p_90,'T',self.T_90,'Water')
        self.s_90 = CP.PropsSI('S','P',self.p_90,'T',self.T_90,'Water')
        self.e_90 = (self.h_90 - self.h_ref) - self.T_ref*(self.s_90 - self.s_ref)
        print("State 9_0 : %f %f %f %f %f %f" % (self.T_90,self.p_90,self.x_90,self.h_90,self.s_90,self.e_90))

        # state 2
        self.h_2 = self.h_1 + CP.PropsSI('V','P',self.p_1,'T',self.T_1,'Water')*(self.p_2 - self.p_1)/self.eta_pump # [J/kg] enthalpy at state 2
        self.T_2 = CP.PropsSI('T','P',self.p_2,'H',self.h_2,'Water')    # [K] temperature at state 2
        self.s_2 = CP.PropsSI('S','P',self.p_2,'H',self.h_2,'Water')    # [J/kg/K] entropy at state 2   
        self.x_2 = CP.PropsSI('Q','P',self.p_2,'H',self.h_2,'Water')    # [-] vapor quality at state 2  
        self.e_2 = (self.h_2 - self.h_ref) - self.T_ref*(self.s_2 - self.s_ref) # [J/kg] exergy at state 2
        print("State 2 : %f %f %f %f %f %f" % (self.T_2,self.p_2,self.x_2,self.h_2,self.s_2,self.e_2))
    


        # States --------------------------------------------------------------
        self.p           = self.p_1, self.p_2,self.p_3,self.p_4,self.p_5,self.p_6I,self.p_6II,self.p_6III,self.p_6IV,self.p_6V, self.p_6VI,self.p_6VII,self.p_6VIII,self.p_6,self.p_7,self.p_7I,self.p_7II,self.p_7III,self.p_7IV,self.p_7V,self.p_7VI,self.p_7VII,self.p_7VIII,self.p_8,self.p_90,self.p_9I,self.p_9II,self.p_9III,self.p_9IV,self.p_9V,self.p_9VI,self.p_9VII, self.p_9VIII # [Pa]     tuple containing the pressure at each state
        self.T           = self.T_1, self.T_2,self.T_3,self.T_4,self.T_5,self.T_6I,self.T_6II,self.T_6III,self.T_6IV,self.T_6V, self.T_6VI,self.T_6VII,self.T_6VIII,self.T_6,self.T_7,self.T_7I,self.T_7II,self.T_7III,self.T_7IV,self.T_7V,self.T_7VI,self.T_7VII,self.T_7VIII,self.T_8,self.T_90,self.T_9I,self.T_9II,self.T_9III,self.T_9IV,self.T_9V,self.T_9VI,self.T_9VII, self.T_9VIII # [K]      temperature at each state
        self.s           = self.s_1, self.s_2,self.s_3,self.s_4,self.s_5,self.s_6I,self.s_6II,self.s_6III,self.s_6IV,self.s_6V, self.s_6VI,self.s_6VII,self.s_6VIII,self.s_6,self.s_7,self.s_7I,self.s_7II,self.s_7III,self.s_7IV,self.s_7V,self.s_7VI,self.s_7VII,self.s_7VIII,self.s_8,self.s_90,self.s_9I,self.s_9II,self.s_9III,self.s_9IV,self.s_9V,self.s_9VI,self.s_9VII, self.s_9VIII # [J/kg/K] entropy at each state
        self.h           = self.h_1, self.h_2,self.h_3,self.h_4,self.h_5,self.h_6I,self.h_6II,self.h_6III,self.h_6IV,self.h_6V, self.h_6VI,self.h_6VII,self.h_6VIII,self.h_6,self.h_7,self.h_7I,self.h_7II,self.h_7III,self.h_7IV,self.h_7V,self.h_7VI,self.h_7VII,self.h_7VIII,self.h_8,self.h_90,self.h_9I,self.h_9II,self.h_9III,self.h_9IV,self.h_9V,self.h_9VI,self.h_9VII, self.h_9VIII # [J/kg]   enthalpy at each state
        self.x           = self.x_1, self.x_2,self.x_3,self.x_4,self.x_5,self.x_6I,self.x_6II,self.x_6III,self.x_6IV,self.x_6V, self.x_6VI,self.x_6VII,self.x_6VIII,self.x_6,self.x_7,self.x_7I,self.x_7II,self.x_7III,self.x_7IV,self.x_7V,self.x_7VI,self.x_7VII,self.x_7VIII,self.x_8,self.x_90,self.x_9I,self.x_9II,self.x_9III,self.x_9IV,self.x_9V,self.x_9VI,self.x_9VII, self.x_9VIII # [-]      vapor quality at each state
        self.e           = self.e_1, self.e_2,self.e_3,self.e_4,self.e_5,self.e_6I,self.e_6II,self.e_6III,self.e_6IV,self.e_6V, self.e_6VI,self.e_6VII,self.e_6VIII,self.e_6,self.e_7,self.e_7I,self.e_7II,self.e_7III,self.e_7IV,self.e_7V,self.e_7VI,self.e_7VII,self.e_7VIII,self.e_8,self.e_90,self.e_9I,self.e_9II,self.e_9III,self.e_9IV,self.e_9V,self.e_9VI,self.e_9VII, self.e_9VIII # [J/kg]   exergy at each state (use ref conditions)
        self.DAT         = self.p,self.T,self.s,self.h,self.x,self.e



        # Flow rate fraction 6I, 6II, 6III ----------------------------------------------
        #
        # X_6I (h_6I-h_7I) + (X_6II + X_6III) (h_7II-h_7I) = (1 + X_6I + X_6II + X_6III) (h_9I - h_8)
        # X_6II (h_6II-h_7II) + X_6III (h_7III-h_7II) = (1 + X_6I + X_6II + X_6III) (h_9II - h_9I)
        # X_6III (h_6III-h_7III) = (1 + X_6I + X_6II + X_6III) (h_9III - h_9II)
        #
        # Ax = B
        
        A11 = self.h_6I - self.h_7I - self.h_9I + self.h_8
        A12 = self.h_7II - self.h_7I - self.h_9I + self.h_8
        A13 = self.h_7II - self.h_7I - self.h_9I + self.h_8
        A21 = self.h_9I - self.h_9II
        A22 = self.h_6II - self.h_7II - self.h_9II + self.h_9I
        A23 = self.h_7III - self.h_7II - self.h_9II + self.h_9I
        A31 = self.h_9II - self.h_9III
        A32 = self.h_9II - self.h_9III
        A33 = self.h_6III - self.h_7III - self.h_9III + self.h_9II
        B1 = self.h_9I - self.h_8
        B2 = self.h_9II - self.h_9I
        B3 = self.h_9III - self.h_9II

        A = np.array([[A11, A12, A13],
                      [A21, A22, A23],  
                      [A31, A32, A33]])
        
        B = np.array([B1, B2, B3])
        X = np.linalg.solve(A, B)
        self.X_6I = X[0]
        self.X_6II = X[1]
        self.X_6III = X[2]

        #
        #
        #
        #
        #



        # Mass flow rates -----------------------------------------------------
        self.dotm_v      =  0
        self.dotm_tot    =  0
        self.MASSFLOW    = self.dotm_v,self.dotm_tot
        #      o dotm_v         [kg/s]  mass flow rate of water at the condenser
        #      o dotm_tot       [kg/s]  total mass flow rate
        # Efficiencies --------------------------------------------------------
        self.eta_cyclen  = 0
        self.eta_cyclex  = 0
        self.eta_condex  = 0
        self.eta_rotex   = 0
        self.ETA         = self.eta_cyclen,self.eta_cyclex,self.eta_condex,self.eta_rotex
        # -> see text book pp. 53-94
        #      o eta_cyclen    [-]      cycle energy efficiency
        #      o eta_cyclex    [-]      cycle exergy efficiency
        #      o eta_condex    [-]      condenser exergy efficiency
        #      o eta_rotex     [-]      pumps and turbines exergy efficiency
        # Energy losses -------------------------------------------------------
        self.loss_mec    = 0
        self.loss_conden = 0
        self.DATEN       = self.loss_mec,self.loss_conden
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_conden   [W]      condenser energy losses
        # Exergy losses -------------------------------------------------------
        self.loss_rotex  = 0
        self.loss_condex = 0
        self.DATEX       = self.loss_mec,self.loss_rotex,self.loss_condex
        #      o loss_mec      [W]      mechanical energy losses
        #      o loss_rotex    [W]      pumps and turbines exergy losses
        #      o loss_condex   [W]      condenser exergy losses
        self.dotm_6I    = 0
        self.dotm_6II   = 0
        self.dotm_6III  = 0
        self.dotm_6IV   = 0
        self.dotm_6V    = 0
        self.dotm_6VI   = 0
        self.dotm_6VII  = 0
        self.dotm_6VIII = 0
        self.XMASSFLOW  = self.dotm_6I,self.dotm_6II,self.dotm_6III,self.dotm_6IV,self.dotm_6V,self.dotm_6VI,self.dotm_6VII,self.dotm_6VIII        
        # -> massflow rate in each feedheater (w.r.t. fig. 2.33, pp. 91)
        #      o dotm_6I          [kg/s]
        #      o dotm_6II         [kg/s]
        #      o dotm_6...        [kg/s]    
        #      o dotm_6VIII       [kg/s]
        # Energy and Exergy pie charts ----------------------------------------
        #if self.display: self.FIG = self.fig_pie_en,self.fig_pie_ex, self.fig_Ts, self.fig_hs
        #      o fig_pie_en: pie chart of energy losses
        #      o fig_pie_ex: pie chart of exergy losses
        #      o fig_Ts_diagram: T-s diagram of the ST cycle
        #      o fig_hs_diagram: h-s diagram of the ST cycle 