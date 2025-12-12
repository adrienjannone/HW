from scipy import optimize
import CoolProp.CoolProp as CP
import numpy as np
from CO2_battery import CO2_battery
import pyomo.environ as pyo
from pyomo.opt import SolverFactory


# Goal 
#   optimize the water inlet temperature to minimize CO2 mass flow rate
#   optimize the CO2 storage temperature to minimize CO2 mass flow rate
#   optimize the dome temperature to minimize CO2 mass flow rate


# From the slides
p_storage_co2_liquid = 70 * 1e5     # Pa
T_storage_water = 22 + 273.15       # K
T_storage_TES = 450 + 273.15        # K
Pe = 20e6                           # W

# From the Paper
T_amb = 15 + 273.15                 # K
p_dome = 1e5                        # Pa
k_dome = 0.005                      # Pressure drop coefficient dome (inlet, outlet)
k_TES = 4e-2                        # Pressure drop coefficient TES (inlet, outlet)
eta_comp = 0.85                     # Compressor adiabatic efficiency (isentropic)
eta_turb = 0.9                      # Turbine adiabatic efficiency (isentropic)
eta_pump = 0.85                     # Pump adiabatic efficiency (isentropic)
eta_mec = 0.99                      # Mechanical efficiency (shaft)
eta_elec = 0.985                    # Electrical efficiency (motor/generator)
pinch_TES0 = 2                      # K
pinch_TES = 7.5                     # K
pinch_PCHX = 5                      # K
p_amb = 1e5                         # Pa
inputs = Pe
params = {
    'p_storage_co2_liquid': p_storage_co2_liquid,
    'T_storage_water': T_storage_water,
    'T_storage_TES': T_storage_TES,
    'T_amb': T_amb,
    'p_amb': p_amb,
    'k_dome': k_dome,
    'k_TES': k_TES,
    'eta_comp': eta_comp,
    'eta_turb': eta_turb,
    'eta_pump': eta_pump,
    'eta_mec': eta_mec,
    'eta_elec': eta_elec,
    'pinch_TES0': pinch_TES0,
    'pinch_TES': pinch_TES,
    'pinch_PCHX': pinch_PCHX,
    'fluid': 'CO2',
}


"""
Play with the TES input and ouput temperatures to minimize the CO2 mass flow rate
"""
def objective(x):
    T_storage_water_opt = x[0]
    T_storage_TES_opt = x[1]
    params_opt = params.copy()
    params_opt['T_storage_water'] = T_storage_water_opt
    params_opt['T_storage_TES'] = T_storage_TES_opt
    dome_opt = CO2_battery(inputs, params_opt, False)
    dome_opt.evaluate()
    return dome_opt.m_dot_CO2

x0 = [T_storage_water, T_storage_TES]
bounds = [(273.15 + 10, 273.15 + 40), (273.15 + 400, 273.15 + 600)]
result = optimize.minimize(objective, x0, bounds=bounds, method='L-BFGS-B')

print("Optimal water storage temperature (K): ", result.x[0])
print("Optimal TES storage temperature (K): ", result.x[1])
print("Minimum CO2 mass flow rate (kg/s): ", result.fun)