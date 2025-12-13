from scipy import optimize
import CoolProp.CoolProp as CP
import numpy as np
from CO2_battery import CO2_battery
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import matplotlib.pyplot as plt



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
    #'eta_turb': eta_turb,
    'eta_pump': eta_pump,
    'eta_mec': eta_mec,
    'eta_elec': eta_elec,
    'pinch_TES0': pinch_TES0,
    'pinch_TES': pinch_TES,
    'pinch_PCHX': pinch_PCHX,
    'fluid': 'CO2',
}


#dome = CO2_battery(inputs, params, False)
#p, T, h, s, x, e, m_dot_CO2, m_dot_TS0, m_dot_TSE1, m_dot_TSE2, m_dot_TSE3, m_dot_TSE4, eta_rotex, eta_transex_TES0, eta_transex_TES1, eta_transex_TES2, eta_transex_TES3, eta_transex_TES4, eta_transex_PCHX, loss_rotex, loss_evaporator, loss_TES1, loss_TES2, loss_TES3, loss_TES4, loss_PCHX, loss_mec, loss_elec, total_energy = dome.evaluate()


#ANALYSE PARAMÉTRIQUE SUR ETA_TURB 

eta_turb_min = 0.90
eta_turb_max = 0.95
num_steps = 11 
eta_turb_range = np.linspace(eta_turb_min, eta_turb_max, num_steps)


results = {'eta_turb': [],'total_energy': [],'eta_rotex': [], 'm_dot_CO2': []}

for eta_turb_val in eta_turb_range:
    
    params_iter = params.copy()
    params_iter['eta_turb'] = eta_turb_val
    dome = CO2_battery(inputs, params_iter, False)
    results_tuple = dome.evaluate()
    m_dot_CO2 = results_tuple[6]
    eta_rotex = results_tuple[12]
    total_energy = results_tuple[-1] 
    results['eta_turb'].append(eta_turb_val)
    results['total_energy'].append(total_energy)
    results['eta_rotex'].append(eta_rotex * 100) 
    results['m_dot_CO2'].append(m_dot_CO2)



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

ax1.plot(results['eta_turb'], results['total_energy'], 'r-o', label='Total energy (W)')
ax1.set_xlabel('Efficiency of the turbine, $\eta_{turb}$ [-]')
ax1.set_ylabel('Total energy used [W]')
ax1.set_title('Sensibility of the totan energy used with $\eta_{turb}$')
ax1.grid(True)



ax2.plot(results['eta_turb'], results['eta_rotex'], 'g-^', label='Efficiency Rotex (%)')
ax2.set_xlabel('Efficiency of the turbine,, $\eta_{turb}$ [-]')
ax2.set_ylabel('Efficiency Rotex (Cycle), $\eta_{rotex}$ [%]')
ax2.set_title('$\eta_{rotex}$ en fonction de $\eta_{turb}$')
ax2.grid(True)

plt.tight_layout()
plt.show()