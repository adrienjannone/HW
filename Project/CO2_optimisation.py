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
    'eta_turb': eta_turb,
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

ETATURB = False
OPTITES = True

"""
OPTIMIZATION OF TES STORAGE TEMPERATURE AND PRESSURE
Variables:
- T_w_hot : Water inlet temperature to TES0
- T_TES1_in : CO2 inlet temperature to TES1
- T_TES2_in : CO2 inlet temperature to TES2
- T_TES3_in : CO2 inlet temperature to TES3
- T_TES4_in : CO2 inlet temperature to TES4
- T_TES1_out : CO2 outlet temperature from TES1
- T_TES2_out : CO2 outlet temperature from TES2
- T_TES3_out : CO2 outlet temperature from TES3
- T_TES4_out : CO2 outlet temperature from TES4
- p_TES1 : CO2 pressure in TES1
- p_TES2 : CO2 pressure in TES2
- p_TES3 : CO2 pressure in TES3
- p_TES4 : CO2 pressure in TES4
"""

if OPTITES:
    def obj_function(x):
        dome = CO2_battery(inputs, params, False)
        dome.TES_Parameters(x)
        results = dome.evaluate()

        mass_flow_CO2 = results[6]
        loss_TES0 = results[18]
        loss_TES1 = results[19]
        loss_TES2 = results[20]
        loss_TES3 = results[21]
        loss_TES4 = results[22]
        TD2 = dome.T_D2
        TD3 = dome.T_D3
        TD4 = dome.T_D4
        TD5 = dome.T_D5

        # Gros coup de burin si la physique n'est pas respectée
        if mass_flow_CO2 > 100: return np.inf
        if loss_TES0 < 0: return np.inf
        if loss_TES1 < 0: return np.inf
        if loss_TES2 < 0: return np.inf
        if loss_TES3 < 0: return np.inf
        if loss_TES4 < 0: return np.inf
        if TD2 > x[3]: return np.inf
        if TD3 > x[4]: return np.inf
        if TD4 > x[5]: return np.inf
        if TD5 > x[6]: return np.inf

        loss  = loss_TES0 + loss_TES1 + loss_TES2 + loss_TES3 + loss_TES4
        return loss

    # x: [T_TES1_in, T_TES2_in, T_TES3_in, T_TES1_out, T_TES2_out, T_TES3_out, T_TES4_out, p_TES1, p_TES2, p_TES3, p_TES4, TS0_in]
    x0 = [42+273.15, 97+273.15, 300+273.15,
          25+273.15, 43+273.15, 98+273.15, 301+273.15,
          1.1e5, 1.1e5, 1.1e5, 1.1e5,
          25+273.15] 
    
    # Bounds for the variables
    bounds = [(30+273.15, 60+273.15),   # T_TES1_in
              (80+273.15, 99+273.15), # T_TES2_in
              (280+273.15, 320+273.15), # T_TES3_in
              (15+273.15, 40+273.15),   # T_TES1_out
              (30+273.15, 60+273.15),   # T_TES2_out
              (80+273.15, 120+273.15),  # T_TES3_out
              (300+273.15, 400+273.15), # T_TES4_out
              (1e5, 5e5),                # p_TES1
              (1e5, 5e5),                # p_TES2
              (1e5, 5e5),                # p_TES3
              (1e5, 5e5),                # p_TES4
              (20+273.15, 40+273.15)]    # TS0_in
    
    # Perform the optimization
    result = optimize.minimize(obj_function, x0, bounds=bounds, method='SLSQP')
    optimized_x = result.x

    dome = CO2_battery(inputs, params, True)
    dome.TES_Parameters(optimized_x)
    res = dome.evaluate()


    print("Optimized Variables:")
    dome.print_results()




#ANALYSE PARAMÉTRIQUE SUR ETA_TURB 
if ETATURB:
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

