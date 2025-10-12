from gas_turbine_codes.gas_turbine_group_09 import gas_turbine
import numpy as np
import matplotlib.pyplot as plt


p_1, T_1 = 1e+5, 288.15 # [Pa], [K]
r_c = 18
T_3 = 1673.15 # [K]
eta_pi = 0.9 # [-]
k_mec = 0.015 # [-]
k_cc = 0.95
P_e = 230e+6

inputs = p_1,T_1,P_e
params =  {'T_3':       T_3,
           'r_c':       r_c,
           'eta_pi_c':  eta_pi,
           'eta_pi_t':  eta_pi,
           'k_cc':      k_cc,
           'k_mec':     k_mec,
           'air':       ['N2','O2'],
           'air_prop':  [0.79,0.21],
           'alkane':    [1,4]
           }

parametric_study_T3 = True
parametric_study_r = True


if parametric_study_T3:
    # Parametric study on T3
    T3_range = np.linspace(1000,3000,100)
    COMBUSTION = np.zeros((len(T3_range),3))
    MASSFLOW = np.zeros((len(T3_range),3))
    ETA = np.zeros((len(T3_range),7))
    DATEN = np.zeros((len(T3_range),2))
    DATEX = np.zeros((len(T3_range),4))

    for i in range(len(T3_range)):
        print("Temperature : ", T3_range[i])
        params['T_3'] = T3_range[i]
        GT = gas_turbine(inputs,params, False)
        GT.evaluate()
        COMBUSTION[i] = [GT.COMBUSTION[0], GT.COMBUSTION[1], GT.COMBUSTION[2]]
        MASSFLOW[i] = [GT.MASSFLOW[0], GT.MASSFLOW[1], GT.MASSFLOW[2]]
        ETA[i] = [GT.ETA[0], GT.ETA[1], GT.ETA[2], GT.ETA[3], GT.ETA[4], GT.ETA[5], GT.ETA[6]]
        DATEN[i] = [GT.DATEN[0], GT.DATEN[1]]
        DATEX[i] = [GT.DATEX[0], GT.DATEX[1], GT.DATEX[2], GT.DATEX[3]]

    params["T_3"] = 1673.15

    # save data in file
    np.savetxt("DATA_T3/COMBUSTION.csv", COMBUSTION, delimiter=",")
    np.savetxt("DATA_T3/MASSFLOW.csv", MASSFLOW, delimiter=",")
    np.savetxt("DATA_T3/ETA.csv", ETA, delimiter=",")
    np.savetxt("DATA_T3/DATEN.csv", DATEN, delimiter=",")
    np.savetxt("DATA_T3/DATEX.csv", DATEX, delimiter=",")


if parametric_study_r:
    # Parametric study on T3
    r_range = np.linspace(3,100,100)
    COMBUSTION = np.zeros((len(r_range),3))
    MASSFLOW = np.zeros((len(r_range),3))
    ETA = np.zeros((len(r_range),7))
    DATEN = np.zeros((len(r_range),2))
    DATEX = np.zeros((len(r_range),4))

    for i in range(len(r_range)):
        print("ratio : ", r_range[i])
        params['r_c'] = r_range[i]
        GT = gas_turbine(inputs,params, False)
        GT.evaluate()
        COMBUSTION[i] = [GT.COMBUSTION[0], GT.COMBUSTION[1], GT.COMBUSTION[2]]
        MASSFLOW[i] = [GT.MASSFLOW[0], GT.MASSFLOW[1], GT.MASSFLOW[2]]
        ETA[i] = [GT.ETA[0], GT.ETA[1], GT.ETA[2], GT.ETA[3], GT.ETA[4], GT.ETA[5], GT.ETA[6]]
        DATEN[i] = [GT.DATEN[0], GT.DATEN[1]]
        DATEX[i] = [GT.DATEX[0], GT.DATEX[1], GT.DATEX[2], GT.DATEX[3]]

    # save data in file
    np.savetxt("DATA_r/COMBUSTION.csv", COMBUSTION, delimiter=",")
    np.savetxt("DATA_r/MASSFLOW.csv", MASSFLOW, delimiter=",")
    np.savetxt("DATA_r/ETA.csv", ETA, delimiter=",")
    np.savetxt("DATA_r/DATEN.csv", DATEN, delimiter=",")
    np.savetxt("DATA_r/DATEX.csv", DATEX, delimiter=",")
