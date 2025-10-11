import numpy as np
import matplotlib.pyplot as plt

# T3 impact study

T3_range = np.linspace(1000,3000,100)

# # IMPACT ON MASSFLOW RATES
# dotm_fuel = np.loadtxt("DATA_T3/MASSFLOW.csv", delimiter=",")[:,1]
# dotm_air = np.loadtxt("DATA_T3/MASSFLOW.csv", delimiter=",")[:,0]
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T3_range, dotm_fuel, 'g-')
# ax2.plot(T3_range, dotm_air, 'b-')
# ax1.set_xlabel('T3 (K)')
# ax1.set_ylabel('Fuel mass flow (kg/s)', color='g')
# ax2.set_ylabel('Air mass flow (kg/s)', color='b')
# plt.title('Impact of T3 on mass flow rates')
# ax1.grid()
# plt.savefig('figures/T3_mass_flow_rates.png')

# # IMPACT ON EFFICIENCIES
# eta_cyclen = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,0]
# eta_meca = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,2]
# eta_toten = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,1]
# eta_totex = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,3]
# eta_cyclex = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,4]
# eta_rotex = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,5]
# eta_combex = np.loadtxt("DATA_T3/ETA.csv", delimiter=",")[:,6]

# plt.figure()
# plt.plot(T3_range, eta_cyclen, 'g-', label='Cycle efficiency')
# plt.plot(T3_range, eta_toten, 'r-', label='Energy efficiency')
# plt.plot(T3_range, eta_totex, 'b-', label='Exergy efficiency')
# plt.plot(T3_range, eta_meca, 'y-', label='Mechanical efficiency')
# plt.plot(T3_range, eta_cyclex, 'c-', label='Cycle exergy efficiency')
# plt.plot(T3_range, eta_rotex, 'm-', label='Rotational exergy efficiency')
# plt.plot(T3_range, eta_combex, 'k-', label='Combustion exergy efficiency')
# plt.xlabel('T3 (K)')
# plt.ylabel('Efficiency')
# plt.title('Impact of T3 on efficiencies')
# plt.legend()
# plt.grid()
# plt.savefig('figures/T3_efficiencies.png')
# plt.show()

# # IMPACT ON EN LOSSES
# Loss_mec = np.loadtxt("DATA_T3/DATEN.csv", delimiter=",")[:,0]
# Loss_echen = np.loadtxt("DATA_T3/DATEN.csv", delimiter=",")[:,1]
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T3_range, Loss_mec, 'g-', label='Mechanical losses')
# ax2.plot(T3_range, Loss_echen, 'b-', label='Energy losses')
# ax1.set_xlabel('T3 (K)')
# ax1.set_ylabel('Mechanical losses (W)', color='g')
# ax2.set_ylabel('Energy losses (W)', color='b')
# plt.title('Impact of T3 on losses')
# ax1.grid()
# fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
# plt.savefig('figures/T3_losses_en.png')
# plt.show()

# # IMPACT ON EX LOSSES
# L_mec = np.loadtxt("DATA_T3/DATEX.csv", delimiter=",")[:,0]
# L_rotex = np.loadtxt("DATA_T3/DATEX.csv", delimiter=",")[:,1]
# L_combex = np.loadtxt("DATA_T3/DATEX.csv", delimiter=",")[:,2]
# L_echex = np.loadtxt("DATA_T3/DATEX.csv", delimiter=",")[:,3]
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(T3_range, L_mec, 'g-', label='Mechanical exergy losses')
# ax1.plot(T3_range, L_rotex, 'b-', label='Rotational exergy losses')
# ax2.plot(T3_range, L_combex, 'r-', label='Combustion exergy losses')
# ax2.plot(T3_range, L_echex, 'y-', label='Chimney exergy losses')
# ax1.set_xlabel('T3 (K)')
# ax1.set_ylabel('Mechanical  and Rotational exergy losses (W)', color='g')
# ax2.set_ylabel('Combustion and Chimney exergy losses (W)', color='b')
# plt.title('Impact of T3 on exergy losses')
# ax1.grid()
# fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
# plt.savefig('figures/T3_losses_ex.png')
# plt.show()





####################################################################################################################################################

# r impact study

r_range = np.linspace(3,100,100)

# # IMPACT ON MASSFLOW RATES
# dotm_fuel = np.loadtxt("DATA_r/MASSFLOW.csv", delimiter=",")[:,1]
# dotm_air = np.loadtxt("DATA_r/MASSFLOW.csv", delimiter=",")[:,0]
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(r_range, dotm_fuel, 'g-')
# ax2.plot(r_range, dotm_air, 'b-')
# ax1.set_xlabel('Compression ratio')
# ax1.set_ylabel('Fuel mass flow (kg/s)', color='g')
# ax2.set_ylabel('Air mass flow (kg/s)', color='b')
# plt.title('Impact of compression ratio on mass flow rates')
# ax1.grid()
# plt.savefig('figures/r_mass_flow_rates.png')
# plt.show()

# # IMPACT ON EFFICIENCIES
# eta_cyclen = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,0]
# eta_meca = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,2]
# eta_toten = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,1]
# eta_totex = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,3]
# eta_cyclex = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,4]
# eta_rotex = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,5]
# eta_combex = np.loadtxt("DATA_r/ETA.csv", delimiter=",")[:,6]

# plt.figure()
# plt.plot(r_range, eta_cyclen, 'g-', label='Cycle efficiency')
# plt.plot(r_range, eta_toten, 'r-', label='Energy efficiency')
# plt.plot(r_range, eta_totex, 'b-', label='Exergy efficiency')
# plt.plot(r_range, eta_meca, 'y-', label='Mechanical efficiency')
# plt.plot(r_range, eta_cyclex, 'c-', label='Cycle exergy efficiency')
# plt.plot(r_range, eta_rotex, 'm-', label='Rotational exergy efficiency')
# plt.plot(r_range, eta_combex, 'k-', label='Combustion exergy efficiency')

# plt.xlabel('Compression ratio')
# plt.ylabel('Efficiency')
# plt.title('Impact of compression ratio on efficiencies')
# plt.legend()
# plt.grid()
# plt.savefig('figures/r_efficiencies.png')
# plt.show()

# # IMPACT ON EN LOSSES
# Loss_mec = np.loadtxt("DATA_r/DATEN.csv", delimiter=",")[:,0]
# Loss_echen = np.loadtxt("DATA_r/DATEN.csv", delimiter=",")[:,1]

# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(r_range, Loss_mec, 'g-', label='Mechanical losses')
# ax2.plot(r_range, Loss_echen, 'b-', label='Energy losses')
# ax1.set_xlabel('Compression ratio')
# ax1.set_ylabel('Mechanical losses (W)', color='g')
# ax2.set_ylabel('Energy losses (W)', color='b')
# plt.title('Impact of compression ratio on losses')
# ax1.grid()
# fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
# plt.savefig('figures/r_losses_en.png')
# plt.show()

# # IMPACT ON EX LOSSES
# L_mec = np.loadtxt("DATA_r/DATEX.csv", delimiter=",")[:,0]
# L_rotex = np.loadtxt("DATA_r/DATEX.csv", delimiter=",")[:,1]
# L_combex = np.loadtxt("DATA_r/DATEX.csv", delimiter=",")[:,2]
# L_echex = np.loadtxt("DATA_r/DATEX.csv", delimiter=",")[:,3]
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(r_range, L_mec, 'g-', label='Mechanical exergy losses')
# ax1.plot(r_range, L_rotex, 'b-', label='Rotational exergy losses')
# ax2.plot(r_range, L_combex, 'r-', label='Combustion exergy losses')
# ax2.plot(r_range, L_echex, 'y-', label='Chimney exergy losses')
# ax1.set_xlabel('Compression ratio')
# ax1.set_ylabel('Mechanical  and Rotational exergy losses (W)', color='g')
# ax2.set_ylabel('Combustion and Chimney exergy losses (W)', color='b')
# plt.title('Impact of compression ratio on exergy losses')
# ax1.grid()
# fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
# plt.savefig('figures/r_losses_ex.png')
# plt.show()