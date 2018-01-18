import numpy as np
import matplotlib.pyplot as plt
PLOT_FIXED_CHARGE = True
PLOT_BAND = True
PLOT_TRAN = False
PLOT_CHARGE = True
gridx = np.genfromtxt('gridx.out')
gridy = np.genfromtxt("gridy.out")
## fixed charge
if (PLOT_FIXED_CHARGE):
	fn_fixed_charge = 'fixed_charge.out'
	fixed_charge = np.genfromtxt(fn_fixed_charge)
	fixed_charge = np.reshape(fixed_charge, (-1, 11))
	print(fixed_charge.shape)
	plt.plot(gridy, fixed_charge[:,5])
	plt.show()
	
## band diagram
if (PLOT_BAND):
	fn_band_diag = 'datiout/Phi0.0.out'
	band_diag = np.genfromtxt(fn_band_diag)
	band_diag = np.reshape(band_diag, (-1, 11))
	print(band_diag.shape)
	plt.plot(gridy, band_diag[:,5])
	plt.plot(gridy, band_diag[:,5]-0.93)
	plt.show()

if (PLOT_TRAN):
	fn_tran = 'datiout/T0.2.out'
	energy_tran = np.genfromtxt(fn_tran,delimiter=' ')
	energy_tran = energy_tran[energy_tran[:,1] > 0.0, :]
	print(energy_tran.shape)
	plt.plot(energy_tran[:,0], energy_tran[:,1])
	plt.show()

if (PLOT_CHARGE):
	fn_charge = 'datiout/ncar0.0.out'
	charge_density = np.genfromtxt(fn_charge)
	charge_density = np.reshape(charge_density, (-1, 11))
	plt.semilogy(gridy, np.abs(charge_density[:,5]))
	plt.show()