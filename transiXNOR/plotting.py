import numpy as np
import matplotlib.pyplot as plt
PLOT_FIXED_CHARGE = False
T_WINDOW = True
PLOT_BAND = False
PLOT_TRAN = False
PLOT_CHARGE = False
PLOT_CURRENT = True
model_path = './D4'

Eg = 0.252
gridx = np.genfromtxt(model_path + '/gridx.out')
Nx = gridx.shape[0]
gridy = np.genfromtxt(model_path + '/gridy.out')

Vds=0.2; Vbg=0.0; Vtg=0.0
## fixed charge
if (PLOT_FIXED_CHARGE):
	fn_fixed_charge = model_path + '/fixed_charge.npy'
	fixed_charge = np.load(fn_fixed_charge)
	fixed_charge = np.reshape(fixed_charge, (-1, Nx))
	print(fixed_charge.shape)
	plt.plot(gridy, fixed_charge[:,Nx/2])
	plt.show()
	
## band diagram
if (T_WINDOW or PLOT_BAND):
	fn_band_diag = model_path + '/data/phi_%s_%s_%s.npy'%(Vds, Vbg, Vtg)
	band_diag = np.load(fn_band_diag)
	band_diag = np.reshape(band_diag, (-1, Nx))
	Ec = band_diag[:,Nx/2]+Eg/2
	Ev = band_diag[:,Nx/2]-Eg/2
	print('D_Ev: %s' % (Ev[0] - 0))
	print('D_Ec: %s' % (-Vds - Ec[-1]))
	print('Tunneling Window: %s' % (Ev[0] - Ec[-1]))
	# print(band_diag.shape)
	if (PLOT_BAND):
		plt.plot(gridy, Ec)
		plt.plot(gridy, Ev)
		plt.plot(gridy, -np.ones_like(gridy)*Vds, 'r--')
		plt.plot(gridy, np.zeros_like(gridy), 'r--')
		plt.show()

if (PLOT_TRAN):
	fn_tran = model_path + '/data/T_%s_%s_%s.npy'%(Vds, Vbg, Vtg)
	energy_tran = np.load(fn_tran,delimiter=' ')
	energy_tran = energy_tran[energy_tran[:,1] > 0.0, :]
	print(energy_tran.shape)
	plt.semilogy(energy_tran[:,0], energy_tran[:,1])
	plt.show()

if (PLOT_CHARGE):
	fn_charge = model_path + '/data/ncar_%s_%s_%s.npy'%(Vds, Vbg, Vtg)
	charge_density = np.load(fn_charge)
	charge_density = np.reshape(charge_density, (-1, Nx))
	plt.semilogy(gridy, np.abs(charge_density[:,Nx/2]))
	plt.show()

if (PLOT_CURRENT):
	cur = np.load(model_path + '/current_20.npy')
	print(cur)