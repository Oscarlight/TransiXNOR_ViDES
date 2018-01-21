import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import rcParams
from matplotlib.ticker import AutoMinorLocator
import glob
from NanoTCAD_ViDES import *

FIGURE_SIZE = (5, 9)
# FONT_SIZE = 26
LINE_WIDTH = 1
MAJOR_LABEL_SIZE = 22
MINOR_LABEL_SIZE = 0
rcParams['figure.figsize'] = FIGURE_SIZE
rcParams['axes.linewidth'] = LINE_WIDTH
fig, ax = plt.subplots()

COMPUTE_CURRENT_FROM_T = False
PLOT_FIXED_CHARGE      = False
T_WINDOW               = True
PLOT_BAND              = False
PLOT_TRAN              = False
PLOT_CHARGE            = False
PLOT_CURRENT           = True
PLOT_CURRENT_SPECTRUM  = True
model_path             = './D5'

Eg = 0.252
vt=kboltz*300/q;
gridx = np.genfromtxt(model_path + '/gridx.out')
Nx = gridx.shape[0]
gridy = np.genfromtxt(model_path + '/gridy.out')
Ny = gridy.shape[0]

Vds=0.2; Vbg=0.2; Vtg=0.0
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
	fn_band_diag = model_path + '/data/phi_%.2f_%.2f_%.2f.npy'%(Vds, Vbg, Vtg)
	band_diag = np.load(fn_band_diag)
	band_diag = np.reshape(band_diag, (-1, Nx))
	Ec = band_diag[:,Nx/2]+Eg/2
	Ev = band_diag[:,Nx/2]-Eg/2
	print('D_Ev: %s' % (Ev[0] - 0))
	print('D_Ec: %s' % (-Vds - Ec[-1]))
	print('Tunneling Window: %s' % (Ev[0] - Ec[-1]))
	print('Source-size barrier: %s' % (Ec[Ny/2] - Ev[0]))
	print('Drain-size barrier:  %s' % (Ec[-1] - Ev[Ny/2]))
	# print(band_diag.shape)
	if (PLOT_BAND):
		plt.plot(gridy, Ec)
		plt.plot(gridy, Ev)
		plt.plot(gridy, -np.ones_like(gridy)*Vds, 'r--')
		plt.plot(gridy, np.zeros_like(gridy), 'r--')
		plt.show()

if (PLOT_TRAN or PLOT_CURRENT_SPECTRUM):
	fn_tran = model_path + '/data/T_%.2f_%.2f_%.2f.npy'%(Vds, Vbg, Vtg)
	energy_tran = np.load(fn_tran)
	energy_tran = energy_tran[energy_tran[:,1] > 0.0, :]
	E = energy_tran[:, 0]
	T = energy_tran[:, 1]
	if (PLOT_CURRENT_SPECTRUM):
		jE = np.abs(2*q*q/(2*pi*hbar)*T*(Fermi((E-0)/vt)-Fermi((E-Vds)/vt)))
		plt.semilogy(E, jE, linewidth=2, color='k')
		plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
		plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
		ax.xaxis.set_minor_locator(AutoMinorLocator())
		plt.savefig('current_spectrum.png')
	if (PLOT_TRAN):
		print(energy_tran.shape)
		plt.semilogy(energy_tran[:,0], energy_tran[:,1])
		plt.show()

if (PLOT_CHARGE):
	fn_charge = model_path + '/data/ncar_%.2f_%.2f_%.2f.npy'%(Vds, Vbg, Vtg)
	charge_density = np.load(fn_charge)
	charge_density = np.reshape(charge_density, (-1, Nx))
	plt.semilogy(gridy, np.abs(charge_density[:,Nx/2]))
	plt.show()

if (PLOT_CURRENT):
	cur = np.load(model_path + '/current_20.npy')
	print(cur)

if (COMPUTE_CURRENT_FROM_T):
	vdsmin=0.0; vdsmax=0.2; vdsN=21;
	vbgmin=0.0; vbgmax=0.2; vbgN=21;
	vtgmin=0.0; vtgmax=0.2; vtgN=21;
	vds_cur = []
	for vds in np.linspace(vdsmin, vdsmax, vdsN):
		vbg_cur = []
		for vbg in np.linspace(vbgmin, vbgmax, vbgN):
			vtg_cur = []
			for vtg in np.linspace(vtgmin, vtgmax, vtgN):
				tran = np.load(model_path + '/data/T_%.2f_%.2f_%.2f.npy'%(vds, vbg, vtg))
				E = tran[:, 0]
				T = tran[:, 1]
				dE = 1e-3
				vtg_cur.append(
					sum(2*q*q/(2*pi*hbar)*T*(Fermi((E)/vt)-Fermi((E-vds)/vt))*dE)
				)
			print('length of vtg_cur = %s' % len(vtg_cur))
			vbg_cur.append(vtg_cur)
		vds_cur.append(vbg_cur)
	np.save(model_path+'/current', np.array(vds_cur))

