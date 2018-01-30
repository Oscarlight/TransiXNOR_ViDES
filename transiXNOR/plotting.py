import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pylab import rcParams
from matplotlib.ticker import AutoMinorLocator
import glob
from NanoTCAD_ViDES import *
import argparse

# -------------------------------------------------
# Because the 0-thickness 2D material, apply the 
# following assumption:
# 	  (Vds, Vbg, Vtg) <--> (Vds, 0, Vtg + Vbg) 
# therefore, we compute:
# 	  current_vds.npy --> (vds, 0, vtg:0 -> 0.4)
# -------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--vtg", default=0.0, type=float)
parser.add_argument("--vbg", default=0.2, type=float)
parser.add_argument("--vds", default=0.2, type=float)
parser.add_argument("model", default='./D6', type=str)
args = parser.parse_args()

FIGURE_SIZE = (5, 6)
LINE_WIDTH = 1
MAJOR_LABEL_SIZE = 22
MINOR_LABEL_SIZE = 0
rcParams['figure.figsize'] = FIGURE_SIZE
rcParams['axes.linewidth'] = LINE_WIDTH
# rcParams['figure.autolayout'] = True
fig, ax = plt.subplots()

COMBINE_CURRENT_VIA_SYMMETRY = False
COMPUTE_CURRENT_FROM_T = False
PLOT_FIXED_CHARGE      = False
T_WINDOW               = True
PLOT_BAND              = False
PLOT_TRAN              = False
PLOT_CHARGE            = False
PLOT_CURRENT           = False
PLOT_CURRENT_SPECTRUM  = False
PLOT_FAMILY_CURVES     = False
PRINT_CURRENT_ONLY     = False
model_path             = args.model

Eg = 0.252
vt=kboltz*300/q;
gridx = np.genfromtxt(model_path + '/gridx.out')
Nx = gridx.shape[0]
gridy = np.genfromtxt(model_path + '/gridy.out')
Ny = gridy.shape[0]

Vds=args.vds; Vbg=args.vbg; Vtg=args.vtg
voltage = '%.2f_%.2f_%.2f' % (Vds, Vbg, Vtg)
## fixed charge
if (PLOT_FIXED_CHARGE):
	fn_fixed_charge = model_path + '/fixed_charge.npy'
	fixed_charge = np.load(fn_fixed_charge)
	fixed_charge = np.reshape(fixed_charge, (-1, Nx))
	plt.plot(gridy, fixed_charge[:,Nx/2])
	plt.show()
	
## band diagram
if (T_WINDOW or PLOT_BAND):
	fn_band_diag = model_path + '/data/phi_' + voltage + '.npy'
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
		plt.plot(gridy, Ec, linewidth=2, color='k')
		plt.plot(gridy, Ev, linewidth=2, color='k')
		plt.plot(gridy[-20:], -np.ones_like(gridy[-20:])*Vds,
			linewidth=2, color='#e74c3c', linestyle='-')
		plt.plot(gridy[:20], np.zeros_like(gridy[:20]),
			linewidth=2, color='#3498db', linestyle='-')
		plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
		plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
		ax.xaxis.set_minor_locator(AutoMinorLocator())
		ax.yaxis.set_minor_locator(AutoMinorLocator())
		plt.ylim([-0.6, 0.4])
		plt.xlim([-2, 40])
		plt.savefig(model_path+'/plots/band_' + voltage + '.pdf',
			bbox_inches='tight', transparent=True)
		plt.clf()

if (PLOT_TRAN or PLOT_CURRENT_SPECTRUM):
	fn_tran = model_path + '/data/T_' + voltage + '.npy'
	energy_tran = np.load(fn_tran)
	energy_tran = energy_tran[energy_tran[:,1] > 0.0, :]
	E = energy_tran[:, 0]
	T = energy_tran[:, 1]
	if (PLOT_CURRENT_SPECTRUM):
		# Fig. 1
		jE = np.abs(2*q*q/(2*pi*hbar)*T*(Fermi((E-0)/vt)-Fermi((E-Vds)/vt)))
		plt.semilogx(jE, -E, linewidth=2, color='k')
		plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
		plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
		# ax.xaxis.set_minor_locator(AutoMinorLocator())
		ax.yaxis.set_minor_locator(AutoMinorLocator())
		plt.ylim([-0.6, 0.4])
		plt.xlim([1e-3, 5e3])
		plt.savefig(model_path+'/plots/current_spectrum_' + voltage + '.pdf', 
			bbox_inches='tight', transparent=True)
		plt.clf()
	if (PLOT_TRAN):
		print(energy_tran.shape)
		plt.semilogy(energy_tran[:,0], energy_tran[:,1])
		plt.show()

if (PLOT_CHARGE):
	fn_charge = model_path + '/data/ncar_' + voltage + '.npy'
	charge_density = np.load(fn_charge)
	charge_density = np.reshape(charge_density, (-1, Nx))
	plt.semilogy(gridy, np.abs(charge_density[:,Nx/2]))
	plt.show()

if (PLOT_CURRENT or PRINT_CURRENT_ONLY):
	# Fig. 1
	cur = np.abs(np.load(model_path + '/current_20.npy'))
	if (PRINT_CURRENT_ONLY):
		print(cur); quit()
	vtg_array = np.linspace(0.0, 0.2, 21)
	# vbg = 0. and vtg from 0 to 0.2 <--> vbg = 0 and vtg from 0. to 0.2
	plt.semilogy(vtg_array, cur[0,0,:21], linewidth=2, color='k')
	plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
	plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	plt.ylim([1e-2, 1e2])
	plt.savefig(model_path+'/plots/current_vds_0.2_vbg_0.0.pdf',
		bbox_inches='tight', transparent=True)
	plt.clf()
	# vbg = 0.1 and vtg from 0 to 0.2 <--> vbg = 0 and vtg from 0.1 to 0.3
	plt.semilogy(vtg_array, cur[0,0,10:31], linewidth=2, color='k')
	plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
	plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	plt.ylim([1e-2, 1e2])
	plt.savefig(model_path+'/plots/current_vds_0.2_vbg_0.1.pdf',
		bbox_inches='tight', transparent=True)
	plt.clf()
	# vbg = 0.2 and vtg from 0 to 0.2 <--> vbg = 0 and vtg from 0.2 to 0.4
	plt.semilogy(vtg_array, cur[0,0,20:41], linewidth=2, color='k')
	plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
	plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	plt.ylim([1e-2, 1e2])
	plt.savefig(model_path+'/plots/current_vds_0.2_vbg_0.2.pdf',
		bbox_inches='tight', transparent=True)
	plt.clf()	

if (COMBINE_CURRENT_VIA_SYMMETRY):
	vdsmin=0.01; vdsmax=0.2; vdsN=20;
	vbgmin=0.0; vbgmax=0.2; vbgN=21;
	vtgmin=0.0; vtgmax=0.2; vtgN=21;
	vds_cur = []
	print('Start combine all current together via symmetry from current_*.npy')
	for vds in np.linspace(vdsmin, vdsmax, vdsN):
		cur_array = np.abs(np.load(model_path + '/current_' + str(int(vds*100)) + '.npy'))
		vbg_cur = []
		for vbg in np.linspace(vbgmin, vbgmax, vbgN):
			vtg_cur = []
			for vtg in np.linspace(vtgmin, vtgmax, vtgN):
				cur = cur_array[0, 0, int(vbg*100) + int(vtg*100)]
				vtg_cur.append(cur)
			vbg_cur.append(vtg_cur)
		vds_cur.append(vbg_cur)
	cur_map = np.array(vds_cur)
	cur_map = np.concatenate((np.zeros_like(cur_map[0:1,:,:]), cur_map), axis=0)
	print(cur_map.shape)
	np.save(model_path+'/current', cur_map)

if (PLOT_FAMILY_CURVES):
	cur = np.abs(np.load(model_path + '/current.npy'))
	vtg_list = [0.0, 0.05, 0.1, 0.15, 0.2]
	vds_array = np.linspace(0, 0.2, 21)
	plt.clf()
	# Vbg = 0
	for vtg in vtg_list:
		plt.plot(vds_array, cur[:, 0, int(vtg*100)], 
			linewidth=2, color='k')
	plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
	plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	plt.ylim([-1, 85])
	plt.savefig(model_path+'/plots/family_curve_vbg_0.pdf',
		bbox_inches='tight', transparent=True)
	plt.clf()	
	# Vbg = 0.2
	for vtg in vtg_list:
		# cur[7, 20, int(vtg*100)] = (cur[8, 20, int(vtg*100)] + cur[6, 20, int(vtg*100)])/2
		plt.plot(vds_array, cur[:, 20, int(vtg*100)], 
			linewidth=2, color='k')
	plt.tick_params(axis='both', which='major', length=10, labelsize=MAJOR_LABEL_SIZE)
	plt.tick_params(axis='both', which='minor', length=5, labelsize=MINOR_LABEL_SIZE)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	plt.ylim([-1, 85])
	plt.savefig(model_path+'/plots/family_curve_vbg_1.pdf',
		bbox_inches='tight', transparent=True)
	plt.clf()

if (COMPUTE_CURRENT_FROM_T):
	vdsmin=0.01; vdsmax=0.2; vdsN=20;
	vbgmin=0.0; vbgmax=0.2; vbgN=21;
	vtgmin=0.0; vtgmax=0.2; vtgN=21;
	vds_cur = []
	print('Start computing the current from the T files...')
	for vds in np.linspace(vdsmin, vdsmax, vdsN):
		vbg_cur = []
		for vbg in np.linspace(vbgmin, vbgmax, vbgN):
			vtg_cur = []
			for vtg in np.linspace(vtgmin, vtgmax, vtgN):
				tmp_volt = '%.2f_%.2f_%.2f' % (vds, vbg, vtg)
				tran = np.load(model_path + '/data/T_' + tmp_volt + '.npy')
				E = tran[:, 0]
				T = tran[:, 1]
				dE = 1e-3
				vtg_cur.append(
					sum(2*q*q/(2*pi*hbar)*T*(Fermi((E)/vt)-Fermi((E-vds)/vt))*dE)
				)
			# print('length of vtg_cur = %s' % len(vtg_cur))
			vbg_cur.append(vtg_cur)
		vds_cur.append(vbg_cur)
	np.save(model_path+'/current', np.array(vds_cur))




