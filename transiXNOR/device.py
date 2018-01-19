from NanoTCAD_ViDES import * 
import sys
import os
from module_TMD import *
import pickle
from math import sqrt
import argparse

rank = 0
OVERWRITE=True
model_path = './D3'

parser = argparse.ArgumentParser()
parser.add_argument("--vdsmin", default=0.2, type=float)
parser.add_argument("--vdsmax", default=0.2, type=float)
parser.add_argument("--vdsN",   default=1,   type=int)
args = parser.parse_args()

if not os.path.exists(model_path):
    os.makedirs(model_path)
    os.makedirs(model_path + '/data')
elif(not OVERWRITE):
    print("model path: %s exists. Quit!" % model_path)
    quit()
# I create the grid
# gate direction mesh from -2.0 to 2.0
# from -2.0 to 0: mesh size from 1 to 0.05
# from 0 to 2.0: mesh size from 0.05 to 1
# TODO: It seems certain grid can lead to segmentation faults.
xg=nonuniformgrid(
    array([-1.0,0.5,
              0,0.1,
            1.0,0.5]))
## MoTe2 (http://pubs.acs.org/doi/pdf/10.1021/acs.jpcc.5b02950)
# semi = {
#     'me': 0.65,
#     'mh': 0.64,
#     'Eg': 0.93,
#     'acc': 0.36,
#     'relative_EA' : 0.2
# }
## Bi2Se3 (Qin's Paper)
semi = {
    'me': 0.124,
    'mh': 2.23,
    'Eg': 0.252,
    'acc': 0.414/sqrt(3), # ref: http://iopscience.iop.org/article/10.1088/1367-2630/12/6/065013/meta 
     # e.g. 0.2 in MoS2: the distance between
     # Mo and S, the lattice constant = acc * sqrt(3)
    'relative_EA': 0.09,      # relative to workfunction of Gr, 
                              # e.g. 0.2 for MoS2
    'fraction_source': 0.004, # p-dope
    'fraction_drain': -0.0015, # n-dope
}
FLAKE=TMD(semi,30.0,"n");

acc=FLAKE.acc;
kF=2*pi/(3*sqrt(3)*acc);
kymax=pi/FLAKE.delta;
Nky=32.0;
dk=kymax/Nky;
FLAKE.kmax=pi/FLAKE.delta;
FLAKE.kmin=0;
FLAKE.dk=dk;

FLAKE.dE=0.001
grid=grid2D(xg,FLAKE.y,FLAKE.x,FLAKE.y);

with open(model_path+"/material.p", "wb") as f:
    pickle.dump(semi,f)

savetxt(model_path+"/gridx.out",grid.gridx)
savetxt(model_path+"/gridy.out",grid.gridy)

# I take care of the solid
Cox = 25/1.1; # 1.1 nm HfO2 for er=25
Csemi = 100/0.7;
er_equ = Cox * Csemi / (Cox + Csemi)
Oxide1=region("hex",grid.xmin,0,grid.ymin,grid.ymax)
Oxide1.eps=er_equ; # !!!

Oxide2=region("hex",0,grid.xmax,grid.ymin,grid.ymax)
Oxide2.eps=er_equ; # !!!

top_gate=gate("hex",grid.xmax,grid.xmax,10.0,20.0);
bottom_gate=gate("hex",grid.xmin,grid.xmin,10.0,20.0);


p=interface2D(grid,Oxide1,Oxide2,top_gate,bottom_gate);

# molar fraction
dope_reservoir(grid,p,FLAKE,
    semi['fraction_source'],
    array([-1,1,grid.ymin,10.0]));
dope_reservoir(grid,p,FLAKE,
    semi['fraction_drain'],
    array([-1,1,20.0,grid.ymax]));

savetxt(model_path+"/er.out", p.eps)
savetxt(model_path+"/fixed_charge.out", p.fixed_charge)

# ------------------------------------------#
p.underel=0.01; # ?????????

Vtgmax=0.2;
Vtgmin=0.0;
VtgN=21;

Vbgmax=0.2;
Vbgmin=0.0;
VbgN=21;

Vdsmax=args.vdsmax;
Vdsmin=args.vdsmin;
VdsN=args.vdsN;

vds_cur = []
for vds in np.linspace(Vdsmin, Vdsmax, VdsN):
    FLAKE.mu1=0.0
    FLAKE.mu2=vds
    vbg_cur = []
    for vbg in np.linspace(Vbgmin, Vbgmax, VbgN):
        vtg_cur = []
        for vtg in np.linspace(Vtgmin, Vtgmax, VtgN):
            print('>>> Vds=%s, Vbg=%s, Vtg=%s' % (vds, vbg, vtg))
            bottom_gate.Ef=vbg 
            set_gate(p,bottom_gate)
            top_gate.Ef=vtg; 
            set_gate(p,top_gate)
            p.normpoisson=1e-1;
            p.normd=1e-3; # 5e-3;
            solve_self_consistent(grid,p,FLAKE);
            vtg_cur.append(FLAKE.current());
            # I save the output files
            if (rank==0):
                np.save(model_path+"/data/phi_%s_%s_%s" % (vds, vbg, vtg),
                    p.Phi)
                np.save(model_path+"/data/ncar_%s_%s_%s" % (vds, vbg, vtg),
                    p.free_charge);
                np.save(model_path+"/data/T_%s_%s_%s" % (vds, vbg, vtg),
                    transpose([FLAKE.E,FLAKE.T]));
        vbg_cur.append(vtg_cur)
    vds_cur.append(vbg_cur)

np.save(model_path+'/current_%s'%int(Vdsmin*100), 
    np.array(vds_cur))
