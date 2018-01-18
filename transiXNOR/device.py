from NanoTCAD_ViDES import * 
import sys 
from module_TMD import *

rank = 0
  
# I create the grid
# gate direction mesh from -2.0 to 2.0
# from -2.0 to 0: mesh size from 1 to 0.05
# from 0 to 2.0: mesh size from 0.05 to 1
xg=nonuniformgrid(
    array([-0.4,0.2,
              0,0.1,
            0.4,0.2]))
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
    # 'mh': 0.124,
    'Eg': 0.252,
    'acc': 0.3,  # ?????????
    'relative_EA' : 0.2
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
savetxt("gridx.out",grid.gridx)
savetxt("gridy.out",grid.gridy)

# I take care of the solid
Oxide1=region("hex",grid.xmin,0,grid.ymin,grid.ymax)
Oxide1.eps=3.9;

Oxide2=region("hex",0,grid.xmax,grid.ymin,grid.ymax)
Oxide2.eps=3.9;

top_gate=gate("hex",grid.xmax,grid.xmax,10.0,20.0);
bottom_gate=gate("hex",grid.xmin,grid.xmin,10.0,20.0);


p=interface2D(grid,Oxide1,Oxide2,top_gate,bottom_gate);

# molar fraction
fraction_source=0.002 # p-doped
fraction_drain=-0.002 # n-doped
dope_reservoir(grid,p,FLAKE,fraction_source,array([-1,1,grid.ymin,10.0]));
dope_reservoir(grid,p,FLAKE,fraction_drain,array([-1,1,20.0,grid.ymax]));

savetxt("er.out", p.eps)
savetxt("fixed_charge.out", p.fixed_charge)

# ------------------------------------------#
p.underel=0.01; # ?????

Vtgmax=0.2;
Vtgmin=0.0;
VtgN=2;

Vbgmax=0.0;
Vbgmin=0.2;
VbgN=1;

Vdsmax=0.2;
Vdsmin=0.2;
VdsN=1;

vds_cur = []
for vds in np.linspace(Vdsmin, Vdsmax, VdsN):
    FLAKE.mu1=0.0
    FLAKE.mu2=vds
    vbg_cur = []
    for vbg in np.linspace(Vbgmin, Vbgmax, VbgN):
        vtg_cur = []
        for vtg in np.linspace(Vtgmin, Vtgmax, VtgN):
            bottom_gate.Ef=vbg 
            set_gate(p,bottom_gate)
            top_gate.Ef=vtg; 
            set_gate(p,top_gate)
            p.normpoisson=1e-1;
            p.normd=5e-2; # 5e-3;
            solve_self_consistent(grid,p,FLAKE);
            vtg_cur.append(FLAKE.current());
            # I save the output files
            if (rank==0):
                savetxt("./datiout/phi_%s_%s_%s.out" % (vds, vbg, vtg),
                    p.Phi)
                savetxt("./datiout/ncar_%s_%s_%s.out" % (vds, vbg, vtg),
                    p.free_charge);
                savetxt("./datiout/T_%s_%s_%s.out" % (vds, vbg, vtg),
                    transpose([FLAKE.E,FLAKE.T]));
        vbg_cur.append(vtg_cur)
    vds_cur.append(vbg_cur)
np.save('current.npy', np.array(vds_cur))