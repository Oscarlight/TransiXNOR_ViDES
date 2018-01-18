from NanoTCAD_ViDES import * 
import sys 
from module_TMD import *

rank = 0
  
# I create the grid
# gate direction mesh from -2.0 to 2.0
# from -2.0 to 0: mesh size from 1 to 0.05
# from 0 to 2.0: mesh size from 0.05 to 1
xg=nonuniformgrid(
    array([-2.0,1,
              0,0.5,
            2.0,1]))
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
    'acc': 0.3,
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
fraction_source=0.001 # p-doped
fraction_drain=-0.001 # n-doped
dope_reservoir(grid,p,FLAKE,fraction_source,array([-1,1,grid.ymin,10.0]));
dope_reservoir(grid,p,FLAKE,fraction_drain,array([-1,1,20.0,grid.ymax]));

# solve_init(grid,p,FLAKE);

Vgmin=0.0;
Vgmax=0.0;
Vgstep=0.05;

Np=int(abs(Vgmin-Vgmax)/Vgstep)+1;
vg=zeros(Np);
current=zeros(Np);
p.underel=0.1;

counter=0;
Vgs=Vgmin;
FLAKE.mu1=0.0
# Vds = mu2 - mu1
FLAKE.mu2=0.3


savetxt("er.out", p.eps)
savetxt("fixed_charge.out", p.fixed_charge)

while (Vgs<=Vgmax):
    bottom_gate.Ef=Vgs 
    set_gate(p,bottom_gate)
    top_gate.Ef=Vgs; 
    set_gate(p,top_gate)
    p.normpoisson=1e-1;
    p.normd=5e-2; # 5e-3;
    solve_self_consistent(grid,p,FLAKE);
    vg[counter]=Vgs;
    current[counter]=FLAKE.current();
    # I save the output files
    if (rank==0):
        savetxt("./datiout/Phi%s.out" %Vgs,
            p.Phi)
        savetxt("./datiout/ncar%s.out" %Vgs,p.free_charge);
        savetxt("./datiout/T%s.out" %Vgs,
            transpose([FLAKE.E,FLAKE.T]));
    counter=counter+1;
    Vgs=Vgs+Vgstep;
savetxt("./datiout/idvgs.out",transpose([vg,current]));