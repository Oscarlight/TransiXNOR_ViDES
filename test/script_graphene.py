#!/usr/bin/python
from NanoTCAD_ViDES import *
from pylab import *
# atoms per slice (or cell)
atoms=1;
slices=8;
thop=-2.7;
k=30
#4*pi/(3*sqrt(3)*0.144);
delta=0.144*sqrt(3);

h=zeros((16,3),dtype=complex);
h[0][0]=1;
for i in range(1,9):
    h[i][0]=i
    h[i][1]=i

h[9][0]=1;
h[9][1]=2;
h[9][2]=thop;

h[10][0]=2;
h[10][1]=3;
h[10][2]=thop+thop*exp(k*delta*1j);

h[11][0]=3;
h[11][1]=4;
h[11][2]=thop;

h[12][0]=4;
h[12][1]=5;
h[12][2]=thop+thop*exp(-k*delta*1j);

h[13][0]=5;
h[13][1]=6;
h[13][2]=thop;

h[14][0]=6;
h[14][1]=7;
h[14][2]=thop+thop*exp(k*delta*1j);

h[15][0]=7;
h[15][1]=8;
h[15][2]=thop;

H = Hamiltonian(atoms, slices)
H.Eupper = 10
H.Elower = -10
H.H = h
H.mu1=H.mu2=0.3
H.dE=1e-2;
H.charge_T()

plot(H.E,H.T)
show()
#plot(H.charge);
#hold
#plot(GNR2.charge,'o');
#show()

a=[H.E,H.T]

savetxt("T",transpose(a));
