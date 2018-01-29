# Copyright (c) 2018-2019, Mingda Li
# All rights reserved.

from NanoTCAD_ViDES import *
VERBAL=False

class Bi2Se3: 
    def __init__(self,semi,L):
        self.me=semi['me'];   # electron effective mass
        self.mh=semi['mh'];   # hole effective mass
        self.Egap=semi['Eg']; # bandgap
        self.delta=semi['lattice_constant'];
                              # the lattice constant
        self.acc=self.delta/sqrt(3); 
                              # also used in GNR_atoms_coordinates
        self.BC_MX2=semi['relative_EA'];
                              # relative to workfunction of Gr, 
                              # e.g. 0.2 for MoS2
        self.BV_MX2=self.BC_MX2-self.Egap;
        # Derived material parameters:
        self.coeff_Ec = hbar**2 /(2*self.me*m0*q)
        self.coeff_Ev = hbar**2 /(2*self.mh*m0*q)
        self.vf = 6.21e5                    # m/s (2nm Bi2Se3)
        #
        self.deg=1;
        self.n=1;
        ymin=-20;
        ymax=L+20;
        self.Nc=int(4*(floor((floor(L/self.acc)-1)/3)));
        self.Phi=zeros(self.Nc); # electrostatic potential 
        self.Ei=zeros(self.Nc);  # mid-gap potential
        self.Eupper=1000.0;  # upper limit for the energy
        self.Elower=-1000.0; # lower limit for the energy
        self.kmax=pi/self.delta;
        self.kmin=0;
        self.dk=0.1;
        self.dE=1e-3;
        self.thop=-2.59; # Assuming the contact is graphene
        self.mt = self.me * self.mh / (self.me + self.mh)
        self.thop_elec=-sqrt(2*q*self.Egap/(3*(self.acc*1e-9*sqrt(3))**2*m0*self.mt))*hbar/q;
        # print(self.thop_elec)
        # self.thop_elec = hbar * self.vf * (self.kmax * 10^9) / q;
        # print(self.thop_elec)
        self.eta=1e-5;
        self.mu1=0.0;
        self.mu2=0.0;
        self.Temp=300;
        self.E=zeros(NEmax);
        self.T=zeros(NEmax);
        self.charge=zeros(self.Nc);
        self.rank=0;
        self.atoms_coordinates();
        self.gap()
        self.T2D="no"
        self.ymin=ymin;
        self.ymax=ymax;
        self.atoms_coordinates();

    def atoms_coordinates(self):
        GNR_atoms_coordinates(self);
        self.y=array(self.z);
        self.x=zeros(size(self.y));
        return;

    def gap(self):
        return self.Egap;
    def charge_T(self):
        # Number of slices and atoms
        slices=self.Nc;
        atoms=1;
        # I define the vector of the k-wave vector
        kvect=arange(self.kmin,self.kmax,self.dk)
        # I start defining the Hamiltonian for the graphene flake
        h=zeros((2*slices,3),dtype=complex);
        h[0][0]=1;
        for i in range(1,slices+1):
            h[i][0]=i
            h[i][1]=i

        kk=1;
        for ii in range(slices+1,2*slices):
            h[ii][0]=kk;
            h[ii][1]=kk+1;
            kk=kk+1;

        # I then compute the charge and the T for each energy and k and perform the integral
        i=0;
        k=self.kmin;
        H = Hamiltonian(atoms, slices)
        if (self.T2D=="yes"):
            EE=arange(self.Elower,self.Eupper,self.dE);
            kvect=arange(self.kmin,self.kmax+self.dk,self.dk);
            X,Y=meshgrid(EE,kvect);
            Z=zeros((size(EE),size(kvect)))
        while (k<=(self.kmax+self.dk*0.5)):
            if (self.rank==0 and VERBAL): 
                print("----------------------------------")
                print(("    kx range: [%s,%s] ") %(self.kmin,self.kmax));
                print(("    iteration %s ") %i);
                print("----------------------------------")

            # I fill the Hamiltonian for the actual wavevector k in the cycle
            # k [1/m]
            # E0, Egap [eV]
            # coeff_Ec [eV*m^2] <-- hbar [ m^2*kg/s], me [kg]
            # vf: [m/s]
            h[:slices+1:2,2]  = self.BC_MX2 - self.coeff_Ec * k * k * 1e18;
            h[0][2] = 0
            h[1:slices+1:2,2] = self.BV_MX2 + self.coeff_Ev * k * k * 1e18;
            # self.thop_elec = hbar * self.vf * (1j*k + self.kmax) * 1e9 / q;
            self.thop_elec = - hbar * self.vf * k * 1e9 * 0.1 / q;
            h[slices+1::2,2]  = self.thop_elec;
            h[slices+2::4,2]  = self.thop_elec # +self.thop_elec*exp(k*self.delta*1j);
            h[slices+4::4,2]  = self.thop_elec # +self.thop_elec*exp(-k*self.delta*1j);

            H.Eupper = self.Eupper;
            H.Elower = self.Elower;
            H.rank=self.rank;
            H.H = h
            H.dE=self.dE;
            H.Phi=self.Phi;
            H.Ei=-self.Phi;
            H.eta=self.eta;
            H.mu1=self.mu1;
            H.mu2=self.mu2;
            H.Egap=self.gap();
            
            # I then compute T and the charge for the actual kx
            H.charge_T()

            # I sum up all the contribution
            if (i==0):
                self.E=H.E;
                # the factor 2 is because I integrate over kx>0
                self.T=self.deg*H.T*(2*self.dk/(2*pi));
                self.charge=self.deg*H.charge*(2*self.dk/(2*pi));
            else:
                # the factor 2 is because I integrate over kx>0
                self.T=self.T+self.deg*H.T*(2*self.dk/(2*pi));
                self.charge=self.charge+self.deg*H.charge*(2*self.dk/(2*pi));

            if (self.T2D=="yes"):
                print(size(Z[:,i]),size(H.T),size(EE));
                Z[:,i]=H.T[:size(EE)];
            k=k+self.dk
            i=i+1;

        if (self.T2D=="yes"):
            plt.imshow(Z, interpolation='bilinear', cmap=cm.gray,
                       origin='lower', extent=[self.kmin,self.kmax,self.Elower,self.Eupper])
            show()

        del H;
        self.E=array(self.E);
        self.T=array(self.T)*1e9;
        #        self.charge=array(self.charge)*1e9;
        temporary=zeros(size(self.charge));
        #Let's give back the average charge within 2 near atoms
        for ii in range(0,self.Nc,2):
            temporary[ii]=(self.charge[ii]+self.charge[ii+1])*0.5*1e9;
            temporary[ii+1]=(self.charge[ii]+self.charge[ii+1])*0.5*1e9;
        self.charge=temporary
        #        savetxt("ncar.temp",self.charge)
        #    exit(0);
        del kvect,h;
        return;

    def current(self):
        vt=kboltz*self.Temp/q;
        E=array(self.E);
        T=array(self.T);
        arg=2*q*q/(2*pi*hbar)*T*(Fermi((E-self.mu1)/vt)-Fermi((E-self.mu2)/vt))*self.dE
        return sum(arg);
