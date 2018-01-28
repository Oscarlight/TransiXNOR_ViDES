import numpy as np
import matplotlib.pyplot as plt

hbar = 1.05459e-34 # m2kg/s
q = 1.6e-19
m0 = 9.1095e-31
mc = 0.124*m0
mv = 2.23*m0
# Effective mass approximation
Ec_curve = hbar**2 /(2*mc*q)
Ev_curve = hbar**2 /(2*mv*q)
Eg = 0.252
Ec = 0.125
Ev = Ec - Eg
E0 = (Ec + Ev)/2
# Only consider the parabolic term
D = -(Ec_curve - Ev_curve)/2
B = -(Ec_curve + Ev_curve)/2
vf = 6.21e5        # m/s


def bands(k):
	a = E0 - D*k*k    # eV
	b = Eg/2 - B*k*k  # eV
	c = np.sqrt(hbar**2*vf*k*k/q + b*b)
	return (a + c, a - c)

kmax = 3.14/0.414
k = np.linspace(0, kmax, 100) # 1/nm
k = k * 1e9  # 1/m
Ec_k, Ev_k = bands(k)

Ec_eff = Ec + hbar**2 * k*k/(2*mc*q)
Ev_eff = Ev - hbar**2 * k*k/(2*mv*q)
plt.plot(k, Ec_k, 'r')
plt.plot(k, Ev_k, 'r')
plt.plot(k, Ec_eff, 'b--')
plt.plot(k, Ev_eff, 'b--')
plt.show()

