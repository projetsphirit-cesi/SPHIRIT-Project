import numpy as np

rho0 = 1000.0
mu = 1.0  # Viscosité dynamique (Pa.s)
g_acc = 9.81
g = np.array([0.0, -g_acc])
H_ref = 1.0
c0    = 10.0 * np.sqrt(g_acc * H_ref)
gamma = 7.0
B = (rho0 * c0**2) / gamma
periodic = False
domain_size = 1.0