import numpy as np

def get_density(p,B,gamma,rho0,c0):
    p_clamped = np.clip(p, -0.9*B, 5.0*rho0*c0**2)
    return rho0 * np.power(np.maximum((p_clamped / B) + 1.0, 1e-6),
                       1.0 / gamma)
