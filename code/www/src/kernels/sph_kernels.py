import numpy as np

def kernel_cubic(r, h):
    q = r / h
    sigma = 10.0 / (7.0 * np.pi * h**2)
    return np.where(q <= 1.0, sigma * (1 - 1.5*q**2 + 0.75*q**3), 
           np.where(q <= 2.0, sigma * 0.25 * (2 - q)**3, 0.0))

def kernel_grad(rij, r, h):
    q = r / h
    sigma_g = 10.0 / (7.0 * np.pi * h**3)
    dw_r = np.where(q <= 1.0, sigma_g * (-3.0*q + 2.25*q**2),
           np.where(q <= 2.0, -sigma_g * 0.75 * (2.0 - q)**2, 0.0))
    return (dw_r / (r + 1e-12)) * rij, dw_r

def get_pressure(rho,B,rho0,gamma):
    # Équation de Tait (sans clipping agressif pour garder la cohésion)
    return B * (np.power(rho / rho0, gamma) - 1.0)