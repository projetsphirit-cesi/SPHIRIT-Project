import numpy as np

def minmod(a, b):
    if a*b <= 0:
        return 0
    else:
        return np.sign(a) * min(abs(a), abs(b))    

def riemann_solver(rho_l, u_l, p_l, rho_r, u_r, p_r, c0):
    rho_avg = 0.5 * (rho_l + rho_r)
    #if rho_avg < 1e-6: return 0.5 * (uL + uR), 0.5 * (pL + pR)
    beta = 0.5  
    p_star = 0.5 * (p_l + p_r) - beta * rho_avg * c0 * (u_r - u_l)
    u_star = 0.5 * (u_l + u_r) - beta * (p_r - p_l) / (rho_avg * c0+ 1e-12)
    #p_star = max(0.0, p_star)
    return u_star, p_star