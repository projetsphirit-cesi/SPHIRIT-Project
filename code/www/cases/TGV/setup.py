#setup testcase TGV_SPH — Taylor-Green Vortex

import numpy as np, os
import constant.physical_properties as _phys_mod
import constant.sph_properties as _sph_mod

def init(control, physical, sph, ale, numerical):
    FLUID = sph['FLUID']
    domain_size     = 1.0    
    rho0  = 1.0      
    U0    = 1.0     
    nu    = 0.01    
    mu    = rho0 * nu  
    gamma = 7.0
    c0    = 10.0 * U0  
    B     = rho0 * c0**2 / gamma
    dx    = sph.get('dx', 0.02)      
    h_sph = 1.2 * dx
    CFL   = 0.01      # ← L211
    g     = np.array([0.0, 0.0])    
 
    # C4 — Override complet du physical{} — Solver lit depuis phys{}
    _phys_mod.rho0     = rho0;  _phys_mod.c0         = c0
    _phys_mod.B        = B;     _phys_mod.mu         = mu
    _phys_mod.g        = g;     _phys_mod.g_acc      = 9.81
    _phys_mod.gamma    = gamma ; 
    _phys_mod.periodic = True; _phys_mod.domain_size = domain_size
    _sph_mod.CFL       = CFL

    def get_tgv_solution(pos, t, L, U0, nu):     
        k     = 2 * np.pi / L  
        decay = np.exp(-2 * k**2 * nu * t)     
        u =  U0 * np.sin(k*pos[:,0]) * np.cos(k*pos[:,1]) * decay      
        v = -U0 * np.cos(k*pos[:,0]) * np.sin(k*pos[:,1]) * decay      
        p = (rho0*U0**2/4.0) * (np.cos(2*k*pos[:,0]) + np.cos(2*k*pos[:,1])) * decay**2      # ← L322
        return np.column_stack((u, v)), p      
    x = np.arange(dx/2, domain_size, dx)      # ← L2403
    X, Y = np.meshgrid(x, x)
    pos   = np.column_stack((X.ravel(), Y.ravel()))     
    N     = len(pos)
 
    vel, p_init = get_tgv_solution(pos, 0.0, domain_size, U0, nu)     
    rho   = rho0 * (1.0 + p_init / (rho0 * c0**2))
    vol   = np.full(N, dx**2)
    m     = rho * vol
    pres  = p_init.copy()
    types = np.full(N, FLUID, dtype=int)
    wn    = np.zeros((N, 2))   # pas de murs
 
    dt = 0.01 * h_sph / c0      # ← L2416
    print(f'TGV : {N} particules  domain_size={domain_size}  dx={dx}  c0={c0}  periodic=True')
 
    return {'pos':pos, 'vel':vel, 'rho':rho, 'pres':pres, 'types':types,
            'vol':vol, 'm':m, 'wall_normals':wn,
            'num_neighbors':np.zeros(N,dtype=int),
            'shepard_coeff':np.zeros(N),
            'v_mesh':np.zeros((N,2)), 'accel':np.zeros((N,2)),
            'drho':np.zeros(N), 't':0.0, 'step':0, 'dt':dt,
            'finaltime':control['finaltime'],
            'periodic':True, 'domain_size':domain_size,  # utilisé par Euler.py ligne 26
            'U0':U0, 'nu':nu}   # pour comparaison analytique en fin de run
