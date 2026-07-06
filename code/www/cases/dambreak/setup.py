#setup testcase DamBreak_SPH — rupture de barrage

import numpy as np
from src.physics.get_density import get_density
import constant.physical_properties as _phys_mod
import constant.sph_properties as _sph_mod

def init(control,physical, sph, ale, numerical):
    
    # ── Lire les paramètres depuis les dicts ─────────────────
    dx       = sph['dx']
    h        = sph['h']
    n_layers = sph['n_layers']
    FLUID    = sph['FLUID']
    WALL     = sph['WALL']
    rho0     = physical['rho0']
    g_acc    = physical['g_acc']
    gamma    = physical['gamma']
    rp       = numerical['random_pos']
    
    # ── Paramètres du cas DamBreak ─────────────────────────
    H_column = 2  #2 m à terme
    L_column = 1  # 1 m à terme
    c0 = 10.0 * np.sqrt(g_acc * H_column)
    B = (rho0 * c0**2) / gamma
    noise_amp = 0.05 * dx  #scalar for random pos [0.001; 0.15]
    physical['c0'] = c0
    physical['B'] = B
    
    _phys_mod.c0 = c0
    _phys_mod.B  = B
    _sph_mod.h   = h
    _sph_mod.dx  = dx
    
    # ── Générateur de murs ─────────────────────────────────
    box_w, box_h, layers = 4.0, 4.0, 5  #4 par 4 à terme
    # On s'assure que les murs finissent juste avant 0
    x_wall = np.arange(-layers*dx, box_w + layers*dx, dx)
    y_wall = np.arange(-layers*dx, box_h, dx)
    x_grid, y_grid = np.meshgrid(x_wall, y_wall)
    p_grid = np.vstack([x_grid.ravel(), y_grid.ravel()]).T

    # Masque pour ne garder que le "U" de la boite
    mask_wall = (p_grid[:,1] < 0) | (p_grid[:,0] < 0) | (p_grid[:,0] > box_w - dx)
    pos_wall = p_grid[mask_wall]
    types_wall = np.full(len(pos_wall), WALL)

    # Calcul des normales géométriques (Unitaires)
    wall_normals = np.zeros((len(pos_wall), 2))
    for i, p in enumerate(pos_wall):
    # Mur Gauche (x < 0) -> Normale vers la droite (1, 0)
        if p[0] < 0:
            wall_normals[i] = [1.0, 0.0]
        # Mur Droit (x > box_w - dx) -> Normale vers la gauche (-1, 0)
        elif p[0] > box_w - 1.5*dx:
            wall_normals[i] = [-1.0, 0.0]
        # Fond (y < 0) -> Normale vers le haut (0, 1)
        if p[1] < 0:
            # Gestion des coins : on additionne et normalise
            wall_normals[i] += [0.0, 1.0]
    
        # Normalisation pour avoir un vecteur unitaire
        norm = np.linalg.norm(wall_normals[i])
        if norm > 1e-9:
            wall_normals[i] /= norm
            
    # ── Génération de fluide ────────────────────────────────────────
    xf, yf = np.meshgrid(np.arange(0, L_column, dx), 
                     np.arange(0, H_column, dx))
    pos_fluid = np.vstack([xf.ravel(), yf.ravel()]).T
    if rp:
        pos_fluid += (np.random.rand(*pos_fluid.shape) - 0.5) * noise_amp
    
    # ── Assenblage ───────────────────────────────────────────────
    pos   = np.vstack([pos_wall, pos_fluid])
    types = np.concatenate([np.full(len(pos_wall),WALL), np.full(len(pos_fluid),FLUID)])
    n_p   = len(pos)
    wn    = np.zeros((n_p, 2))
    wn[:len(pos_wall)] = wall_normals
    vel   = np.zeros((n_p, 2))
    vol   = np.full(n_p, dx**2)
    rho   = np.full(n_p, rho0)
    pres  = np.zeros(n_p)
    m     = np.zeros(n_p)
    mf = (types==FLUID)
    p_hydro = rho0*g_acc*np.maximum(0., H_column-pos[mf,1])
    rho[mf]  = get_density(p_hydro, B, gamma, rho0, c0)
    pres[mf] = p_hydro
    m[mf]    = rho[mf]*vol[mf]
    rho[~mf] = rho0; m[~mf] = rho0*dx**2
    dt = 0.01*h/c0

    return {'pos':pos,'vel':vel,'rho':rho,'pres':pres,'types':types,
            'vol':vol,'m':m,'wall_normals':wn,
            'num_neighbors':np.zeros(n_p,dtype=int),
            'shepard_coeff':np.zeros(n_p),
            'v_mesh':np.zeros((n_p,2)),'accel':np.zeros((n_p,2)),
            'drho':np.zeros(n_p),'t':0.0,'step':0,'dt':dt,
            'finaltime':control['finaltime']}

