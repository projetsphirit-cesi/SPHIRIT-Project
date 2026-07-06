#setup testcase Funnel_SPH — entonnoir

import numpy as np
from src.physics.get_density import get_density
from src.mesh.primitives import add_position_noise
from src.mesh.funnel_mesh import generate_funnel_wall, get_funnel_radius ,compute_wall_normals
from scipy.spatial import cKDTree

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
    periodic = False
    domain_size        = 1.0
    
    # ── Paramètres du cas funnel ─────────────────────────
    H_column = 0.04  #0.0821 m à terme
    # --- Paramètres de l'entonnoir ---
    H_funnel = H_column
    R_top = 0.8
    R_bottom = 0.15
    #L_column = 1  # 1 m à terme
    c0 = 20.0 * np.sqrt(g_acc * 0.0621)
    noise_amp = 0.05 * dx  #scalar for random pos [0.001; 0.15]
    B = (rho0 * c0**2) / gamma
    physical['c0'] = c0
    physical['B'] = B
    
    # ── Générateur de murs ─────────────────────────────────
    pos_wall = generate_funnel_wall(dx, n_layers,H_funnel)
            
    # ── Génération de fluide ────────────────────────────────────────
    # --- Génération du Fluide ---
    # --- Génération du Fluide (Adaptée à la forme) ---
    pos_fluid_list = []

    # On définit les hauteurs (y) de chaque rangée de particules fluides
    y_fluid_levels = np.arange(dx, H_column, dx)

    for y in y_fluid_levels:
        # Rayon maximal disponible à cette hauteur (avec une petite marge dx)
        R_max_fluid = get_funnel_radius(y) - dx
    
        # Nombre de particules pouvant tenir sur cette ligne
        # On centre la ligne sur x = 0
        n_particles_row = int(2 * R_max_fluid / dx)
    
        if n_particles_row > 0:
            # On génère les x de manière symétrique et régulière
            x_row = np.linspace(-R_max_fluid, R_max_fluid, n_particles_row)
            for x in x_row:
                pos_fluid_list.append([x, y])

    pos_fluid = np.array(pos_fluid_list)
    pos_fluid = np.asarray(pos_fluid).reshape(-1,2)
    if rp :    
        pos_fluid = add_position_noise(pos_fluid, noise_amp)
    
    # ── Assenblage ───────────────────────────────────────────────
        # --- Après avoir défini pos_wall et pos_fluid ---
    pos = np.vstack([pos_wall, pos_fluid])
    types = np.concatenate([np.full(len(pos_wall), WALL), np.full(len(pos_fluid), FLUID)])
    n_p = len(pos)
    # --- Juste après la génération de pos_wall ---
    wall_normals = np.zeros((len(pos),2))
    wall_normals[types==WALL] = compute_wall_normals(pos_wall)

    # Initialisation des tableaux à la bonne taille
    vel = np.zeros((n_p, 2))
    rho = np.full(n_p, rho0)
    pres = np.zeros(n_p)
    vol = np.full(n_p, dx**2)
    m = np.zeros(n_p)

    m = rho * vol
    # --- Définition des masques initiaux ---
    mask_f = (types == FLUID)
    mask_w = (types == WALL)

    # 1. Calcul du volume SPH réel (1/somme des noyaux)
    if periodic: 
        tree = cKDTree(pos, boxsize=domain_size)
    else :
        tree = cKDTree(pos)
    # Utilisation d'une requête globale plus rapide
    vol[:] = dx**2
    # 2. Pression hydrostatique théorique
    # On utilise la hauteur maximale du fluide pour le calcul (H_column)
    y_fluid = pos[mask_f, 1]
    p_hydro = rho0 * g_acc * np.maximum(0.0, (H_column - y_fluid))

    # 3. Application des valeurs initiales
    rho[mask_f] = get_density(p_hydro, B, gamma, rho0, c0)
    pres[mask_f] = p_hydro
    m[mask_f] = rho[mask_f] * vol[mask_f]

    # Pour les murs : densité de base et masse calculée sur le volume SPH
    rho[mask_w] = rho0
    pres[mask_w] = 0.0
    m[mask_w] = rho0 * vol[mask_w]
#e_int[mask_f] = pres[mask_f] / ((gamma - 1.0) * rho[mask_f])
#E[mask_f] = e_int[mask_f] + 0.5 * np.sum(vel[mask_f]**2, axis=1)
#e = np.zeros(n_p)
#e[mask_f] = pres[mask_f] / ((gamma - 1) * rho[mask_f])
    t, step, dt = 0.0, 0, 0.01 * h / c0


    vel[types == FLUID] = 0
    vel[types == WALL]  = 0
    rho[types == WALL]  = rho0
    #d_m = np.zeros(n_p)
    drho = np.zeros(n_p)

    accel = np.zeros((n_p, 2))
    num_neighbors = np.zeros(n_p, dtype=int)
    shepard_coeff = np.zeros(n_p)
    
    return {'pos':pos,
            'vel':vel,
            'rho':rho,
            'pres':pres,
            'types':types,
            'vol':vol,
            'm':m,
            'wall_normals':wall_normals,
            'num_neighbors':np.zeros(n_p,dtype=int),
            'shepard_coeff':np.zeros(n_p),
            'v_mesh':np.zeros((n_p,2)),'accel':np.zeros((n_p,2)),
            'drho':np.zeros(n_p),'t':0.0,'step':0,'dt':dt,
            'finaltime':control['finaltime']}

