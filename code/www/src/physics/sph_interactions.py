import numpy as np
import system.control_dict as _ctrl
import constant.physical_properties as _phys
import constant.sph_properties as _sph
from scipy.spatial import cKDTree
from src.kernels.sph_kernels import kernel_cubic, kernel_grad
from src.physics.riemann_solver import riemann_solver
from src.discretization.ale_mesh import periodic_rij
from constant.numerical_properties import visco_on, car_spring_k, car_spring_damping
from system.sph_schemes import mode_p_refinement, shepard_corr

CAR = 2
BODY = 3

def compute_shepard_coefficients(pairs, pos, vol, h,
                                 shepard_coeff,
                                 num_neighbors,
                                 kernel_cubic,periodic,domain_size):
    """
    Calcule les coefficients de normalisation de Shepard
    et le nombre de voisins SPH.
    """

    h2 = 2.0 * h

    for i, j in pairs:
        if periodic:
            rij = periodic_rij(pos[i], pos[j],domain_size)
        else :
            rij = pos[i]-pos[j]# vector from j to i
        
        r = np.linalg.norm(rij)
        
        # Skip interactions hors support
        if r < 1e-12 or r > h2:
            continue

        w_ij = kernel_cubic(r, h)

        # --- contribution symétrique ---
        shepard_coeff[i] += w_ij * vol[j]
        shepard_coeff[j] += w_ij * vol[i]

        num_neighbors[i] += 1
        num_neighbors[j] += 1

def compute_fluid_fluid_interaction(i, j, pos, vel, v_mesh, rho, pres, vol, m, 
                                    shepard_coeff, accel, drho, 
                                    h, rho0, c0, mu, visco_on,grad_rho, grad_p, grad_vel,shepard,periodic,domain_size):
    """
    Calcule l'interaction entre deux particules de fluide (i, j) 
    selon le formalisme SPH-ALE avec solveur de Riemann.
    """
    if periodic:
            rij = periodic_rij(pos[i], pos[j], domain_size)
    else :
            rij = pos[i]-pos[j]# vector from j to i
        
    r = np.linalg.norm(rij)
    
    if r < 1e-12 or r > 2*h:
        return

    # 1. Géométrie de l'interface
    grad_w, grad_w_mag = kernel_grad(rij, r, h)  #vector from i to j              
    n_ij = rij / r
    
    
    
    # 2. Formalisme ALE : Vitesse de transport (Vitesse de l'interface)
    # On utilise la vitesse vmesh (déplacement des particules) pour définir l'interface
    v_interface = 0.5 * np.dot(v_mesh[i] + v_mesh[j], n_ij)
    
    # 3. Solveur de Riemann (États de surface)
    # n_ij pointe de j vers i, donc j est l'état Gauche (L) et i l'état Droit (R)
    u_l   = np.dot(vel[j], n_ij) 
    u_r   = np.dot(vel[i], n_ij)
    rho_l = rho[j]
    rho_r = rho[i]
    p_l   = pres[j]
    p_r   = pres[i]
    if mode_p_refinement ==1:
        mid_dist = 0.5 * rij        
        # --- RECONSTRUCTION MUSCL ---
        # Etat Gauche (i) et Droit (j)
        rho_l +=  np.dot(grad_rho[j], mid_dist)
        rho_r += - np.dot(grad_rho[i], mid_dist)
        p_l   +=  np.dot(grad_p[j], mid_dist)
        p_r   +=  - np.dot(grad_p[i], mid_dist)
        u_l  = np.dot(vel[j] + grad_vel[j] @ mid_dist, n_ij) 
        u_r  = np.dot(vel[i] - grad_vel[i] @ mid_dist, n_ij) 
        
    # Appel au solveur de Riemann (attention à bien passer j en premier, puis i)
    u_star, p_star = riemann_solver(rho_l, u_l, p_l, rho_r, u_r, p_r, c0)
    
    # Vitesse relative (Vitesse matérielle Riemann - Vitesse de transport ALE)
    u_relative = u_star - v_interface 
    p_max = 2*rho0*c0**2
    #p_star = np.clip(p_star,-0.5*B,p_max)
    #rho_star = get_density(p_star)        
    rho_star = 0.5*(rho[j]+rho[i])
    # 4. Calcul du Vecteur Surface (A_vec)
    # On peut inclure une correction de Shepard ici si nécessaire (s_i, s_j)
    A_vec = vol[i] * vol[j] * grad_w           
    A_mag = abs(np.dot(A_vec, n_ij))     #positive value bc n_ij from j to i and A_vec from i to j      
    if shepard_corr:
        A_mag /= (max(shepard[i], 0.8)*max(shepard[j], 0.8) +1e-12 )       
    # 5. Flux de Masse (Équation de continuité)
    # mass_flux = rho * (u_matériel - v_transport) * Surface
    mass_flux = rho_star * u_relative * A_mag           
    
    drho[i] += 2.0 * mass_flux / vol[i]
    drho[j] -= 2.0 * mass_flux / vol[j]

    # 6. Flux de Quantité de Mouvement (Équation de Momentum)
    # F_mom = (Flux de masse * Vitesse de transport) + Force de pression
    f_mom = (rho_star * u_star * u_relative + p_star) * A_mag * n_ij #from j to i
    
    accel[i] += 2.0 * f_mom / m[i]  
    accel[j] -= 2.0 * f_mom / m[j]          

    # 7. Viscosité Physique (Morris et al.)
    if visco_on:
        v_diff = vel[i] - vel[j]
        # On utilise une formulation stable pour les petits r
        visc_term = (2.0 * mu * r * grad_w_mag) / (rho[i] * rho[j] * (r**2 + 0.01*h**2))
        f_visc_ij = visc_term * v_diff
        
        accel[i] += 2.0 * m[j] * f_visc_ij / m[i]  
        accel[j] -= 2.0 * m[i] * f_visc_ij / m[j]


def compute_fluid_wall_interaction(i, j, pos, vel, types, rho, pres, vol, m, 
                                   shepard_coeff, wall_normals_full, accel, 
                                   h, mu, gamma, visco_on, FLUID,shepard,periodic, grad_rho, grad_p, grad_vel):
    """
    Calcule l'interaction entre une particule de fluide et une particule de mur
    en utilisant un solveur de Riemann (Miroir No-Slip) et le formalisme SPH-ALE.
    """
    
    # 1. Identification Fluide / Mur
    idx_f, idx_w = (i, j) if types[i] == FLUID else (j, i)
    
    
    rij = pos[idx_w] - pos[idx_f]
    r = np.linalg.norm(rij)
    
    if r < 1e-12 or r > 2*h:
        return

    # 2. Géométrie et Normale du mur
    n_fw = wall_normals_full[idx_w] #from w to f
    # Sécurité orientation : la normale doit pointer vers le fluide
    if np.dot(n_fw, pos[idx_f] - pos[idx_w]) < 0:
        n_fw = -n_fw
        
    grad_w, grad_w_mag_fw = kernel_grad(rij, r, h)
    
    # 3. Correction de surface (Shepard)
    # On utilise le coefficient de la particule fluide pour compenser la troncature
    
    # Vecteur surface efficace
    A_vec = (vol[idx_f] * vol[idx_w] * grad_w) 
    # Amplitude de la surface projetée sur la normale
    A_mag = abs(np.dot(A_vec, n_fw))
    
    if shepard_corr:
        A_mag /= (max(shepard[idx_f], 0.8) +1e-12 )  
    # 4. États de Riemann (Miroir No-Slip)
    u_l = np.dot(vel[idx_f], n_fw)
    p_l = pres[idx_f]
    rho_l = rho[idx_f]
    cL = _phys.c0 * (rho_l / _phys.rho0)**((gamma - 1.0)/2.0)
    
    # Pour un mur fixe No-Slip, la vitesse cible à l'interface est 0
    u_star =  0.0 
    # Stabilisation de Riemann pour la pression au mur
    #p_star = pL - rhoL * cL * (u_star - uL)
    
    #u_star, p_star = riemann_solver(rho[idx_f], uL, pres[idx_f], rho[idx_f], -uL, pres[idx_f], c0)
    #u_star, p_star = riemann_solver(rhoL, uL, pL, rho[i], uR, pres[i], c0) #test 29/03
    # Vitesse relative ALE (ici v_interface_wall = 0 car le mur ne bouge pas)
    u_rel_w = u_star - 0.0 
    rho_star = rho_l
    # Pression à l'interface p* (Solveur de Riemann acoustique exact pour un mur)
    # Si le fluide tape le mur (v_impact > 0), le mur crée une surpression acoustique.
    # S'il s'éloigne (v_impact < 0), on garde simplement la pression locale.
    v_impact = np.dot(vel[idx_f], -n_fw)
    p_star = p_l + rho_l * _phys.c0 * max(v_impact, 0.0)
    # 5. Flux de Quantité de Mouvement
    # Momentum = (Pression + Terme advectif ALE) * Surface
    # Note: u_star * u_rel_w est nul pour un mur fixe, seul p_star reste.
    momentum_flux_w = (rho_star * u_star * u_rel_w + p_star) * A_mag * n_fw

    # 6. Application des Accélérations    
    accel[idx_f] += 2.0 * momentum_flux_w / m[idx_f] 

    # 7. Viscosité au mur (Condition No-Slip physique)
    if visco_on:
        # La vitesse relative au mur est simplement vel[idx_f] car vel[wall]=0
        visc_w = (2.0 * mu * r * grad_w_mag_fw) / (rho[idx_f] * rho[idx_w] * (r**2 + 0.01*h**2))
        f_visc_wall = visc_w * vel[idx_f]
        
        accel[idx_f] -= 2.0 * m[idx_w] * f_visc_wall / m[idx_f]        
    
    if (_ctrl.test_case == 'car' and types[idx_w] == CAR) or \
    (_ctrl.test_case == 'bodyF' and types[idx_w] == BODY):
        accel[idx_w] -= 2.0 * momentum_flux_w / m[idx_w]

def apply_density_filter(pos, m, rho, vol, types, h, periodic,domain_size):
    """
    Réinitialisation du champ de densité (Filtre de Shepard d'ordre 0)
    pour stabiliser la pression et corriger la troncature en surface libre.
    """
    
    if periodic: 
        tree = cKDTree(pos, boxsize=domain_size)
    else :
        tree = cKDTree(pos)
        
    n_p = len(pos)
    new_rho = np.copy(rho)
    
    # On ne filtre que les particules fluides
    mask_f = (types == 0) # FLUID
    
    for i in np.where(mask_f)[0]:
        idx = tree.query_ball_point(pos[i], 2*h)
        if len(idx) == 0: continue
        if periodic:
            rij = periodic_rij(pos[i], pos[idx])
        else :
            rij = pos[i]-pos[idx]
        
        r = np.linalg.norm(rij, axis=1)
       
        W = kernel_cubic(r, h)
        
        # Formule de Shepard pour la densité
        num = np.sum(m[idx] * W)
        den = np.sum((m[idx] / rho[idx]) * W)
        
        if den > 1e-12:
            new_rho[i] = num / den

    return new_rho  


def update_car_rigid_body(
        pos,
        vel,
        accel,
        m,
        types,
        dt,
        car_vel,
        car_omega,
        car_mass,
        car_I,
        g):
    """
    Mise à jour du corps rigide 'voiture' à partir
    des forces SPH accumulées dans accel.
    """

    # Particules appartenant à la voiture
    mask_car = (types == CAR)

    if not np.any(mask_car):
        return car_vel, car_omega

    # ----------------------------
    # Centre de masse
    # ----------------------------
    car_center = np.mean(pos[mask_car], axis=0)

    # ----------------------------
    # Forces totales
    # ----------------------------
    forces_indiv = m[mask_car, None] * (accel[mask_car] + g)
    total_force_car = np.sum(forces_indiv, axis=0)

    # ----------------------------
    # Moments
    # ----------------------------
    rel_pos = pos[mask_car] - car_center

    total_torque_car = np.sum(
        rel_pos[:, 0] * forces_indiv[:, 1]
        - rel_pos[:, 1] * forces_indiv[:, 0]
    )

    # ----------------------------
    # Accélérations rigides
    # ----------------------------
    acc_rigid = total_force_car / car_mass
    alpha_rigid = total_torque_car / car_I

    # ----------------------------
    # Intégration vitesses rigides
    # ----------------------------
    car_vel += acc_rigid * dt
    car_omega += alpha_rigid * dt

    # ----------------------------
    # Mise à jour du centre
    # ----------------------------
    car_center += car_vel * dt

    # ----------------------------
    # Rotation
    # ----------------------------
    d_theta = car_omega * dt

    rot_mat = np.array([
        [np.cos(d_theta), -np.sin(d_theta)],
        [np.sin(d_theta),  np.cos(d_theta)]
    ])

    new_rel_pos = rel_pos @ rot_mat.T

    # ----------------------------
    # Update positions
    # ----------------------------
    pos[mask_car] = car_center + new_rel_pos

    # ----------------------------
    # Update vitesses particules
    # ----------------------------
    vel[mask_car] = (
        car_vel
        + np.vstack([
            -car_omega * new_rel_pos[:, 1],
             car_omega * new_rel_pos[:, 0]
        ]).T
    )

    return car_vel, car_omega

def compute_derivatives(        
        pos, vel, rho, pres,
        vol, m,
        accel, drho,
        types,
        h, dt,
        shepard_coeff,
        wall_normals,
        v_mesh, tree, pairs,
        periodic=False,
        domain_size=1.0):  

    grad_rho = grad_p = grad_vel = None  # calculés uniquement si mode_p_refinement==1
    
    # Compute interactions (except for mesh analysis : testCase TGVmesh)
    if _ctrl.test_case != 'TGVmesh':
        for i, j in pairs:
            if types[i] == _sph.FLUID and types[j] == _sph.FLUID:

                compute_fluid_fluid_interaction(
                    i, j,
                    pos, vel, v_mesh,
                    rho, pres, vol, m,
                    shepard_coeff,
                    accel, drho,
                    h, _phys.rho0, _phys.c0, _phys.mu, visco_on,
                    grad_rho, grad_p, grad_vel,
                    shepard_coeff,periodic, domain_size
                )

            elif (types[i] == _sph.FLUID and types[j] == _sph.WALL) or \
                 (types[i] == _sph.WALL and types[j] == _sph.FLUID):

                compute_fluid_wall_interaction(
                    i, j,
                    pos, vel, types,
                    rho, pres, vol, m,
                    shepard_coeff,
                    wall_normals,
                    accel,
                    h, _phys.mu, _phys.gamma, visco_on,
                    _sph.FLUID,
                    shepard_coeff,
                    periodic,
                    grad_rho, grad_p, grad_vel
                )

            if _ctrl.test_case == 'car':

                if (types[i] == _sph.FLUID and types[j] == CAR) or \
                   (types[j] == _sph.FLUID and types[i] == CAR):

                    compute_fluid_wall_interaction(
                    i, j,
                    pos, vel, types,
                    rho, pres, vol, m,
                    shepard_coeff,
                    wall_normals,
                    accel,
                    h, _phys.mu, _phys.gamma, visco_on,
                    _sph.FLUID,
                    shepard_coeff,
                    periodic,
                    grad_rho, grad_p, grad_vel
                )

                # Interaction CAR - WALL (Collision simple)
                if (types[i] == CAR and types[j] == _sph.WALL) or \
                   (types[j] == CAR and types[i] == _sph.WALL):

                    idx_c, idx_w = (i, j) if types[i] == CAR else (j, i)

                    # collision
                    rij = pos[idx_c] - pos[idx_w]
                    r = np.linalg.norm(rij)

                    if r > 1e-12:
                        nij = rij / r
                        R_car = _sph.dx
                        R_wall = _sph.dx

                        penetration = (R_car + R_wall) - r

                        if penetration > 0:
                            vrel = vel[idx_c] - vel[idx_w]
                            vn = np.dot(vrel, nij)

                            force = car_spring_k * penetration * nij - car_spring_damping * vn * nij
                            accel[idx_c] += force / m[idx_c]

    # -------------------- chantier ----------------------------------------
    return drho, accel