import numpy as np
import constant.ale_properties as _ale_props
import constant.physical_properties as _phys
import constant.numerical_properties as _num_props
from src.kernels.sph_kernels import kernel_cubic, kernel_grad

def periodic_rij(xi,xj,domain_size):

    r_ij=xi-xj
    r_ij-=domain_size*np.round(r_ij/domain_size)

    return r_ij


def compute_v_mesh_full(pos, vel_f, vol_f, h, dt,  mode_ale,shepard, tree, pairs,periodic,domain_size):
    n_p = len(pos)

    # Calcul v_mesh selon le mode
    if mode_ale == 'lagrangian':
        #v_mesh = vel_f.copy()
        return vel_f.copy(), np.zeros(n_p), np.zeros(n_p), np.zeros((n_p,2)), np.zeros((n_p,2))
    if mode_ale == 'eulerian':
        return (np.zeros_like(vel_f),   # v_mesh : grille fixe
                np.zeros(n_p),           # S_G_local
                np.zeros(n_p),           # S_dist
                np.zeros((n_p, 2)),      # grad_S_G
                np.zeros((n_p, 2)))  
        
    if mode_ale == 'wass':
        v_mesh = vel_f.copy() 
        n_p = len(pos)
        s_g_local, s_g_mean = compute_geometric_entropy(pos, vol_f, h,shepard)        
        grad_s_g = np.zeros_like(pos)    
        error_pu = shepard - 1.0
        for i, j in pairs:
            if periodic:
                rij = periodic_rij(pos[i], pos[j],domain_size)
            else :
                rij = pos[i]-pos[j]
            r = np.linalg.norm(rij)
            if r < 1e-12: continue
            gW, _ = kernel_grad(rij, r, h)
            f_mag = 2.0 * (error_pu[i] * vol_f[j] + error_pu[j] * vol_f[i])
            grad_s_g[i] += f_mag * gW
            grad_s_g[j] -= f_mag * gW       
        grad_s_dist, s_dist = compute_distribution_entropy(pos, vol_f, h,shepard,tree,pairs,periodic,domain_size)    
    # Normalisation (éviter div0)
    #grad_norm = np.linalg.norm(grad_S_G + grad_S_dist, axis=1)
    #grad_norm[grad_norm < 1e-12] = 1.0
        beta_geom = _ale_props.beta_geom
        sub_mode_ale  = _ale_props.mode_ale
        if sub_mode_ale == 'grad': #wass or grad 
            grad_total = beta_geom * grad_s_g 
            norm = np.linalg.norm(grad_total, axis=1)
            max_norm = np.max(norm) + 1e-12
            grad_total /= max_norm
            v_shift =  -_phys.c0 * grad_total
           # momentum = np.sum(m[:,None]*v_shift,axis=0)/np.sum(m)
           #v_shift -= momentum
        elif sub_mode_ale == 'wass': #wass or grad 
            grad_s_dist = compute_wasserstein_shift(pos, vol_f, h, tree, pairs,periodic,domain_size)    
            # Calcul de la norme pour chaque particule
            # norm shape: (n,)
            norm_shift = np.linalg.norm(grad_s_dist, axis=1, keepdims=True)    
            # Direction unitaire (avec epsilon pour la stabilité)
            direction = grad_s_dist / (norm_shift + 1e-12)    
            # Vitesse de shift basée sur c0 
            max_v = np.max(np.linalg.norm(vel_f, axis=1)) if len(vel_f) > 0 else 0.0
            v_shift = -beta_geom * max_v * direction  
            # Application : on déplace la "grille" (mesh) vers le centroïde
        else:
            raise ValueError(
                f"sub_mode_ale inconnu : '{sub_mode_ale}'. "
        "Valeurs acceptées : 'grad', 'wass'."
    )
        v_mesh += v_shift
        
    return v_mesh, s_g_local, s_dist, grad_s_g, grad_s_dist
    

# --------------------------
# Entropie géométrique
# --------------------------
def compute_geometric_entropy(pos, vol, h,shepard,):          
    s_g_local = (shepard - 1)**2   
    return s_g_local, np.mean(s_g_local)
# --------------------------
# Entropie de distribution compute_distribution_entropy(pos, vol_f, h,shepard,tree,pairs)
# --------------------------
def compute_distribution_entropy(pos, vol, h, shepard, tree, pairs, periodic, domain_size):
    n_p             = len(pos)
    s_dist          = np.zeros(n_p)
    grad_s_dist     = np.zeros_like(pos)

    num_centroid    = np.zeros((n_p, 2))
    den_centroid    = np.zeros(n_p)

    for i, j in pairs:
        rij = periodic_rij(pos[i], pos[j], domain_size) if periodic else pos[i] - pos[j]
        r = np.linalg.norm(rij)
        if r < 1e-12 or r > 2 * h:
            continue
        W = kernel_cubic(r, h)
        # contribution symétrique pondérée par le volume
        num_centroid[i] += pos[j] * W * vol[j]
        den_centroid[i] += W * vol[j]
        num_centroid[j] += pos[i] * W * vol[i]
        den_centroid[j] += W * vol[i]

    mask = den_centroid > 1e-12
    centroid = np.where(mask[:, None], num_centroid / np.where(mask[:, None], den_centroid[:, None], 1.0), pos)

    grad_s_dist[mask] = pos[mask] - centroid[mask]          # direction Lloyd
    s_dist[mask]      = np.sum(grad_s_dist[mask]**2, axis=1)  # |x_i - c_i|²

    return grad_s_dist, s_dist

def compute_wasserstein_shift(pos, vol, h, tree, pairs,periodic,domain_size):   
    n = len(pos)
    # On initialise des accumulateurs pour le centre de gravité
    num_centroid = np.zeros((n, 2))
    den_centroid = np.zeros(n)    
    # On parcourt les paires pour accumuler les poids
    for i, j in pairs:
        if periodic:
            rij = periodic_rij(pos[i], pos[j],domain_size)
        else :
            rij = pos[i]-pos[j]
        r = np.linalg.norm(rij)
        
        if r < 1e-12 or r > 2*h: continue 
        
        W = kernel_cubic(r, h)
        weight = W * vol[j]
        
        # Accumulation pour la particule i
        num_centroid[i] += pos[j] * weight
        den_centroid[i] += weight
        
        # Accumulation pour la particule j (symétrie)
        # On utilise vol[i] pour la contribution de i sur j
        weight_rev = W * vol[i]
        num_centroid[j] += pos[i] * weight_rev
        den_centroid[j] += weight_rev

    # Calcul final du shift : (Centroïde - Position Actuelle)
    shift = np.zeros_like(pos)
    # On ne calcule que là où il y a des voisins pour éviter div par 0
    mask = den_centroid > 1e-12
    
    # Formule : Shift = (Somme(pos_j * W_ij * V_j) / Somme(W_ij * V_j)) - pos_i
    centroid = num_centroid[mask] / den_centroid[mask][:, None]
    shift[mask] = centroid - pos[mask]
    
    return shift
   

def compute_mls_gradients(pos, vol, field, h, tree,periodic,domain_size):
    """
    Calcule le gradient d'un champ (scalaire ou vecteur) via Moving Least Squares.
    field: array de taille (n,) ou (n, 2)
    """
    n_p = len(pos)
    dim = pos.shape[1]
    
    # Si field est scalaire (n,), on le reshape en (n,1) pour généraliser
    f_is_scalar = (field.ndim == 1)
    f = field[:, None] if f_is_scalar else field
    n_components = f.shape[1]
    
    gradients = np.zeros((n_p, n_components, dim))
    
    for i in range(n_p):       
            
        idx= tree.query_ball_point(pos[i], 2*h)
        if periodic:
            rij = periodic_rij(pos[i], pos[idx],domain_size)
        else :
            rij = pos[i]-pos[idx]
            
        r = np.linalg.norm(rij, axis=1)
        W = kernel_cubic(r, h)
        
        # Matrice de Moment M
        # M = sum( rij * rij^T * W * V )
        M = np.zeros((dim, dim))
        rhs = np.zeros((n_components, dim))
        
        for k in range(len(idx)):
            if r[k] < 1e-12: continue
            weight = W[k] * vol[idx[k]]
            
            # Produit extérieur pour la matrice de moment
            M += weight * np.outer(rij[k], rij[k])
            
            # Terme source pour le gradient
            df = f[idx[k]] - f[i]
            rhs += weight * np.outer(df, rij[k])
            
        # Résolution du système M * grad = rhs
        try:
            # On utilise pseudo-inverse pour la stabilité
            inv_M = np.linalg.pinv(M)
            gradients[i] = rhs @ inv_M
        except np.linalg.LinAlgError:
            gradients[i] = 0.0

    return gradients[:, 0, :] if f_is_scalar else gradients

def compute_energy_functionals(pos, vol, h, shepard, pairs,periodic,domain_size):
    """
    Calcule F1 (Wasserstein-2), F_S (Shepard defect), F2 = F1 + alpha*F_S.
    Retourne également le gradient complet pour la descente.
    """
    n = len(pos)
    
    # --- Terme F_S : défaut de Shepard ---
    F_S = 0.5 * np.sum(vol * (shepard - 1.0)**2)
    
    # --- Terme F_1 : Wasserstein-2 centroïdal ---
    F_1 = 0.0
    num_c = np.zeros((n, 2))  # numérateur centroïde
    den_c = np.zeros(n)       # = shepard (vérification)
    variance_dist = np.zeros(n)  # variance locale des |x_ij|^2

    for i, j in pairs:
        
        if periodic:
            rij = periodic_rij(pos[i], pos[j],domain_size)
        else :
            rij = pos[i]-pos[j]
            
        r = np.linalg.norm(rij)
        if r < 1e-12 or r > 2*h: continue
        
        W = kernel_cubic(r, h)
        r2 = r * r
        
        # Contributions centroïde pour i (utilisation de j)
        num_c[i] += pos[j] * W * vol[j]
        den_c[i] += W * vol[j]
        variance_dist[i] += W * vol[j] * r2
        
        # Contributions centroïde pour j (utilisation de i)
        num_c[j] += pos[i] * W * vol[i]
        den_c[j] += W * vol[i]
        variance_dist[j] += W * vol[i] * r2

    # Centroïdes et shift de domain_sizeloyd
    mask = den_c > 1e-12
    centroid = np.where(mask[:, None], num_c / np.where(mask[:, None], den_c[:, None], 1), pos)
    lloyd_shift = pos - centroid  # = x_i - c_i
    
    # F_1 = (1/2) * sum_i (variance_dist[i] / sigma_i)
    F_1 = 0.5 * np.sum(
        np.where(mask, variance_dist / np.where(mask, den_c, 1), 0.0)
    )
    
    return F_1, F_S, lloyd_shift, centroid, variance_dist

def compute_grad_shepard_corrected(pos, vol, h, shepard, pairs):
    """
    Gradient de F_S = (1/2) * sum_i V_i * (sigma_i - 1)^2
    CORRIGÉ : inclut le facteur V_i manquant.
    grad_FS[i] = V_i * sum_j V_j * [(sigma_i-1) + (sigma_j-1)] * grad_W_ij
    """
    n = len(pos)
    grad_FS = np.zeros_like(pos)
    error = shepard - 1.0  # sigma_i - 1
    
    for i, j in pairs:
        rij = pos[i] - pos[j]
        r = np.linalg.norm(rij)
        if r < 1e-12 or r > 2*h: continue
        
        gW, _ = kernel_grad(rij, r, h)
        
        # Facteur symétrique + facteur V_i * V_j (CORRECTION vs code actuel)
        coeff_i = vol[i] * vol[j] * (error[i] + error[j])
        
        grad_FS[i] += coeff_i * gW
        grad_FS[j] -= coeff_i * gW  # antisymétrie du grad W
    
    return grad_FS



def compute_grad_wasserstein_full(pos, vol, h, shepard, pairs, lloyd_shift, variance_dist):
    """
    Gradient de F_1 incluant :
      - Terme Lloyd : x_i - c_i
      - Correction d'ordre 2 (terme en nabla W * |x_ij|^2)
    """
    n = len(pos)
    grad_W2 = lloyd_shift.copy()  # terme dominant : x_i - c_i
    lloyd_sq = np.sum(lloyd_shift**2, axis=1)  # |x_i - c_i|^2
    
    mask_sigma = shepard > 1e-12
    inv_sigma = np.where(mask_sigma, 1.0 / np.where(mask_sigma, shepard, 1), 0.0)
    
    for i, j in pairs:
        rij = pos[i] - pos[j]
        r = np.linalg.norm(rij)
        if r < 1e-12 or r > 2*h: continue
        
        gW, _ = kernel_grad(rij, r, h)
        r2 = r * r
        
        # Correction ordre 2 pour i : (1/2*sigma_i) * V_j * nabla W_ij * (|x_ij|^2 - |c_i - x_i|^2)
        corr2_i = 0.5 * inv_sigma[i] * vol[j] * (r2 - lloyd_sq[i]) * gW
        grad_W2[i] += corr2_i
        
        # Correction ordre 2 pour j
        corr2_j = 0.5 * inv_sigma[j] * vol[i] * (r2 - lloyd_sq[j]) * (-gW)
        grad_W2[j] += corr2_j
    
    return grad_W2

def compute_v_mesh_coupled(pos, vel, vol, h, dt, shepard, tree, pairs,
                           alpha=1.0, beta_geom=0.001, mode_ale='coupled'):
    """
    Calcule la vitesse de grille ALE basée sur la descente de gradient de F2.
    
    Modes :
      'lloyd'   : descente F1 seul (terme dominant Lloyd)
      'wass'    : descente F1 complet (avec correction ordre 2)  
      'shepard' : descente F_S seul (correction Shepard)
      'coupled' : descente F2 = F1 + alpha*F_S  ← NOUVEAU
    """
    n = len(pos)
    v_mesh = vel.copy()
    
    # 1. Calcul des fonctionnelles et termes intermédiaires
    F1, FS, lloyd_shift, centroid, variance_dist = compute_energy_functionals(
        pos, vol, h, shepard, pairs
    )
    F2 = F1 + alpha * FS
    
    # 2. Gradient selon le mode
    if mode_ale in ('lloyd', 'wass', 'coupled'):
        grad_W2 = compute_grad_wasserstein_full(
            pos, vol, h, shepard, pairs, lloyd_shift, variance_dist
        )
    else:
        grad_W2 = np.zeros_like(pos)
    
    if mode_ale in ('shepard', 'coupled'):
        grad_FS = compute_grad_shepard_corrected(pos, vol, h, shepard, pairs)
    else:
        grad_FS = np.zeros_like(pos)
    
    # 3. Gradient total de F2
    grad_total = grad_W2 + alpha * grad_FS
    
    # 4. Normalisation locale (preserves spatial information)
    # ||x_i - c_i|| est bornée par 2h ; normaliser par h donne un shift adimensionnel
    norms = np.linalg.norm(grad_total, axis=1, keepdims=True)
    h_inv_norm = norms / h  # adimensionnel
    safe_norm = np.maximum(h_inv_norm, 1e-10)
    direction = grad_total / safe_norm  # unité : longueur
    
    # 5. Vitesse de shift proportionnelle à c0 et beta_geom
    v_shift = -beta_geom * _phys.c0 * (direction / h)  # dimension : m/s
    
    # 6. Application (seulement pour les particules bien supportées)
    mask_interior = shepard > _num_props.shepard_interior  # éviter les zones mal définies
    v_shift[~mask_interior] *= 0.0
    v_mesh += v_shift
    
    # 7. Conservation de la quantité de mouvement globale du shift
    # (optionnel : recentrage pour éviter la dérive)
    # momentum = np.sum(vol[:, None] * v_shift, axis=0) / np.sum(vol)
    # v_shift -= momentum
    
    return v_mesh, F1, FS, F2, grad_W2, grad_FS, lloyd_shift