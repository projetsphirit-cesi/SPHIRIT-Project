import numpy as np

def apply_periodic_bc(pos, domain_size):
    '''Conditions périodiques — uniquement pour TGV.'''
    pos[:, 0] = pos[:, 0] % domain_size
    pos[:, 1] = pos[:, 1] % domain_size


def update_jet_buffers(
        test_case,
        pos, vel, types,
        rho, pres, vol, m,
        num_neighbors, shepard_coeff,
        bx, dx,
        y_inlet_top, y_inlet_bottom,
        v_jet, rho0,
        BUFFER_TOP, BUFFER_BOTTOM, FLUID,
        DOUBLE_JET=True,
        random_pos=False,
        add_position_noise=None
    ):
    """
    Gestion des buffers d'injection pour le test 'jet'.
    - conversion buffer -> fluide
    - régénération automatique des couches buffer
    """

    if test_case != 'jet':
        return pos, vel, types, rho, pres, vol, m, num_neighbors, shepard_coeff

    # =====================================================
    # 1. Conversion buffer -> fluide
    # =====================================================
    mask_top = (types == BUFFER_TOP)
    mask_bot = (types == BUFFER_BOTTOM)

    to_fluid_top = mask_top & (pos[:, 1] < y_inlet_top)
    to_fluid_bot = mask_bot & (pos[:, 1] > y_inlet_bottom)

    types[to_fluid_top | to_fluid_bot] = FLUID

    # =====================================================
    # 2. Régénération buffers
    # =====================================================
    to_add_pos = []
    to_add_vel = []
    to_add_types = []

    tol = 1e-7

    # ---------- TOP JET ----------
    if np.any(types == BUFFER_TOP):

        y_max_top = np.max(pos[types == BUFFER_TOP, 1])

        if y_max_top < (y_inlet_top + 3 * dx + tol):

            y_new = y_max_top + dx
            new_layer = np.vstack([bx, np.full_like(bx, y_new)]).T

            if random_pos and add_position_noise is not None:
                new_layer = add_position_noise(new_layer, 0.15 * dx)

            to_add_pos.append(new_layer)
            to_add_vel.append(np.tile([0, -v_jet], (len(bx), 1)))
            to_add_types.append(np.full(len(bx), BUFFER_TOP))

    # ---------- BOTTOM JET ----------
    if DOUBLE_JET and np.any(types == BUFFER_BOTTOM):

        y_min_bot = np.min(pos[types == BUFFER_BOTTOM, 1])

        if y_min_bot > (y_inlet_bottom - 3 * dx - tol):

            y_new = y_min_bot - dx
            new_layer = np.vstack([bx, np.full_like(bx, y_new)]).T

            if random_pos and add_position_noise is not None:
                new_layer = add_position_noise(new_layer, 0.15 * dx)

            to_add_pos.append(new_layer)
            to_add_vel.append(np.tile([0, v_jet], (len(bx), 1)))
            to_add_types.append(np.full(len(bx), BUFFER_BOTTOM))

    # =====================================================
    # 3. Ajout des nouvelles particules
    # =====================================================
    if to_add_pos:

        new_p = np.vstack(to_add_pos)
        new_v = np.vstack(to_add_vel)
        new_t = np.concatenate(to_add_types)

        n_add = len(new_p)

        pos = np.vstack([pos, new_p])
        vel = np.vstack([vel, new_v])
        types = np.concatenate([types, new_t])

        rho = np.concatenate([rho, np.full(n_add, rho0)])
        pres = np.concatenate([pres, np.zeros(n_add)])
        vol = np.concatenate([vol, np.full(n_add, dx**2)])
        m = np.concatenate([m, np.full(n_add, rho0 * dx**2)])

        num_neighbors = np.concatenate(
            [num_neighbors, np.zeros(n_add, dtype=int)]
        )

        shepard_coeff = np.concatenate(
            [shepard_coeff, np.zeros(n_add)]
        )

    return (
        pos, vel, types,
        rho, pres, vol, m,
        num_neighbors, shepard_coeff
    ) 