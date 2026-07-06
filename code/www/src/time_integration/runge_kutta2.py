import numpy as np
import constant.ale_properties as _ale_props
from src.physics.sph_interactions import compute_derivatives
from src.discretization.ale_mesh import compute_v_mesh_full
from src.kernels.sph_kernels import get_pressure


def rk2_step(pos, vel, rho, pres, vol, m,
             types, wall_normals, shepard_coeff,
             tree, pairs, mask_f, dt, g,
             B, rho0, gamma, h,
             periodic=False, domain_size=1.0):
    n = len(pos)
    mesh_mode = _ale_props.mesh_mode

    # ── Prédicteur ──────────────────────────────────────────
    accel1 = np.zeros((n, 2)); drho1 = np.zeros(n)
    vmesh1, _, _, _, _ = compute_v_mesh_full(
        pos, vel, vol, h, dt, mesh_mode,
        shepard_coeff, tree, pairs, periodic, domain_size)
    compute_derivatives(pos, vel, rho, pres, vol, m, accel1, drho1,
                         types, h, dt, shepard_coeff, wall_normals,
                         vmesh1, tree, pairs, periodic, domain_size)

    rho_star = rho.copy(); vel_star = vel.copy(); pos_star = pos.copy()
    vol_star = vol.copy()
    rho_star[mask_f] += drho1[mask_f] * dt
    vel_star[mask_f] += (accel1[mask_f] + g) * dt
    pos_star[mask_f] += vmesh1[mask_f] * dt
    vol_star[mask_f] = m[mask_f] / rho_star[mask_f]
    pres_star = pres.copy()
    pres_star[mask_f] = get_pressure(rho_star[mask_f], B, rho0, gamma)

    # ── Correcteur ──────────────────────────────────────────
    accel2 = np.zeros((n, 2)); drho2 = np.zeros(n)
    vmesh2, _, _, _, _ = compute_v_mesh_full(
        pos_star, vel_star, vol_star, h, dt, mesh_mode,
        shepard_coeff, tree, pairs, periodic, domain_size)
    compute_derivatives(pos_star, vel_star, rho_star, pres_star, vol_star, m,
                         accel2, drho2, types, h, dt, shepard_coeff,
                         wall_normals, vmesh2, tree, pairs, periodic, domain_size)

    # ── Moyenne des deux pentes ─────────────────────────────
    rho[mask_f] += 0.5 * (drho1[mask_f] + drho2[mask_f]) * dt
    vol[mask_f] = m[mask_f] / rho[mask_f]
    vel[mask_f] += 0.5 * (accel1[mask_f] + accel2[mask_f]) * dt + g * dt
    pres[mask_f] = get_pressure(rho[mask_f], B, rho0, gamma)
    pos[mask_f] += 0.5 * (vmesh1[mask_f] + vmesh2[mask_f]) * dt

    if periodic:
        pos[:, 0] = pos[:, 0] % domain_size
        pos[:, 1] = pos[:, 1] % domain_size