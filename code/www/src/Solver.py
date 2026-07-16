import numpy as np
import os
import constant.ale_properties as _ale_props
from scipy.spatial import cKDTree
from src.kernels.sph_kernels import kernel_cubic
from src.physics.sph_interactions import compute_shepard_coefficients, compute_derivatives
from src.physics.get_density import get_pressure
from src.discretization.ale_mesh import compute_v_mesh_full
from src.time_integration.euler import euler_step
from src.time_integration.runge_kutta2 import rk2_step
from io_.vtk import write_vtk


def run(state,ctrl,phys,sph_params):
    # Lire les variables depuis les dicts
    pos=state['pos']; vel=state['vel']; rho=state['rho']
    pres=state['pres']; vol=state['vol']; m=state['m']
    types=state['types']; wall_normals=state['wall_normals']
    t=state['t']; step=state['step']; dt=state['dt']
    h=sph_params['h']; CFL=sph_params['CFL']
    B=phys['B']; rho0=phys['rho0']; c0=phys['c0']
    gamma=phys['gamma']; g=phys['g']; mu=phys['mu']
    FLUID=sph_params['FLUID']; WALL=sph_params['WALL']
    finaltime=ctrl['finaltime']; save=ctrl['save']
    vtk_dir=ctrl['vtk_dir']; mode_rk=ctrl.get('mode_rk',1)
    mesh_mode = _ale_props.mesh_mode
    os.makedirs(vtk_dir, exist_ok=True)

    while t < finaltime:
        n_current = len(pos) # La taille peut changer à chaque itération    


        periodic = phys.get('periodic', False)
        domain_size = phys.get('domain_size', 1.0)
        tree = cKDTree(pos, boxsize=domain_size) if periodic else cKDTree(pos)
        pairs = list(tree.query_pairs(2*h))
        mask_f = (types == FLUID)

        v_mesh=np.zeros((n_current,2)); accel=np.zeros((n_current,2))
        drho=np.zeros(n_current); shepard=np.zeros(n_current); num_nb=np.zeros(n_current,dtype=int)
        shepard[:] = kernel_cubic(0., h) * vol[:]  # self-contribution
        compute_shepard_coefficients(pairs, pos, vol, h, shepard, num_nb, kernel_cubic, periodic,domain_size)

        # Notebook lignes 2790-2791 : vitesse de maillage
        v_mesh,s_g,s_d,grad_s_g,grad_s_d = compute_v_mesh_full(pos, vel, vol, h, dt, mesh_mode, shepard, tree, pairs, periodic,domain_size)

        compute_derivatives(pos, vel, rho, pres, vol, m, accel, drho,
                    types, h, dt, shepard, wall_normals, v_mesh, tree, pairs,periodic,domain_size)

        if mode_rk == 1:
            euler_step(pos, vel, rho, pres, vol, m, accel, drho,
                       mask_f, dt, g, B, rho0, gamma, shepard, v_mesh, periodic, domain_size)
        elif mode_rk == 2:
            rk2_step(pos, vel, rho, pres, vol, m,
                     types, wall_normals, shepard,
                     tree, pairs, mask_f, dt, g,
                     B, rho0, gamma, h, periodic, domain_size)

        if np.any(mask_f):
            max_v = np.max(np.linalg.norm(vel[mask_f],axis=1))
            dt = CFL*h / (c0 + max_v + 1e-12)

        if step % save == 0:
            fname = os.path.join(vtk_dir, f'test_case{step:05d}.vtk')
            write_vtk(fname, pos, vel, v_mesh-vel, pres, types, rho,
                      vol, m, num_nb, shepard, wall_normals,
                      s_g, s_d, grad_s_g, grad_s_d,
                      drho=drho, accel=accel, t=t, step=step, dt=dt)
            print(f'Step {step} | t={t:.4f}s | dt={dt:.2e} | N={int(np.sum(mask_f))}')

        t += dt; step += 1
    state['t']    = t
    state['step'] = step
    state['dt']   = dt