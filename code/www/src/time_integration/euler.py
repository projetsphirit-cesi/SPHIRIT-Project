import numpy as np
import constant.ale_properties as _ale_props
import constant.numerical_properties as _num_props
import system.control_dict as _ctrl
from src.kernels.sph_kernels import get_pressure

def euler_step(pos, vel, rho, pres, vol, m, accel, drho,
               mask_f, dt, g,
               B, rho0, gamma,
               shepard_coeff, v_mesh,
               periodic=False, domain_size=1.0):
        mesh_mode = _ale_props.mesh_mode
        if not np.any(mask_f): return   
        rho[mask_f] += drho[mask_f] * dt
        vol[mask_f] = m[mask_f] / rho[mask_f]
        vel[mask_f] += (accel[mask_f] + g) * dt
        pres[mask_f] = get_pressure(rho[mask_f],B,rho0,gamma)  
        if mesh_mode == 'wass':
            pos[mask_f] += np.where(
                shepard_coeff[mask_f][:, None] > _num_props.shepard_threshold,
                v_mesh[mask_f],
                vel[mask_f]
                ) * dt
        else :
            pos[mask_f] += v_mesh[mask_f]* dt


        if periodic:
            if _ctrl.test_case == 'TGV':
                pos[:, 1] = pos[:, 1] % domain_size  #only TGV testcase
            else:
                pos[:, 0] = pos[:, 0] % domain_size