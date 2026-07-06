import numpy as np
import re

def write_vtk(filename, pos, vel,v_mesh, pres, types, rho, vol, m, neighbors, shepard, wall_norms,
              s_g_local, s_dist,grad_s_g, grad_s_dist,
              # ── Nouvelles données pour le restart ──
              drho=None, accel=None, t=0.0, step=0, dt=0.0):
    
    n = len(pos)
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"SPHIRIT_RESTART t={t:.10e} step={step} dt={dt:.10e}\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        # Points
        f.write(f"POINTS {n} float\n")
        for p in pos:
            f.write(f"{p[0]:.6f} {p[1]:.6f} 0.0\n")

        # Vertices (obligatoire en POLYDATA)
        f.write(f"\nVERTICES {n} {2*n}\n")
        for i in range(n):
            f.write(f"1 {i}\n")

        # Données associées aux points
        f.write(f"\nPOINT_DATA {n}\n")

        f.write("SCALARS Pressure float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for p in pres:
            f.write(f"{float(p):.6f}\n")

        f.write("\nSCALARS Density float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for r in rho:
            f.write(f"{float(r):.6f}\n")

        f.write("\nSCALARS Volume float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for v in vol:
            f.write(f"{float(v):.6f}\n")

        f.write("\nSCALARS Type int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for tp in types:
            f.write(f"{int(tp)}\n")

        f.write("\nVECTORS Velocity float\n")
        for v in vel:
            f.write(f"{v[0]:.6f} {v[1]:.6f} 0.0\n")
            
        f.write("\nVECTORS VelocityMeshSPH float\n")
        for v in v_mesh:
            f.write(f"{v[0]:.6f} {v[1]:.6f} 0.0\n")
   
        f.write("\nSCALARS Neighbors int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for nb in neighbors:
            f.write(f"{int(nb)}\n")

        f.write("\nSCALARS Shepard float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for s in shepard:
            f.write(f"{float(s):.6f}\n")
        
        f.write("\nSCALARS Masse float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for s in m:
            f.write(f"{float(s):.6f}\n") 
            
        f.write("\nVECTORS WallNormals float\n")
        for vn in wall_norms:
            f.write(f"{vn[0]:.6f} {vn[1]:.6f} 0.0\n") 
            
        f.write("\nSCALARS S_G float 1\nLOOKUP_TABLE default\n")
        for s in s_g_local:
            f.write(f"{float(s):.6f}\n")
            
        f.write("\nSCALARS S_dist float 1\nLOOKUP_TABLE default\n")
        for i in range(n):
            f.write(f"{s_dist[i]:.6f}\n")
            
        f.write("\nVECTORS grad_S_G float\n")
        for v in grad_s_g:
            f.write(f"{v[0]:.6f} {v[1]:.6f} 0.0\n")
            
        f.write("\nVECTORS grad_S_dist float\n")
        for v in grad_s_dist:
            f.write(f"{v[0]:.6f} {v[1]:.6f} 0.0\n")

        # ── Champs restart (drho et accel) ─────────────────────────
        # drho : taux de variation de densité [kg/m³/s], nécessaire pour RK2
        if drho is not None:
            f.write("\nSCALARS drho float 1\nLOOKUP_TABLE default\n")
            for v in drho:
                f.write(f"{float(v):.10e}\n")

        # accel : accélération [m/s²], vecteur 2D, nécessaire pour RK2
        if accel is not None:
            f.write("\nVECTORS Acceleration float\n")
            for a in accel:
                f.write(f"{a[0]:.10e} {a[1]:.10e} 0.0\n")

def write_vtk_wf(filename, pos, vel, v_mesh, pres, types, rho, vol, m,
                     num_neighbors, shepard, wall_normals,
                     body_cm, body_theta, body_vel, body_omega):
        """VTK enrichi avec les champs spécifiques WaveFloat."""
        n = len(pos)
        with open(filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\nSPH_WAVEFLAT\nASCII\nDATASET POLYDATA\n")
            f.write(f"POINTS {n} float\n")
            for p in pos:
                f.write(f"{p[0]:.6f} {p[1]:.6f} 0.0\n")
            f.write(f"\nVERTICES {n} {2*n}\n")
            for i in range(n):
                f.write(f"1 {i}\n")
            f.write(f"\nPOINT_DATA {n}\n")
            # Champs scalaires
            for name, arr in [('Pressure', pres), ('Density', rho),
                               ('Volume', vol), ('Shepard', shepard), ('Mass', m)]:
                f.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n")
                for v in arr:
                    f.write(f"{float(v):.6f}\n")
            f.write("SCALARS Type int 1\nLOOKUP_TABLE default\n")
            for tp in types:
                f.write(f"{int(tp)}\n")
            f.write("SCALARS Neighbors int 1\nLOOKUP_TABLE default\n")
            for nb in num_neighbors:
                f.write(f"{int(nb)}\n")
            # Champs vectoriels
            for name, arr in [('Velocity', vel), ('VelocityMesh', v_mesh),
                               ('WallNormals', wall_normals)]:
                f.write(f"VECTORS {name} float\n")
                for v in arr:
                    f.write(f"{v[0]:.6f} {v[1]:.6f} 0.0\n")
                    

# ── Extraction des blocs de données ─────────────────────────────
def extract_scalar(lines, name,n , dtype=float):

    for idx, line in enumerate(lines):
        if line.startswith(f'SCALARS {name}'):
            # Sauter la ligne LOOKUP_TABLE
            data_start = idx + 2
            values = [dtype(lines[data_start + k].strip())
                        for k in range(n)]
            return np.array(values)
    return None

def extract_vector(lines, name,n):

    for idx, line in enumerate(lines):
        if line.startswith(f'VECTORS {name}'):
            data_start = idx + 1
            values = []
            for k in range(n):
                parts = lines[data_start + k].strip().split()
                values.append([float(parts[0]), float(parts[1])])
            return np.array(values)
    return None

def read_vtk(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    # ── Ligne 2 : extraction de t, step, dt ─────────────────────────
    header = lines[1].strip()
    t    = float(re.search(r't=([\d.eE+\-]+)',    header).group(1))
    step = int(  re.search(r'step=(\d+)',          header).group(1))
    dt   = float(re.search(r'dt=([\d.eE+\-]+)',   header).group(1))

    # ── Lecture du nombre de points ──────────────────────────────────
    n = None
    for line in lines:
        if line.startswith('POINTS'):
            n = int(line.split()[1])
            break
    if n is None:
        raise ValueError("Nombre de points introuvable dans le fichier VTK.")


    # ── Positions : bloc POINTS ──────────────────────────────────────
    pos = []
    for idx, line in enumerate(lines):
        if line.startswith('POINTS'):
            for k in range(n):
                parts = lines[idx + 1 + k].strip().split()
                pos.append([float(parts[0]), float(parts[1])])
            break
    pos = np.array(pos)

    # ── Champs physiques ─────────────────────────────────────────────
    pres         = extract_scalar(lines, 'Pressure',n)
    rho          = extract_scalar(lines, 'Density',n)
    vol          = extract_scalar(lines, 'Volume',n)
    types        = extract_scalar(lines, 'Type',n, dtype=int)
    m            = extract_scalar(lines, 'Masse',n)
    shepard      = extract_scalar(lines, 'Shepard',n)
    num_neighbors= extract_scalar(lines, 'Neighbors',n, dtype=int)
    drho         = extract_scalar(lines, 'drho',n)         # None si fichier ancien
    vel          = extract_vector(lines, 'Velocity',n)
    wall_normals = extract_vector(lines, 'WallNormals',n)
    accel        = extract_vector(lines, 'Acceleration',n)  # None si fichier ancien

    # ── Valeurs par défaut si champs restart absents ─────────────────
    # (ex: fichiers VTK générés avant l'ajout de drho et accel)
    if drho is None:
        drho = np.zeros(n)
    if accel is None:
        accel = np.zeros((n, 2))

    return {
        'pos':          pos,
        'vel':          vel,
        'rho':          rho,
        'pres':         pres,
        'vol':          vol,
        'm':            m,
        'types':        types,
        'wall_normals': wall_normals,
        'num_neighbors':num_neighbors if num_neighbors is not None else np.zeros(n, dtype=int),
        'shepard_coeff':shepard      if shepard       is not None else np.zeros(n),
        'v_mesh':       np.zeros((n, 2)),
        'accel':        accel,
        'drho':         drho,
        't':            t,
        'step':         step,
        'dt':           dt,
    }