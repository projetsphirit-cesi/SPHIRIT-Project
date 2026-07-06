import numpy as np
from src.mesh.primitives import smoothstep

def get_funnel_radius(y):
    y1 = 0.020
    y2 = 0.0321
    if y < y1:
        return 0.003
    elif y < y2:
        s = smoothstep(y1,y2,y)
        return (1-s)*0.003 + s*0.025
    else:
        return 0.025

def get_funnel_derivative(y):
    y1=0.020
    y2=0.0321
    if y<y1 or y>y2:
        return 0.0
    t=(y-y1)/(y2-y1)
    ds=6*t*(1-t)/(y2-y1)
    return (0.025-0.003)*ds

def generate_funnel_wall(dx, n_layers, H_funnel):
    pos_wall=[]
    s=0.0
    y=0.0
    while y<=H_funnel:
        r=get_funnel_radius(y)
        drdy=get_funnel_derivative(y)
        t=np.array([drdy,1.0])
        t/=np.linalg.norm(t)
        n=np.array([t[1],-t[0]])
        for layer in range(n_layers):
            offset=layer*dx
            xw=r+n[0]*offset
            yw=y+n[1]*offset
            pos_wall.append([ xw, yw ])
            pos_wall.append([ -xw, yw ])
        # avance réelle sur la surface
        ds=dx
        dy=ds/np.sqrt(1+drdy**2)
        y+=dy
    return np.array(pos_wall)



def compute_wall_normals(pos_wall):
    normals = np.zeros_like(pos_wall)
    for i,(x,y) in enumerate(pos_wall):
        # projection radiale vers la surface
        r_surface = get_funnel_radius(y)
        # dérivée vraie locale
        drdy = get_funnel_derivative(y)
        t = np.array([drdy,1.0])
        t /= np.linalg.norm(t)
        n = np.array([t[1],-t[0]])
        if x < 0:
            n[0] *= -1
        normals[i] = n
    return -normals