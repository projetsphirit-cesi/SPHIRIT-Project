import numpy as np

def fill_polygon(polygon,dx):

    polygon = np.array(polygon)

    xmin, ymin = polygon.min(axis=0)
    xmax, ymax = polygon.max(axis=0)

    pts = []

    xs = np.arange(xmin, xmax, dx)
    ys = np.arange(ymin, ymax, dx)

    for x in xs:
        for y in ys:
            if point_in_polygon(x, y, polygon):
                pts.append([x, y])

    return pts


# --------------------------------------------------
# AJOUT ROUE (CERCLE)
# --------------------------------------------------
def fill_circle(cx, cy, radius,dx):

    pts = []

    xs = np.arange(cx-radius, cx+radius, dx)
    ys = np.arange(cy-radius, cy+radius, dx)

    for x in xs:
        for y in ys:
            if (x-cx)**2 + (y-cy)**2 <= radius**2 :
                pts.append([x, y])

    return pts

def point_in_polygon(x, y, polygon):
    inside = False
    n = len(polygon)

    x1, y1 = polygon[0]
    for i in range(n + 1):
        x2, y2 = polygon[i % n]

        if y > min(y1, y2):
            if y <= max(y1, y2):
                if x <= max(x1, x2):
                    if y1 != y2:
                        xinters = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                    if x1 == x2 or x <= xinters:
                        inside = not inside

        x1, y1 = x2, y2

    return inside

def distance_to_segment(p, a, b):
    ap = p - a
    ab = b - a
    t = np.dot(ap, ab) / np.dot(ab, ab)
    t = np.clip(t, 0, 1)
    closest = a + t * ab
    return np.linalg.norm(p - closest), closest

def add_position_noise(points, amplitude):
    r = amplitude * np.sqrt(np.random.rand(len(points)))
    theta = 2*np.pi*np.random.rand(len(points))
    noise = np.column_stack((r*np.cos(theta), r*np.sin(theta)))
    return points + noise

def smoothstep(a,b,y):
    t=np.clip((y-a)/(b-a),0,1)
    return t*t*(3-2*t)