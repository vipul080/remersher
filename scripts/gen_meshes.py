#!/usr/bin/env python3
"""Generate deterministic OBJ test meshes for the ground-truth corpus."""
import math, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "meshes")

def write_obj(path, verts, faces):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for v in verts:
            f.write("v %.6f %.6f %.6f\n" % v)
        for fa in faces:
            f.write("f " + " ".join(str(i + 1) for i in fa) + "\n")
    return path

def uv_sphere(seg_u, seg_v, radius=1.0, displ=None, center=(0,0,0)):
    verts, faces = [], []
    for iv in range(seg_v + 1):
        v = iv / seg_v
        theta = v * math.pi
        for iu in range(seg_u + 1):
            u = iu / seg_u
            phi = u * 2 * math.pi
            r = radius
            if displ:
                r *= displ(u, v)
            x = r * math.sin(theta) * math.cos(phi)
            y = r * math.cos(theta)
            z = r * math.sin(theta) * math.sin(phi)
            verts.append((x + center[0], y + center[1], z + center[2]))
    for iv in range(seg_v):
        for iu in range(seg_u):
            a = iv * (seg_u + 1) + iu
            b = a + 1
            c = a + seg_u + 1
            d = c + 1
            if iv < seg_v - 1:
                faces.append((a, b, d))
                faces.append((a, d, c))
            else:  # pole fan
                faces.append((a, d, b))
    return verts, faces

def torus(R=1.0, r=0.4, su=32, sv=16):
    verts, faces = [], []
    for iv in range(sv):
        v = iv / sv * 2 * math.pi
        for iu in range(su):
            u = iu / su * 2 * math.pi
            x = (R + r * math.cos(v)) * math.cos(u)
            y = r * math.sin(v)
            z = (R + r * math.cos(v)) * math.sin(u)
            verts.append((x, y, z))
    for iv in range(sv):
        for iu in range(su):
            a = iv * su + iu
            b = a + 1 if iu < su - 1 else iv * su
            c = (iv + 1) % sv * su + iu
            d = (iv + 1) % sv * su + (iu + 1 if iu < su - 1 else 0)
            faces.append((a, b, d, c))  # quads as input
    return verts, faces

def cube():
    # 8 corners, 6 quad faces
    c = 1.0
    verts = [(-c,-c,-c),(c,-c,-c),(c,c,-c),(-c,c,-c),(-c,-c,c),(c,-c,c),(c,c,c),(-c,c,c)]
    faces = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(2,3,7,6),(1,2,6,5),(3,0,4,7)]
    return verts, faces

def ngon_cube():
    # cube where top face is a 5-gon (extra vertex mid-edge)
    v, f = cube()
    # split edge (5,6) -> new vertex m at midpoint
    m = ((v[5][0]+v[6][0])/2, (v[5][1]+v[6][1])/2, (v[5][2]+v[6][2])/2)
    v.append(m); mi = len(v) - 1
    # top face (4,5,6,7) -> (4,5,mi,6,7) 5-gon
    f[1] = (4, 5, mi, 6, 7)
    return v, f

def stair():
    # L-profile extrusion: sharp edges, flat patches
    prof = [(0,0),(1,0),(1,0.5),(0.5,0.5),(0.5,1),(0,1)]
    n = len(prof)
    depth = 1.0
    verts, faces = [], []
    for z in (0.0, depth):
        for (x, y) in prof:
            verts.append((x, y, z))
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, j + n, i + n))       # side walls
    faces.append(tuple(range(n)))                 # bottom
    faces.append(tuple(range(n, 2 * n)))          # top
    return verts, faces

def plane_with_hole(grid=8, hole=2):
    # flat plane with a square hole: open borders (outer + inner)
    verts, faces = [], []
    for j in range(grid + 1):
        for i in range(grid + 1):
            verts.append((i, j, 0.0))
    def vid(i, j):
        return j * (grid + 1) + i
    h0 = (grid - hole) // 2
    h1 = h0 + hole
    for j in range(grid):
        for i in range(grid):
            if h0 <= i < h1 and h0 <= j < h1:
                continue  # hole
            a = vid(i, j); b = vid(i + 1, j); c = vid(i + 1, j + 1); d = vid(i, j + 1)
            faces.append((a, b, c, d))
    return verts, faces

def gear(teeth=24, radius=1.0, inner=0.45, depth=0.3):
    # extruded gear polygon with center hole: hard surface + borders
    n = teeth * 2
    verts, faces = [], []
    pts = []
    for i in range(n):
        ang = i / n * 2 * math.pi
        r = radius if i % 2 == 0 else radius * 0.8
        pts.append((r * math.cos(ang), r * math.sin(ang)))
    for z in (0.0, depth):
        for (x, y) in pts:
            verts.append((x, y, z))
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, j + n, i + n))
    faces.append(tuple(range(n)))         # top ngon (24-gon input)
    faces.append(tuple(range(n, 2 * n)))
    return verts, faces

def bumpy_sphere():
    def displ(u, v):
        return 1.0 + 0.12 * math.sin(5 * u * 2 * math.pi) * math.cos(3 * v * math.pi) \
                   + 0.07 * math.sin(11 * u * 2 * math.pi + 1.7) * math.sin(7 * v * math.pi)
    return uv_sphere(64, 32, 1.0, displ=displ)

def two_spheres():
    v1, f1 = uv_sphere(32, 16, 1.0, center=(-0.8, 0, 0))
    v2, f2 = uv_sphere(32, 16, 1.0, center=(0.8, 0, 0))
    off = len(v1)
    return v1 + v2, f1 + [(a + off, b + off, c + off) for (a, b, c) in f2]

def main():
    meshes = {
        "sphere_4160": lambda: uv_sphere(64, 32),
        "sphere_32768": lambda: uv_sphere(256, 128),
        "cube": cube,
        "ngon_cube": ngon_cube,
        "stair": stair,
        "torus": lambda: torus(1.0, 0.4, 32, 16),
        "plane_hole": plane_with_hole,
        "gear": gear,
        "bumpy": bumpy_sphere,
        "twospheres": two_spheres,
    }
    for name, gen in meshes.items():
        verts, faces = gen()
        path = write_obj(os.path.join(OUT, name + ".obj"), verts, faces)
        tris = sum(max(0, len(f) - 2) for f in faces)
        print(f"{name}: {len(verts)}v {len(faces)}f ({tris} tris)")

if __name__ == "__main__":
    main()
