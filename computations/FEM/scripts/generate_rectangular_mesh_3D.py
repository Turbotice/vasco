#%%
# USER PARAMETERS (edit these lines as you wish)
# distances in cm :
plate_length_x = 8      # total extent in x-direction
plate_width_y = 4       # total extent in y-direction
plate_thickness_z = 0.4   # total extent in z-direction

nx_cells = 20              # number of elements along x
ny_cells = 10              # number of elements along y
nz_cells = 5              # number of elements along z

region_id = 1             # region id for hexahedra

# --- Script starts here ---
import numpy as np

# generate coordinate arrays
x_nodes = np.linspace(0.0, plate_length_x, nx_cells+1)
y_nodes = np.linspace(0.0, plate_width_y, ny_cells+1)
z_nodes = np.linspace(0.0, plate_thickness_z, nz_cells+1)

# build node list (x,y,z) ordered with x fastest, then y, then z
nodes = []
for kz, z in enumerate(z_nodes):
    for ky, y in enumerate(y_nodes):
        for kx, x in enumerate(x_nodes):
            nodes.append((x, y, z))
nodes = np.array(nodes)
num_nodes = nodes.shape[0]
print(f"Generated nodes: {num_nodes} ({len(x_nodes)} x {len(y_nodes)} y {len(z_nodes)} z)")

# helper to compute node index (1-based) given (ix,iy,iz)
def nid(ix,iy,iz):
    nxn = len(x_nodes)
    nyn = len(y_nodes)
    return iz*(nyn*nxn) + iy*nxn + ix + 1

# build hexahedral elements (8-node) with MEDIT ordering:
hexes = []
for iz in range(nz_cells):
    for iy in range(ny_cells):
        for ix in range(nx_cells):
            n0 = nid(ix  , iy  , iz)
            n1 = nid(ix+1, iy  , iz)
            n2 = nid(ix+1, iy+1, iz)
            n3 = nid(ix  , iy+1, iz)
            n4 = nid(ix  , iy  , iz+1)
            n5 = nid(ix+1, iy  , iz+1)
            n6 = nid(ix+1, iy+1, iz+1)
            n7 = nid(ix  , iy+1, iz+1)
            hexes.append((n0,n1,n2,n3,n4,n5,n6,n7, region_id))
num_hex = len(hexes)
print(f"Generated hexahedra: {num_hex}")

# write .mesh file
disk = 'D:'
out_path = f"{disk}/FEM/meshes/3d/3drect_plate.mesh"
with open(out_path, "w") as f:
    f.write("MeshVersionFormatted 2\n")
    f.write("Dimension 3\n\n")
    f.write("Vertices\n")
    f.write(f"{num_nodes}\n")
    for (x,y,z) in nodes:
        f.write(f"{x:.12e} {y:.12e} {z:.12e} 0\n")
    f.write("\n")
    f.write("Hexahedra\n")
    f.write(f"{num_hex}\n")
    for h in hexes:
        f.write(" ".join(str(int(v)) for v in h) + "\n")
    f.write("\nEnd\n")

print(f"Saved structured hexahedral mesh to: {out_path}")
print("\nPreview (first 12 vertices):")
for i,(x,y,z) in enumerate(nodes[:12],1):
    print(f"{i:3d}: {x:.3f}, {y:.3f}, {z:.3f}")
print("\nPreview (first 8 hexa elements):")
for i,h in enumerate(hexes[:8],1):
    print(f"{i:3d}: {h[:-1]}")
# Provide link for download
out_path

# %%
import pyvista as pv
import numpy as np

# assume `nodes` and `hexes` are already defined
points = nodes
cells = []
celltypes = []

for h in hexes:
    n = [idx-1 for idx in h[:-1]]  # convert to 0-based indices
    cells.append([8] + n)          # prepend number of points
    celltypes.append(pv.CellType.HEXAHEDRON)

cells = np.hstack(cells)
celltypes = np.array(celltypes)

grid = pv.UnstructuredGrid(cells, celltypes, points)

# --- Plot with axes ---
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color="lightblue", opacity=0.8)

# 1) Small interactive orientation widget in the corner
plotter.show_axes()

# 2) Global axes with tick values
plotter.show_bounds(
    grid='back',        # show a box around the geometry
    location='outer',   # put axes/ticks outside the box
    ticks='both',       # show ticks on all axes
    xlabel='X (length)',
    ylabel='Y (width)',
    zlabel='Z (thickness)'
)

plotter.show()

# %%
