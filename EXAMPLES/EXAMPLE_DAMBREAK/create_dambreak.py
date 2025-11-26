#!/usr/bin/env python3
"""
create_dambreak_static.py

Create an IMEX/SfloW2D input and restart for a simple 2‑D dam‑break problem.

This script sets up a rectangular channel domain with a block of material
initially occupying part of the channel.  The block represents a vertical
dam of thickness ``bed_length`` and height ``bed_height`` at the left
end of the channel.  The remainder of the channel is initially empty.

Compared with ``create_example_Roche_static.py`` this script

* Accepts channel length and width explicitly rather than min/max
  coordinates.  The domain always starts at ``x=0`` and ``y=0``.
* Places a rectangular bed (dam) spanning the full channel width from
  ``x=0`` to ``x=bed_length``.  Its height above the bed is ``bed_height``.
* Writes ``topography_dem.asc`` (a flat DEM), ``example_2D_0000.q_2d``
  (restart file) and ``IMEX_SfloW2D.inp`` (input file) in the working
  directory.
* Optionally produces a PNG of the initial bed for quick inspection.

Usage::

    python create_dambreak_static.py NX NY X_LENGTH Y_WIDTH BED_LENGTH BED_HEIGHT ALFAS T PORE_PRESSURE_FRACT PLOT_FLAG

Arguments:

``NX``, ``NY``
    Integer grid resolution in the x and y directions.
``X_LENGTH``
    Length of the computational domain in metres (channel length).
``Y_WIDTH``
    Width of the channel in metres.
``BED_LENGTH``
    Length of the initial dam (bed) at the left end of the channel in metres.
``BED_HEIGHT``
    Height of the dam above the bed in metres.
``ALFAS``
    Initial solid volume fraction (0—1).  The mixture density will
    be computed from this value, assuming a solid density of 2500 kg/m³
    and an ambient air density derived from ``PRES`` and ``T``.
``T``
    Ambient temperature in Kelvin (for computing air density).
``PORE_PRESSURE_FRACT``
    Fraction of hydrostatic pore pressure to initialize (0—1).
``PLOT_FLAG``
    Either ``true`` or ``false``.  If ``true``, a PNG named
    ``initial_bed.png`` will be created showing the initial bed shape.

Example::

    # A 3 m long channel, 0.15 m wide, with a 0.2 m long, 0.4 m high dam at the left end.
    python create_dambreak_static.py 300 15 3.0 0.15 0.2 0.4 0.6 300 0.5 true

The script relies on an ``IMEX_SfloW2D.template`` file in the current
directory.  It replaces placeholders ``runname``, ``restartfile``,
``x_min``, ``y_min``, ``nx_cells``, ``ny_cells`` and ``dx`` as in the
original example script.  The ``runname`` is set to ``dambreak2D``.

Author: ECP_BREARD
"""

import sys
import numpy as np
from pathlib import Path

def usage():
    print("Usage: python create_dambreak_static.py NX NY X_LENGTH Y_WIDTH "
          "BED_LENGTH BED_HEIGHT ALFAS T PORE_PRESSURE_FRACT PLOT_FLAG")
    sys.exit(1)

if len(sys.argv) != 11:
    usage()

def as_int(x): return int(x)
def as_float(x): return float(x)

nx_cells = as_int(sys.argv[1])
ny_cells = as_int(sys.argv[2])
x_length = as_float(sys.argv[3])
y_width = as_float(sys.argv[4])
bed_length = as_float(sys.argv[5])
bed_height = as_float(sys.argv[6])
alfas = as_float(sys.argv[7])
T = as_float(sys.argv[8])
pore_pressure_fract = as_float(sys.argv[9])
plot_flag = sys.argv[10].lower() == 'true'

# Physical constants (consistent with the Roche example)
n_solid = 1
rho_s = 2500.0              # kg/m^3, solid density
SP_HEAT_S = 1617.0
SP_HEAT_A = 998.0
SP_GAS_CONST_A = 287.051
PRES = 101300.0
grav = 9.81

# Derived mixture properties
rho_a = PRES / (T * SP_GAS_CONST_A)          # ambient air density
rho_m = alfas * rho_s + (1.0 - alfas) * rho_a
xs = alfas * rho_s / rho_m
SP_HEAT_MIX = xs * SP_HEAT_S + (1.0 - xs) * SP_HEAT_A

# Reduced gravity
red_grav = grav * (rho_m - rho_a) / rho_m
print(f"density difference = {rho_m - rho_a:.2f} kg/m³")
print(f"reduction factor = {(rho_m - rho_a)/rho_m:.4f}")
print(f"reduced gravity = {red_grav:.4f} m/s²")

# Domain extents (always start at x=0, y=0)
x_min = 0.0
y_min = 0.0
x_max = x_length
y_max = y_width

# Compute grid spacing
dx = (x_max - x_min) / float(nx_cells)
dy = (y_max - y_min) / float(ny_cells)

# Cell centres
nx_points = nx_cells + 1
ny_points = ny_cells + 1
x_cent = np.linspace(x_min + 0.5 * dx, x_max - 0.5 * dx, nx_cells)
y_cent = np.linspace(y_min + 0.5 * dy, y_max - 0.5 * dy, ny_cells)
Xc, Yc = np.meshgrid(x_cent, y_cent)

# Initialize bed (dam) height array
H = np.where((Xc >= x_min) & (Xc <= bed_length), bed_height, 0.0)

# Initial velocities (dam at rest)
U = np.zeros_like(H)
V = np.zeros_like(H)

# Write flat topography (DEM).  We use dx as cellsize so the header is consistent
Z = np.zeros((ny_points, nx_points))
header = (
    f"ncols     {nx_points}\n"
    f"nrows    {ny_points}\n"
    f"xllcorner {x_min - 0.5 * dx}\n"
    f"yllcorner {y_min - 0.5 * dx}\n"
    f"cellsize {dx}\n"
    "NODATA_value -9999\n"
)
with open('topography_dem.asc', 'w') as f:
    np.savetxt(f, Z, header=header, fmt='%1.12f', comments='')

# Write restart file (example_2D_0000.q_2d)
init_file = 'example_2D_0000.q_2d'
for j in range(ny_cells):
    q0 = np.zeros((7 + n_solid, nx_cells))
    # x coordinate
    q0[0, :] = x_cent
    # y coordinate
    q0[1, :] = y_cent[j]
    print(f"Row {j+1}/{ny_cells}, y={y_cent[j]:.4f} m, Column 1,x={x_cent[j]:.4f} m, H={np.max(H[j,0]):.4f} m")
    # mass per cell: rho_m * H
    q0[2, :] = rho_m * H[j, :]
    # x momentum: rho_m * H * U
    q0[3, :] = rho_m * H[j, :] * U[j, :]
    # y momentum: rho_m * H * V
    q0[4, :] = rho_m * H[j, :] * V[j, :]
    # internal energy: rho_m * H * cp_mix * T
    q0[5, :] = rho_m * H[j, :] * SP_HEAT_MIX * T
    # solid mass fraction: rho_s * H * alfas / n_solid
    q0[6, :] = rho_s * H[j, :] * alfas / n_solid
    # excess pore pressure: rho_m * g * H (Pa)
    P_excess = rho_m * red_grav * H[j, :]
    #print(f"pore pressure: {(P_excess+PRES)/1000} Pa")
    print(f"excess pore pressure: {(P_excess[0])/1000} kPa")
    # conservative variable for pore pressure: (P_excess) * (rho_m * H) * pore_pressure_fract
    q0[-1, :] = P_excess * q0[2, :] * pore_pressure_fract
    # Append to file
    mode = 'w+' if j == 0 else 'a'
    with open(init_file, mode) as f:
        np.savetxt(f, q0.T, fmt='%19.12e')
    with open(init_file, 'a') as f:
        f.write(' \n')

# Prepare IMEX input by replacing placeholders in the template
tmpl_path = Path('IMEX_SfloW2D.template')
if not tmpl_path.exists():
    sys.exit('Error: IMEX_SfloW2D.template not found in current directory.')
tmpl = tmpl_path.read_text(encoding='utf-8', errors='ignore')
filedata = tmpl
filedata = filedata.replace('runname', 'dambreak2D')
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_min', f"{x_min}")
filedata = filedata.replace('y_min', f"{y_min}")
filedata = filedata.replace('nx_cells', str(nx_cells))
filedata = filedata.replace('ny_cells', str(ny_cells))
filedata = filedata.replace('dx', f"{dx}")
Path('IMEX_SfloW2D.inp').write_text(filedata, encoding='utf-8')

# Print basic info
print(f"Created dam-break setup with NX={nx_cells}, NY={ny_cells}")
print(f"Domain length {x_length} m, width {y_width} m")
print(f"Dam length {bed_length} m, height {bed_height} m")
print(f"Pore pressure fraction: {pore_pressure_fract}")
print(f"Saved topography_dem.asc, {init_file}, IMEX_SfloW2D.inp")

# Optional plot of the initial bed
if plot_flag:
    import matplotlib.pyplot as plt
    Xp, Yp = np.meshgrid(x_cent, y_cent)
    plt.figure(figsize=(6, 3))
    pcm = plt.pcolormesh(Xp, Yp, H, shading='auto', cmap='viridis')
    plt.colorbar(pcm, label='H (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Initial dam height distribution')
    plt.tight_layout()
    plt.savefig('initial_bed.png', dpi=150)
    plt.close()