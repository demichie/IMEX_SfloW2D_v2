----------------------------------------------------------------------------

CLAUDIA ELIJAS-PARRA - regularised mu(I) & dynamic permeability & compaction effect

Simulation example that accounts for the effect of compaction and diffusion occurring simultaneously in thin fluidised granular flows. The conditions are defined in the following namelists:

&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL = 11,
 MU_S = 0.55D0,
 MU_2 = 0.73D0, 
 MUI_INF = 0.05D0,
 I_0 = 0.279D0,
 /

These parameters correspond to the regularised mu(I) rheology (Barker & Gray, 2017). If MUI_INF = 0, the implemented mu(I) will correspond to the unregularized mu(I).

&PORE_PRESSURE_PARAMETERS
 hydraulic_permeability = 1.10D-11,
 GAS_LOSS_FLAG = T,
 ALPHA_TRANS = 0.2,
 F_INHIBIT_MODE = 'DYNAMIC',
 DYNAMIC_PERMEABILITY_FLAG = T,
 /

- F_INHIBIT_MODE defines the way we describe the transition from completely unhindered to fully inhibited gas escape as the solid fraction (alpha) approaches the maximum solid packing (alpha_max):
 'OFF' (always = 1, even if alpha>alpha_max),
 'STEP' (step function that switches at alpha_max),
 'STATIC' (a single smooth function independent of height, for which you provide ALPHA_TRANS, the value of alpha at which the transition from 1 to 0 starts),
 'DYNAMIC' (smooth function that depends on height so you DON'T have to provide ALPHA_TRANS). 
- GAS_LOSS_FLAG must be T, otherwise you'll lose pressure but not gas, and it will yield unphysical results. 
- if DYNAMIC_PERMEABILITY_FLAG = T, the provided  hydraulic_permeability  = 5.85e-12 will be ignored and permeability is calculated dynamically at each time step and cell depending on the alphas and Sauter diameter. 

It uses the set up of a dam break experiment, explained below. 

----------------------------------------------------------------------------

ECP BREARD — Dam-break example (channel)

This sets up a 2-D dam-break in a long rectangular channel with no-flux side walls (lateral) and no-flux at the right boundary, runs IMEX to produce a NetCDF, and post-processes the centreline flow height to PNGs and a movie (equal x–y scale, grey fill under the surface, blue outline).

Scripts used: create_dambreak_static.py (setup + IMEX_SfloW2D.inp), height_summary_nc.py (centreline height → PNGs + MP4/GIF), clean_folder.py (optional tidy), IMEX_SfloW2D.template (input template).

Setup script arguments (create_dambreak.py): 
NX NY X_LENGTH Y_WIDTH BED_LENGTH BED_HEIGHT ALFAS T PORE_PRESSURE_FRACT PLOT_FLAG

Where:
- NX, NY = grid cells in x and y directions
- X_LENGTH (m) = channel length (total x-length)
- Y_WIDTH (m) = channel width (y-width)
- BED_LENGTH (m) = initial dam thickness in x
- BED_HEIGHT (m) = initial dam height in z
- ALFAS = initial solid volume fraction (0-1)
- T (K) = ambient temperature
- PORE_PRESSURE_FRACT = fraction of hydrostatic pore pressure (0-1)
- PLOT_FLAG = true or false for PNG output

----------------------------------------------------------------------------

USAGE EXAMPLE:

1. CREATE DAM BREAK RESTART FILE AND INPUT FILE

>> python create_dambreak.py 300 15 3.0 0.15 0.20 0.40 0.58 300 1.0 true

Writes: IMEX_SfloW2D.inp, topography_dem.asc, example_2D_0000.q_2d (and a preview PNG).


2. LINK SOLVER TO THIS FOLDER (once)

>> ln -s ../../bin/IMEX_SfloW2D .


3. RUN THE SOLVER WITH THE INPUT FILE PRODUCED BY THE SETUP

>> ./IMEX_SfloW2D IMEX_SfloW2D.inp

Produces: dambreak2D.nc (plus snapshots).


4. CONVERT TO NC FILE

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

>> python p2d_to_netCDF4.py dambreak2D.bak


5. FIND THE CENTERLINE HEIGHT ALONG Y = center_fraction (0=bottom, 0.5=mid, 1=top). ADD --animate FOR A MOVIE

>> python height_summary_nc.py dambreak2D.nc --center_fraction 0.5 --animate

Outputs: frames/dambreak2D_hcenter_####.png and dambreak2D_height_movie.mp4 (or GIF if ffmpeg not available).

NOTE (MP4 support): if needed install ffmpeg bindings and system ffmpeg
>> pip install -U "imageio[ffmpeg]" imageio-ffmpeg
>> sudo apt-get update && sudo apt-get install -y ffmpeg


6. CLEAN FOLDER (dry-run first; then actually delete non-essentials)

>> python clean_folder.py

>> python clean_folder.py --apply

----------------------------------------------------------------------------
