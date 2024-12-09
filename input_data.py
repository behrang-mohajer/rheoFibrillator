######## Use this file to input the attached .py code.	########
######## There is no need to modify any other library.  ########
######## Run either main_loop, Allrun, or the jupyter.  ########

# Address to .foam file:
filename			        = 	'GitHub_Die_swell.foam'
directory                   =   './'    #it directory if read NOT from the current address

# List of the fields that must be read from teh openFoam Case:
stress_varibale             =   'tauMF' # The one that corresponds to the stress tensor in the OpenFOAM case; it is usually 'tau'
list_of_imp_fileds          =   ['U', 'p', 'tau', 'T', 'grad(U)']   # required feilds
particles_initation_patch 	= 	'inlet' # Name of the patch on openfoam case from which the particles are released

# The time velue to be set as active. This code reads only one time step
active_time_step            =   120     # must be a the same as the directory to be read

stretching_axis             =   'x'     # available options: 'x', 'y', or 'z'
# streching axis. The packages assumes all ellipsods strech along one of the coordinate axes.

volume_release              =   'cube'  # Volume release shape; either 'cube' or 'cylinder'
# Also, if cylindrical selected, the cylinder will have its height along input_data.stretching axis

# Number of particles releasing  at a time
n_points                    =    20     # at each release time n_points new particles/Ellipsoids are released.
d_average                   =    0.02   # average particle diameter [m]

# With seed  = int(time.time()) in particle_init.py, we set a log-normal distribution
mean  = 0
sigma_dist = 1                          # A log-normal parameter; I assigned 1 out of blue
particle_dimater_variation  =    0.005  # max variation in d_average

# setting Cylindrical or Cubical distribution volume
distribution_radius         =   0.1 # If cylinder was selected

X_min = -35.0;        X_max =   -34.0     # \Delta X release
Y_min = 0;            Y_max =   0.97
Z_min = 0;            Z_max =   0.97

translate_vector=[0, 0, 0]              # This vector translates all the initial particle points.

# rheoFibrillator simluation details    t_start = 0s cannot be changed
delta_t_release             =   1       # [s] time between the particle realeas 
t_stop_release              =   2.5 * delta_t_release     # time after which no particle will be released.
delta_t                     =   1       # [s] time step for PyVista simulation
t_end                       =   20      # [s] end of simulation
time_out_limit              =   500     # seconds
write_step                  =   4       # After how many time steps it rights (NOT [s])
coelescence_threshold       =   0       # (psi_coal.)particles with surface distance less than ``` coelescence_threshold * minor diameter ``` coalesce: a relative measure
psi_break_up                =   1.2     # breakup time extension factor
coalece_check_range         =   1.5     # How far should the code check coalesce. if 1, every droplet is checked against n_points. if 2, twice, etc. The higher, the more compuation. 

# Export all data as CSV in the given locations
report_locations            =   [-30, -20, -10, -0.002, 0.1, 5, 10, 20, 24]  # [m] along the ```stretching_axis```
tolerance_from_locations    =   1       # Saves a csv if any droplets reaches any of the locations above.

# Animation output
animation_output            =   True    # saves an animation using `imageio` library
animation_name              =   'particle motion.gif'
frames_directory_name       =   'frames'#creates a directory and saves all frames
fps                         =   20      # Frames per second
keep_all_frames             =   True    #if false, it'll delete all frames from the directory
screenshot_initiation       =   False
animation_zoom              =   3       # zoom to the center

# Visualizations:
screenshot_initiation       =   True
dynamic_renderer            =   True     # If true, it tries to show dynamic renderer in jupyter notebook which takes RAM. opions: static", "client", "server", "trame", "html", "none"
default_color               =   'gray'
breakup_color               =   'blue'    
coalesce_color              =   'red'
breakup_and_coalesced_color =   'blueviolet'
include_OF_background       =   True

# spunbonding conditions
T_die                       =   210     # [degC]

# Material properties
eta_s_m                     =   16      # [Pa.s] the matrix solvent viscosity
eta_p_m                     =   3       # [Pa.s] the matrix polymeric viscosity
rho_d_ref                   =   1125    # [kg/m3] at the reference Temperature
rho_m_ref                   =   925     # [kg/m3] at the reference Temperature
eta_d                       =   15      # [Pa.s] droplet viscosity
p                           =   eta_d / eta_s_m #viscosity ratio

# Technical details for further coding, debugging, and scientific developments:
remove_droplets_with_warning=   True     #It keeps any droplet/particle/Ellipsoid that has been marked with warning due to numerical errors.
numPy_percission            =   5         # percission control of the whole code


def sigma(T):                           # [N/m] Interfacial tension
    """It returns the interfacial tension as afunction of tempearture:
    float T -> float sigma"""
    return 1.0168e-3                       # [N/m] see my article

def rho(T):
    """It returns the densities as a function of Temperature.
    float T -> float rho_d, float rho_m"""
    pass

