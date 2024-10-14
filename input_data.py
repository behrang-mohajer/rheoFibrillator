######## Use this file to input the attached .py code.	########

# Address to .foam file:
filename			        = 	'channel.foam'
directory                   =   './'      #it directory if read NOT from the current address

# List of the fields that must be read from teh openFoam Case:
list_of_imp_fileds          =   ['U', 'p', 'tau', 'T', 'grad(U)']

# the time velue to be set as active. This code reads only one time step
# amond the ones availbale in the OpenFOAM directory. Make sure it matches
active_time_step            =   30      # must be a the same as the directory

# Name of the patch on openfoam case from which the particles are released
particles_initation_patch 	= 	'inlet'

# Number of particles releasing  at a time
n_points                    =   3       # at each release time
d_average                   =   0.08    #average particle diameter [m]

# PyVista simluation details
delta_t_release             =   4       # [s] time between the particle realeas 
delta_t                     =   0.1     # [s] time step for PyVista simulation
t_end                       =   5      # [s] end of simulation
time_out_limit              =   30     # seconds

# Animation output
animation_output            =   True    # saves an animation using `imageio` library
animation_name              =   'particle motion.gif'
fps                         =   20      # Frames per second

# spunbonding conditions
R_die                       =   2e-4    # [m]
T_die                       =   210     # [degC]
screw_rpm                   =   100     # [rpm] or 14 Hz
drafter_freq                =   30      # [Hz]
feed_rate_SB                =   4/10800 # [kg/s] or 4/3 [kg/h] or 4 [Hz]
n_die_holes                 =   90      # -
mass_rate_m                 =   feed_rate_SB / n_die_holes  # [kg/s]

# Material properties
eta_s_m                     =   0.8     # [Pa.s] the matrix solvent viscosity
eta_p_m                     =   0.2     # [Pa.s] the matrix polymeric viscosity
rho_d_ref                   =   900     # [kg/m3] at the reference Temperature
rho_m_ref                   =   900     # [kg/m3] at the reference Temperature

eta_d                       =   1.0     # [Pa.s] droplet viscosity
p                           =   eta_d / eta_s_m #viscosity ratio
def sigma(T):                           # [N/m] Interfacial tension
    """It returns the interfacial tension as afunction of tempearture:
    float T -> float sigma"""
    return 2 / (1.1*T_die - T)              # out of blue XXXXXXXXXXXXXXX

def rho(T):
    """It returns the densities as a function of time.
    float T -> float rho_d, float rho_m"""
    pass

