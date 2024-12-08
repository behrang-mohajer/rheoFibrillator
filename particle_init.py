##### This file (class) initiates circular particles and return their locations and radii in NumPy arrays. ####
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import input_data
import Ellipsoid
import sys
import time

n_points = input_data.n_points

# Setting a log-normal distribution
seed  = int(time.time())                # Setting a log-normal distribution
mean  = input_data.mean
sigma = input_data.sigma_dist    # out of blue
particle_average_diamter = input_data.d_average    
particle_dimater_variation = input_data.particle_dimater_variation

# Cylindrical distribution
radius = input_data.distribution_radius       # distribution radius / radius of die capillary hole - 5%
#height is the same az Z below.

# this vector translates all the initial particle points:
translate_vector= input_data.translate_vector


def init_cylindrical_particles(n_points = n_points, show_and_save_fig=False, plot_size_sacle = 10, sigma = sigma, mean = mean):
    """This function inputs the particles number and seed to initiate their distribution on the release patch.
    It provides a log-normal distribution based on He et al.'s suggetsion.
    
    The mass distribution is still under development.
    
    int, float, float, float  -> pv.PolyData()
    """


    # log-normal on r and uniform on theta
    r       = radius * np.random.lognormal(mean, sigma, n_points) 
        # bounding the radii with a mathematical trick:
    r[r > radius] = r[ r > radius] / ( 1.1 * r[ r > radius]) * radius

    
    theta   = np.random.uniform(0, 2*np.pi, n_points) #* 2 * np.pi()

    # Cubical distribution. If cylindrical selected, the cylinder will have its height along input_data.stretching axis 
    
    # Finding the heigh axis for the cylinder
    if      input_data.stretching_axis == 'x':  
        X = np.random.uniform( input_data.X_min,  input_data.X_max, n_points)
        particles_locations = np.column_stack( (X, r * np.cos(theta), r * np.sin(theta) ) )    + translate_vector
        
    elif    input_data.stretching_axis == 'y':  
        Y = np.random.uniform( input_data.Y_min,  input_data.Y_max, n_points)
        particles_locations = np.column_stack( (r * np.sin(theta), Y, r * np.cos(theta) ) )    + translate_vector

    elif    input_data.stretching_axis == 'z':   
        Z = np.random.uniform( input_data.Z_min,  input_data.Z_max, n_points)
        particles_locations = np.column_stack( ( r * np.cos(theta), r * np.sin(theta), Z ) )   + translate_vector

    else:   sys.exit("\n The stretching axis was not properly defined in input_data.py\nTerminating program.")

    # Correcting Y and Z > 0 always:
    particles_locations[:, 1:] = np.abs(particles_locations[:, 1:])

    particles_radii     = np.random.uniform( particle_average_diamter - particle_dimater_variation,
                                            particle_average_diamter + particle_dimater_variation,
                                            n_points )
    
    # ensure are particles have positive radii:
    particles_radii = np.abs(particles_radii)
    particles_radii[ particles_radii==0 ] = 1e-6    #numerical reasons
    
    particles_mesh = pv.PolyData(particles_locations) 
    particles_mesh['radii'] = particles_radii
    particles_mesh['mass']  = 4 /3 * particles_radii ** 3 * np.pi * input_data.rho_d_ref

    if show_and_save_fig:
        plt.figure()
        ax=plt.subplot(111, polar=True)
        ax.scatter(theta, r, s=particles_radii * 2 * plot_size_sacle)
        plt.title("The inital distribution of the particles' centers")
        plt.savefig("The inital distribution of the particles.png", dpi=300)
        plt.show()

    return particles_mesh

def init_cubical_particles(n_points = n_points, show_and_save_fig=False, plot_size_sacle = 10):
    """This function releases pre-defined particles within a defined withn a cube of boundaries the same as x, y, and z.
    It provides a log-normal distribution based on He et al.'s suggetsion.
    
    The mass distribution is still under development.
    
    int, float, float, float  -> pv.PolyData()
    """

    # Cubical distribution. If cylindrical selected, the cylinder will have its height along input_data.stretching axis 
    X = np.random.uniform( input_data.X_min,  input_data.X_max, n_points)
    Y = np.random.uniform( input_data.Y_min,  input_data.Y_max, n_points)
    Z = np.random.uniform( input_data.Z_min,  input_data.Z_max, n_points)

    particles_locations = np.column_stack( (X, Y, Z ) )    + translate_vector
    particles_radii     = np.random.uniform( particle_average_diamter - particle_dimater_variation,
                                            particle_average_diamter + particle_dimater_variation,
                                            n_points )
    # ensure are particles have positive radii:
    particles_radii = np.abs(particles_radii)
    particles_radii[ particles_radii==0 ] = 1e-6    #numerical reasons
    
    particles_mesh = pv.PolyData(particles_locations) 
    particles_mesh['radii'] = particles_radii
    particles_mesh['mass']  = 4 /3 * particles_radii ** 3 * np.pi * input_data.rho_d_ref

    if show_and_save_fig:
        plt.figure()
        ax=plt.subplot(111, polar=False)
        ax.scatter(X, Y, Z, s=particles_radii * 2 * plot_size_sacle)
        plt.title("The inital distribution of the particles' centers")
        plt.savefig("The inital distribution of the particles.png", dpi=300)
        plt.show()

    return particles_mesh

def release_new_particles(list_of_partilces, t=0):
    """It calls one of the functions above depeding on the choice listed in ```input_data.py```, adds more particles at the same patch, and updates
     the input list with more Ellipoids instances. The Ellipsoid ID will be the continuation of the 
     size of the input list. Also, this does not change the seed. FYI, time rounds to prevent numeric errors.
     
     list[Ellipoid] ->  - """

    if input_data.volume_release == 'cube':
        new_particles = init_cubical_particles()
    elif input_data.volume_release == 'cylinder':
        new_particles = init_cylindrical_particles()
    else: sys.exit(f"\nError: shape {input_data.volume_release} in input_data.py is not defined.")

    for i in range( new_particles.n_points) :
        list_of_partilces.append( Ellipsoid.Ellipsoid(new_particles.points[i], len( list_of_partilces) + 1, new_particles['radii'][i], t= np.round(t, 3)) )
