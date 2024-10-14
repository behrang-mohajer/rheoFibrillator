##### This file (class) initiates circular particles and return their locations and radii in NumPy arrays. ####
import numpy as np
import matplotlib.pyplot as plt
import time
import pyvista as pv

import input_data
import particle_init
import Ellipsoid

n_points = input_data.n_points
radius = 1          # distribution radius
particle_average_diamter = input_data.d_average    
particle_dimater_variation = 0.05

# Setting a log-normal distribution
seed  = int(time.time())
mean  = 0
sigma = 1    # out of blue

# this vector translates all the initial particle points:
translate_vector=[0, 0, 0]

def mass_distribution (particle_mesh, d_mass_ration):
    pass


def init_particles(n_points = n_points, seed = seed, show_and_save_fig=False, plot_size_sacle = 10, sigma = sigma, mean = mean):
    """This function inputs the particles number and seed to initiate their distribution on the release patch.
    It provides a log-normal distribution based on He et al.'s suggetsion.
    
    The mass distribution is still under development.
    
    int, float, float, float  -> pv.PolyData()
    """


    # log-normal on r and uniform on theta
    r       = radius * np.random.lognormal(mean, sigma, n_points) 
        # bounding the radii with a mathematical trick:
    r[r > radius] = r[ r > radius] / ( 1 + r[ r > radius])

    theta   = np.random.uniform(0, 2*np.pi, n_points) #* 2 * np.pi()


    particles_locations = np.column_stack( ( np.zeros( n_points), 
                                            r * np.cos(theta), 
                                            r * np.sin(theta) ) )    + translate_vector
    particles_radii     = np.random.uniform( particle_average_diamter - particle_dimater_variation,
                                            particle_average_diamter + particle_dimater_variation,
                                            n_points )
    
    particles_mesh = pv.PolyData(particles_locations) 
    particles_mesh['radii'] = particles_radii
    particles_mesh['mass']  = np.ones( particles_mesh.n_points) * 0.1       # to be devveloped

    if show_and_save_fig:
        plt.figure()
        ax=plt.subplot(111, polar=True)
        ax.scatter(theta, r, s=particles_radii * 2 * plot_size_sacle)
        plt.title("The inital distribution of the particles' centers")
        plt.savefig("The inital distribution of the particles.png", dpi=300)
        plt.show()

    return particles_mesh



def release_new_particles(list_of_partilces):
    """It calls particle_init.init_particles, adds more particles at the same patch, and updates
     the input list with more Ellipoids instances. The Ellipsoid ID will be the continuation of the 
     size of the input list. Also, this does not change the seed.
     
     list[Ellipoid] ->  - """

    new_particles = particle_init.init_particles( input_data.n_points)

    for i in range( new_particles.n_points) :
        list_of_partilces.append( Ellipsoid.Ellipsoid(new_particles.points[i], len( list_of_partilces) + i, new_particles['radii'][i]) )


