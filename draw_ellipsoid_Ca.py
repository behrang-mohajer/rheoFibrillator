import pyvista as pv
import numpy as np

L = 2; B = 0.5; N = 200
# Create a sphere
sphere = pv.Sphere(radius=1.0, theta_resolution=100, phi_resolution=100)

# Scale the sphere to create an ellipsoid
a, b, c = 2.0, 0.5, 0.5  # Semi-axes lengths
ellipsoid = sphere.scale([a, b, c])

# Compute lines
theta = np.linspace(0, 2 * np.pi, N)

# plot a circle on yz
z = B * np.cos( theta )
y = B * np.sin( theta )
yz = pv.Spline( np.column_stack( ( np.zeros(N) , y, z)) )

# plot the ellipse on xy
x = L * np.cos( theta)
xy = pv.Spline( np.column_stack( (x, y, np.zeros(N)) ))

# Plot the ellipsoid with curvature mapped to the surface
plotter = pv.Plotter()

# Making a line
line_orange = np.array( [[0, 0.7, 0,], [2, 0.7, 0,]] ) 
line_black  = np.array( [[0, 0, 0], [0, 0, 0.5], [0, 0, 0], [0, 0.5, 0]])

# Add the ellipsoid to the plotter
plotter.add_mesh(ellipsoid, opacity=0.4)

#plotting the lines:
plotter.add_mesh(yz, line_width=4, color='black')
plotter.add_mesh(xy, line_width=4, color='orange')

# adding arrows
plotter.add_arrows( np.array([1, 0, 0.5]) , np.array([ 1, 0, 0]), mag=0.8, color='orange' )
plotter.add_arrows( np.array([0, 0, 0]),np.array([ 1, 0, 0]) , mag=0.8, color='black')

# add a line
plotter.add_lines(line_orange,  color='orange', connected=False, width=2, label=r" Major diameter, L, along x") 
plotter.add_lines(line_black,   color='black' , connected=False, width=2, label=r"Minor diameters, B, along y and z") 
plotter.add_axes()
#plotter.add_legend(loc='lower right', size=(0.4, 0.4) )

plotter.save_graphic("Ca surfaces.svg")

# Show the plot
plotter.show()