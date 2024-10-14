import pyvista  as pv
import numpy    as np

import input_data

import input_data

class Ellipsoid: # Ellipsoid with a long x axis
    """This is a class of one single Ellipsoid; there is no time data, and the inputs are
        the location [x, y, z]
        the original radius,
        the time elapased since the sphirical formation,
        the particle mass,
        the droplet ID,
        and existance of any warning for this droplet.
        """
    def __init__(self, X, ID=0, r_0=input_data.d_average/2 ,t=0, mass=0.1, U=np.zeros(3), N_11=0, shear=0, T=input_data.T_die, warning = False):
        # constructor:
        self.X = X; self.U =U;  self.N_11 = N_11; self.t = t; self.shear  = shear
        self.r_0 = float(r_0);  self.mass= mass;  self.T = T; self.warning= warning; self.ID =ID
        self.location = {}      # A dictionary saving location (value) over time (key).
        
        # saving breakup and coalescenc history:
        self.breakup = False;    self.coalesce=False

        self.p = input_data.p       # viscosity ratio
        self.L = self.r_0           # Ellipsoid diameter
        self.B = self.r_0
        self.acceleration = 0       # particle accelaration after calculating the force balance

        # calculating parameters:
        self.B = ( self.r_0 ** 3 / self.L ) ** 0.5        # continuity of incompressible particle
        self.sigma = input_data.sigma(self.T)               # Surface interface d and m

        # Critical Ca for elongation
        if -0.01< np.log( self.p) -0.00645 < 0.01:                                    # checking the risk of division by zero
            self.Ca_cr_e = np.exp( 0.0331 * ( np.log( self.p) - 0.5 ) ** 2  - 0.699)    # Peters et al.
            self.warning = True
        else:                                                                               # Grace correlation
            self.Ca_cr_e = np.exp( -0.064853 - 0.2442 * np.log( self.p) + 0.02221 * np.log( self.p) ** 2 - 0.00056 / (np.log( self.p) - 0.00645) )

        # Critical Ca for shear
        if self.p < 4.01:                                         # from Fortelny 2019 review 
            self.Ca_cr_s     = np.exp( -0.506 - 0.0995 * np.log (self.p)  + 0.124 * np.log( self.p ) ** 2 - 0.115/( np.log (self.p) - np.log(4.08) ) )
        else:
            self.warning = True
            self.Ca_cr_s = 1e-4                              # for numerical reasons

    def update(self, time, OpenFOAM_mesh):
        """It updates everything of a cell using the methods within this class."""

        # location update -> update velocity -> Ca
        self.update_location_time(time).\
            update_fields_from_OF(OpenFOAM_mesh).\
                strain_rate(OpenFOAM_mesh, save_data_on_OF_mesh=True).\
                    calcuale_Ca_numbers()


    def calcuale_Ca_numbers(self):
        """This function returNs Ca_e and Ca_s.
        It calculates the surface tension as a function of T with input_data.py"""
        self.Ca_e = 2 *  self.N_11  * self.B / self.sigma
        self.Ca_s =      self.shear * self.B * self.L / self.sigma / (self.B + self.L)

        return self

    def update_location_time (self, simulation_time):
        """This function moves all the particles in the particle_mesh of PolyData format one step forward.
        It uses X = at^2 / 2 + Ut + X0 and assumes the particle mesh contains the ['mass'] attribute.
        It moves the time step and creates a new location to calculate the deformation later.
        It sores the particle past locations by updating `self.location` dictionary."""

        self.location[np.round(simulation_time, 5)] = self.X
        self.t += input_data.delta_t
        self.X = 0.5 * self.acceleration * input_data.delta_t **2 + input_data.delta_t * self.U + self.X 
        
        return self
    
    def update_fields_from_OF(self, OpenFOAM_mesh):
        '''It goes into the active time step directory, list all the avaible data (no directory),
         interpoltes them all at the location self.X, and updates the other attributes.
         This uses a buffer mesh (named pdata of PolyData) to speed up the interpolation.'''
        
        # Getting the cell centers of the OpenFOAM case:
        foam_cell_centers = OpenFOAM_mesh.cell_centers().points

        # Creating a buffer mesh for the purpose of interpolation
        pdata = pv.PolyData(foam_cell_centers)

        # An empty dictionary to hold the field arrays
        hold = {}

        for field in OpenFOAM_mesh.array_names:
            if field in input_data.list_of_imp_fileds:      # this saves some computation

                # Interpoltaing N_11 and shear:
                if 'tau' in field:
                    pdata['N_11']   = OpenFOAM_mesh['tau'][:, 0] - OpenFOAM_mesh['tau'][:, 3]
                    pdata['shear']  = OpenFOAM_mesh['tau'][:, 1]

                    self.N_11  = pv.PolyData(self.X).interpolate(pdata)['N_11'][0]      # [0] must be added because interpolate returns a numpy array
                    self.shear = pv.PolyData(self.X).interpolate(pdata)['shear'][0]

                # Interpolating other fields:
                else:
                    pdata[field]    = OpenFOAM_mesh.cell_centers()[field]
                    hold[field]     = pv.PolyData(self.X).interpolate(pdata)[field]
        
        self.U = hold['U'][0]
        self.grad_U = hold['grad(U)'][0]

        if 'T' in OpenFOAM_mesh.array_names:
            self.T = hold['T'][0]
        else: self.T = input_data.T_die     # this is used while developing the code

        return self

    def update_deformation(self):                   ### under development
        """This function calculates the reduced Ca numbers, i.e., either Ca_e^* or Ca_s^* (only the larger one matters).
        Next finds the defomration type via Huneault model, and calculates the deformation given the delta_t [s].
        It also updates the droplet location and time."""
        Ca_star = max( self.Ca_s/self.Ca_cr_s,  self.Ca_e/self.Ca_cr_e)
        #self.t_b = -( 19 * self.p)

        # Huneault conditions
        #if 0.1< Ca_star < 1: elliptical Deformation to be completed

        #if Ca_star > 4:        #Affine creep motion

        return self
    
    def strain_rate(self, OpenFOAM_mesh, save_data_on_OF_mesh = True):
        """This updates the attributes values of elongational and shear strain rates.
        Importantly, it saves the data in the OpenFOAM_mesh as well.
        
        pv.OpenFOAM_mesh ->  """

        self.grad_U
        # find the deformation matrix
        D = 0.5 * ( self.grad_U + np.transpose(self.grad_U )  )
        D = D.reshape(3,3)
        
        # Calculating shear and elongational strain rates to the Ellipsoid:
        self.shear_rate = np.sqrt(2 * (D[0, 1]**2 + D[0, 2]**2 + D[1, 2]**2))
        self.elongation_rate = np.sqrt(D[0, 0]**2 + D[1, 1]**2 + D[2, 2]**2)

        # Adding shear and elongational strain rates to the OPenFOAM mesh
        if save_data_on_OF_mesh:
            OpenFOAM_mesh.point_data['shear_rate'] = self.shear_rate
            OpenFOAM_mesh.point_data['elongation_rate'] = self.shear_rate

        return self
    
    def breakup(self, list_of_ALL_Elliposids):    ### under development
        """This function reads the time, and if it is time to break, makes it work.
        The parent Ellipsod becomes sphirical, gets halved, and moves to X+B.
        Another sphirical Ellipsoid instance will be created at X-B and retuenred, so that, it can assign a new ID.
        Note, once we have a new born all IDs bigger than the mother must shift to one higher. This way `coalesce1
        method remains functional with the IDs order. Hence, `list_of_ALL_Elliposids` updates here."""

        pass

    def coalesce(self, list_of_ALL_Elliposids):   ###under development
        """This function assesses the surface-to-surface distance to all NEARBY Ellipsoids in the inputted list.
        So, first the list members with ID + (0, self.ID + input_data.n_points]. The self.ID will optimize this process,
        as self will be checked only with members with higher IDs.
        If there exist a touch the list member will send all its mass to self.
        Lastly, the the ID number of the corresponding lits member will be returnd.
        Because it'll need to be removed from an external `list_of_ALL_Elliposids`.
        So, this list gets updated, automatically, as well."""

        pass

    def __str__(self):
        if self.warning:
            return "check the warning list: \n\t Peters correalation for Ca_cr_e was adopted. \n\t p > 4, and Ca_cr for shear cann't be used."
        return f"An Ellipsoid of L={np.round(self.L, 3)}, B={np.round(self.B, 3)}, Cr_cr={np.round(self.Ca_cr_e, 3)}"

    # plot
    def show(self): 
        """ provides a simple plot of the instance.
        Of note, It creates a new pv.Ellipsoid and trasnlates it."""

        print("Main diameter, L =", np.round(self.L, 3), ", and other diameter, B =", np.round(self.B, 3) )
        ellipse = pv.ParametricEllipsoid(self.L, self.B, self.B)
        ellipse.translate( (self.X[0], self.X[1], self.X[2]) , inplace=True)
        ellipse.plot(color='red') if self.warning else ellipse.plot(color='lightgreen')

    def add_to_Plotter(self, plotter, opacity = 1.0, default_color='darkgray', breakup_color='darkred', coalesce_color='lightpink'):
        """It inputs a pyvista.Plotter() object and addes an ellipse to it; however, 
        does not plot anything yet. Of note, It creates a new pv.Ellipsoid and trasnlates it.
        It plots"""
        ellipse = pv.ParametricEllipsoid(self.L, self.B, self.B)
        ellipse.translate( (self.X[0], self.X[1], self.X[2]) , inplace=True)

        # plotting with the proper color:
        if self.breakup:    plotter.add_mesh(ellipse, opacity=opacity, color=breakup_color )
        elif self.coalesce: plotter.add_mesh(ellipse, opacity=opacity, color=coalesce_color)
        else:               plotter.add_mesh(ellipse, opacity=opacity, color=default_color )

    def add_path_to_Plotter(self, plotter, path_opacity = 1, default_color='darkgray', breakup_color='darkred', coalesce_color='lightpink'):
        """This creates a SPline of the location of the Ellipse instance and addes it to the input plotter."""
      
        path = np.zeros( ( len(self.location.keys() )  , 3) )              # initializing the trajectory

        i = 0
        for time in self.location.keys():
            path[i, 0] = self.location[time][0]     #X
            path[i, 1] = self.location[time][1]     #Y
            path[i, 2] = self.location[time][2]     #Z
            i+=1
    
        track = pv.Spline(path)                 # creating the track with pv.Spline

        # plotting with the proper color:
        if self.breakup:    plotter.add_mesh(track, opacity=path_opacity, color=breakup_color  )
        elif self.coalesce: plotter.add_mesh(track, opacity=path_opacity, ccolor=coalesce_color)
        else:               plotter.add_mesh(track, opacity=path_opacity, color= default_color )

ellip1= Ellipsoid([3,3,10], 1, 1.3 )
print(ellip1.sigma)