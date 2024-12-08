import pyvista  as pv
import numpy    as np

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
    def __init__(self, X, ID=0, r_0=input_data.d_average/2 ,t=0.0, U=np.zeros(3), N_11=10, shear=10, T=input_data.T_die, warning = False):
        # constructor:
        self.X = X; self.U =U;  self.N_11 = N_11; self.t = t; self.shear  = shear; self.Peters = False  # To mark the use of Peters eq. if so
        self.r_0 = float(r_0);  self.T = T; self.warning= warning; self.ID =ID; self.fibrillation_index = {self.t: 1} # fibrillation_index, D = (L-B)/(L+B)
        self.location = {self.t: self.X}      # A dictionary saving location (value) over time (key).
        
        # saving breakup and coalescenc history:
        self.breakup = False;    self.coalesce=False;  self.last_breakup = 0.0; self.strain = 0.0   # the last one is under development

        self.p = input_data.p       # viscosity ratio
        self.L = self.r_0           # Ellipsoid diameter
        self.B = self.r_0
        self.acceleration = 0       # particle accelaration after calculating the force balance
        self.mass= 4 * np.pi / 3 * np.power(self.r_0, 3) * input_data.rho_d_ref

        # Selecting between Grace et al. and Peters et al. correlations for Ca_cr. the former is more accurate but not applicable aroound a certain range
            # Critical Ca for elongation
        if -0.02< np.log( self.p) -0.00645 < 0.02:                                    # checking the risk of division by zero
            self.Ca_cr_e = np.exp( 0.0331 * ( np.log( self.p) - 0.5 ) ** 2  - 0.699)    # Peters et al.
            self.Peters  = True
        else:                                                                               # Grace correlation
            self.Ca_cr_e = np.exp( -0.064853 - 0.2442 * np.log( self.p) + 0.02221 * np.log( self.p) ** 2 - 0.00056 / (np.log( self.p) - 0.00645) )

        # Critical Ca for shear
        if self.p < 4.01:                                    # from Fortelny 2019 review 
            self.Ca_cr_s     = np.exp( -0.506 - 0.0995 * np.log (self.p)  + 0.124 * np.log( self.p ) ** 2 - 0.115/( np.log (self.p) - np.log(4.08) ) )
        else:
            self.Ca_cr_s = 1e-4                              # for numerical reasons

    def update(self, time, OpenFOAM_mesh):
        """It updates everything of a cell using the methods within this class."""

        # location update -> update velocity -> Ca
        self.update_location_time(time).\
            update_fields_from_OF(OpenFOAM_mesh).\
                strain_rate(OpenFOAM_mesh, save_data_on_OF_mesh=True).\
                    calcuale_Ca_numbers().deform()


    def calcuale_Ca_numbers(self):                              #------------ under development the 2nd line gives division by 0
        """This function returNs Ca_e and Ca_s.
        It calculates the surface tension as a function of T with input_data.py"""
        self.Ca_e = 2 *     self.N_11  * self.B / self.sigma
        self.Ca_s = np.abs( self.shear * self.B * self.L / self.sigma / (self.B + self.L)   )

        # Ca*  based on constant p throughout the domian
        self.Ca_star_e = self.Ca_e / self.Ca_cr_e
        self.Ca_star_s = self.Ca_s / self.Ca_cr_s

        return self

    def update_location_time (self, simulation_time):
        """This function moves all the particles in the particle_mesh of PolyData format one step forward.
        It uses X = at^2 / 2 + Ut + X0 and assumes the particle mesh contains the ['mass'] attribute.
        It moves the time step and creates a new location to calculate the deformation later.
        It sores the particle past locations by updating `self.location` dictionary."""

        self.location[np.round(simulation_time, 5)] = self.X
        self.t = np.round(self.t + input_data.delta_t, 5)
        self.X = 0.5 * self.acceleration * input_data.delta_t **2 + input_data.delta_t * self.U + self.X 

        # checking for errors:
        if self.t % 5 == 0 and self.ID % 5 == 0:
            if np.isnan(self.X).any(): 
                print(f"Warning: location error for particle {self.ID} at {self.X}")
                self.warning = True

        # deformation index saving at every 2 steps and at the end.
        if np.round( self.t / input_data.delta_t) % 2 == 0 or np.round( input_data.t_end - self.t ) < 2 * input_data.delta_t:
            self.fibrillation_index[ np.round(simulation_time, 3)] = ( 2* self.B) / (self.L + self.B) 
        
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

        important_fields= input_data.list_of_imp_fileds + ['N_11', 'shear']

        for field in OpenFOAM_mesh.array_names:
            if field in important_fields:
                pdata[field]    = OpenFOAM_mesh.cell_centers()[field]
                hold[field]     = pv.PolyData(self.X).interpolate(pdata)[field]
        
        self.U = hold['U'][0]
        self.grad_U = hold['grad(U)'][0]
        self.N_11 = hold['N_11']
        self.shear= hold['shear']

        if 'T' in OpenFOAM_mesh.array_names:
            self.T = hold['T'][0]
        else: self.T = input_data.T_die     # this is used while developing the code

        self.sigma = input_data.sigma(self.T)               # Surface interface d and m

        return self

    def deform(self):                   ### under development
        """This function calculates the reduced Ca numbers, i.e., either Ca_e^* or Ca_s^* (only the larger one matters).
        Next finds the defomration type via Huneault model, and calculates the deformation given the delta_t [s].
        It also updates the droplet location and time."""
        self.Ca_star = max( self.Ca_s/self.Ca_cr_s,  self.Ca_e/self.Ca_cr_e)

        # Huneault conditions
        if self.Ca_star < 0.1:                # no deformation
            self.L = self.r_0; self.B = self.r_0

        elif 0.1< self.Ca_star < 4:          #elliptical Deformation to be completed
            # Max DI by Cox approved by Padsalgikar
            DI = max(self.Ca_e, self.Ca_s) /2 *   (19 * self.p + 16) / (16 * self.p + 16) / ( ( 19 * self.p * max(self.Ca_e, self.Ca_s)  / 40 ) **2 + 1 ) ** 0.5
            
            # Eq. 10 in my article
            L_max = self.r_0 * np.power( (1+ DI) / (1-DI) , 2/3)
            if L_max < self.L: self.L = L_max
            self.B = np.power( np.power(self.r_0, 3) / self.L, 0.5)

        else:        #if Ca_star > 4:        #Affine creep motion
            t_init_particle  = min( self.location.keys() )         # np.min did not work on this dictionary
            self.strain += max(self.shear_rate,  self.elongation_rate) * input_data.delta_t
        
            self.L +=  2 * self.r_0 * self.elongation_rate * np.exp( self.elongation_rate * ( self.t - t_init_particle ) )  * input_data.delta_t
            self.B = np.power( np.power(self.r_0, 3) / self.L, 0.5)

        if (self.L < 0 or self.B < 0): 
            self.warning = True
            print(f"Dimension error with particle {self.ID}, r_0: {self.r_0}, L: {self.L}")

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
    
    def check_breakup(self):    ### under development
        """This function reads the self.time, and if it is time to break, enforces it.
        The parent Ellipsod becomes sphirical, gets halved, and moves to X-3L.
        Another sphirical Ellipsoid instance will be created at X+3L and returned with self.ID+1.
        Note, once we have a new born all IDs bigger than the mother must shift to one higher in the main loop.
        This way `coalesce` method remains functional with the IDs order."""

        if self.Ca_star > 0.1:
            # breakup equation by Taylor ------------### under development & only elongational
            t_b_e = - 19 * input_data.p * self.Ca_e / 20 / self.elongation_rate * np.log( 1 -  1 / 3 / self.Ca_e * ( 19 * input_data.p +19)  / ( 19 * input_data.p + 16) )
            t_b_s = - 19 * input_data.p * self.Ca_s / 20 / self.shear_rate      * np.log( 1 -  1 / 4 / self.Ca_s * ( 19 * input_data.p +19)  / ( 19 * input_data.p + 16) )
            # Pick the smaller = >
            t_b = min(t_b_e, t_b_s) * input_data.psi_break_up

                # break up occuring:
            if self.t - self.last_breakup > t_b:
                self.breakup = True
                
                # printing for developing purposes: ------------------------------------ under development
                print('At t: ', self.t, 'Elliposid ', self.ID, 'broke at X: ', self.X)
                #print(self.t - input_data.delta_t)
                
                x_previous_t = self.location[ np.round(self.t - input_data.delta_t  , 3)  ]
                x_daughter1  = self.X - (self.X - x_previous_t) / np.linalg.norm(self.X - x_previous_t)   * self.L * 3      # 3 for numerical reasons
                x_daughter2  = self.X + (self.X - x_previous_t) / np.linalg.norm(self.X - x_previous_t)   * self.L * 3

                # Reset to a new sphere
                self.X = x_daughter1
                self.last_breakup = self.t
                self.mass /= 2; self.strain = 0.0
                self.r_0  /= 2 ** (1/3); self.L = self.r_0; self.B = self.r_0

                # create the proceeding daughter Ellipsoid:
                new_ellipsoid = Ellipsoid(x_daughter2, r_0=self.r_0, ID=self.ID +1,  T = self.T, t=self.t)
                new_ellipsoid.last_breakup = self.t

                return new_ellipsoid
        return False

    def __str__(self):
        if self.Peters: print("Peters equation was used for numerical reasons.")
        if self.warning:
            return "Check the warning list: Peters correalation for Ca_cr_e, Odd dimensions, etc."
        return f"An Ellipsoid of L={self.L}, B={self.B}, t= {self.t}, Cr_cr={np.round(  max(self.Ca_cr_e, self.self.Ca_cr_s) , 3)}"

    # plot
    def show(self): 
        """ provides a simple plot of the instance.
        Of note, It creates a new pv.Ellipsoid and trasnlates it."""

        print("Main diameter, L =", np.round(self.L, 3), ", and other diameter, B =", np.round(self.B, 3) )
        ellipse = pv.ParametricEllipsoid(self.L, self.B, self.B)
        ellipse.translate( (self.X[0], self.X[1], self.X[2]) , inplace=True)
        ellipse.plot(color='red') if self.warning else ellipse.plot(color='lightgreen')

    def add_to_Plotter(self, plotter, opacity = 1.0, default_color='darkgray'):
        """It inputs a pyvista.Plotter() object and addes an ellipse to it; however, 
        does not plot anything yet. Of note, It creates a new pv.Ellipsoid and trasnlates it.
        It plots"""
        ellipse = pv.ParametricEllipsoid(self.L, self.B, self.B)
        ellipse.translate( (self.X[0], self.X[1], self.X[2]) , inplace=True)

        # plotting with the proper color:
        if self.breakup:
            if self.coalesce: plotter.add_mesh(ellipse, opacity=opacity, color=input_data.breakup_and_coalesced_color)
            else: plotter.add_mesh(ellipse, opacity=opacity, color=input_data.breakup_color )
        elif self.coalesce: plotter.add_mesh(ellipse, opacity=opacity, color=input_data.coalesce_color)
        else:               plotter.add_mesh(ellipse, opacity=opacity, color=input_data.default_color )

    def add_path_to_Plotter(self, plotter, path_opacity = 1, default_color='darkgray', breakup_color='red', coalesce_color='green'):
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
        if self.breakup:    plotter.add_mesh(track, opacity=path_opacity, color = breakup_color  )
        elif self.coalesce: plotter.add_mesh(track, opacity=path_opacity, color = coalesce_color )
        else:               plotter.add_mesh(track, opacity=path_opacity, color = default_color  )