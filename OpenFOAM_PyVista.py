# This library  contains all function at the interface of OpenFOAM and Pyvista.
import pyvista as pv
import numpy as np
import os
import imageio
import sys
import csv

import input_data
np.set_printoptions(precision=input_data.numPy_percission)

def read_and_prepare_OpenFOAM(active_time_step=input_data.active_time_step, print_output=True):
    """This function inputs the .foam address and creats a pv mesh.
    Next, it activates the given time step and prints the mesh data if set True.
    If also sets N_11, shear, and U_magnitude within the mesh, 
    assuming the elongation is along XX and shear along XY."""

    # Reading
    reader =  pv.POpenFOAMReader(input_data.filename) 
    mesh = reader.read()
    reader.set_active_time_value( input_data.active_time_step )

    reader.cell_to_point_creation = True  # Need point data for streamlines
    mesh = reader.read()
    boundaries = mesh["boundary"]
    internal_mesh = mesh["internalMesh"]
    print(boundaries)
    reader.set_active_time_value(active_time_step)

    # Adding more data to the mesh:
    internal_mesh.cell_data['U_magnitude'] = np.linalg.norm(internal_mesh['U'], axis=1)

    # Adding elongational and shear stress based on the choince in input_data.stretching_axis:
    tau = internal_mesh[input_data.stress_varibale]

    #Elongational stress difference along the input_data.stretching_axis
    if   input_data.stretching_axis=='x':
        internal_mesh.cell_data['N_11'] =  tau[:, 0] - np.sqrt( tau[:, 1] ** 2 + tau[:, 5] ** 2)
    elif input_data.stretching_axis=='y':
        internal_mesh.cell_data['N_11'] =  tau[:, 1] - np.sqrt( tau[:, 0] ** 2 + tau[:, 5] ** 2)
    elif input_data.stretching_axis=='z':
        internal_mesh.cell_data['N_11'] =  tau[:, 5] - np.sqrt( tau[:, 0] ** 2 + tau[:, 1] ** 2)
    else: sys.exit("wrong input in input_data.stretching_axis")

    # Shear - reading the max
    internal_mesh.cell_data['shear'] = np.maximum(   tau[:, 3], tau[:, 4], tau[:, 5]  )
    
    #printing
    if print_output:
        print("Availbale fields at the active time step: ", internal_mesh.array_names)
        print(f"All patch names: {reader.patch_array_names}")
        print(f"All patch status: {reader.all_patch_arrays_status}")
        print(f"Boundaries patches: {boundaries.keys()} \n")
        print(f"Available Time Values: {reader.time_values} \nSelected: ", active_time_step)
        print("Stress tensor (tau) exists and is in shape of ", tau.shape)
        print("\nN_11 exists along ", input_data.stretching_axis, " and is in shape of", internal_mesh['N_11'].shape)
        print("Maxium shear at any point exists andis in shape of ", internal_mesh['N_11'].shape)
        print("\nCell Data:")
        print(internal_mesh.cell_data)
        print("\nPoint Data:")
        print(internal_mesh.point_data)
        print('max(tau_XX)= ', np.max( tau[:, 0]), '\t[Pa.s]' )
        print('max(tau_YY)= ', np.max( tau[:, 1]), '\t[Pa.s]'  )
        print('max(tau_ZZ)= ', np.max( tau[:, 2]), '\t[Pa.s]'  )
        print('max(tau_XY)= ', np.max( tau[:, 3]), '\t[Pa.s]'  )
        print('max(tau_XZ)= ', np.max( tau[:, 4]), '\t[Pa.s]'  )
        print('max(tau_ZZ)= ', np.max( tau[:, 5]), '\t[Pa.s]'  )

    return internal_mesh, boundaries


def plot_all_Ellipoids_in_OF (OpenFOAM_mesh, list_of_Ellipsoids, plotter = False, save_figure = False, scalars = 'U', 
                              show_plotter = False, file_name = 'screenshot', OF_opacity= 0.3, Ellipsoids_opacity= 1,
                              OF_color = False, background_OF = True):
    """This functions plots the whole OpenFOAM U_magnitude and shows the Ellipsoid particles within it at their most recent location.
    If the scalars are defined, they must be a scalar in the OpenFOAM mesh read by PyVista, e.g.,
    if scalars='U_magnitude' is selected, OpenFOAM_mesh['U_magnitude'] must be available.
    If a Plotter object is inputted, it will add all Ellipsoids to it.

    pv.mesh, list[Ellipsods], pv.Plotter() -> -   """

    if plotter:
        for i in range( len(list_of_Ellipsoids) ):
            list_of_Ellipsoids[i].add_to_Plotter(plotter, Ellipsoids_opacity)
    
    else:
        plotter = pv.Plotter()
        for i in range( len(list_of_Ellipsoids) ):
            list_of_Ellipsoids[i].add_to_Plotter(plotter, Ellipsoids_opacity)

    # ploting the OpenFOAM (OF)
    if background_OF:
        if OF_color:    plotter.add_mesh(OpenFOAM_mesh, color=OF_color, opacity=OF_opacity)
        else:           plotter.add_mesh(OpenFOAM_mesh, scalars=OpenFOAM_mesh[scalars], opacity=OF_opacity)

    # Output managing
    if save_figure :    plotter.screenshot(file_name)
    if show_plotter:    plotter.show()

def delete_Ellipsoid_from_list(Ellipsoid, list_of_particles):
    """This function removes the input ellipsoid from the list and renews the IDs.
    It is used only in method ```delete_ellipsoids_with_warning``` and ```check_coalesce```.

    Elliposid, list -> -"""

    for i in list_of_particles:
        if Ellipsoid.ID == i.ID:

            #shifting all next IDs down:
            for j in range(list_of_particles.index(i)+1, len(list_of_particles) ):
                list_of_particles[j].ID -= 1

            #deleting i-th Ellipsoid
            list_of_particles.remove(i)  

def delete_ellipsoids_with_warning(list_of_all_particles):
    """This functions removes all ellipsoids that have been marked with a warning, e.g., the ones with invalid location.
    list -> -"""

    i = 0
    while i < len(list_of_all_particles):
        if list_of_all_particles[i].warning:
            delete_Ellipsoid_from_list( list_of_all_particles[i], list_of_all_particles)
        i+= 1

def check_coalesce(list_of_partilces, check_range= int(input_data.n_points * int(input_data.coalece_check_range) ) ):   ###under development
    """This function assesses the surface-to-surface distance to all NEARBY Ellipsoids in the range of self to self+check_range.
    So, first the list members with ID + (0, self.ID + check_range]. The self.ID is used in this process, as self will be checked
    only with members with higher IDs. If there exist a touch the list member will send all its mass to self.
    Lastly, the the ID number of the corresponding lits member will be deleted within this method from external `list_of_ALL_Elliposids`.
    So, this list gets updated, automatically, as well.

    self, list(Ellipsoids), int -> - """
    
    i = 0; Number_of_coalescences = 0
    while i < len( list_of_partilces):
        #print(list_of_partilces[i].X)

        #checking all particles from 0 to -check_range:
        if len(list_of_partilces) - i <= check_range:

            # Another loop for the one-to-one coalesce check; ```for```` loop is not applicable as the ```len(list_of_partilces)``` changes.
            j=1
            while j < min(check_range, len(list_of_partilces) - i):

                # Eucleadian_distance between Ellipsoid centers:
                D = np.sqrt( (list_of_partilces[i].X[0] - list_of_partilces[i+j].X[0]) ** 2\
                                             + (list_of_partilces[i].X[1] - list_of_partilces[i+j].X[1]) ** 2\
                                             + (list_of_partilces[i].X[2] - list_of_partilces[i+j].X[2]) ** 2 )
                # Saving some computation using L:
                if D < list_of_partilces[i].L + list_of_partilces[i+j].L:
                    
                    # Getting the streching direction:
                    if      input_data.stretching_axis == 'x': axis =0
                    elif    input_data.stretching_axis == 'y': axis =1
                    elif    input_data.stretching_axis == 'z': axis =2
                    else:   sys.exit("\n Wrong input for ```stretching_axis``` in ```input_data.py```")

                    # Setting up a temp coordinate system X and Y (from x, y, and z)
                    delta_X = np.abs( list_of_partilces[i].X[axis] - list_of_partilces[i+j].X[axis]  )
                    delta_Y = 0 

                    for k in range(3):
                        if i != axis:
                            delta_Y += (list_of_partilces[i].X[k] - list_of_partilces[i+j].X[k] ) ** 2
                    delta_Y = np.sqrt(delta_Y)

                        # Angle between vector L (the streching direction)
                    theta = np.arctan( delta_Y / delta_X )

                    # Surface to surface distance:
                    d = D - \
                        np.sqrt(  ( list_of_partilces[i].L   * np.cos(theta) ) ** 2   +  ( list_of_partilces[i].B   * np.sin(theta) ) ** 2  )\
                       -np.sqrt(  ( list_of_partilces[i+j].L * np.cos(theta) ) ** 2   +  ( list_of_partilces[i+j].B * np.sin(theta) ) ** 2  )

                    #print(f"d between {list_of_partilces[i].ID} and {list_of_partilces[j+i].ID} = ", d)
                    #Coalescence check. For numerical resones no coalescence occurs for freshly released variables.
                    if d <= input_data.coelescence_threshold  *  min( list_of_partilces[i].B, list_of_partilces[i+j].B ) and list_of_partilces[i].t > 3 * input_data.delta_t and list_of_partilces[i+j].t > 3 * input_data.delta_t:

                        # partilce i gets updated and j removed:
                        list_of_partilces[i].coalesce= True
                        print(list_of_partilces[i].ID," and ",list_of_partilces[i+j].ID," coalesced.")

                        #Update: mass, location, r_0, L, B, delete i+j
                            # weight averaged location:
                        list_of_partilces[i].X = (list_of_partilces[i].X        * list_of_partilces[i].mass\
                                                + list_of_partilces[i+j].X      * list_of_partilces[i+j].mass)\
                                                /(list_of_partilces[i].mass     + list_of_partilces[i+j].mass)
                            # r_0:
                        r_0 = ( np.power(list_of_partilces[i].r_0 , 3) + np.power(list_of_partilces[i+j].r_0 , 3) ) ** (1/3)
                        
                        # debugging -------------------
                        L_temp = r_0 * max( list_of_partilces[i].L  / list_of_partilces[i].r_0, list_of_partilces[i+j].L / list_of_partilces[i+j].r_0)
                        B_temp = (r_0 ** 3 / L_temp ) ** 0.5

                            # L
                        list_of_partilces[i].L = L_temp
                            # B
                        list_of_partilces[i].B = B_temp

                        list_of_partilces[i].r_0 = r_0

                            # mass
                        list_of_partilces[i].mass += list_of_partilces[i+j].mass
                        #print(f"after coalese, Ellipsoid { list_of_partilces[i].ID} has B: {list_of_partilces[i].B}")

                        # Deleting [i+j]-th particle
                        delete_Ellipsoid_from_list( list_of_partilces[i+j], list_of_partilces)

                        Number_of_coalescences +=1
                # inner loop counter    
                j+= 1    

        # Loop iteration continues after updating the length of list_of_partilces
        i+= 1

                        #-----------------
    return Number_of_coalescences

def add_pathes_to_Plotter (OpenFOAM_mesh, list_of_Ellipsoids, plotter, save_figure = False, scalars = 'U', 
                           file_name = 'all_pathes', OF_opacity= 0.2, Ellipsoids_opacity= 1, OF_color = False):
    """This plots all tracking pathes using the `Ellipsoid.add_path_to_Plotter`.
    
    pv.mesh, list[Ellipsods], pv.Plotter() -> -   """

    for i in range( len(list_of_Ellipsoids) ):
        list_of_Ellipsoids[i].add_path_to_Plotter(plotter, Ellipsoids_opacity)

    # ploting the OpenFOAM (OF)
    if OF_color:    plotter.add_mesh(OpenFOAM_mesh, color=OF_color, opacity=OF_opacity)
    else:           plotter.add_mesh(OpenFOAM_mesh, scalars=OpenFOAM_mesh[scalars], opacity=OF_opacity)

    if save_figure: plotter.screenshot(file_name)  

def start_frames_folder (folder_name = input_data.frames_directory_name):
    """This function handles creates a dierctory to save all frames for scienetific development.
    
    str, bool ->  - """

    if os.path.exists(folder_name):
        for filename in os.listdir(folder_name):
            file_path = os.path.join(folder_name, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                for subfile in os.listdir(file_path):
                    os.remove(os.path.join(file_path, subfile))  # Remove files inside subdirectory
                os.rmdir(file_path)  # Remove empty subdirectory
        os.rmdir(folder_name)  # Remove the main directory

    os.makedirs(folder_name)  # Recreate the folder

def create_anim(file_name =input_data.animation_name, frame_folder_name = input_data.frames_directory_name):
    """This function inputs a folder name and creates an animation with all the sorted png files in it as frames.
    
    folder name -> -"""
    
    # Get a list of all PNG files in the frames directory, sorted by time (assuming the filenames are str(t)+'.png')
    frame_files = sorted([os.path.join(frame_folder_name, file) 
                        for file in os.listdir(input_data.frames_directory_name) 
                        if file.endswith('.png')])

    # Create the animation from the files without loading them all into memory
    with imageio.get_writer(file_name, fps=input_data.fps) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)  # Read each frame
            writer.append_data(image)  # Write each frame to the animation

def delete_frames_folder(condition = input_data.keep_all_frames):
    """this fundtion deletes the frame folder if it as indicated so in ```input_data.py``` 
    bool -> - """

    if condition: pass  # no need to delete
    else: 
        folder_name = input_data.frames_directory_name
        
        if not input_data.keep_all_frames and os.path.exists(folder_name):
            for filename in os.listdir(folder_name):
                file_path = os.path.join(folder_name, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.remove(file_path)  # Remove files or symbolic links
                elif os.path.isdir(file_path):
                    os.rmdir(file_path)  # Remove empty subdirectories
            os.rmdir(folder_name)  # Finally, remove the empty 'frames' folder

def save_vtu_snapshot (list_of_particles, t):
    """This function inputs all Ellipsoids in the simulation within a list and saves them as stl as t.vtp.
    
    list, float ->  VTP file"""
    
    # Create a list to save blocks of mesh
    wrting_mesh = pv.MultiBlock()

    for i in list_of_particles: 
        # drawing Ellipsoids with low quality, i.e., decimation
        
        ellipse = pv.ParametricEllipsoid(i.L, i.B, i.B).translate( (i.X[0], i.X[1], i.X[2]) , inplace=True).decimate(0.3)
        #if i.breakup: 
        ellipse.color='red'
        wrting_mesh.append(ellipse) 

    # merging them
    wrting_mesh = pv.merge( [j for j in wrting_mesh])

    # write
    wrting_mesh.save(str(  np.round(t, 4) )+".vtp")
    print("A snapshopt mesh of ", np.round(t, 3), "was saved.")

    del wrting_mesh

import csv
import numpy as np

def save_properties_at_location(locations, list_of_particles, tolerance):
    """This function checks if any particle is very close to one of the locations and saves data in a CSV file.
    
    list(float), list(Ellipsoids) -> csv file"""

    # pick the axis
    if input_data.stretching_axis == 'x':   i = 0
    elif input_data.stretching_axis == 'y': i = 1
    elif input_data.stretching_axis == 'z': i = 2


    # Iterate through a copy of the locations to avoid modification issues during iteration
    for location in locations[:]:
        # Check if any particle's X[i] is close to the location
        if any(abs(particle.X[i] - location) < tolerance for particle in list_of_particles):
            #print(f"Removing {location} from {locations}")
            locations.remove(location)

            # Prepare data for CSV
            csv_data = []
            for particle in list_of_particles:
                csv_data.append([
                    particle.X[0], particle.X[1], particle.X[2], particle.t, particle.r_0,
                    particle.L, particle.B, particle.fibrillation_index[max(particle.fibrillation_index.keys())],
                    particle.breakup, particle.coalesce, particle.last_breakup, particle.N_11, 
                    particle.shear, particle.Ca_e, particle.Ca_s, particle.Ca_star_e, particle.Ca_star_s
                ])

            # Save to CSV
            filename = f'at_{location}_properties.csv'
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([
                    'X', 'Y', 'Z', 't', 'r_0', 'L', 'B', 'fibrillation_index',
                    'breakup', 'coalesce', 'last_breakup', 'N_11', 'shear',
                    'Ca_e', 'Ca_s', 'Ca_star_e', 'Ca_star_s'
                ])
                writer.writerows(csv_data)

            print(f"Data saved to {filename}")