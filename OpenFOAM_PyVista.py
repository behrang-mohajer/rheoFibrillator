# This library  contains all function at the interface of OpenFOAM and Pyvista.
import pyvista as pv
import numpy as np
import os

import input_data

def read_and_prepare_OpenFOAM(filename, active_time_step=input_data.active_time_step, print_output=True, return_boundaries=False):
    """This function inputs the .foam address and creats a pv mesh.
    Next, it activates the given time step and prints the mesh data if set True.
    If also sets N_11, shear, and U_magnitude within the mesh, 
    assuming the elongation is along XX and shear along XY."""

    # Reading
    reader =  pv.POpenFOAMReader(input_data.filename) 
    mesh = reader.read()
    internal_mesh = mesh["internalMesh"]   
    boundaries = mesh["boundary"]
    reader.set_active_time_value(active_time_step)

    # Adding more data to the mesh:
    internal_mesh.cell_data['U_magnitude'] = np.linalg.norm(internal_mesh['U'], axis=1)
    tau = internal_mesh['tau']
    tau_XX = tau[:, 0]
    tau_YY = tau[:, 3]
    internal_mesh.cell_data['N_11']  = tau_XX - tau_YY
    internal_mesh.cell_data['shear'] = tau[:, 1]          # tau_XY


    #printing
    if print_output:
        print("Availbale fields at the active time step: ", mesh.array_names)
        print(f"All patch names: {reader.patch_array_names}")
        print(f"All patch status: {reader.all_patch_arrays_status}")
        print(f"Boundaries patches: {boundaries.keys()} \n")
        print(f"Available Time Values: {reader.time_values} \nSelected: ", active_time_step)
        print("Stress tensor (tau) exists and is in shape of ", internal_mesh['tau'].shape)
        print("N_11 exists and is in shape of ", internal_mesh.cell_data['N_11'].shape)
        print("Cell Data:")
        print(internal_mesh.cell_data)
        print("\nPoint Data:")
        print(internal_mesh.point_data)
 
    if return_boundaries:
        return internal_mesh, boundaries
    return internal_mesh

def plot_all_Ellipoids_in_OF (OpenFOAM_mesh, list_of_Ellipsoids, plotter = False, save_figure = False, scalars = 'U', 
                              show_plotter = False, file_name = 'screenshot', OF_opacity= 0.3, Ellipsoids_opacity= 1):
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
    
    plotter.add_mesh(OpenFOAM_mesh, scalars=OpenFOAM_mesh[scalars], opacity=OF_opacity)

    if save_figure :    plotter.screenshot(file_name)
    if show_plotter:    plotter.show()

def add_pathes_to_Plotter (OpenFOAM_mesh, list_of_Ellipsoids, plotter, save_figure = False, scalars = 'U', file_name = 'all_pathes', OF_opacity= 0.3, Ellipsoids_opacity= 1):
    """This plots all tracking pathes using the `Ellipsoid.add_path_to_Plotter`.
    
    pv.mesh, list[Ellipsods], pv.Plotter() -> -   """

    for i in range( len(list_of_Ellipsoids) ):
        list_of_Ellipsoids[i].add_path_to_Plotter(plotter, Ellipsoids_opacity)

    plotter.add_mesh(OpenFOAM_mesh, scalars=OpenFOAM_mesh[scalars], opacity=OF_opacity)
    #plotter.show()

    if save_figure: plotter.screenshot(file_name)  

