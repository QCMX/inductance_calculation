mlsi_path = 'C:\\Users\\jgr\\Documents\\3D_MLSI_VK'
data_folder = 'data/'

max_points_per_polygon = 32 # empirical estimation, to warn user if exceeded

import gdspy
import subprocess
import os
import numpy as np
import time


def convert_file(filename, cell_name='TOP', lmbd=0.1, conds=None, ah=1, ahb=0.25, find_terminals=False, verbose=True, res=0.1):
    '''
    Convert a GDS file to an MLSCS file.
    Parameters:
    -----------
    filename : str
        The name of the GDS file to convert.
    cell_name : str
        The name of the cell to extract from the GDS file.
    lmbd : float
        The London penetration depth (in um) to use in the MLSCS file. Default is 0.1 um.
    nc : int
        The number of conductors to use in the MLSCS file.  Default is 2.
    conds : list of list of float
        The conductors to use in the MLSCS file. Each conductor is a list of two floats (in um): bottom and top.
        Default is None, putting everyone to [0, 0.1].
    ah: float
        The mesh size parameter for the MLSCS file. Default is 1.0
    ahb: float
        The mesh size parameter for the MLSCS file. Default is 0.25. No idea what it does.
    find_terminals : bool
        Whether to find terminals in the GDS file. Terminals are assumed to be paths on layer 99.
        Default is False.
    res: float
        The resolution for the gds to mlscs conversion. Default is 0.1 um.
    '''
    filename_noext = filename.split('.')[0]
    if verbose:
        print("Converting {} to {}.mlscs...".format(filename, filename_noext))
    lib = gdspy.GdsLibrary()
    test = lib.read_gds(data_folder + filename)
    cell_names = list(test.cells.keys())
    if cell_name not in cell_names:
        raise ValueError(f"Cell '{cell_name}' not found in GDS file.")
    test = test.cells[cell_name]
    polygons = test.get_polygons(by_spec=True)
    if find_terminals:
        paths = test.get_paths()
        paths = [path for path in paths if path.layers[0] == 99]
        paths = np.array([path.points for path in paths])
        paths = np.round(paths/res, 0)*res
        n_terms = len(paths)
        if verbose:
            print(n_terms, "terminals found.")
    keys = polygons.keys()
    layerlist = list(keys)
    for i, layer in enumerate(layerlist):
        if layer[0] == 99:
            layerlist.remove(layer)
    nc = len(layerlist)
    if verbose:
        print(nc, "layers found. Creating one conductor per layer.")

    if conds is None:
        conds = [[0, 0.1] for _ in range(nc)]

    header = write_header(filename_noext, lmbd=lmbd, nc=nc, conds=conds, ah=ah, ahb=ahb)
    with open(f"{data_folder}{filename_noext}.mlscs", "w") as f:
        f.writelines(header)
    
    add_to_file = ""
    # Add terminals if there are any, assuming they are in pairs
    # Make sure it works!
    if find_terminals:
        if n_terms % 2 == 0:
            for i in range(0, n_terms//2):
                add_to_file += f"tp {2*i+1}->{2*i+2}\n"
    i_term = 1

    max_points = 0
    # Add polygons
    for i_cond, layer in enumerate(layerlist):
        poly_layer = ""
        polygons_layer = polygons[layer]
        for polygon in polygons_layer:            
            points = np.array(polygon)
            points = np.round(points/res, 0)*res
            if len(points) > max_points:
                max_points = len(points)
            if verbose:
                print('Writing polygon with {} points on layer {}.'.format(len(points), layer))
            if not(check_poly_orientation(points)):
                points = change_poly_orientation(points)
                if verbose:
                    print("Changed polygon orientation to clockwise.")
            paths_to_add = []
            for path in paths:
                fp = find_path(path, points)
                if fp >= 0:
                    if verbose:
                        print("Found terminal path in polygon on layer {}.".format(layer))
                    paths_to_add.append(fp)
            n_terms_in_poly = len(paths_to_add)
            paths_to_add = sorted(paths_to_add)


            for i, point in enumerate(points):
                if i < len(points)-1:
                    if not(point[0] == points[i+1][0] and point[1] == points[i+1][1]):
                        poly_layer += f"ell {i_cond} 0 {point[0]} {point[1]} {points[i+1][0]} {points[i+1][1]}\n"
                else:
                    poly_layer += f"ell {i_cond} 0 {point[0]} {point[1]} {points[0][0]} {points[0][1]}\n"
                if n_terms_in_poly > 0:
                    if i in paths_to_add:
                        if verbose:
                            print("{} is a terminal path, adding terminal {}.".format(i, i_term))
                        poly_layer = poly_layer[:-1] + f"  t {i_term}\n"
                        i_term += 1
                
            add_to_file += poly_layer
    with open(f"{data_folder}{filename_noext}.mlscs", "a") as f:
        f.writelines(add_to_file)   
    if verbose:
        print(f"Conversion finished. File saved as {data_folder}{filename_noext}.mlscs")
        


def check_poly_orientation(polygon):
    '''
    Check the orientation of a polygon.
    Parameters:
    -----------
    polygon : np.ndarray
        The polygon to check.
    Returns:
    --------
    str
        'cw' if the polygon is clockwise, 'ccw' if the polygon is counter-clockwise.
    '''
    area = 0
    n = len(polygon)
    for i in range(n):
        j = (i + 1) % n
        area += polygon[i][0] * polygon[j][1]
        area -= polygon[j][0] * polygon[i][1]
    if area > 0:
        return True
    else:
        return False

def change_poly_orientation(polygon):
    '''
    Change the orientation of a polygon.
    Parameters:
    -----------
    polygon : np.ndarray
        The polygon to change.
    Returns:
    --------
    np.ndarray
        The polygon with the orientation changed.
    '''
    return polygon[::-1]

def write_header(name, nc, lmbd, conds, ah, ahb):
    assert len(conds) == nc, "Number of conductors must match 'nc'"
    header = f'''cc {name}
pb=0
out=1
usetriangle
ah={ah}
ahb={ahb}
avg=off
manualambda
nc={nc}
lmbd={lmbd}
'''
    for i, cond in enumerate(conds):
        header += f"cond {i} " + " ".join(map(str, cond)) + "\n"
    return header

def find_path(path, polygon):
    '''
    Find the path in the polygon.
    Parameters:
    -----------
    path : np.ndarray
        The path to find.
    polygon : np.ndarray
        The polygon to search.
    '''
    for i in range(len(polygon)-1):
        if np.array_equal(path[0], polygon[i]) and np.array_equal(path[1], polygon[i+1]):
            return i
        elif np.array_equal(path[1], polygon[i]) and np.array_equal(path[0], polygon[i+1]):
            return i
    if np.array_equal(path[0], polygon[-1]) and np.array_equal(path[1], polygon[0]):
        return len(polygon)-1
    elif np.array_equal(path[1], polygon[-1]) and np.array_equal(path[0], polygon[0]):
        return len(polygon)-1
    return -1


def check_mlscs(name):
    '''
    Run 3DMLSI. Nice to check that the .mlscs file is correct.
    Parameters:
    -----------
    name : str
        The name of the MLSCS file (without extension).
    Returns:
    --------
    '''
    print("Running 3DMLSI...")
    result = subprocess.call([mlsi_path + '\\wpm.exe', data_folder + name + '.mlscs'])

def run_upm(name, verbose=True):
    if verbose:
        print("Running UPM...")
    result = subprocess.call([mlsi_path + '\\upm.exe', data_folder + name + '.mlscs'])
    if result == 0:
        if verbose:
            print("UPM finished successfully. Output saved in {}.upm".format(name))
        return result
    else:
        if verbose:
            print("UPM finished with errors.")
        return result

def run_mlw(name, verbose=True):
    if verbose:
        print("Running MLW...")
    result = subprocess.call([mlsi_path + '\\mlw.exe', data_folder + name + '.upm'])
    if result == 0:
        if verbose:
            print("MLW finished successfully. Output saved in {}.out".format(name))
        return result
    else:
        if verbose:
            print("MLW finished with errors.")
        return result

def simulate_inductance(name, verbose=True, recalc=False):
    '''
    Run the MLW tool. Creates the .upm file if it does not exist.
    Save the results in a .out file.
    Parameters:
    -----------
    name : str
        The name of the MLSCS file (without extension). 
    verbose : bool
        Whether to print progress messages. Default is True.
    recalc : bool
        Whether to recalculate the inductance even if the .out file already exists. Default is False.
    Returns:
    --------
    '''
    files = os.listdir(data_folder)
    if recalc:
        if verbose:
            print('Recalculation requested...')
        run_upm(name, verbose=verbose)
    else:
        if name + '.upm' in files:
            if verbose:
                print('upm file found, proceeding to MLW...')
        else:
            if verbose:
                print('upm file not found, running UPM first...')
            run_upm(name, verbose=verbose)
    result = run_mlw(name, verbose=verbose)
    return result

def get_inductance_matrix(name):
    '''
    Get the inductance matrix from the MLSCS output file.
    Parameters:
    -----------
    name : str
        The name of the MLSCS file (without extension).
    Returns:
    --------
    dict
        A dictionary with the inductance matrix. The keys are tuples (i, j) and the values are the inductance values (in pH).

    '''
    with open(data_folder + name + '.out', 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'Inductance matrix' in line:
            start = lines.index(line) + 2
            break
    inductance = []
    for line in lines[start:]:
        if line.strip() == '':
            break
        inductance.append([float(x) for x in line.split()])
    inductance = np.array(inductance)[:, :3]
    inductance_dict = {}
    for i in range(len(inductance)):
        inductance_dict[(int(inductance[i, 0]), int(inductance[i, 1]))] = inductance[i, 2]
    return inductance_dict

def reset_sample(name, verbose=True):
    '''
    Remove all files related to a sample.
    Parameters:
    -----------
    name : str
        The name of the sample (without extension).
    Returns:
    --------
    '''
    files = os.listdir(data_folder)
    for file in files:
        if file.startswith(name + '.') and file.split('.')[-1] not in ['gds', 'GDS']:
            os.remove(data_folder + file)
            if verbose:
                print('Removed file:', file)
    return


def move_conductor(mlscsname, conductor, dx=0, dy=0, verbose=False):
    '''
    Move a conductor directly in the MLSCS file.
    Create a new MLSCS file with the moved conductor.
    Parameters:
    -----------
    mlscsname : str
        The name of the MLSCS file (without extension).
    conductor : int
        The conductor to move (0-indexed).
    dx : float
        The distance to move in the x direction (in um). Default is 0.
    dy : float
        The distance to move in the y direction (in um). Default is 0.
    verbose : bool
        Whether to print progress messages. Default is False.
    Returns:
    --------
    '''
    with open(data_folder + mlscsname + '.mlscs', 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('ell ' + str(conductor)):
            elts = line.split(' ')
            if dx != 0:
                elts[3] = str(float(elts[3]) + dx)
                elts[5] = str(float(elts[5]) + dx)
            if dy != 0:
                elts[4] = str(float(elts[4]) + dy)
                elts[6] = str(float(elts[6]) + dy)
            newline = ' '.join(elts)
            lines[i] = newline
    if verbose:
        print('Moved conductor %d by (%.1f, %.1f) µm' % (conductor, dx, dy))
    with open(data_folder + mlscsname + '_moved.mlscs', 'w') as f:
        f.writelines(lines)
    return

def move_conductor_and_simulate(mlscsname, conductor, dx=0, dy=0, verbose=True, recalc=False, del_files=False):
    '''
    Move a conductor directly in the MLSCS file and run the simulation.
    Create a new MLSCS file with the moved conductor.
    Parameters:
    -----------
    mlscsname : str
        The name of the MLSCS file (without extension).
    conductor : int
        The conductor to move (0-indexed).
    dx : float
        The distance to move in the x direction (in um). Default is 0.
    dy : float
        The distance to move in the y direction (in um). Default is 0.
    verbose : bool
        Whether to print progress messages. Default is True.
    recalc : bool
        Whether to recalculate the inductance even if the .out file already exists. Default is False.
    del_files : bool
        Whether to delete the generated files after the simulation. Default is False.
    Returns:
    --------
    dict
        A dictionary with the inductance matrix. The keys are tuples (i, j) and the values are the inductance values (in pH).
    '''
    move_conductor(mlscsname, conductor, dx=dx, dy=dy, verbose=verbose)
    newname = mlscsname + '_moved'
    result = simulate_inductance(newname, verbose=verbose, recalc=recalc)
    if result != 0:
        print("Error in simulation. Returning None.")
        with open(data_folder + newname + '.upm.mlw.log', 'r') as f:
            lines = f.readlines()
            print("Last 10 lines of the log file:")
            for line in lines[-10:]:
                print(line.strip())
        return None
    inductance = get_inductance_matrix(newname)
    if verbose:
        print("Inductance matrix (in pH):")
        for (i, j), value in inductance.items():
            print(f"({i}, {j}): {value}")
    if del_files:
        reset_sample(newname, verbose=verbose)
    return inductance


def get_inductance_from_gds(name, check_conv=False, cell_name='TOP', lmbd=0.1, conds=None, ah=1, ahb=0.25, verbose=True, res=0.1, recalc=False):
    '''
    Convert a GDS file to an MLSCS file and run 3DMLSI software to get inductances.
    Parameters:
    -----------
    name : str
        The name of the GDS file to convert (without .gds). Must be in the data directory.
    check_conv : bool
        Whether to run 3DMLSI to check the conversion. Default is False.
    cell_name : str
        The name of the cell to extract from the GDS file. Default is 'TOP'.
    lmbd : float
        The London penetration depth (in um) to use in the MLSCS file. Default is 0.1 um.
    conds : list of list of float
        The conductors to use in the MLSCS file. Each conductor is a list of two floats (in um): bottom and top.
        You can add a third float to that list, defining a London length for that conductor.
        Default is None, putting everyone to [0, 0.1].
    ah: float
        The mesh size parameter for the MLSCS file. Default is 1.0
    ahb: float
        The mesh size parameter for the MLSCS file. Default is 0.25. No idea what it does.
    verbose : bool
        Whether to print progress messages. Default is True.
    res: float
        The resolution for the gds to mlscs conversion. Default is 0.1 um.
    recalc : bool
        Whether to recalculate the inductance even if the .out file already exists. Default is False.
    Returns:
    --------
    dict
        A dictionary with the inductance matrix. The keys are tuples (i, j) and the values are the inductance values (in pH).
    '''
    t0 = time.time()
    convert_file(name + '.gds', find_terminals=True, cell_name=cell_name, lmbd=lmbd, conds=conds, ah=ah, ahb=ahb, verbose=verbose, res=res)
    t_conv = time.time() - t0
    if check_conv:
        check_mlscs(name)
    result = simulate_inductance(name, verbose=verbose, recalc=recalc)
    if result != 0:
        print("Error in simulation. Returning None.")
        with open(data_folder + name + '.upm.mlw.log', 'r') as f:
            lines = f.readlines()
            print("Last 10 lines of the log file:")
            for line in lines[-10:]:
                print(line.strip())
        return None
    t_sim = time.time() - t0 - t_conv
    if verbose:
        print(f"Conversion took {t_conv:.2f} s, simulation took {t_sim:.2f} s.")
    inductance = get_inductance_matrix(name)
    if verbose:
        print("Inductance matrix (in pH):")
        for (i, j), value in inductance.items():
            print(f"({i}, {j}): {value}")
    return inductance

