import os
from sys import exit
from openmm import *
from openmm.app import *
from openmm.unit import *

# checks whether the parameters set by the user in BXD.py are valid.
# returns nothing if all parameters are valid, calls sys.exit() and prints error message otherwise.
def CheckAllParamaters(params):
    if os.path.exists(params['pdb_name']) == False:
        print("Error with pdb_name parameter, file not found")
        exit()

    available_platforms=['CPU','Cuda','OpenCL']

    if params['platform_choice'] not in available_platforms:
        print("Error with platform_choice parameter: ,must be one of the following:")
        print(available_platforms)
        exit()

    IsPositiveInt(params['n_timesteps'],'n_timesteps')
    IsPositiveInt(params['dcd_write_frequency'],'dcd_write_frequency')
    IsPositiveInt(params['standard_output_write_frequency'],'standard_output_write_frequency')
    IsPositiveInt(params['max_minimisation_steps'],'max_minimisation_steps')

    supported_CVs=['rog','distance']

    if params['CV_type'] not in supported_CVs:
        print("Error with CV_type parameter: supported CV_types are as follows:")
        print(supported_CVs)
        exit()

    # get  list of all atom indexes and all atom names that are in pdb file to help us check that atom selection is valid
    pdb = PDBFile(params['pdb_name'])
    topology = pdb.topology
    all_atom_indexes=[atom.index for atom in topology.atoms()]
    all_atom_types=[atom.name for atom in topology.atoms()]

    # Checking the atom selection and bounds parameters is more coplicated so seperate functions are called to do this.
    # Function to check atom selection depends on the choice of CV

    if params['CV_type'] == 'rog':
        CheckAtomSelectionForRoG(params['atom_selection'],all_atom_indexes,all_atom_types)

    # if we are using the distance CV then call a function to check the atom selection is valid for this CV.
    elif params['CV_type'] == 'distance':
        CheckAtomSelectionForDistance(params['atom_selection'],all_atom_indexes)

    CheckBounds(params['bounds'])

    if params['initial_direction'] != 'up' and params['initial_direction'] != 'down':
        print("Error with initial_direction parameter, intitial_direction mut be either 'up' or 'down'")
        exit()

    IsPositiveInt(params['bxd_write_frequency'],'bxd_write_frequency')

    return

# checks if a parameter is a positive integer. Input is n, the value of a parameter, and the name of the parameter for printing the error message
# returns nothing if parameter is a positive integer, calls sys.exit() if not.
def IsPositiveInt(n,param_name):
    if isinstance(n, str) == True or int(n) != n or n < 1:
        print('Error with', param_name, ' parameter: must be positive integer.')
        exit()
    
    return

# checks whether valid atom selections have been made if the distance CV has been chosen.
# atom_selection is an array of the atoms selected.
# all_atom_indexes is an array of the indexes of every atom present in the openMM system.
# OpenMM gives the first atom index 0 and many pdb files start at 1
# returns nothing if the test is passed, calls sys.exit() otherwise
def CheckAtomSelectionForDistance(atom_selection,all_atom_indexes):
    if len(atom_selection)!=2:
        print("Error with atom_selection parameter: when using distance CV, atom_selection must be [x,y] where x and y are atom indexes")
        exit()
 
    for atom in atom_selection:
        try:
            atom = int(atom)

        except ValueError:
            print('Error with atom_selection parameter: Atom indexes must be integers')
            exit()

        if atom < 0:
            print('Error with atom_selection parameter: Atom indexes cannot be < 0')
            exit()

        if atom not in all_atom_indexes:
            print('Error with atom_selection parameter: Atom with index', atom, 'not found in pdb.')
            print('If first atom in pdb file has index 1, remember to subract 1 from your selection as OpenMM atom indexes start at 0.')
            exit()

    if atom_selection[0] == atom_selection[1]:
        print('Error with atom_selection parameter: Atom indexes cannot be the same')
        exit()

    return

# checks if atom selection is valid if the RoG CV has been chosen
# returns nothing is selction is valid, calls sys.exit() if not.
# openMM gives the first atom index 0 and many pdb files start at 1 
def CheckAtomSelectionForRoG(atom_selection,all_atom_indexes,all_atom_types):
    atom_selection_words=['all','backbone']

    if atom_selection in atom_selection_words:
        return

    # if atom selection is done by keyword, check it is valid
    if isinstance(atom_selection,str) == True:
        if atom_selection not in atom_selection_words:
            print('Error with atom_selection parameter: supported single-word selections are as follows:')
            print(atom_selection_words)
            exit()

    # if atom selection not by keyword, check if it is by atom index or atom type
    # if atom selection is all integers, assume user wants to select by index
    # if atom selection is all strings, assume user wants to select by type
    # if neither, call sys.exit() as selection is not valid
    all_ints=True
    all_strings=True

    for item in atom_selection:
        if isinstance(item,str) == False:
            all_strings=False

            if int(item) != item:
                all_ints = False

            if item < 0:
                all_ints=False
        
        else:
            all_ints=False

    if all_ints == True:
        for idx in atom_selection:
            if idx not in all_atom_indexes:
                print('Error with atom_selection parameter: Atom with index', idx, 'not found in pdb.')
                print('If first atom in pdb file has index 1, remember to subract 1 from your selection as OpenMM atom indexes start at 0.')
                exit()

    if all_strings == True:
        for name in atom_selection:
            if name not in all_atom_types:
                print('Error with atom_selection parameter: Atom with name', name, 'not found in pdb.')
                exit()

    if all_ints == False and all_strings == False:
        print("Error with atom_selection parameter: all atoms must be ints if selecting by index, or all strings if selecting by atom names")
        exit()

    return

# checks that boundaries are valid. Returns nothing if they are, calls sys.exit() if not
def CheckBounds(bounds):
    if len(bounds)<2:
        print('Error: at least two boundaries must be provided')
        exit()

    for bound in bounds:
        try:
            bound = float(bound)
        
        except ValueError:
            print('Error with boundaries: all boundaries must be floats or ints')
            exit()

    if sorted(bounds) != bounds:
        print('Error with boundaries: boundaries must be in ascending order')
        exit()

    if len(set(bounds)) != len(bounds):
        print('Error with boundaries: repeated boundary found. All boundaries must be unique.')
        exit()

    return

# if RoG is selected as the CV, find out what kind of atom selection has been chosen
# possible selections are single word, atom index or atom type
def GetSelectionInfoForRoG(atom_selection):
    if isinstance(atom_selection,str) == True:
        selection_type='single_word'
        return selection_type
    
    all_ints=True
    all_strings=True

    for item in atom_selection:

        if isinstance(item,str) == False:
            all_strings=False

            if int(item) != item:
                all_ints = False

            if item < 0:
                all_ints=False
        
        else:
            all_ints=False

    if all_ints == True:
        selection_type='by_index'
    
    if all_strings == True:
        selection_type='by_type'

    return selection_type
