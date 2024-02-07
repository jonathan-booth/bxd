from openmm import *
from openmm.app import *
from openmm.unit import *
from sys import exit
from CheckParameters import CheckAllParamaters, GetSelectionInfoForRoG

# Parameters for the user to set ------------------------------------------------

# I/O 
pdb_name='6-GLY.pdb'
output_name='test'

# Platform (can be CPU, CUDA or OpenCL)
platform_choice='CPU'

# Simulation params
n_timesteps=100000000
dcd_write_frequency=1000
standard_output_write_frequency=10000
max_minimisation_steps=2000
forcefield = ForceField('charmm36.xml', 'implicit/GBn2.xml')

# BXD params, currently CV = rog or distance
CV_type='rog' 
atom_selection='backbone'
bounds=[0.31,0.32,0.33,0.34,0.35,0.4,0.45,0.5,0.55,0.575,0.6,0.61,0.62,0.63]
threshold=200
initial_direction='up'
bxd_write_frequency=200

# End of parameter section ----------------------------------------------------

CV_type=CV_type.lower()
initial_direction=initial_direction.lower()

# Parameter sense checks
params_to_check={'pdb_name':pdb_name,
                'platform_choice':platform_choice,
                 'n_timesteps':n_timesteps,
                 'dcd_write_frequency':dcd_write_frequency,
                 'standard_output_write_frequency':standard_output_write_frequency,
                 'max_minimisation_steps':max_minimisation_steps,
                 'CV_type':CV_type,
                 'atom_selection':atom_selection,
                 'bounds':bounds,
                 'initial_direction':initial_direction,
                 'bxd_write_frequency':bxd_write_frequency
                 }

CheckAllParamaters(params_to_check)

# set the upper and lower limits of the bounded area of phase space 
# this is so we know if and when the trajectory enters the bounded region and BXD should begin.
lower_limit=bounds[0]
upper_limit=bounds[-1]

# OpenMM's custom integrator needs the variables to be numeric rather than strings
if initial_direction == 'down':
    initial_direction=-1

elif initial_direction == 'up':
    initial_direction=1

# an OpenMM custom integrator to implement the BXD algorithm.
# integration is done with the velocity verlet algorithm with velocity reflection if a boundary is hit
# an Anderson thermostat is used
# OpenMM custom integrators don't allow integration on the loops so they are labelled with letters in the comments
def BXD_anderson_vv(temperature=298*kelvin, collision_rate=91.0/picoseconds, timestep=1.0*femtoseconds):

    integrator = CustomIntegrator(timestep)
    integrator.setIntegrationForceGroups({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}) # force 31 is reserved for the CV
    integrator.addPerDofVariable('x_old',0) # old positions, for integration
    integrator.addPerDofVariable('v_old',0) # old velocities, for integration
    integrator.addPerDofVariable('f_old',0) # old forces, for integration
    integrator.addGlobalVariable('hit_count',0) # number of collisions against boundary in direction of travel
    integrator.addGlobalVariable('BXD_started',0) # has the trajectory moved into the bounded region of phase space, 0 if false 1 if true
    integrator.addGlobalVariable('CV',0) # the value of collective variable
    integrator.addGlobalVariable('collision',0) # have a boundary been crossed on the current timestep, 0 if false 1 if true
    integrator.addGlobalVariable('boundary_hit',0) # which boundary was hit, 1 if upper , -1 if lower
    integrator.addGlobalVariable('threshold',threshold) # how many collisions must happen before the trajectory moves into the next box
    integrator.addGlobalVariable('direction',initial_direction) # is the trajectory going up or down the CV to begin with
    integrator.addGlobalVariable('first_time_in_box',0) # has the trajectory just crossed into a new box? 1 if true, -1 if false.
    integrator.addGlobalVariable('lower_bound',-1) # value of lower boundary of current box
    integrator.addGlobalVariable('upper_bound',-1) # value of upper boundary of current box
    integrator.addGlobalVariable('upper_bound_index',0) # position of current upper bound in list of boundaries, for tracking how far along the list the trajectory is

    # make an array of the boundaries, called bounds, from the boundaries parameter set by the user
    boundaries=Discrete1DFunction(bounds)
    integrator.addTabulatedFunction("B",boundaries)

    integrator.addGlobalVariable('lower_limit',lower_limit)
    integrator.addGlobalVariable('upper_limit',upper_limit)

    # vv part of integration
    integrator.addUpdateContextState()
    integrator.addComputePerDof('x_old','x')
    integrator.addComputePerDof('v_old','v')
    integrator.addComputePerDof('f_old','f31')

    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.addComputePerDof("x", "x+dt*v")
    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.addUpdateContextState()

    integrator.addComputeGlobal('CV','energy31^0.5')

    # check if trajectory has moved into one of the boxes so BXD can start
    integrator.beginIfBlock('BXD_started = 0') # A open
    integrator.beginIfBlock('CV < upper_limit') # B open
    integrator.beginIfBlock('CV > lower_limit') # C open

    # if BXD hasn't started but the trajectory is in the boudned region of space, find the boudnaries of the current box
    # also set the BXD_started variable to 1
    integrator.beginWhileBlock('B(upper_bound_index)<CV') # P open
    integrator.addComputeGlobal('upper_bound_index','upper_bound_index+1')
    integrator.endBlock() # P close

    integrator.addComputeGlobal('lower_bound','B(upper_bound_index-1)')
    integrator.addComputeGlobal('upper_bound','B(upper_bound_index)')
    integrator.addComputeGlobal('BXD_started','1')

    integrator.endBlock() # C close
    integrator.endBlock() # B close
    integrator.endBlock() # A close

    integrator.addComputeGlobal('collision','0')

    # if BXD has started, check if this integration step took the trajectory accross a boundary
    integrator.beginIfBlock('BXD_started = 1') # D
    integrator.addComputeGlobal('collision','0')

    # assume that this is not the first timestep in the box
    # if it is then this is corrected later
    # if this is the first timestep for which BXD is started then this is set to -1 even though this is incorrect
    # first_time_in_box is used later to decide, after a collision is detected, whether to reflect the velocities or allow the trajectory to pass into the next box
    integrator.addComputeGlobal('first_time_in_box','-1')

    # check to see if the trajectory has gone outside a box by seeing
    # this is done by checking if the value of the CV has gone above the current upper bound or below the current lower bound
    # if either is true, set collision variable to 1 and chqange boundary_hit accordingly
    integrator.beginIfBlock('CV > upper_bound') # E open
    integrator.addComputeGlobal('collision','1')
    integrator.addComputeGlobal('boundary_hit','1')
    integrator.endBlock() # E close

    integrator.beginIfBlock('CV < lower_bound') # F open
    integrator.addComputeGlobal('collision','1')
    integrator.addComputeGlobal('boundary_hit','-1')
    integrator.endBlock() # F close

    # if a boundary was hit then count this as a collision for comparison to the threshold value
    # when counting the number of hits we only care if the hit was in the direction of travel
    # e.g. if going up through the CV, only count the hit if it was on the upper boundary
    # reflection is not done yet, but collision is set to 1 and first_time_in_box is set to -1 so we will now to do this later
    integrator.beginIfBlock('collision = 1') # G open

    integrator.beginIfBlock('boundary_hit=-1') # H open
    integrator.beginIfBlock('direction=-1') # I open
    integrator.addComputeGlobal('hit_count','hit_count+1')
    integrator.endBlock() # I close
    integrator.endBlock() # H close

    integrator.beginIfBlock('boundary_hit=1') # J open
    integrator.beginIfBlock('direction=1') #  K open
    integrator.addComputeGlobal('hit_count','hit_count+1')
    integrator.endBlock() # K close
    integrator.endBlock() # J close

    # if the number of hits on the relevent boundary has reached the threshold then allow the trajectory to pass into the next box
    # the trajectory passes through because first_time_in_box is set to 0 so reflection is not carried out later.
    integrator.beginIfBlock('hit_count=threshold') # L open

    integrator.addComputeGlobal('hit_count','0')
    integrator.addComputeGlobal('first_time_in_box','1')

    # if the system is moving up through the CV then set the new boundaries to be those of the next box up
    integrator.beginIfBlock('direction=1') # M open
    integrator.addComputeGlobal('upper_bound_index','upper_bound_index+1')
    integrator.addComputeGlobal('upper_bound','B(upper_bound_index)')
    integrator.addComputeGlobal('lower_bound','B(upper_bound_index-1)')
    integrator.endBlock() # M close

    # if the system is moving down through the CV then set the new boundaries to be those of the next box down
    integrator.beginIfBlock('direction=-1') # N open
    integrator.addComputeGlobal('upper_bound_index','upper_bound_index-1')
    integrator.addComputeGlobal('upper_bound','B(upper_bound_index)')
    integrator.addComputeGlobal('lower_bound','B(upper_bound_index-1)')
    integrator.endBlock() # N close

    # if the trajectory has just moved into the highest box then reverse the direction so it starts going down
    integrator.beginIfBlock('upper_bound=upper_limit') # Q open
    integrator.addComputeGlobal('direction','-1')
    integrator.endBlock() # Q close

    # if the trajectory has just moved into the lowest box then reverse the direction so it starts going up
    integrator.beginIfBlock('lower_bound=lower_limit') # P open
    integrator.addComputeGlobal('direction','1')
    integrator.endBlock() # P close

    integrator.endBlock() # L close

    # if a collision has happened (we are still inside the if collision == 1 loop) then check if the trajectory should be reflected
    # if the trajectory has just entered a new box then don't refect it
    # reflect if the trajectory has not just entered a new box i.e. hit_count < threshold so it must stay in its current box for now
    integrator.beginIfBlock('first_time_in_box=-1') # O open

    # for the previous timestep before the trajectory hit the boundary get the projection of the velocity along the CV
    integrator.addComputeSum('V_dot_F','v_old*f_old')
    integrator.addComputeSum('F_dot_F','f_old*f_old')
    integrator.addGlobalVariable('scalar_product',0)
    integrator.addComputeGlobal('scalar_product','V_dot_F/F_dot_F')
    integrator.addPerDofVariable('V_proj_F',0)
    integrator.addComputePerDof('V_proj_F','scalar_product*f_old')

    # change the positions back to what they were in the timestep before hte boundary was hit
    integrator.addComputePerDof('x','x_old')
    # reflect the velocities (of the timestep before the boundary was hit) with respect to the CV
    integrator.addComputePerDof('v','v_old-V_proj_F-V_proj_F')

    integrator.endBlock() # O close
    integrator.endBlock() # G close
    integrator.endBlock() # D close

    return integrator


pdb = PDBFile(pdb_name)
nonbondedMethod = CutoffNonPeriodic
nonbondedCutoff = 1*nanometers
constraints = None # constraints are none because the timestep is small enough to capture C-H bond vibration

platform = Platform.getPlatformByName(platform_choice)
T = 300.0

dcdReporter = DCDReporter(output_name+'.dcd', dcd_write_frequency, append=False)
print('Building system...')

topology = pdb.topology
positions = pdb.positions
system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints)

# define the CV
# CV must be calculated as a 'force' for OpenMM custom integrator purposes
if CV_type == 'rog':
    selection_type=GetSelectionInfoForRoG(atom_selection)

    if selection_type == 'by_type':
        chosen_atoms=[]
        chosen_masses=[]

        for atom in topology.atoms():
            if atom.name in atom_selection:
                chosen_atoms.append(atom.index)
                chosen_masses.append(system.getParticleMass(atom.index).value_in_unit(amu))

    elif selection_type == 'by_index':
        chosen_atoms=[]
        chosen_masses=[]

        for atom in topology.atoms():
            if atom.index in atom_selection:
                chosen_atoms.append(atom.index)
                chosen_masses.append(  system.getParticleMass(atom.index).value_in_unit(amu))

    elif selection_type == 'single_word':
        if atom_selection == 'all':
            chosen_atoms=[]
            chosen_masses=[]

            for atom in topology.atoms():
                chosen_atoms.append(atom.index)
                chosen_masses.append(  system.getParticleMass(atom.index).value_in_unit(amu))

        if atom_selection == 'backbone':
            atom_names=['C','CA','N','O']
            chosen_atoms=[]
            chosen_masses=[]

            for atom in topology.atoms():
                if atom.name in atom_names:
                    chosen_atoms.append(atom.index)
                    chosen_masses.append(  system.getParticleMass(atom.index).value_in_unit(amu))

    # RoG is the mean mass weighted distance between each selected atom and the centre of mass of the selected atoms
    RoG = CustomCentroidBondForce(2,"(mass*(distance(g1,g2)^2))/N")

    N_param=sum(chosen_masses)
    RoG.addGlobalParameter("N",N_param)
    RoG.addPerBondParameter("mass")

    # add group that is all the selection of atoms for the calculation of the centre of mass
    RoG.addGroup(chosen_atoms)

    # add group of each atom
    for atom_index in chosen_atoms:
        RoG.addGroup([atom_index])

    # add bond for each atom and whole molecule centre of mass
    for i in range(len(chosen_atoms)):
        RoG.addBond([0,i+1],[chosen_masses[i]])

    RoG.setForceGroup(31)

    system.addForce(RoG)

# this force is just the distance between the two atoms
if CV_type == 'distance':

    dist = CustomCentroidBondForce(2,"distance(g1,g2)")

    dist.addGroup([atom_selection[0]])
    dist.addGroup([atom_selection[1]])

    dist.addBond([0,1])
    dist.setForceGroup(31)
    system.addForce(dist)

system.addForce(AndersenThermostat(T*kelvin, 1/picosecond))

integrator = BXD_anderson_vv()
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)
simulation.reporters.append(dcdReporter)
simulation.minimizeEnergy(maxIterations=max_minimisation_steps)
simulation.currentStep = 0
simulation.reporters.append(StateDataReporter(output_name+'.out',standard_output_write_frequency, step=True, totalEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('output.pdb', 1000000))

for i in range(n_timesteps):

    CV=integrator.getGlobalVariableByName('CV')

    # write some output to the .bxd if a collision happens
    if integrator.getGlobalVariableByName('collision') == 1:
        if integrator.getGlobalVariableByName('boundary_hit') == 1:
            active_boundary=integrator.getGlobalVariableByName('upper_bound')
        elif integrator.getGlobalVariableByName('boundary_hit')==-1:
            active_boundary=integrator.getGlobalVariableByName('lower_bound')

        active_boundary=str(round(active_boundary,3))

        with open(output_name+'.bxd','a') as bxd_file:
            bxd_file.write(str(i) + ' ' + str(CV) + ' ' + active_boundary + ' ' + str(round(integrator.getGlobalVariableByName('lower_bound'),3)) + ' ' + str (round(integrator.getGlobalVariableByName('upper_bound'),3)) + ' ' + str(integrator.getGlobalVariableByName('first_time_in_box'))+ "\n")
        bxd_file.close()

    # if a collision has not happened but it is time to write based on the bxd_write_frequency
    if i % bxd_write_frequency == 0  :
        with open(output_name+'.bxd','a') as bxd_file:
            bxd_file.write(str(i) + ' '+ str(CV) + "\n")
        bxd_file.close()

    simulation.step(1)

# Write file with final simulation state

state = simulation.context.getState(getPositions=True)
with open(output_name+'.pdb', mode="w") as file:
    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)
    print("done")
