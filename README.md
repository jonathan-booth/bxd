# BXD

## Introduction
This is an implementation of Boxed Molecular Dynamics in OpenMM,[1] an open-source simulation package. BXD is an enhanced sampling method for Molecular Dynamics (MD) simulation. A Collective Variable (CV) is defined which is a single dimensional metric to describe a process of interest. BXD works by placing reflecting boundaries in phase space along the CV which allow the trajectory to explore the free energy landscape along the CV. This allows simulation of processes which are far beyond the timescale of conventional MD simulations, as well as a calculation of the free energy along the CV. See refs [2] to [5] for a more detailed explanation. BXD is implemented as an input file in the OpenMM simulation script, with a custom integrator object carrying out the BXD part of the simulation.

## Installation
1) Ensure you have a recent version of Python installed
https://www.python.org/

2) Install Anaconda or Miniconda 
https://conda.io/projects/conda/en/latest/user-guide/install/index.html

3) Create a conda environment with a chosen name, for example: 
conda create -n BXD

4) Activate the environment: 
conda activate BXD

5) Install OpenMM: 
conda install -c conda-forge openmm

6) Install tqdm:
conda install tqdm

7) Clone this repository via either SSH or HTTPS
git clone $REPOSITORY_PATH

## Usage
1) cd into the bxd directory created by cloning this repository.

2) Place your pdb file in the bxd directory. You can modify BXD.py so that PDBx/mmCIF files can be used, change the line pdb = PDBFile(pdb_name) to pdb = PDBxFile(pdb_name)
3) Alternatively, place the pdb file in any location and pass the path to BXD.py (instructions later).

4) Open BXD.py and set the parameters at the start of the file. There is some skill to setting these parameters. See Appendix A of Jonathan Booth's thesis[6] for details of how to choose the best values, or contact me via Github. The parameters are as follows:

### pdb_name (string or path)
The name of your pdb file. If you have not placed your pdb file in the bxd directory then set this as the path to wherever it is located.

### output_name (string)
The name root of the output files. For example, setting this to 'test' will give test.bxd, test.dcd, test.out etc.

### platform_choice ('CPU','CUDA' or 'OpenCL')
The platform for the simulation. OpenMM supports CPU as well as cuda or openCL for parallelisation. See OpenMM documentation for details:
http://docs.openmm.org/latest/developerguide/03_writing_plugins.html#creating-new-platforms

### n_timesteps (int)
The numper of timestpes (in femtoseconds) for the simulation.

### dcd_write_frequency (int)
The interval in timesteps for writing trajectory snapshots to the dcd file.

### standard_output_write_frequency (int)
The interval in timesteps for writing the system energy and temperature to the .out file. The properties written to this file can be altered by changing the following line in BXD.py:
simulation.reporters.append(StateDataReporter(output_name+'.out',standard_output_write_frequency, step=True, totalEnergy=True, temperature=True))

### max_minimisation_steps (int)
The maximumm number of steps to attempt to minimise the system before simulation starts. To alter or disable minimisation, change the following line in BXD.py: 
simulation.minimizeEnergy(maxIterations=max_minimisation_steps)

### forcefield (string)
The choice of forcefield(s) and whether to use implicit solvent. See OpenMM documentation for more details:
https://github.com/openmm/openmmforcefields

### CV_type ('rog' or 'distance')
What to use for the collective variable. Currently there are two choices:
A) rog - Radius of gyration, the mean mass weighted distances between each selected atom and the centre of mass of the selected atoms
B) dist - distance between two atoms

### atom_selection (string, array of ints or array of strings)
The collection of atoms to be used for calculating the CV.

For distance, this is two atoms and must be in the form of an array containing the indexes of the two atoms. Note that these indexes are those of the OpenMM system not in the pdb file. PDB files often start at index 0 whereas OpenMM starts at index 0. If the pdb file starts at index 1 then subtract 1 from the pdb file indexes of the atoms you want to select. For example, if the pdf file starts at index 1, and you wish to use atoms 1 and 12, set the following:
atom_selection = [0,11]

For rog, there are three ways of selecting atoms:

1) By keyword. Setting the parameter to be a string allows some preset selections. Currently these are 'all' for all atoms, or 'backbone' for the backbone of a peptide, which may be more useful if studying conformational dynamics. For example:
atom_selection = 'backbone'

2) By index. Pass in an array of the atom indexes to be used. Note the above warning from the distance section on the diference in atom numbering convention between pdb files and OpenMM. For example:
atom_selection = [0,4,26,57,66]

3) By type. Pass in an array of strings of atom types to be used, typically from the third column of a pdb file. For example, to use all Nitrogen and terminal Oxygen atoms:
atom_selection = ['N','OXT']

### bounds (array of floats)
An array of the location of boundaries in phase space. Must all be unique and in ascending order. Units are nanometres. For example:
bounds = [0.3,0.4,0.5]

### threshold (int)
How many collisions should occur before the trajectory is allowed through into the next box. Collisions must be against the correct boundary to be counted. For example if the trajectory is moving upwards through the boxes then only collisions against the upper boundary of a box will be counted. 

### initial_direction ('up or 'down')
Should the trajectory move upwards or downwards through the boxes at the start of the simulation. This is not very important, as the trajectory will cycle through the boxes in both directions. I often set this as the direction in which free energy is increasing, to increase the chances of dropping into a different local minimum to that of the starting structure.

### bxd_write_frequency (int)
The interval in timesteps for writing the value of the CV to the bxd file. Note that, regardless of this value, output is written to the bxd file whenever a collision occurs.

4) Run the script:
python BXD.py

5) To calculate the free energy, use the following command:
python GetFreeEnergy.py [bxd_file] [tau] [T] [lower_limit] [upper_limit]

bxd_file is the .bxd file generated by the simulation

tau is the decorrelation time where collisions that occur less than tau timesteps after a previous collision are dropped from the calculation of the rate constants

T is the temperature of the simulation

lower_limit and upper_limit are the values of the CV between which the user wants to calculate the free energy.

## Project status, future work and contributions
This is a work in progress and there may be bugs or omissions. GetFreeEnergy.py could be made more sophisticated, for example it could provide the box to box rate constants or allow for higher resolution free energies. More reaction coordinates could be added to BXD.py and experiments could be carried out to find out how BXD behaves with constraints which may allow for a longer timestep.

Feel free to contact me via github for any questions, issues, ideas or offers of contribution.

## Authors and acknowledgment
Created by Jonathan Booth at the STFC Hartree Centre, with help from the OpenMM Setup input generator (for the non-BXD parts of BXD.py). The author is grateful for the help provided by Peter Eastman at Stanford University via the OpenMM Forums, from Dmitrii Shalashilin at the University of Leeds, and from Dmitry Makhov at the University at Leeds.

## License
This is provided with an MIT licence, see LICENCE.txt for details

## References
[1]https://openmm.org/

[2]https://pubs.acs.org/doi/10.1021/jp9074898

[3]https://pubs.rsc.org/en/content/articlelanding/2016/fd/c6fd00138f

[4]https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/syst.201900024

[5]https://pubs.acs.org/doi/10.1021/ct200011e

[6]https://etheses.whiterose.ac.uk/13101/
