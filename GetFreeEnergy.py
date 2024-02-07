# this script calcualtes the free energy along the reaction coordinate from the .bxd file

import sys
from math import log
from tqdm import tqdm

# box object makes it easier to keep track of collisions on boundaries as a trajectory can cycle through the CV multiple times
class box():
    def __init__(self,lower_bound,upper_bound):
        self.lower_bound=lower_bound
        self.upper_bound=upper_bound
        self.PT_lower=[] # passage times between collisions on lower boundary of box
        self.PT_upper=[] # passage times between collisions on upper boundary of box
        self.k_down=0 # rate constant for passage to box below
        self.k_up=0 # rate constant for passage to box above

    # calculate rate constants for passage to neighbouring boxes
    def get_ks(self):
        self.k_down=len(self.PT_lower)/sum(self.PT_lower)
        self.k_up=len(self.PT_upper)/sum(self.PT_upper)

    def addPTLower(self,time):
        self.PT_lower.append(time)

    def addPTUpper(self,time):
        self.PT_upper.append(time)

# get the boundarieso of the box of a particular line in the .bxd file
def GetCurrentBounds(CV,bounds):
    for i in range(len(bounds)-1):
        if bounds[i] <= CV and bounds[i+1] > CV:
            return bounds[i],bounds[i+1]
            
bxd_file=sys.argv[1]
threshold = int(sys.argv[2]) # the decorrelation threshold
T = float(sys.argv[3]) # temperature of the simulation
lower_limit=float(sys.argv[4]) # lowest value of the CV the user wants in the free energy profile
upper_limit=float(sys.argv[5]) # highest value of the CV the user wants in the free energy profile

kB = 1.380649e-23
NA = 6.02214076e23

raw_lines=[]

with open(bxd_file) as file:
    for line in file:
        raw_lines.append(line.rstrip())
file.close()
        
bound_values=[]
collision_lines=[]

for line in raw_lines:
    frags=line.split()
    if len(frags)==6: # this indicates that a collision happened
        bound_hit = frags[2]
        lower_bound=float(frags[3])
        upper_bound=float(frags[4])

        if lower_bound >= lower_limit and upper_bound <= upper_limit:
            collision_lines.append(line) # fill an array of lines where a collision happened
        
            if bound_hit not in bound_values: 
                bound_values.append(bound_hit) # fill an array of the boundaries involved in the run that the user wants the free energy for
                
bound_values=sorted(bound_values)
boxes=[]

for i in range(len(bound_values)-1):
    boxes.append(box(bound_values[i],bound_values[i+1]))

# get the boundaries of the box the trajectory is in at the start of the .bxd file
frags_0=collision_lines[0].split()
b=float(frags_0[3])
lower_bound=str(round(b,3))

# remove trailing zeros
if lower_bound[-1]=='0':
    lower_bound=lower_bound[:-1]
if lower_bound[-1]=='0':
    lower_bound=lower_bound[:-1]      

# timestep value of first collision
time_0=int(frags_0[0])

for box in boxes:
    if box.lower_bound==lower_bound:
        current_box=box
        break

# loop through collision lines and add the collision time to the relevent boundary of the box the trajectory was in when that collision happened
for i in tqdm(range(1,len(collision_lines))):
    frags=collision_lines[i].split()
    b=float(frags[3])
    lower_bound=str(round(b,3))
    if lower_bound[-1]=='0':
        lower_bound=lower_bound[:-1]
    if lower_bound[-1]=='0':
        lower_bound=lower_bound[:-1]

    boundary_hit = frags[2]

    new_box=frags[5]
    time_1=int(frags[0])

    if new_box=='1.0':
        for box in boxes:
           if box.lower_bound==lower_bound:
               current_box=box
               break

    if new_box == '-1.0':
        PT=time_1-time_0
        if PT > threshold:
            if boundary_hit == current_box.lower_bound:
                current_box.addPTLower(PT)
            elif boundary_hit == current_box.upper_bound:
                current_box.addPTUpper(PT)

    time_0=time_1

for box in boxes:
    box.get_ks()

# get equillibrium constants of passage between each box and its neighbours
K_list=[]
G_list=[0]
for i in range(1,len(boxes)):
    K=boxes[i].k_down/boxes[i].k_up
    K_list.append(K)
    
# calculate free energy of each box from the equillibrium constants
for i in range(len(K_list)):
    G=(-1*kB*T*NA*log(K_list[i]))/1000.
    G_list.append(G)
    
output_lines=[]
final_G=[]

# add up the box to box free energies to get the cumulative free energy profile
G_total=0
for i in range(len(G_list)):
    G_total+=G_list[i]
    final_G.append(G_total)

# adjust everything so that the lowest free energy is 0
min_G=min(final_G)
for i in range(len(final_G)):
    final_G[i]-=min_G

# write the free energy profile to an output file
for i in range(1,len(bound_values)):
    CV=(float(bound_values[i])+float(bound_values[i-1]))/2.
    output_lines.append(str(CV) + ' ' + str(final_G[i-1]))
    
with open('free_energy.txt','w') as output_file:
    for line in output_lines:
        output_file.write(line + "\n")

output_file.close()

print("Done, output written to free_energy.txt in units of kJ mol-1.")


