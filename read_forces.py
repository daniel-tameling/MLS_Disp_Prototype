import math;


forces_ref= {}
f = file("./input/Ref4000atoms.input")
linenumber = 1
for l in f:
    if linenumber >= 10:
        linesplit = l.split()
        forces_ref[linesplit[0]] = [ (float(linesplit[8])), (float(linesplit[9])), (float(linesplit[10]))]
    linenumber += 1

num_atoms = len(forces_ref)

forces= {}
f = file("output_initial_forces.txt")
linenumber = 1
count = 0
for l in f:
    if linenumber >= 2 and count < num_atoms:
        linesplit = l.split()
        forces[linesplit[0]] = [ (float(linesplit[1])), (float(linesplit[2])), (float(linesplit[3]))]
        count += 1
    linenumber += 1

errsum = 0.
fsum = 0.
for i in forces:
    errsum += math.sqrt( (forces[i][0] - forces_ref[i][0])**2 + (forces[i][1] - forces_ref[i][1])**2 + (forces[i][2] - forces_ref[i][2])**2 ) / math.sqrt( (forces_ref[i][0])**2 + (forces_ref[i][1])**2 + (forces_ref[i][2])**2 )
    
print errsum/num_atoms 
