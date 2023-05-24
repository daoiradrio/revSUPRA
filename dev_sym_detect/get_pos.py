import os
import numpy as np


n_atoms         = 0
element_numbers = []
coords          = []
center          = np.array([0.0, 0.0, 0.0])
pos             = np.array([])


file = os.sys.argv[1]
with open(file, "r") as infile:
    data = infile.readlines()
    for i, line in enumerate(data):
        line = line.split()
        if i == 0:
            n_atoms = int(line[0])
        else:
            element_numbers.append(int(line[0]))
            new_coords = [
                float(line[1]),
                float(line[2]),
                float(line[3])
            ]
            coords.append(np.array(new_coords))
    coords = np.array(coords)


print(element_numbers)
print(coords)


for coord in coords:
    center = center + coord
center = center / float(n_atoms)


print(center)


pos = coords[0] - center
pos = coords[0] + 1.5*pos


print(pos)
