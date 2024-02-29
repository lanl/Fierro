#!/usr/bin/python3

import os
import os.path
import sys
import math

executable = "./../../build-fierro-openmp/bin/fierro-parallel-explicit"
sim_input = "Solver-Inputs/SGH_Sedov_12x12x12.yaml"

# Run simulation
os.system(executable + ' ' + sim_input)


# Functions for reading results from vtk file
def extract_vector_data(filename, keyword):
    data = []
    found_keyword = False

    with open(filename, 'r') as file:
        for line in file:
            if found_keyword:
                if line.strip():  # If the line is not blank
                    # Split each line into fields.
                    fields = line.split()
                    xx = float(fields[0])
                    yy = float(fields[1])
                    zz = float(fields[2])
                    # Group and append to results.
                    vec = [xx, yy, zz]
                    data.append(vec)
                else:
                    break  # Exit the loop when encountering a blank line
            elif keyword in line:
                found_keyword = True

    return data

def extract_scalar_data(filename, keyword):
    data = []
    found_keyword = False

    with open(filename, 'r') as file:
        for line in file:
            if found_keyword:
                if line.strip():  # If the line is not blank
                    fields = line.split()
                    if fields[0] != "LOOKUP_TABLE":
                        data.append(float(fields[0]))
                else:
                    break  # Exit the loop when encountering a blank line
            elif keyword in line:
                found_keyword = True
                # file.next() # skip first line "LOOKUP_TABLE"
 
    return data



def percent_difference_scalars(array1, array2):
    if len(array1) != len(array2):
        raise ValueError("Arrays must have the same length")

    percent_diff = []
    for i in range(len(array1)):
        diff = array2[i] - array1[i]
        percent_diff.append((diff / array1[i]) * 100 if array1[i] != 0 else 0)

    return percent_diff

def percent_difference_vectors(array1, array2):
    if len(array1) != len(array2):
        raise ValueError("Arrays must have the same length")

    percent_diff = []
    for i in range(len(array1)):
        if len(array1[i]) != 3 or len(array2[i]) != 3:
            raise ValueError("Subarrays must have length 3")
        
        diff = [array2[i][j] - array1[i][j] for j in range(3)]
        percent_diff.append([(diff[j] / array1[i][j]) * 100 if array1[i][j] != 0 else 0 for j in range(3)])

    return percent_diff

def magnitude(array):
    mag = math.sqrt(sum(x**2 for x in array))
    return mag


# Read ground truth results for sedov
GT_filename = "standard-results/SGH/Sedov_12x12x12/vtk/data/VTK0.vtk"

velocity_keyword = "VECTORS velocity float"
position_keyword = "POINTS 2197 float"
SIE_keyword = "SCALARS SIE float 1"
density_keyword = "SCALARS element_density float 1"

GT_positions = extract_vector_data(GT_filename, position_keyword)
GT_velocities = extract_vector_data(GT_filename, velocity_keyword)
GT_SIE = extract_scalar_data(GT_filename, SIE_keyword)
GT_densities = extract_scalar_data(GT_filename, density_keyword)


# Read simulation results
results_filename = "vtk/data/VTK0.vtk"

results_positions = extract_vector_data(results_filename, position_keyword)
results_velocities = extract_vector_data(results_filename, velocity_keyword)
results_SIE = extract_scalar_data(results_filename, SIE_keyword)
results_densities = extract_scalar_data(results_filename, density_keyword)


density_diff = percent_difference_scalars(GT_densities, results_densities)

print("Density difference: ")
print(magnitude(density_diff))

