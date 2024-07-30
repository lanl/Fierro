#!/usr/bin/python3

import os
import os.path
import sys
import math
import glob



# Extract data from txt file
def extract_state_data(filename):
    data = []

    # Ensure the file exists before proceeding
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print("The file was not found at the specified location.")
        sys.exit(1)

    # Skip the first two lines assuming it's a header
    lines = lines[2:]
    # Parse the simulation data
    for line in lines:
        values = line.rstrip().split('\t')
        values = [float(val) for val in values]
        data.append(values)
 
    return data



builds = ["openmp"]
solvers = ["Fierro"]


executables = []
tests = []

# Add paths to executable
for i in range(len(builds)):
    for j in range(len(solvers)):
        executables.append("../build-SGH-"+builds[i]+"/bin/"+solvers[j])
        if not os.path.exists(executables[i]):
            raise ValueError("Executable not found in "+executables[i]+" directory")


# # Check that each executable exists
# for i in range(len(solvers)):
#     for j in range(len(builds)):
#         if not os.path.exists(executables[i]):
#             raise ValueError(solvers[i]+" executable not found in build-fierro-openmp directory")
\


# Add names of each test
tests = ["Sedov", "Sod_X", "Sod_Y", "Sod_Z"]
# tests = ["Sedov"]
inputs = []

# Get inputs for tests
for i in range(len(tests)):
    inputs.append("standard_inputs/"+tests[i]+".yaml")

# print(inputs)

standard_results = []
# Get paths to standard results
for i in range(len(tests)):
    pattern = "standard_results/"+tests[i]+"/state/mat_pt_state*"
    file = glob.glob(pattern)
    standard_results.append(file[0])

# print(standard_results)

# Run each tests
for i in range(len(executables)):
    for j in range(len(tests)):
        # Call Fierro with YAML inputs
        print("Running "+tests[j])
        os.system(executables[i] + ' ' + inputs[j])

        # Compare to standard results
        pattern = "state/mat_pt_state*"
        file = glob.glob(pattern)
        file_path = file[-1]

        result_data = extract_state_data(file_path)
        standard_data = extract_state_data(standard_results[j])

        calc_density = [row[5] for row in result_data]
        true_density = [row[5] for row in standard_data]

        for k in range(len(true_density)):
            print(calc_density[k] - true_density[k])

        os.system('rm -rf  state' )


# for i in range(len(solvers)):
#     tmp1 = []
#     tmp2 = []
#     for j in range(len(tests[i])):
#         tmp1.append("Solver-Inputs/"+solvers[i]+"/"+tests[i][j]+".yaml")
#         tmp2.append("standard-results/"+solvers[i]+"/"+tests[i][j]+"/vtk/data/VTK0.vtk")
#     inputs.append(tmp1)
#     standard_results.append(tmp2)

# # Extract vector valued data from vtk output file
# def extract_vector_data(filename, keyword):
#     data = []
#     found_keyword = False

#     with open(filename, 'r') as file:
#         for line in file:
#             if found_keyword:
#                 if line.strip():  # If the line is not blank
#                     # Split each line into fields.
#                     fields = line.split()
#                     xx = float(fields[0])
#                     yy = float(fields[1])
#                     zz = float(fields[2])
#                     # Group and append to results.
#                     vec = [xx, yy, zz]
#                     data.append(vec)
#                 else:
#                     break  # Exit the loop when encountering a blank line
#             elif keyword in line:
#                 found_keyword = True

#     return data

# # Extract scalar valued data from vtk output file
# def extract_scalar_data(filename, keyword):
#     data = []
#     found_keyword = False

#     with open(filename, 'r') as file:
#         for line in file:
#             if found_keyword:
#                 if line.strip():  # If the line is not blank
#                     fields = line.split()
#                     if fields[0] != "LOOKUP_TABLE":
#                         data.append(float(fields[0]))
#                 else:
#                     break  # Exit the loop when encountering a blank line
#             elif keyword in line:
#                 found_keyword = True
#                 # file.next() # skip first line "LOOKUP_TABLE"
 
#     return data

# # Calculate the percent difference between two arrays of scalars
# def percent_difference_scalars(array1, array2):
#     if len(array1) != len(array2):
#         raise ValueError("Arrays must have the same length")

#     percent_diff = []
#     for i in range(len(array1)):
#         diff = array2[i] - array1[i]
#         percent_diff.append((diff / array1[i]) * 100 if array1[i] != 0 else 0)

#     return percent_diff

# # Calculate the percent difference between two arrays of vectors
# def percent_difference_vectors(array1, array2):
#     if len(array1) != len(array2):
#         raise ValueError("Arrays must have the same length")

#     percent_diff = []
#     for i in range(len(array1)):
#         if len(array1[i]) != 3 or len(array2[i]) != 3:
#             raise ValueError("Subarrays must have length 3")
        
#         diff = [array2[i][j] - array1[i][j] for j in range(3)]
#         percent_diff.append([(diff[j] / array1[i][j]) * 100 if array1[i][j] != 0 else 0 for j in range(3)])

#     percent_diff_mag = []
#     for i in range(len(array1)):
#         percent_diff_mag.append(magnitude(percent_diff[i]))

#     return percent_diff_mag

# # Calculate the magnitude of a vector
# def magnitude(array):
#     mag = math.sqrt(sum(x**2 for x in array))
#     return mag

# # Run each test
# for i in range(len(solvers)):
#     for j in range(len(tests[i])):
        
#         # Run simulation
#         print("Running "+tests[i][j])
#         os.system(executables[i] + ' ' + inputs[i][j])

#         GT_positions = extract_vector_data(standard_results[i][j], position_keyword)

#         # Read simulation results
#         results_filename = "vtk/data/VTK0.vtk"

#         if os.path.exists(results_filename):
#             print("Simulation Finished")
#         else:
#             print("Simulation did not finish")
#             raise ValueError("Simulation did not finish")

#         results_positions = extract_vector_data(results_filename, position_keyword)
#         position_diff = percent_difference_vectors(GT_positions, results_positions)


#         for k in range(len(position_diff)):
#             if position_diff[k] >= 1.0e-6:
#                 raise ValueError(" ****************** ERROR: Position difference out of range for "+tests[i][j]+" problem ****************** ")

#         print("Removing simulation outputs")
#         os.system('rm -rf  vtk' )