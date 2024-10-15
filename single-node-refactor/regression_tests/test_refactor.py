#!/usr/bin/python3

import os
import os.path
import sys
import math
import glob

# Builds being tested
builds = ["openmp"]

# Name(s) of the solver being used
solvers = ["Fierro"]

# Add names of each test
tests = ["Sedov", "Sod_X", "Sod_Y", "Sod_Z", "Sedov_Erosion","Sedov_Read_Ensight", "Sedov_rz_polar", "Abaqus_read"]

# Extract data from txt file
def extract_state_data(filename):
    data = []

    # Ensure the file exists before proceeding
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print("The file was not found at the specified location.")
        sys.exit(1)

    # Skip the first line to get table headers
    lines = lines[1:]
    # First, replace '#' with nothing ('') to remove it
    cleaned_header = lines[0].replace('#', '')
    # Then, split the cleaned string into a list of words
    headers = cleaned_header.split()

    # Skip one more line to get data
    lines = lines[1:]
    # Parse the simulation data
    for line in lines:
        values = line.rstrip().split('\t')
        values = [float(val) for val in values]
        data.append(values)
    return data, headers

# Grab paths to executable
executables = []
for i in range(len(builds)):
    for j in range(len(solvers)):
        executables.append("../build-SGH-"+builds[i]+"/bin/"+solvers[j])
        if not os.path.exists(executables[i]):
            raise ValueError("Executable not found in "+executables[i]+" directory")

# Get paths to inputs
inputs = []
for i in range(len(tests)):
    inputs.append("standard_inputs/"+tests[i]+".yaml")

# Grab paths to standard results
standard_results = []
for i in range(len(tests)):
    pattern = "standard_results/"+tests[i]+"/state/mat_pt_state*"
    file = glob.glob(pattern)
    standard_results.append(file[0])

# Run each tests
for i in range(len(executables)):
    for j in range(len(tests)):
        # Call Fierro with YAML inputs
        print("Running "+tests[j])
        os.system(executables[i] + ' ' + inputs[j])

        # Compare to standard results
        pattern = "state/mat_pt_state*"
        files = glob.glob(pattern)
        
        try: fileIndex = list(file.split("/")[-1] for file in files).index(standard_results[j].split("/")[-1])
        except ValueError:
            raise ValueError("State file not found for "+tests[j]+" test")


        file_path = files[fileIndex]

        print(file_path)
        print(standard_results[j])

        result_data, header1 = extract_state_data(file_path)
        standard_data, header2 = extract_state_data(standard_results[j])

        for k in range(len(result_data[0])):
            calc = [row[k] for row in result_data]
            true = [row[k] for row in standard_data]



            for l in range(len(calc)):
                diff = calc[l] - true[l]
                # print(diff)

                if abs(diff) > 1E-11:
                    print(diff)
                    print(calc[l])
                    print(true[l])
                    raise Exception("Results do not match for "+header1[k]+" for the "+tests[j]+" test!")

        # Remove simulated state dump
        os.system('rm -rf  state' )

        print("\n******************************")
        print("\nPASSED : "+tests[j]+"\n")
        print("******************************\n")