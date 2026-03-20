#!/usr/bin/python3

import os
import os.path
import sys
import math
import glob
import argparse

# Builds being tested
builds = ["cuda"]

# Name(s) of the solver being used
solvers = ["Fierro"]

solver_path = "../build/app/Fierro"

# Add names of each test
tests = ["TaylorAnvil", "TaylorAnvil_rz", "Compaction", "Compaction_rz", \
         "Sedov", "Sod_X", "Sod_Y", "Sod_Z", "Sedov_Erosion", \
        "Sedov_Read_Ensight", "Sedov_rz_polar", "Abaqus_read", \
        "Pressure_bc_box","vtu_read", \
        "lin_vol_frac_two_mat", "Bending-3D-plate", "Vel_bc_box", \
        "slanted_block_bounce", "slanted_impact", "SGTM_cooling_cube", \
        "sie_expansion_test", "confined_preload", "unconfined_preload",\
        "edge_flat_test", "billiards", "3by3_stack", "cylinder_contact", "fracture_mode_1",
        "fracture_mode_2", "fracture_reorientation"]
#,"SGTM_cooling_cube" currently broken

# Parse command line arguments
parser = argparse.ArgumentParser(description='Run regression tests.')
parser.add_argument('test_name', nargs='?', help='Name of the specific test to run')
args = parser.parse_args()

if args.test_name:
    if args.test_name in tests:
        tests = [args.test_name]
    else:
        print(f"Error: Test '{args.test_name}' not found.")
        sys.exit(1)

# Extract data from txt file
def extract_state_data(filename):
    data = []

    # Ensure the file exists before proceeding
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {filename} was not found at the specified location.")

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
        executables.append(solver_path)
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
    if not file:
        print(f"Warning: Standard result file not found for {tests[i]}")
        standard_results.append("") # Placeholder
    else:
        standard_results.append(file[0])

failed_tests = []
passed_tests = []

# Run each tests
for i in range(len(executables)):
    for j in range(len(tests)):
        test_name = tests[j]
        try:
            # Call Fierro with YAML inputs
            print("Running "+test_name)
            os.system(executables[i] + ' ' + inputs[j])
    
            # Compare to standard results
            if not standard_results[j]:
                raise FileNotFoundError(f"Standard result file missing for {test_name}")

            pattern = "state/mat_pt_state*"
            files = glob.glob(pattern)
            
            try: 
                fileIndex = list(file.split("/")[-1] for file in files).index(standard_results[j].split("/")[-1])
            except ValueError:
                raise ValueError("State file not found for "+test_name+" test")
    
    
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
    
                    if abs(diff) > 1E-8:
                        print(f"{'Calculated Result:':<20} {calc[l]:.10e}")
                        print(f"{'Expected Result:':<20} {true[l]:.10e}")
                        print(f"{'Difference:':<20} {diff:.10e}")
                        raise Exception("Results do not match for "+header1[k]+" for the "+test_name+" test!")
            
            print("\n******************************")
            print("\nPASSED : "+test_name+"\n")
            print("******************************\n")
            passed_tests.append(test_name)

        except Exception as e:
            print(f"\nFAILED : {test_name}")
            print(f"Error: {e}\n")
            failed_tests.append(test_name)
        
        finally:
            # Remove simulated state dump
            if os.path.exists('state'):
                os.system('rm -rf state')

print("\n" + "="*30)
print("       TEST SUMMARY")
print("="*30)
print(f"Total Tests Run: {len(passed_tests) + len(failed_tests)}")
print(f"Passed: {len(passed_tests)}")
print(f"Failed: {len(failed_tests)}")

if failed_tests:
    print("\nFailed Tests:")
    for test in failed_tests:
        print(f"  - {test}")
    sys.exit(1)
else:
    print("\nAll tests passed!")
    sys.exit(0)