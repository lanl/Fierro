#!/usr/bin/python3
import sys
import matplotlib
import matplotlib.pyplot as plt

# Check if the file path is provided as a command-line argument
if len(sys.argv) < 2:
    print("Usage: python script.py <path_to_file>")
    sys.exit(1)

file_path = sys.argv[1]

# Ensure the file exists before proceeding
try:
    with open(file_path, 'r') as file:
        lines = file.readlines()
except FileNotFoundError:
    print("The file was not found at the specified location.")
    sys.exit(1)

# Skip the first two lines assuming it's a header
lines = lines[2:]


# Parse the data
data = []
for line in lines:
    values = line.rstrip().split('\t')
    values = [float(val) for val in values]
    data.append(values)


# Extract the radius_3D column for the X-axis
x_values = [row[4] for row in data]
y_values = [row[5] for row in data]



plt.plot(x_values, y_values, 'x')
plt.title(f'Density')
plt.xlabel('Radius 3D')
plt.ylabel(f'Density')
plt.savefig(f"test.png", dpi=300)

