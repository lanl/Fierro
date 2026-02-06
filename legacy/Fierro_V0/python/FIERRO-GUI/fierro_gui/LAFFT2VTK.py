import numpy as np

def los_alamos_to_vtk(input_file, output_file):
    # Read the input file
    data = np.loadtxt(input_file)

    # Extract dimensions
    x_max, y_max, z_max = np.max(data[:, 3:6], axis=0).astype(int) + 1

    # Open the output file
    with open(output_file, 'w') as f:
        # Write VTK header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Los Alamos FFT Data\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {x_max} {y_max} {z_max}\n")

        # Write points
        f.write(f"POINTS {x_max * y_max * z_max} float\n")
        for z in range(z_max):
            for y in range(y_max):
                for x in range(x_max):
                    f.write(f"{x} {y} {z}\n")

        # Write cell data
        f.write(f"CELL_DATA {(x_max-1) * (y_max-1) * (z_max-1)}\n")

        # Write Euler angles
        for i, angle_name in enumerate(["Phi1", "Phi", "Phi2"]):
            f.write(f"SCALARS {angle_name} float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for value in data[:, i]:
                f.write(f"{value}\n")

        # Write Feature ID
        f.write("SCALARS Feature_ID int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for value in data[:, 6]:
            f.write(f"{int(value)}\n")

        # Write Phase ID
        f.write("SCALARS Phase_ID int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for value in data[:, 7]:
            f.write(f"{int(value)}\n")
