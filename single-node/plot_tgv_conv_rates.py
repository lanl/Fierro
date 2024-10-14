import os
import numpy as np
import matplotlib.pyplot as plt
import sys

def exact_solution(x, y, z):
    """Compute the exact solution."""
    U_exact_x = np.sin(np.pi * x) * np.cos(np.pi * y)
    U_exact_y = -np.cos(np.pi * x) * np.sin(np.pi * y)
    U_exact_z = np.zeros_like(x)
    return np.sqrt(U_exact_x**2  + U_exact_y**2 + U_exact_z**2)

def calculate_L1_magnitude_error(U_num_mag, U_exact_mag, dx, dy, dz):
    """Calculate L1 error in the magnitude of the velocity."""
    abs_diff = np.abs(U_num_mag - U_exact_mag) * dx * dy * dz
    abs_exact = np.abs(U_exact_mag) * dx * dy * dz
    L1_error = np.sum(abs_diff)/np.sum(abs_exact)
    return L1_error

def read_data(filename, delimiter):
    data = np.genfromtxt(filename, comments='#', delimiter=delimiter, invalid_raise=False)
    return data

# Read input arguments
k = int(sys.argv[1])
q = int(sys.argv[2])
r = int(sys.argv[3])
l = int(sys.argv[4])
mesh_sizes = [(2**m) for m in range(l, r+1)]
EOC = np.zeros((k-q+1, len(mesh_sizes)-1))

# Prepare for plotting
errors = np.zeros(len(mesh_sizes))
inverse_mesh_sizes = 1.0 / np.array(mesh_sizes)

# Create the output directory if it doesn't exist
output_dir = "TGV_CONV_PLOTS"
os.makedirs(output_dir, exist_ok=True)

# Main loop
index = 0
for i in range(q, k + 1):
    j = i - 1
    
    for idx, mesh_size in enumerate(mesh_sizes):
        mesh_size_str = f"{mesh_size}x{mesh_size}x1"
        data_file = f"TGV_convergence/TGV_Q{i}Q{j}_{mesh_size_str}_nodes.txt"
        mat_pt_data_file = f"TGV_convergence/TGV_Q{i}Q{j}_{mesh_size_str}_matpt.txt"
        current_data = read_data(data_file, delimiter='\t') 
        current_mat_pt_data = read_data(mat_pt_data_file, delimiter='\t')

        print(f"Loaded data from: {data_file}")
    
        x, y, z = current_data[:, 0], current_data[:, 1], current_data[:, 2]
        u_x, u_y, u_z = current_data[:, 3], current_data[:, 4], current_data[:, 5]

        U_num_mag = np.sqrt(u_x**2 +  u_y**2 + u_z**2)
        
        U_exact_mag = exact_solution(x, y, z)
        
        h = min(current_mat_pt_data[:, 10])  # Assuming the mesh size is in column 11 (index 10)
    
        # Calculate the L1 error in the magnitude
        L1_error = calculate_L1_magnitude_error(U_num_mag, U_exact_mag, h, h, h)
        errors[idx] = L1_error
        print(f'L1 Error for {mesh_size}x{mesh_size}x1 elements: {errors[idx]}')
       
    # Plotting the error against the mesh size
    plt.figure()
    plt.loglog(inverse_mesh_sizes, errors, 'o-', linewidth=2)
    plt.xlabel('1/h')
    plt.ylabel('L1 Error')
    plt.grid(True)
    plt.xticks(inverse_mesh_sizes)
    plt.gca().tick_params(axis='both', which='major', labelsize=15)

    p = np.polyfit(np.log(inverse_mesh_sizes), np.log(errors), 1)
    plt.loglog(inverse_mesh_sizes, np.exp(np.polyval(p, np.log(inverse_mesh_sizes))), '--', linewidth=2)
    plt.legend(['Error', f'Fit: slope = {p[0]:.2f}'], loc='upper left', fontsize=15, edgecolor='white')
    plt.tight_layout()

    # Save the plot to the TGV_CONV_PLOTS folder
    plot_filename = f"{output_dir}/convergence_Q{i}Q{j}.png"
    plt.savefig(plot_filename)
    print(f'Plot saved to {plot_filename}')

    # Display the slope of the fit line
    print(f'Best fit slope: {p[0]:.2f}')

    # Compute EOC for this pair
    for idx in range(len(mesh_sizes)-1):
        EOC[i-q, idx] = np.log(errors[idx]/errors[idx+1])/np.log(2)

# Create the output directory if it doesn't exist
eoc_output_dir = "TGV_EOC_TABLE"
os.makedirs(eoc_output_dir, exist_ok=True)

# Write the EOC values to a LaTeX table
latex_filename = f"{eoc_output_dir}/eoc_table_Q{k}_Q{q}.tex"
with open(latex_filename, 'w') as f:
    f.write("\\begin{table}\n")
    f.write("\t\\centering\n")
    f.write("\t\\begin{tabular}{|c|" + "c|" * (k-q+1) + "}\n")
    f.write("\t\t\\hline\n")
    f.write("\t\t$h$ & " + " & ".join([f"$L^1$ EOC Q{i}-Q{i-1}" for i in range(q, k+1)]) + " \\\\\n")
    f.write("\t\t\\hline\\hline\n")
    for idx in range(len(mesh_sizes)):
        h_value = f"$\\frac{{1}}{{{2**(idx + l)}}}$"
        if idx < len(mesh_sizes) - 1:
            eoc_values = " & ".join([f"{EOC[i-q, idx]:.4f}" for i in range(q, k+1)])
        else:
            eoc_values = " & ".join(["-" for _ in range(q, k+1)])
        f.write(f"\t\t{h_value} & {eoc_values}\\\\\n")
        f.write("\t\t\\hline\n")
    f.write("\t\\end{tabular}\n")
    f.write("\t\\caption{Experimental order of convergence with limited artificial viscosity at $t=0.1$.}\n")
    f.write("\t\\label{{tab:eoc_Q{q}_Q{k}}}\n")
    f.write("\\end{table}\n")

print(f"LaTeX table saved to {latex_filename}")

