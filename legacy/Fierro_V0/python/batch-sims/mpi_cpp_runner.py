from mpi4py import MPI
import subprocess
import os

def run_cpp_code(rank):
    """ Runs a C++ program with the rank as an argument. """
    # Call the compiled C++ code using subprocess
    # Assuming the compiled C++ executable is named 'cpp_program'
    result = subprocess.run(['./cpp_program', str(rank)], capture_output=True, text=True)
    
    # Capture the output (assuming C++ program prints the result)
    output = result.stdout.strip()
    
    # Return the result from the C++ program
    return output

def main():
    # Initialize the MPI environment
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # Process ID
    size = comm.Get_size()  # Total number of processes

    # Each rank runs a C++ program with its own rank as an argument
    cpp_output = run_cpp_code(rank)

    # Gather the results from all processes to the root process (rank 0)
    all_outputs = comm.gather(cpp_output, root=0)

    # Root process prints the final results
    if rank == 0:
        print("Results from all ranks:")
        for i, output in enumerate(all_outputs):
            print(f"Rank {i} output: {output}")

if __name__ == "__main__":
    main()
