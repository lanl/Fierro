import vtk
import os

def convert_exodus_to_vtm(input_exo, output_vtm, output_full_mesh_vtu):
    """
    Convert a multi-part Exodus mesh file to a VTM (vtkMultiBlockDataSet) file,
    exporting individual parts as ASCII .vtu files and the full mesh as a single ASCII .vtu file.

    Parameters:
        input_exo (str): Path to the input Exodus file.
        output_vtm (str): Path to the output .vtm file.
        output_full_mesh_vtu (str): Path to the output full mesh .vtu file.
    """
    # Create the Exodus reader
    reader = vtk.vtkExodusIIReader()
    reader.SetFileName(input_exo)
    reader.UpdateInformation()

    print("\n")
    print("    ....Processing")

    # Get number of element blocks
    num_blocks = reader.GetNumberOfElementBlockArrays()
    print(f"    Number of parts in the Exodus file: {num_blocks}")

    # Enable all element blocks for extraction
    for i in range(num_blocks):
        block_name = reader.GetElementBlockArrayName(i)
        reader.SetElementBlockArrayStatus(block_name, 1)

    # Update reader to load all blocks
    reader.Update()

    # Get the full Exodus dataset
    exodus_output = reader.GetOutput()

    # Create a MultiBlock dataset
    multiblock = vtk.vtkMultiBlockDataSet()
    append_filter = vtk.vtkAppendFilter()

    vtmtext = []

    for i in range(num_blocks):
        block_name = reader.GetElementBlockArrayName(i)

        # saving the text for the .vtm file
        text = f"    <DataSet index=\"{i}\" name=\"{block_name}\" file=\"data/{block_name}.vtu\" />"
        vtmtext.append(text)


        # Extract the block using vtkUnstructuredGrid
        block_mesh = reader.GetOutput().GetBlock(0).GetBlock(i)  # Correct way to get blocks

        if isinstance(block_mesh, vtk.vtkUnstructuredGrid) and block_mesh.GetNumberOfCells() > 0:
            # Assign block to the MultiBlock dataset
            multiblock.SetBlock(i, block_mesh)
            multiblock.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), block_name)

            # Write each block as a separate VTU file in ASCII format
            vtu_filename = f"{block_name}.vtu"
            vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
            vtu_writer.SetFileName(vtu_filename)
            vtu_writer.SetInputData(block_mesh)
            vtu_writer.SetDataModeToAscii()  # Ensure ASCII format
            vtu_writer.Write()
            print(f"    Saved mesh block {block_name} as {vtu_filename}")

            # Append to full mesh
            append_filter.AddInputData(block_mesh)

    # Merge all blocks into a full mesh
    append_filter.Update()
    full_mesh = vtk.vtkUnstructuredGrid()
    full_mesh.DeepCopy(append_filter.GetOutput())

    # Clean the merged mesh to remove duplicate points
    clean_filter = vtk.vtkCleanUnstructuredGrid()
    clean_filter.SetInputData(full_mesh)
    clean_filter.Update()

    # Get the cleaned, merged mesh
    full_mesh_cleaned = clean_filter.GetOutput()


    # Write the full mesh as a separate VTU file in ASCII format
    vtu_writer_full_mesh = vtk.vtkXMLUnstructuredGridWriter()
    vtu_writer_full_mesh.SetFileName(output_full_mesh_vtu)
    vtu_writer_full_mesh.SetInputData(full_mesh_cleaned)
    vtu_writer_full_mesh.SetDataModeToAscii()  # Ensure ASCII format
    vtu_writer_full_mesh.Write()
    print(f"    Saved full mesh as {output_full_mesh_vtu}")

    # save the vtm text for the full mesh
    text = f"    <DataSet index=\"{num_blocks}\" name=\"{output_full_mesh_vtu}\" file=\"data/{output_full_mesh_vtu}\" />"
    vtmtext.append(text)
    

    # write vtm file and move the files
    with open(f"{output_vtm}", 'w') as file:
        file.write("<?xml version=\"1.0\"?>\n")
        file.write("<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n")
        file.write("  <vtkMultiBlockDataSet>\n")
        for i in range(num_blocks+1):
            file.write(vtmtext[i]+'\n')
        file.write("  </vtkMultiBlockDataSet>\n")
        file.write("</VTKFile>") 
    
    # make data folder to store vtu files
    directory_path = "data"
    os.makedirs(directory_path, exist_ok=True)
    os.system('mv *.vtu data/.')


def show_help():
    help_text = """
    Usage: python3 exodus_to_ascii.py <exodus_file> <mesh_out_name>

    The script converts the supplied exodus file into the ASCII 
    multi-block VTK format.  Each mesh part will have a unique mesh 
    .vtu file along with a single full mesh .vtu file.  The VTK 
    files can be viewed by opening the .vtm file in Paraview.
    """
    print(help_text)

# Example usage
if __name__ == "__main__":

    import sys
    
    input_exodus_file = None # Initialize the file name
    output_vtm_file = None # Initialize the file name
    output_full_mesh_vtu = "full_mesh.vtu"  # Output full mesh VTU file

    success = 0

    if len(sys.argv) == 2:
        if '--help' in sys.argv:
            show_help()
        else:
            input_exodus_file = sys.argv[1]  # Get the file name
            output_vtm_file = "mesh_output.vtm"   
            convert_exodus_to_vtm(input_exodus_file, output_vtm_file, output_full_mesh_vtu)
            success = 1
    elif len(sys.argv) == 3:
        input_exodus_file = sys.argv[1]  # First file
        output_vtm_file = sys.argv[2]  # Second file
        output_vtm_file += ".vtm"
        convert_exodus_to_vtm(input_exodus_file, output_vtm_file, output_full_mesh_vtu)
        success = 1
    elif len(sys.argv) > 3:
        print("\n")
        print("    Too many input arguments supplied")
        show_help()
    else:
        print("\n")
        print("    Please provide an exodus mesh file \n")
        show_help()

    if (success == 1):
        print("\n")
        print("    ....Finished \n")



