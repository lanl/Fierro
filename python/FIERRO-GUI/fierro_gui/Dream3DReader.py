import simplnx as nx
import itkimageprocessing as nxitk
import orientationanalysis as nxor
import pygraphviz as pgv
import os
import matplotlib.pyplot as plt
import numpy as np

def show_data_structure_heirarchy(data_structure: nx.DataStructure) -> None:
    """
    This method will create an image and then show that image in MatPlotLib
    """

    # This will generate the hierarchy as a GraphViz formatted string that you can
    # print or save to a file
    graphviz_content = data_structure.hierarchy_to_graphviz()

    # Create a graph from the DOT string using pygraphviz
    G = pgv.AGraph(string=graphviz_content)
    temp_file_path = '.graphviz_output.png'

    # Render the graph to a PNG file
    G.draw(temp_file_path, format='png', prog='dot')
  
    # Use Matplotlib to display the generated image
    img = plt.imread(temp_file_path)
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.axis('off')  # Hide axes
    plt.show()

    # Check if the file exists to avoid an error if the file is not found
    if os.path.exists(temp_file_path):
        os.remove(temp_file_path)
        print("File has been deleted successfully.")

# Create a Data Structure
# Create the DataStructure object
def Dream3DReader(fileIn, out, out_dir): 
    data_structure = nx.DataStructure()

    '''
    result = nx.ReadStlFileFilter.execute(
        data_structure=data_structure, 
        face_attribute_matrix_name="FaceData", 
        face_normals_name="FaceNormals", 
        scale_factor=1, scale_output=False, 
        stl_file_path="Cylinder.stl", 
        output_triangle_geometry_path=nx.DataPath("TriangleDataContainer"), 
        vertex_attribute_matrix_name="Vertex Data"
    )

    output_file_path = "output/ImportSTLFile.dream3d"
    result = nx.WriteDREAM3DFilter.execute(
        data_structure=data_structure, 
        export_file_path=output_file_path, 
        write_xdmf_file=True)
    if len(result.errors) != 0:
        print('Errors: {}', result.errors)
        print('Warnings: {}', result.warnings)
    else:
        print("No errors running the WriteDream3D filter")

    if len(result.errors) != 0:
        print('Errors: {}', result.errors)
        print('Warnings: {}', result.warnings)
    else:
        print("No errors running ImportSTLFile filter")
    '''

    # Create a nx.Dream3dImportParameter.ImportData object and set its values
    import_data = nx.Dream3dImportParameter.ImportData()
    # Set the path to the file on the file system
    print ("THIS IS THE FILE: " + fileIn)
    import_data.file_path = "output/" + os.path.basename(os.path.normpath(fileIn))
    # Set the import_data.data_paths value to 'None' which signals to the filter to
    # import EVERY piece of data from the file.
    import_data.data_paths = None

    # Instantiate and execte the filter immediately.
    result = nx.ReadDREAM3DFilter.execute(
        data_structure=data_structure, 
        import_data_object=import_data
    )

    # Check for any execution warnings or errors
    if len(result.errors) != 0:
        print('Errors: {}', result.errors)
        print('Warnings: {}', result.warnings)
    else:
        print("No errors running the ReadDREAM3DFilter filter")

    #geometry: nx.TriangleGeom = data_structure[nx.DataPath("TriangleDataContainer")]
    #vertices_view = geometry.vertices.npview()
    #faces_view = geometry.faces.npview()

    #show_data_structure_heirarchy(data_structure=data_structure)

    result = nx.CreateDataArrayFilter.execute(
        data_structure=data_structure,
        component_count=2,
        data_format="",
        initialization_value_str="10",
        numeric_type_index=nx.NumericType.int32,
        output_array_path=nx.DataPath("TriangleDataContainer/FaceData/FaceLabels"),
        set_tuple_dimensions=False,
        tuple_dimensions = [[5, 4, 3]]
    )

    #geometry: nx.TriangleGeom = data_structure[nx.DataPath("TriangleDataContainer")]
    #vertices_view = geometry.vertices.npview()
    #faces_view = geometry.faces.npview()

    #show_data_structure_heirarchy(data_structure=data_structure)

    #vertices_view = data_structure["TriangleDataContainer/SharedVertexList"].npview()

    stl = nx.WriteStlFileFilter.execute(
        data_structure = data_structure, 
        feature_ids_path = nx.DataPath("TriangleDataContainer/FaceData/FaceLabels"),  
        feature_phases_path = nx.DataPath(["Phases"]), 
        grouping_type_index = 0, 
        input_triangle_geometry_path = nx.DataPath("TriangleDataContainer"), 
        output_stl_directory= out_dir, 
        output_stl_file = out + ".stl", 
        output_stl_prefix = "", 
        part_number_path = nx.DataPath(["Part Number"]))

    print (stl)

    if len(stl.errors) != 0:
        print('Errors: {}', result.errors)
        print('Warnings: {}', result.warnings)
    else:
        print("No errors running the STL filter")

'''
x = nx.CreateDataGroupFilter.execute(
    data_structure=data_structure,
    data_object_path=nx.DataPath('RectGridCoords')
)
if len(x.errors) != 0:
    print('Errors: {}', result.errors)
    print('Warnings: {}', result.warnings)
else:
    print("No errors running the CreateDataGroupFilter filter")


vtk = nx.WriteVtkRectilinearGridFilter.execute(
    data_structure = data_structure,
    input_data_array_paths = [nx.DataPath("RectGridCoords")], 
    input_image_geometry_path = nx.DataPath("Output Geometry"),
    output_file = "output/test.vtk",
    write_binary_file = True
)

if len(vtk.errors) != 0:
    print('Errors: {}', result.errors)
    print('Warnings: {}', result.warnings)
else:
    print("No errors running the VTK filter")
'''