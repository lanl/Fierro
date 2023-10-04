#include "MeshIO.h"
#include "IOUtilities.h"
#include "matar.h"

namespace {
    std::string get_format_string(const Mesh& mesh) {
        switch (mesh.num_dim) {
            case 2:
                return "quad4";
                break;
            case 3:
                return "hexa8";
                break;
            default:
                throw std::runtime_error("Unsupported number of dimensions: " + mesh.num_dim);
        }
    }

    /**
     * Writes the .geo file for a mesh.
     * Currently fixed to write *.00001.geo.
     * 
     * Contains the point coordinates along with the connectivity of the mesh.
    */
    void write_geometry(const char* name, const Mesh& mesh) {
        FILE* out = fopen(name, "w");
        
        fprintf(out, "This is the 1st description line of the EnSight Gold geometry file\n");
        fprintf(out, "This is the 2nd description line of the EnSight Gold geometry file\n");
        fprintf(out, "node id assign\n");
        fprintf(out, "element id assign\n");
        
        fprintf(out, "part\n");
        fprintf(out, "%10d\n",1);
        fprintf(out, "Mesh\n");
        fprintf(out, "coordinates\n");
        fprintf(out, "%10ld\n", mesh.points.dims(0));
        
        for (size_t j = 0; j < 3; j++) {
            for (size_t i = 0; i < mesh.points.dims(0); i++) {
                fprintf(out , "%12.5e\n", j < mesh.points.dims(1) ? mesh.points(i, j) : 0.0);  
            }
        }
        
        fprintf(out, "%s", (get_format_string(mesh) + "\n").c_str());
        fprintf(out, "%10ld\n", mesh.element_point_index.dims(0));
        
        auto& map = MeshIO::_Impl::fea_to_ijk();
        for (size_t i = 0; i < mesh.element_point_index.dims(0); i++) {
            for (size_t j = 0; j < mesh.element_point_index.dims(1); j++)
                fprintf(out, "%10d", mesh.element_point_index(i, map[j]) + 1); // EnSight array indexing starts at 1
            fprintf(out, "\n");
        }
        
        fclose(out);
    }

    /**
     * Writes the .var file for the mesh.
     * Currently just writes dummy data at the node level.
     * TODO: Have this write mesh properties once the mesh properties are set.
    */
    void write_variables(const char* name, const Mesh& mesh) {
        FILE* out = fopen(name, "w");
        fprintf(out, "Per_elem scalar values\n");
        fprintf(out, "part\n");
        fprintf(out, "%10d\n", 1);

        fprintf(out, "%s\n", get_format_string(mesh).c_str());
        
        // TODO: Pull from actual mesh data here.
        for (size_t i = 0; i < mesh.element_point_index.dims(0); i++) {
            fprintf(out, "%12.5e\n", 1.);
        }
        
        fclose(out);
    }

    /**
     * Writes the metadata file for the ensight mesh.
     * Just points data/*.*****.geo and data/*.*****.var
     * 
     * Doesn't support any time series mesh data.
    */
    void write_case(const char* name, const std::string& mesh_name, const Mesh& mesh) {

        FILE* out = fopen(name, "w");
        
        fprintf(out, "FORMAT\n");
        fprintf(out, "type: ensight gold\n");
        fprintf(out, "GEOMETRY\n");
        fprintf(out, "model: data/%s.*****.geo\n", mesh_name.c_str());
        fprintf(out, "VARIABLE\n");
        fprintf(out, "scalar per element: mat_index data/%s.*****.var\n", mesh_name.c_str());
        
        // TODO: Allow for mutliple timestamps
        fprintf(out, "TIME\n");
        fprintf(out, "time set: 1\n");
        fprintf(out, "number of steps: %4d\n", 1);
        fprintf(out, "filename start number: 1\n");
        fprintf(out, "filename increment: 1\n");
        fprintf(out, "time values: \n");
        fprintf(out, "%12.5e\n", 0.);
        
        fclose(out);
    }
}

Mesh MeshIO::read_ensight(std::string filename) {
    throw std::runtime_error("Unimplemented");
}

void MeshIO::write_ensight(std::string filename, const Mesh& mesh) {
    if (mesh.p_order > 1)
        throw std::runtime_error("Cannot write a high order mesh to Ensight format.");

    char name[100];
    std::string fp;
    
    IOUtilities::mkdir("ensight/data");

    // 1 is required because Ensight dump 1 is the t=0 dump
    sprintf(name, "%s.%05d.geo", filename.c_str(), 1);
    fp = std::filesystem::path("ensight") / "data" / name;
    std::cout << fp << std::endl;
    write_geometry(fp.c_str(), mesh);
    
    sprintf(name, "%s.%05d.var", filename.c_str(), 1);
    fp = std::filesystem::path("ensight") / "data" / name;
    std::cout << fp << std::endl;
    write_variables(fp.c_str(), mesh);
    
    sprintf(name, "%s.case", filename.c_str());
    fp = std::filesystem::path("ensight") / name;
    std::cout << fp << std::endl;
    write_case(fp.c_str(), filename, mesh);
}