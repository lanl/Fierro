#include "MeshIO.h"
#include "Mesh.h"
#include "IOUtilities.h"
#include <fstream>
#include <vector>
#include <filesystem>
#include <sstream>
#include "VtkIO.h"

namespace {
    void write_geometry(FILE* out, const Mesh& mesh) {
        fprintf(out, "# vtk DataFile Version 2.0\n");  // part 2
        fprintf(out, "Mesh for Fierro\n");             // part 2
        fprintf(out, "ASCII \n");                      // part 3
        fprintf(out, "DATASET UNSTRUCTURED_GRID\n\n"); // part 4
        
        fprintf(out, "POINTS %ld float\n", mesh.points.dims(0));

        // write all components of the point coordinates
        for (size_t i = 0; i < mesh.points.dims(0); i++){
            fprintf(out, 
                "%f %f %f\n",
                mesh.points(i, 0),
                mesh.points(i, 1),
                mesh.num_dim > 2 ? mesh.points(i, 2) : 0.
            );
        }
        
        /*
        ---------------------------------------------------------------------------
        Write the elems
        ---------------------------------------------------------------------------
        */
        int num_elems = mesh.element_point_index.dims(0);
        fprintf(out, "\n");
        fprintf(out, "CELLS %d %ld\n", num_elems, num_elems + num_elems * mesh.element_point_index.dims(1));  
        // size=all printed values
        
        auto lagrange_map = VtkIO::lagrange_to_ijk(mesh.num_dim, mesh.p_order);
        auto& linear_map = MeshIO::_Impl::fea_to_ijk();

        // write all global point numbers for this elem
        for (size_t i=0; i < num_elems; i++) {
            int num_points_in_elem = mesh.element_point_index.dims(1);
            
            // Print number of points in the element
            fprintf(out, "%d ", num_points_in_elem);
            for (size_t j = 0; j < num_points_in_elem; j++) {
                size_t j_reordered;
                if (VtkIO::is_lagrange_ordered(mesh.element_types(i)))
                    j_reordered = lagrange_map[j];
                else
                    j_reordered = linear_map[j];
                fprintf(out, "%d ", mesh.element_point_index(i, j_reordered));
            }
            fprintf(out, "\n");
        }
        
        fprintf(out, "\n");
        fprintf(out, "CELL_TYPES %d \n", num_elems);
        for (size_t i = 0; i < num_elems; i++) {
            fprintf(out, "%d \n", mesh.element_types(i));
        }
    }

    void write_nodal_variables(FILE* out, const Mesh& mesh) {
        size_t num_points = mesh.points.dims(0);
        fprintf(out, "\n");
        fprintf(out, "POINT_DATA %ld \n", num_points);
        fprintf(out, "SCALARS point_var float 1\n"); // the 1 is number of scalar components [1:4]
        fprintf(out, "LOOKUP_TABLE default\n");

        // TODO: Eventually allow for mesh data here.
        for (size_t i = 0; i < num_points; i++) {
            fprintf(out, "%f\n", (double)i);
        }
    }

    void write_vector_variables(FILE* out, const Mesh& mesh) {
        fprintf(out, "\n");
        fprintf(out, "VECTORS point_vec float\n");
        // TODO: Pull from real mesh data.
        for (size_t i = 0; i < mesh.points.dims(0); i++) {
            double var1=0;
            double var2=1;
            double var3=2;
            fprintf(out, "%f %f %f\n", var1, var2, var3);
        }
    }

    void write_element_variables(FILE* out, const Mesh& mesh) {
        size_t num_elems = mesh.element_point_index.dims(0);
        fprintf(out, "\n");
        fprintf(out, "CELL_DATA %ld \n", num_elems);
        fprintf(out, "SCALARS elem_var float 1\n"); // the 1 is number of scalar components [1:4]
        fprintf(out, "LOOKUP_TABLE default\n");
        for (size_t i = 0; i < num_elems; i++) {
            fprintf(out, "%f\n", (double)i);
        }
        
        fprintf(out, "\n");
        fprintf(out, "SCALARS elem_var2 float 1\n"); // the 1 is number of scalar components [1:4]
        fprintf(out, "LOOKUP_TABLE default\n");
        for (size_t i = 0; i < num_elems; i++) {
            fprintf(out, "%f\n", -(double)i);
        }
    }

    template<typename T>
    mtr::CArray<T> vec_to_array(const std::vector<std::vector<T>>& vec) {
        if (vec.size() == 0)
            return mtr::CArray<T>();
        
        auto array = mtr::CArray<T>(vec.size(), vec[0].size());
        for (size_t i = 0; i < vec.size(); i++) {
            const auto& row = vec[i];
            if (row.size() != array.dims(1))
                throw std::runtime_error("Cannot load ragged arrays.");
            for (size_t j = 0; j < row.size(); j++) {
                array(i, j) = row[j];
            }
        }

        return array;
    }

    template<typename T>
    std::vector<T> parse_line(const std::string& line) {
        std::vector<T> result;
        std::stringstream ss(line);
        while (!ss.eof()) {
            T v;
            ss >> v;
            result.push_back(v);
        }
        return result;
    }

    template<typename T>
    void for_n_lines(std::ifstream& in, const size_t n, T f) {
        std::string line;
        for (size_t i = 0; i < n; i++) {
            std::getline(in, line);
            f(i, line);
        }
    }

    template<typename T>
    std::vector<std::vector<T>> read_2darray(std::ifstream& in, const size_t length) {
        std::vector<std::vector<T>> data;
        for_n_lines(in, length, [&](const size_t i, const std::string& line) {
            data.push_back(parse_line<T>(line));
        });
        return data;
    } 

    template<typename T>
    mtr::CArray<T> read_1darray(std::ifstream& in, const size_t length) {
        mtr::CArray<T> data = mtr::CArray<int>(length);
        for_n_lines(in, length, [&](const size_t i, const std::string& line) {
            data(i) = parse_line<T>(line)[0];
        });
        return data;
    }

    bool read_points(std::ifstream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
        if (header_tokens[0] != "POINTS" || header_tokens[2] != "float")
            return false;
        
        int num_points = stoi(header_tokens[1]);
        if (num_points <= 0)
            return false;

        auto points = read_2darray<double>(in, num_points);
        if (points.size() == 0)
            return false;
        
        mesh.points = vec_to_array(points);
        mesh.num_dim = mesh.points.dims(1);
        return true;
    }

    bool read_cells(std::ifstream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
        if (header_tokens[0] != "CELLS")
            return false;
        
        int num_elems = stoi(header_tokens[1]);
        if (num_elems <= 0)
            return false;

        auto elem_point_index = read_2darray<int>(in, num_elems);
        if (elem_point_index.size() == 0)
            return false;
        
        for (auto& row : elem_point_index)
            row.erase(row.begin()); // First column is points per element.
        
        mesh.element_point_index = vec_to_array(elem_point_index);
        return true;
    }

    bool read_cell_types(std::ifstream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
        if (header_tokens[0] != "CELLS")
            return false;
        
        int num_elems = stoi(header_tokens[1]);
        if (num_elems <= 0)
            return false;
        
        mesh.element_types = read_1darray<int>(in, num_elems);
        return true;
    }
}

Mesh MeshIO::read_vtk(std::string filename, bool verbose) {
    std::ifstream in(filename);
    auto mesh = Mesh();

    bool found_points     = false,
         found_cells      = false,
         found_cell_types = false;

    while (!in.eof()) {
        std::string line;
        std::getline(in, line);

        if (line.length() == 0)
            continue;

        std::vector<std::string> tokens = _Impl::split(line, " ");

        if (read_points(in, tokens, mesh))
            found_points = true;
        else if (read_cells(in, tokens, mesh))
            found_cells = true;
        else if (read_cell_types(in, tokens, mesh))
            found_cell_types = true;
    }

    in.close();

    if (!mesh.validate()) {
        throw std::runtime_error("Invalid mesh constructed.");
    }
    
    if (mesh.element_types.dims(0) > 0) {
        if (VtkIO::is_lagrange_ordered(mesh.element_types(0))) {
            _Impl::redorder_columns(
                mesh.element_point_index,
                VtkIO::ijk_to_lagrange(mesh.num_dim, mesh.p_order).data()
            );
        } else {
            _Impl::redorder_columns(
                mesh.element_point_index, 
                _Impl::ijk_to_fea().data()
            );
        }
    }
    
    return mesh;
}

void MeshIO::write_vtk(std::string filename, const Mesh& mesh, bool verbose) {
    IOUtilities::mkdir("vtk");
    auto path = std::filesystem::path("vtk") / (filename + ".vtk");

    FILE* out = fopen(path.c_str(), "w");
    if (verbose)
        std::cout << "Creating file: " << path.string() << std::endl;
        
    write_geometry(out, mesh);
    write_nodal_variables(out, mesh);
    write_vector_variables(out, mesh);
    write_element_variables(out, mesh);
    fclose(out);
}