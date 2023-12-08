#include "MeshIO.h"
#include "Mesh.h"
#include "IOUtilities.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <sstream>
#include <cmath>
#include "VtkIO.h"

namespace {
    void write_geometry(std::ostream& out, const Mesh& mesh) {
        out << "# vtk DataFile Version 2.0" << std::endl;  // part 2
        out << "Mesh for Fierro" << std::endl;             // part 2
        out << "ASCII" << std::endl;                       // part 3
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;   // part 4
        out << std::endl;
        
        out << "POINTS " << mesh.points.dims(0) << " float" << std::endl;

        // write all components of the point coordinates
        for (size_t i = 0; i < mesh.points.dims(0); i++){
            for (size_t j = 0; j < mesh.points.dims(1); j++)
                out << mesh.points(i, j) << " ";
            if (mesh.points.dims(1) == 2)
                out << 0;
            out << std::endl;
        }
        
        /*
        ---------------------------------------------------------------------------
        Write the elems
        ---------------------------------------------------------------------------
        */
        int num_elems = mesh.element_point_index.dims(0);
        out << std::endl;
        auto size = num_elems + num_elems * mesh.element_point_index.dims(1);
        out << "CELLS " << num_elems << " " << size << std::endl;
        // size=all printed values
        
        auto lagrange_map = VtkIO::lagrange_to_ijk(mesh.num_dim, mesh.p_order);
        auto& linear_map = MeshIO::_Impl::fea_to_ijk();

        // write all global point numbers for this elem
        for (size_t i = 0; i < num_elems; i++) {
            int num_points_in_elem = mesh.element_point_index.dims(1);
            
            // Print number of points in the element
            out << num_points_in_elem << " ";
            for (size_t j = 0; j < num_points_in_elem; j++) {
                size_t j_reordered;
                if (VtkIO::is_lagrange_ordered(mesh.element_types(i)))
                    j_reordered = lagrange_map[j];
                else
                    j_reordered = linear_map[j];
                out << mesh.element_point_index(i, j_reordered) << " ";
            }
            out << std::endl;
        }
        
        out << std::endl;
        out << "CELL_TYPES " << num_elems << std::endl;
        for (size_t i = 0; i < num_elems; i++) {
            out << mesh.element_types(i) << std::endl;
        }
    }

    void write_nodal_variables(std::ostream& out, const Mesh& mesh) {
        size_t num_points = mesh.points.dims(0);
        out << std::endl;
        out << "POINT_DATA " << num_points << std::endl;
        out << "SCALARS " << "point_var " << "float " << "1" << std::endl; // the 1 is number of scalar components [1:4]
        out << "LOOKUP_TABLE " << "default" << std::endl;

        // TODO: Eventually allow for mesh data here.
        for (size_t i = 0; i < num_points; i++) {
            out << i << std::endl;
        }
    }

    void write_vector_variables(std::ostream& out, const Mesh& mesh) {
        out << std::endl;
        out << "VECTORS point_vec float" << std::endl;
        // TODO: Pull from real mesh data.
        for (size_t i = 0; i < mesh.points.dims(0); i++) {
            double var1=0;
            double var2=1;
            double var3=2;
            out << var1 << " " << var2 << " " << var3 << std::endl;
        }
    }

    void write_element_variables(std::ostream& out, const Mesh& mesh) {
        size_t num_elems = mesh.element_point_index.dims(0);
        out << std::endl;
        out << "CELL_DATA " << num_elems << std::endl;
        out << "SCALARS " << "elem_var " << "float " << "1" << std::endl; // the 1 is number of scalar components [1:4]
        out << "LOOKUP_TABLE " << "default " << std::endl;
        for (size_t i = 0; i < num_elems; i++) {
            out << i << std::endl;
        }
        out << std::endl;
        out << "SCALARS " << "elem_var2 " << "float " << "1" << std::endl; // the 1 is number of scalar components [1:4]
        out << "LOOKUP_TABLE " << "default " << std::endl;
        for (size_t i = 0; i < num_elems; i++) {
            out << -i << std::endl;
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
        T v;
        while (ss >> v) {
            result.push_back(v);
        }
        return result;
    }

    template<typename T>
    void for_n_lines(std::istream& in, const size_t n, T f) {
        std::string line;
        for (size_t i = 0; i < n; i++) {
            std::getline(in, line);
            f(i, line);
        }
    }

    template<typename T>
    std::vector<std::vector<T>> read_2darray(std::istream& in, const size_t length) {
        std::vector<std::vector<T>> data;
        for_n_lines(in, length, [&](const size_t i, const std::string& line) {
            data.push_back(parse_line<T>(line));
        });
        return data;
    } 

    template<typename T>
    mtr::CArray<T> read_1darray(std::istream& in, const size_t length) {
        mtr::CArray<T> data = mtr::CArray<int>(length);
        for_n_lines(in, length, [&](const size_t i, const std::string& line) {
            data(i) = parse_line<T>(line)[0];
        });
        return data;
    }

    bool read_points(std::istream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
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

    bool read_cells(std::istream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
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

    bool read_cell_types(std::istream& in, const std::vector<std::string>& header_tokens, Mesh& mesh) {
        if (header_tokens[0] != "CELL_TYPES")
            return false;
        
        int num_elems = stoi(header_tokens[1]);
        if (num_elems <= 0)
            return false;
        
        mesh.element_types = read_1darray<int>(in, num_elems);
        return true;
    }
}

Mesh MeshIO::read_vtk(std::istream& in, bool verbose) {
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

    if (found_points && found_cells)
        mesh.p_order = std::pow(mesh.element_point_index.dims(1), 1.0 / mesh.num_dim);

    if (!mesh.validate()) {
        throw std::runtime_error("Invalid mesh constructed.");
    }
    
    if (mesh.element_types.dims(0) > 0) {
        if (VtkIO::is_lagrange_ordered(mesh.element_types(0))) {
            _Impl::reorder_columns(
                mesh.element_point_index,
                VtkIO::ijk_to_lagrange(mesh.num_dim, mesh.p_order).data()
            );
        } else {
            _Impl::reorder_columns(
                mesh.element_point_index, 
                _Impl::ijk_to_fea().data()
            );
        }
    }
    return mesh;
}

Mesh MeshIO::read_vtk(std::string filename, bool verbose) {
    std::ifstream in(filename);
    
    auto mesh = read_vtk(in, verbose);

    in.close();
    return mesh;
}

void MeshIO::write_vtk(std::ostream& out, const Mesh& mesh) {
    write_geometry(out, mesh);
    write_nodal_variables(out, mesh);
    write_vector_variables(out, mesh);
    write_element_variables(out, mesh);
}

void MeshIO::write_vtk(std::string filename, const Mesh& mesh, bool verbose) {
    IOUtilities::mkdir("vtk");
    auto path = std::filesystem::path("vtk") / (filename + ".vtk");

    std::ofstream out(path.c_str(), std::ofstream::out);
    if (verbose)
        std::cout << "Creating file: " << path.string() << std::endl;
    write_vtk(out, mesh);
    out.close();
}
