#pragma once
#include "Mesh.h"
#include <ostream>
#include <string>
#include <vector>

namespace {
    inline bool hasEnding(const std::string &fullString, const std::string &ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
        }
    }
}

namespace MeshIO {
    enum class FileType {
        EnSight,
        VTK
    };

    Mesh read_vtk(std::string filename, bool verbose=false);
    Mesh read_vtk(std::istream& in, bool verbose=false);
    void write_vtk(std::string filename, const Mesh& mesh, bool verbose=false);
    void write_vtk(std::ostream& out, const Mesh& mesh);

    Mesh read_ensight(std::string filename, bool verbose=false);
    void write_ensight(std::string filename, const Mesh& mesh, bool verbose=false);

    namespace _Impl {
        inline std::vector<std::string> split(std::string s, std::string delimiter) {
            size_t pos_start = 0, pos_end, delim_len = delimiter.length();
            std::string token;
            std::vector<std::string> res;

            while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
                token = s.substr(pos_start, pos_end - pos_start);
                pos_start = pos_end + delim_len;
                res.push_back(token);
            }

            res.push_back(s.substr(pos_start));
            return res;
        }

        /**
         * \brief Reorder the columns of an array given a bijective mapping `order`.
         * 
         * Given an n x m matrix, `mat`, reorder the columns of `mat` such that 
         * column j -> order[j].
        */
        template<typename T>
        inline void reorder_columns(mtr::CArray<T>& mat, const int* order) {
            std::vector<T> temp;
            temp.resize(mat.dims(1));
            for (size_t i = 0; i < mat.dims(0); i++) {
                for (size_t j = 0; j < temp.size(); j++)
                    temp[order[j]] = mat(i, j);
                
                for (size_t j = 0; j < temp.size(); j++)
                    mat(i, j) = temp[j];
            }
        }
        
    
        inline std::vector<int> construct_inverse(const std::vector<int>& map) {
            std::vector<int> inv;
            inv.resize(map.size());
            for (size_t i = 0; i < inv.size(); i++)
                inv[map[i]] = i;
            return inv;
        }

        static const std::vector<int>& ijk_to_fea() {
            static const std::vector<int> map { 0, 1, 3, 2, 4, 5, 7, 6 };
            return map;
        }
        
        static const std::vector<int>& fea_to_ijk() {
            static const std::vector<int> map = construct_inverse(ijk_to_fea());
            return map;
        }
    }
}