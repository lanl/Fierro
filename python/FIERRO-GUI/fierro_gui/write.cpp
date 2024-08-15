#include <fstream>
#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void c_write(const py::array_t<bool> array, const char* file_name) {
    std::ofstream fout(file_name, std::ios_base::app);
    
    // RGB values:   0,   0,   0 is black
    //             255, 255, 255 is white
    for (py::ssize_t i = 0; i < array.size(); i++){
        fout << (array.data()[i] ? '1' : '0') << ' ';
    }
}

PYBIND11_MODULE(c_writer, m) {
    m.doc() = "";

    m.def("write", &c_write, "Write a numpy array to a file");
}
