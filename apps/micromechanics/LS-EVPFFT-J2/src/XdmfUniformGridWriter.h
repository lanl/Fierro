#pragma once

#include <cstdio>
#include <array>
#include <string>
#include <stdexcept>

struct XdmfUniformGridWriter
{
  int N1;
  int N2;
  int N3;
  FILE *xdmf;

  XdmfUniformGridWriter(int N1_, int N2_, int N3_);
  ~XdmfUniformGridWriter();
  void open_file(const char* xdmf_filename);
  void write_header(const char* gridName);
  void write_attribute(const char* aName, const char* aType, const char* aNumberType,
                       const char* aPrecision, const char* aLocation);
  void write_footer();
  void close_file();

  template <typename T>
  void deduce_attribute_number_type(T *data_ptr, std::string &aNumberType, std::string &aPrecision);
};

template <typename T>
void XdmfUniformGridWriter::deduce_attribute_number_type(T *data_ptr, std::string &aNumberType, std::string &aPrecision)
{
    if constexpr (std::is_same<T, int>::value) {
        aNumberType = "Int";
        aPrecision = "4";
    }
    else if constexpr (std::is_same<T, float>::value) {
        aNumberType = "Float";
        aPrecision = "4";
    }
    else if constexpr (std::is_same<T, double>::value) {
        aNumberType = "Float";
        aPrecision = "8";
    }
    else {
        throw std::runtime_error("attribute_number_type cannot be deduced.\n");
    }
}
