

#ifndef FIERRO_PARSE_YAML_H
#define FIERRO_PARSE_YAML_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>


#include "Yaml.hpp"
#include "material.h"
#include "region.h"

// checks to see if a path exists
bool DoesPathExist(const std::string &s)
{
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
}


// for string delimiter parsing
std::vector<std::string> exact_array_values (std::string s, std::string delimiter);

// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str);

// prints the contents of a parsed yaml file
void print_yaml(Yaml::Node root);


// parse the region text
void parse_regions(Yaml::Node &root, std::vector <reg_fill_t> &region_fills);


// parse the region text
void parse_materials(Yaml::Node &root,
                     std::vector <material_t> &materials,
                     std::vector <std::vector <double>> &eos_global_vars);



#endif // end Header Guard