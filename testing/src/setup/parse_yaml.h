

#ifndef FIERRO_PARSE_YAML_H
#define FIERRO_PARSE_YAML_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>


#include "Yaml.hpp"
// #include "material.h"
// #include "region.h"
// #include "mesh_inputs.h"
// #include "solver_inputs.h"
// #include "output_options.h"

struct solver_input_t;
struct mesh_input_t;
struct reg_fill_t;
struct material_t;
struct output_options_t;
struct boundary_condition_t;


// checks to see if a path exists
static bool DoesPathExist(const std::string &s)
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

// Parse the solver related data
void parse_solver_input(Yaml::Node &root, std::vector <solver_input_t> &solver_input);

// Parse the mesh related data
void parse_mesh_input(Yaml::Node &root, mesh_input_t &mesh_input);

// Parse output options
void parse_output_options(Yaml::Node &root, output_options_t &output_options);

// parse the region text
void parse_regions(Yaml::Node &root, std::vector <reg_fill_t> &region_fills);


// parse the region text
void parse_materials(Yaml::Node &root,
                     std::vector <material_t> &materials,
                     std::vector <std::vector <double>> &eos_global_vars);


// parse the boundary condition text
void parse_bcs(Yaml::Node &root, std::vector <boundary_condition_t> &boundary_conditions);

#endif // end Header Guard