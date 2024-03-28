// //==============================================================================
// //
// // yaml-parsing-test.cpp
// //
// //  Writen by: Nathaniel Morgan, March 21, 2024
// //
// //==============================================================================
// //
// //  if building standalone
// //
// //     g++ --std=c++17 yaml-parsing-test.cpp Yaml.cpp
// //
// //==============================================================================

// #include <iostream>
// #include <sstream>
// #include <fstream>
// #include <string>
// #include <stdio.h>
// #include <math.h>
// #include <sys/stat.h>
// #include <vector>
// #include <variant>
// #include <algorithm>
// #include <map>

// #include "matar.h"
// #include "Yaml.hpp"

// #define PI 3.141592653589793

// using namespace mtr;


// bool VERBOSE = false;


// // save size_t type vars
// void exact_region_info(Yaml::Node &root, size_t reg_id, std::string name, size_t &val){
//     val = root["regions"][reg_id]["fill_volume"][name].As<size_t>();
//     std::cout << "\tval=" << val << std::endl;
// }; 

// // save double type vars
// void exact_region_info(Yaml::Node &root, size_t reg_id, std::string name, double &val){
//     val = root["regions"][reg_id]["fill_volume"][name].As<double>();
//     std::cout << "\tval=" << val << std::endl;
// }; 


// #define populate_fill_volume(str, varname)\
// do { \
//     if ( str.compare(#varname) == 0)  { \
//         std::cout << "\t=========" << std::endl; \
//         exact_region_info(root, reg_id, #varname, region_fills[reg_id].varname); \
//         std::cout << "\t=========" << std::endl; \
//         std::cout << "\t=========" << std::endl; \
//     } \
// } while(0) 
//         //std::cout << "hello = " << xtractstr(varname) << std::endl; \
//         // int val = root["regions"][reg_id]["fill_volume"][#varname].As<decltype(region_fills[reg_id].varname)>(); \
//         // region_fills[reg_id].varname = val; \
//         // std::cout << "\t string ==" << str << std::endl; \
//         //std::cout << "\tval = " << val << std::endl; \

//         // use macro
// //            for(auto given_word:str_region_inps){
// //                //std::string some_word = "material_id";
// //                populate_fill_volume(a_word, material_id);
// //            }



// //==============================================================================
// //   Fierro enums and structs to populate
// //==============================================================================
// namespace region
// {

//     // for tagging boundary faces
//     enum vol_tag
//     {
//         global = 0,     // tag every elements in the mesh
//         box = 1,        // tag all elements inside a box
//         cylinder = 2,   // tag all elements inside a cylinder
//         sphere = 3,     // tag all elements inside a sphere
//         readVoxelFile = 4,       // tag all elements in a voxel mesh input
//         planes = 5,     // tag all elements between two planes
//     };

// } // end of namespace

// namespace init_conds
// {
    
//     // applying initial conditions
//     enum init_velocity_conds
//     {
//         // uniform
//         cartesian = 0,   // cart velocity
//         radial = 1,      // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
//         spherical = 2,   // spherical
    
//         // linear variation
//         radial_linear = 3,     // linear variation from 0,0,0
//         spherical_linear = 4,   // linear variation from 0,0,0
    
//         // vortical initial conditions
//         tg_vortex = 5
//     };
    
// } // end of initial conditions namespace


// // fill instructions (was called mat_fill_t)
// struct reg_fill_t {
    
//     // type
//     region::vol_tag volume; // global, box, sphere, planes, etc.
    
//     // material id
//     size_t material_id;
    
//     // planes
//     double x1;
//     double x2;
//     double y1;
//     double y2;
//     double z1;
//     double z2;
    
//     // radius
//     double radius1;
//     double radius2;

    
//     // initial conditions
//     init_conds::init_velocity_conds velocity;
    
//     // velocity coefficients by component
//     double u,v,w;
    
//     // velocity magnitude for radial velocity initialization
//     double speed;
    
//     double ie;   // exstenive internal energy
//     double sie;  // specific internal energy
//     double den;  // density
// };



// namespace model
// {

//     // strength model types
//     enum strength_tag
//     {
//         none = 0,
//         hypo = 1,     // hypoelastic plastic model
//         hyper = 2,    // hyperelastic plastic model
//     };

// } // end of namespace

// namespace model_init
// {

//     // strength model setup
//     enum strength_setup_tag
//     {
//         input = 0,
//         user_init = 1,
//     };

// } // end of namespace


// //WARNING: placeholder
// void ideal_gas(double pres, double den, double sie){
//     // do nothing
//     std::cout << "hello from ideal_gas! Replace with actual EOS!" << std::endl;
// };

// //WARNING: placeholder
// void elastic_plastic(double stress, double strain){
//     // do nothing
//     std::cout << "hello from elastic_plastic! Replace with actual strength model!" << std::endl;
// }

// // material model parameters
// struct material_t {

//     size_t id;
    
//     // statev(0) = gamma
//     // statev(1) = minimum sound speed
//     // statev(2) = specific heat c_v
//     // statev(3) = ref temperature
//     // statev(4) = ref density
//     // statev(5) = ref specific internal energy
    
//     // eos fcn pointer
//     void (*eos_model)(double, double, double); // WARNING: a placeholder
    
//     // strength fcn pointer
//     void (*strength_model)(double, double, double); // WARNING: a placeholder
    
//     // hypo or hyper elastic plastic model
//     model::strength_tag strength_type;
    
//     // setup the strength model via the input file for via a user_setup
//     model_init::strength_setup_tag strength_setup=model_init::input;
    
//     size_t num_eos_state_vars;
//     size_t num_strength_state_vars;
//     size_t num_eos_global_vars;
//     size_t num_strength_global_vars;
    
//     double q1;    // acoustic coefficient in Riemann solver for compresion
//     double q1ex;  // acoustic coefficient in Riemann solver for expansion
//     double q2;    // linear coefficient in Riemann solver for compression
//     double q2ex;  // linear coefficient in Riemann solver for expansion
// }; // end material_t






// //==============================================================================
// //   Functions
// //==============================================================================


// // checks to see if a path exists
// bool DoesPathExist(const std::string &s)
// {
//     struct stat buffer;
//     return (stat (s.c_str(), &buffer) == 0);
// }


// // for string delimiter parsing
// std::vector<std::string> exact_array_values (std::string s, std::string delimiter);

// // retrieves multiple values between [ ]
// std::vector<double> extract_list(std::string str);

// // prints the contents of a parsed yaml file
// void print_yaml(Yaml::Node root);


// // parse the region text
// void parse_regions(Yaml::Node &root, std::vector <reg_fill_t> &region_fills);


// // parse the region text
// void parse_materials(Yaml::Node &root,
//                      std::vector <material_t> &materials,
//                      std::vector <std::vector <double>> &eos_global_vars);



// //==============================================================================
// //   Valid inputs to Fierro
// //==============================================================================

// // ----------------------------------
// // valid inputs for a material fill
// // 
// //   materials_text_inp["words"]
// //
// std::vector <std::string> str_material_inps
// {
//     "id",
//     "eos_model",
//     "strength_model",
//     "q1",
//     "q2",
//     "q1ex",
//     "q2ex",
//     "eos_global_vars"
// };


// // add the eos models here
// typedef void (*eos_type)(double, double, double);
// std::map <std::string, eos_type> eos_map
// {
//     {"ideal_gas", ideal_gas}
// };

// // add the strength models here
// typedef void (*strength_type)(double, double);
// std::map <std::string, strength_type> strength_map
// {
//     {"elastic_plastic", elastic_plastic}
// };



// // ----------------------------------
// // valid inputs for a material fill
// //
// //   fill_volume_text_inp["words"]
// //
// std::vector <std::string> str_region_inps
// {
//     "type",
//     "material_id",
//     "x1",
//     "x2",
//     "x1",
//     "x2",
//     "y1",
//     "y2",
//     "z1",
//     "z2",
//     "radius1",
//     "radius2",
//     "velocity",
//     "u",
//     "v",
//     "w",
//     "speed",
//     "sie",
//     "ie",
//     "den"
// }; //






// std::map <std::string, region::vol_tag> region_type_map
// {
//     {"global",   region::global},
//     {"sphere",   region::sphere},
//     {"planes",   region::planes},
//     {"cylinder", region::cylinder},
//     {"readVoxelFile", region::readVoxelFile}
// };

// std::map <std::string, init_conds::init_velocity_conds> velocity_type_map
// {
//     {"cartesian",        init_conds::cartesian},
//     {"radial",           init_conds::radial},
//     {"spherical",        init_conds::spherical},
//     {"radial_linear",    init_conds::radial_linear},
//     {"spherical_linear", init_conds::spherical_linear},
//     {"tg_vortex",        init_conds::tg_vortex}
// };







// //==============================================================================
// //    Main
// //==============================================================================

// int main(int argc, char *argv[])
// {
    
//     // check to see of a mesh was supplied when running the code
//     if (argc == 1) {
//         std::cout << "\n\n**********************************\n\n";
//         std::cout << " ERROR:\n";
//         std::cout << " Please supply a YAML input, \n";
//         std::cout << "   ./mesh-builder input.yaml \n\n";
//         std::cout << "**********************************\n\n" << std::endl;
//         return 0;
//     } // end if
    
    
    
//     Yaml::Node root;
//     try
//     {
//         Yaml::Parse(root, argv[1]);
//     }
//     catch (const Yaml::Exception e)
//     {
//         std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
//         return 0;
//     }
    
//     // print the input file
//     print_yaml(root);
    

//     // parse the region yaml text into a vector of region_fills
//     std::vector <reg_fill_t> region_fills;
//     parse_regions(root, region_fills);

        
        
    
//     // parse the material yaml text into a vector of materials
//     std::vector <material_t> materials;
//     std::vector <std::vector <double>> eos_global_vars;
//     parse_materials(root, materials, eos_global_vars);


    
//     std::cout << "Done " << std::endl;
    
//     return 0;
    
// } // end of main



// // modified code from stackover flow for string delimiter parsing
// std::vector<std::string> exact_array_values (std::string s, std::string delimiter) {
//     size_t pos_start = 0, pos_end, delim_len = delimiter.length();
//     std::string token;
//     std::vector<std::string> res;

//     // remove first and last char in the string, which are [ and ]
//     s.erase(s.begin());
//     s.erase(s.end()-1);

//     // now parse the values in the array into a vector
//     while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
//         token = s.substr (pos_start, pos_end - pos_start);
//         pos_start = pos_end + delim_len;
//         res.push_back (token);
//     }

//     res.push_back (s.substr (pos_start));
//     return res;
    
// } // end of split



// // a function to print a yaml file to 6 levels
// void print_yaml(Yaml::Node root){

//     Yaml::Node & layer0_items = root;

//     if (layer0_items.Size()!=0){

//         std::cout << "\n";

//         if(VERBOSE) std::cout << "Layer 0 member size = " << layer0_items.Size() << "\n";

//         for(auto layer0_item = layer0_items.Begin(); layer0_item != layer0_items.End(); layer0_item++)
//         {

//             // print the outlayer variable
//             std::cout << (*layer0_item).first << "\n";
        
//             Yaml::Node & layer1_items = (*layer0_item).second;

//             // layer 1
//             if (layer1_items.Size()!=0){

//                 size_t count_layer_1 = 0;

//                 if(VERBOSE) std::cout << "\t num_items in this layer = " << layer1_items.Size() << "\n";

//                 for(auto layer1_item = layer1_items.Begin(); layer1_item != layer1_items.End(); layer1_item++)
//                 {

//                     std::string text_here1 = (*layer1_item).first;
//                     if(text_here1.size()>0) std::cout << "\t " << text_here1 << "\n";
//                     else std::cout << "\t " << count_layer_1 << "\n";

//                     // layer 2
//                     Yaml::Node & layer2_items = (*layer1_item).second;
                    
//                     if (layer2_items.Size()!=0){

//                         size_t count_layer_2 = 0; 

//                         if(VERBOSE) std::cout << "\t\t num_items in this layer = " << layer2_items.Size() << "\n";

//                         for(auto layer2_item = layer2_items.Begin(); layer2_item !=  layer2_items.End(); layer2_item++)
//                         {

//                             std::string text_here2 = (*layer2_item).first;
//                             if(text_here2.size()>0) std::cout << "\t\t " << text_here2 << std::endl;
//                             else std::cout << "\t\t " << count_layer_2 << "\n";

//                             // layer 3
//                             Yaml::Node & layer3_items = (*layer2_item).second;

//                             if (layer3_items.Size()!=0){    

//                                 size_t count_layer_3 = 0; 

//                                 if(VERBOSE) std::cout << "\t\t\t num_items in this layer = " << layer3_items.Size() << "\n";

//                                 for(auto layer3_item = layer3_items.Begin(); layer3_item !=  layer3_items.End(); layer3_item++)
//                                 {

//                                     std::string text_here3 = (*layer3_item).first;
//                                     if(text_here3.size()>0) std::cout << "\t\t\t " << text_here3 << std::endl;
//                                     else std::cout << "\t\t\t " << count_layer_3 << "\n";

//                                     // layer 4
//                                     Yaml::Node & layer4_items = (*layer3_item).second;
                                    
//                                     if (layer4_items.Size()!=0){ 

//                                         size_t count_layer_4 = 0; 

//                                         if(VERBOSE) std::cout << "\t\t\t\t num_items in layer 4 = " << layer4_items.Size() << "\n";

//                                         for(auto layer4_item = layer4_items.Begin(); layer4_item !=  layer4_items.End(); layer4_item++)
//                                         {
                                        
//                                             std::string text_here4 = (*layer4_item).first;
//                                             if(text_here4.size()>0) std::cout << "\t\t\t\t " << text_here4 <<std::endl; 
//                                             else std::cout << "\t\t\t\t " << count_layer_4 << "\n";
                                            
//                                             // layer 5
//                                             Yaml::Node & layer5_items = (*layer4_item).second;

//                                             if (layer5_items.Size()!=0){    

//                                                 size_t count_layer_5 = 0; 
//                                                 if(VERBOSE) std::cout << "\t\t\t\t\t num_items in layer 5 = " << layer5_items.Size() << "\n";

//                                                 for(auto layer5_item = layer5_items.Begin(); layer5_item !=  layer5_items.End(); layer5_item++)
//                                                 {
                                                
//                                                     std::string text_here5 = (*layer5_item).first;
//                                                     if(text_here5.size()>0) std::cout << "\t\t\t\t\t " << text_here5 << std::endl;
//                                                     else std::cout << "\t\t\t\t\t " << count_layer_5 << "\n";

//                                                     // layer 6
//                                                     Yaml::Node & layer6_items = (*layer5_item).second;

//                                                     if (layer6_items.Size()!=0){  

//                                                         size_t count_layer_6 = 0;  
//                                                         if(VERBOSE) std::cout << "\t\t\t\t\t\t num_items in layer 6 = " << layer6_items.Size() << "\n";

//                                                         for(auto layer6_item = layer6_items.Begin(); layer6_item !=  layer6_items.End(); layer6_item++)
//                                                         {
                                                        
//                                                             std::string text_here6 = (*layer6_item).first;
//                                                             if(text_here6.size()>0) std::cout << "\t\t\t\t\t\t layer 6 = " << text_here6 << std::endl;
//                                                             else std::cout << "\t\t\t\t\t\t " << count_layer_6 << "\n";

//                                                             count_layer_6++;
//                                                         } // end loop over layer 6
//                                                     } // end of layer6

//                                                     count_layer_5++;
//                                                 } // end loop over layer 5
//                                             } // end if layer5 exists

//                                             count_layer_4++;
//                                         } // end loop over layer4
//                                     } // end if layer4 exists

//                                     count_layer_3 ++;
//                                 } // end loop over layer 3 items
//                             } // end if layer 3 exists

//                             count_layer_2 ++;
//                         } // end if loop over layer 2
//                     } // end if layer 2 exists

//                     count_layer_1 ++;
//                 } // end loop over layer 1
//             } // end if layer 1 exists

//         } // end loop over layer 0
//     } // end if layer0 exists

// } // end print yaml function




// // =================================================================================
// //    Fill regions
// // =================================================================================
// void parse_regions(Yaml::Node &root,
//                    std::vector <reg_fill_t> &region_fills){

//     Yaml::Node & region_yaml = root["regions"];
    
//     size_t num_regions = region_yaml.Size();
    
//     region_fills = std::vector <reg_fill_t>(num_regions);
    
//     // loop over the fill regions specified
//     for(int reg_id=0; reg_id<num_regions; reg_id++){
        
        
//         // read the variables names
//         Yaml::Node & inps_yaml = root["regions"][reg_id]["fill_volume"];
        
        
//         // get the material variables names set by the user
//         std::vector <std::string> user_str_region_inps;
        
        
//         // extract words from the input file and validate they are correct
//         for(auto item = inps_yaml.Begin(); item != inps_yaml.End(); item++)
//         {
            
//             std::string var_name = (*item).first;
            
//             // print the variable
//             std::cout << "this is var name = "<< var_name << "\n";
            
//             user_str_region_inps.push_back(var_name);
            
//             // validate input: user_str_region_inps match words in the str_region_inps
//             if (std::find(str_region_inps.begin(), str_region_inps.end(), var_name) == str_region_inps.end())
//             {
//                 std::cout << "ERROR: invalid input: " << var_name << std::endl;
//             } // end if variable exists
            
//         } // end for item in this yaml input
        
        
        
//         // loop over the words in the material input definition
//         for(auto &a_word : user_str_region_inps){
            
//             std::cout << a_word << std::endl;
        
//             Yaml::Node & material_inps_yaml = root["regions"][reg_id]["fill_volume"][a_word];
            

//             // set the values
//             if(a_word.compare("material_id") == 0){
//                 region_fills[reg_id].material_id = 
//                     root["regions"][reg_id]["fill_volume"][a_word].As<int>();
            
//             } // mat_id
//             else if(a_word.compare("den") == 0){

//                 double den = root["regions"][reg_id]["fill_volume"]["den"].As<double>();
                
//                 // check for a valid density else save it
//                 if(den < 0.0){
//                     std::cout << "ERROR: density is negative: " << den << std::endl;
//                 } else {
//                     region_fills[reg_id].den = den;  // NOTE: GPUs will require a RUN({})
//                 }
            
//             } // den
//             else if(a_word.compare("sie") == 0){
                
//                 // specific internal energy

//                 double sie = root["regions"][reg_id]["fill_volume"]["sie"].As<double>();
//                 std::cout << "\tsie = " << sie << std::endl;
                
//                 region_fills[reg_id].sie = sie;
            
//             } // sie
//             else if(a_word.compare("ie") == 0){
                
//                 // extensive internal energy

//                 double ie = root["regions"][reg_id]["fill_volume"]["ie"].As<double>();
//                 std::cout << "\tie = " << ie << std::endl;
                
//                 region_fills[reg_id].ie = ie;
            
//             } // ie
//             else if(a_word.compare("speed") == 0){

//                 double speed = root["regions"][reg_id]["fill_volume"]["speed"].As<double>();
//                 std::cout << "\tspeed = " << speed << std::endl;
                
//                 region_fills[reg_id].speed = speed;
            
//             } // speed
//             else if(a_word.compare("u") == 0){

//                 // x-component of velocity
                
//                 double u = root["regions"][reg_id]["fill_volume"]["u"].As<double>();
//                 std::cout << "\tu = " << u << std::endl;
                
//                 region_fills[reg_id].u = u;
            
//             } // u
//             else if(a_word.compare("v") == 0){
                
//                 // y-component of velocity

//                 double v = root["regions"][reg_id]["fill_volume"]["v"].As<double>();
//                 std::cout << "\tie = " << v << std::endl;
                
//                 region_fills[reg_id].v = v;
            
//             } // v
//             else if(a_word.compare("w") == 0){
                
//                 // z-component of velocity

//                 double w = root["regions"][reg_id]["fill_volume"]["w"].As<double>();
//                 std::cout << "\tw = " << w << std::endl;
                
//                 region_fills[reg_id].w = w;
            
//             } // w
//             else if(a_word.compare("radius1") == 0){
                
//                 // inner radius of sphere/cylinder

//                 double radius1 = root["regions"][reg_id]["fill_volume"]["radius1"].As<double>();
//                 std::cout << "\tradius1 = " << radius1 << std::endl;
                
//                 region_fills[reg_id].radius1 = radius1;
            
//             } // radius1
//             else if(a_word.compare("radius2") == 0){
                
//                 // outer radius of sphere/cylinder

//                 double radius2 = root["regions"][reg_id]["fill_volume"]["radius2"].As<double>();
//                 std::cout << "\tradius2 = " << radius2 << std::endl;
                
//                 region_fills[reg_id].radius2 = radius2;
            
//             } // radius2
//             else if(a_word.compare("x1") == 0){
                
//                 // inner plane

//                 double x1 = root["regions"][reg_id]["fill_volume"]["x1"].As<double>();
//                 std::cout << "\tx1 = " << x1 << std::endl;
                
//                 region_fills[reg_id].x1 = x1;
            
//             } // x1
//             else if(a_word.compare("x2") == 0){
                
//                 // outer plane

//                 double x2 = root["regions"][reg_id]["fill_volume"]["x2"].As<double>();
//                 std::cout << "\tx2 = " << x2 << std::endl;
                
//                 region_fills[reg_id].x2 = x2;
            
//             } // x2
//             else if(a_word.compare("y1") == 0){
                
//                 // inner plane

//                 double y1 = root["regions"][reg_id]["fill_volume"]["y1"].As<double>();
//                 std::cout << "\ty1 = " << y1 << std::endl;
                
//                 region_fills[reg_id].y1 = y1;
            
//             } // y1
//             else if(a_word.compare("y2") == 0){
                
//                 // outer plane

//                 double y2 = root["regions"][reg_id]["fill_volume"]["y2"].As<double>();
//                 std::cout << "\ty2 = " << y2 << std::endl;
                
//                 region_fills[reg_id].y2 = y2;
            
//             } // y2
//             else if(a_word.compare("z1") == 0){
                
//                 // inner plane

//                 double z1 = root["regions"][reg_id]["fill_volume"]["z1"].As<double>();
//                 std::cout << "\tz1 = " << z1 << std::endl;
                
//                 region_fills[reg_id].z1 = z1;
            
//             } // z1
//             else if(a_word.compare("z2") == 0){
                
//                 // outer plane

//                 double z2 = root["regions"][reg_id]["fill_volume"]["z2"].As<double>();
//                 std::cout << "\tz2 = " << z2 << std::endl;
                
//                 region_fills[reg_id].z2 = z2;
            
//             } // z2
//             else if(a_word.compare("type") == 0){

//                 std::string type = root["regions"][reg_id]["fill_volume"]["type"].As<std::string>();
//                 std::cout << "\ttype = " << type << std::endl;
                
//                 // set the volume tag type
//                 if(region_type_map.find(type) != region_type_map.end()){
                    
//                     region_fills[reg_id].volume = region_type_map[type];
//                     std::cout << "\tvolume_fill = " << type << std::endl;
//                     std::cout << region_fills[reg_id].volume << std::endl;
//                 }
//                 else{
//                     std::cout << "ERROR: invalid input: " << type << std::endl;
//                 } // end if
                    
//             } // end volume fill type
//             else if(a_word.compare("velocity") == 0){

//                 std::string type = root["regions"][reg_id]["fill_volume"]["velocity"].As<std::string>();
//                 std::cout << "\tvelocity = " << type << std::endl;
                
                
//                 // set the volume tag type
//                 if(velocity_type_map.find(type) != velocity_type_map.end()){
                    
//                     region_fills[reg_id].velocity = velocity_type_map[type];
//                     std::cout << "\tvelocity_fill = " << type << std::endl;
//                     std::cout << region_fills[reg_id].velocity << std::endl;
//                 }
//                 else{
//                     std::cout << "ERROR: invalid input: " << type << std::endl;
//                 } // end if
                    
//             } // end velocity
//             else if(a_word.compare("origin") == 0){

//                 std::string origin = root["regions"][reg_id]["fill_volume"]["origin"].As<std::string>();
//                 std::cout << "\torigin = " << origin << std::endl;
                
//                 // get the origin numbers, values are words
//                 std::vector<std::string> numbers = exact_array_values(origin, ",");
                
//                 double x1 = std::stod(numbers[0]);
//                 double y1 = std::stod(numbers[1]);
//                 double z1 = std::stod(numbers[2]);
                
//                 std::cout << "\tx1 = " << x1 << std::endl;
//                 std::cout << "\ty1 = " << y1 << std::endl;
//                 std::cout << "\tz1 = " << z1 << std::endl;
                
//                 // storing the origin values as (x1,y1,z1)
//                 region_fills[reg_id].x1 = x1;
//                 region_fills[reg_id].y1 = y1;
//                 region_fills[reg_id].z1 = z1;
                    
//             } // origin
//             else {
//                 std::cout << "ERROR: invalid input: " << a_word << std::endl;
//             }
            
            
//         } // end for words in material
    
//     } // end loop over regions
    
// } // end of function to parse region



// // =================================================================================
// //      material definitions
// // =================================================================================
// void parse_materials(Yaml::Node &root,
//                      std::vector <material_t> &materials,
//                      std::vector <std::vector <double>> &eos_global_vars){


//     Yaml::Node & material_yaml = root["materials"];

//     size_t num_materials = material_yaml.Size();
    
//     materials = std::vector <material_t>(num_materials);
    
//     // allocate room for each material to store eos_global_vars
//     eos_global_vars = std::vector <std::vector <double>>(num_materials);

//     // loop over the materials specified
//     for(int mat_id=0; mat_id<num_materials; mat_id++){
        
        
//         // read the variables names
//         Yaml::Node & inps_yaml = root["materials"][mat_id]["material"];
        
//         size_t num_vars_set = inps_yaml.Size();
        
//         // get the material variables names set by the user
//         std::vector <std::string> user_str_material_inps;
        
        
//         // extract words from the input file and validate they are correct
//         for(auto item = inps_yaml.Begin(); item != inps_yaml.End(); item++)
//         {
            
//             std::string var_name = (*item).first;
            
//             // print the variable
//             std::cout << "this is var name = "<< var_name << "\n";
            
//             user_str_material_inps.push_back(var_name);
            
//             // validate input: user_str_material_inps match words in the str_material_inps
//             if (std::find(str_material_inps.begin(), str_material_inps.end(), var_name) == str_material_inps.end())
//             {
//                 std::cout << "ERROR: invalid input: " << var_name << std::endl;
//             }
            
//         } // end for item in this yaml input
        
        
//         // loop over the words in the material input definition
//         for(auto &a_word : user_str_material_inps){
            
//             std::cout << a_word << std::endl;
        
//             Yaml::Node & material_inps_yaml = root["materials"][mat_id]["material"][a_word];
            
            
//             // set the values in the input for this word
            
            
//             // set the values
//             if(a_word.compare("q1") == 0){
                
//                 // outer plane

//                 double q1 = root["materials"][mat_id]["material"]["q1"].As<double>();
//                 std::cout << "\tq1 = " << q1 << std::endl;
                
//                 materials[mat_id].q1 = q1;
            
//             } // q1
//             else if(a_word.compare("q1ex") == 0){
                
//                 // outer plane

//                 double q1ex = root["materials"][mat_id]["material"]["q1ex"].As<double>();
//                 std::cout << "\tq1ex = " << q1ex << std::endl;
                
//                 materials[mat_id].q1ex = q1ex;
            
//             } // q1ex
//             else if(a_word.compare("q2") == 0){
                
//                 // outer plane

//                 double q2 = root["materials"][mat_id]["material"]["q2"].As<double>();
//                 std::cout << "\tq2 = " << q2 << std::endl;
                
//                 materials[mat_id].q2 = q2;
            
//             } // q1
//             else if(a_word.compare("q2ex") == 0){
                
//                 // outer plane

//                 double q2ex = root["materials"][mat_id]["material"]["q2ex"].As<double>();
//                 std::cout << "\tq2ex = " << q2ex << std::endl;
                
//                 materials[mat_id].q2ex = q2ex;
            
//             } // q1ex
//             else if(a_word.compare("id") == 0){
                  
//                 int m_id = root["materials"][mat_id]["material"]["id"].As<int>();
//                 std::cout << "\tid = " << m_id << std::endl;
            
//                 materials[mat_id].id = m_id;
            
//             } // id
//             else if(a_word.compare("eos_model") == 0){
                  
//                 std::string eos = root["materials"][mat_id]["material"]["eos_model"].As<std::string>();
                
//                 // set the EOS
//                 if(eos_map.find(eos) != eos_map.end()){
//                     materials[mat_id].eos_model = eos_map[eos];
//                     std::cout << "\teos_model = " << eos << std::endl;
//                     materials[mat_id].eos_model(0.,1.,2.);
//                 }
//                 else{
//                     std::cout << "ERROR: invalid input: " << eos << std::endl;
//                 } // end if
            
//             } // id
//             // exact the eos_global_vars
//             else if(a_word.compare("eos_global_vars") == 0){
                
//                 size_t num_global_vars = material_inps_yaml.Size();
                
//                 std::cout << "num global eos vars = " << num_global_vars << std::endl;
                
//                 for(int global_var_id=0; global_var_id<num_global_vars; global_var_id++){
                    
//                     double eos_var = root["materials"][mat_id]["material"]["eos_global_vars"][global_var_id].As<double>();
                    
//                     eos_global_vars[mat_id].push_back(eos_var);
                    
//                     std::cout << "\t var = " << eos_var<< std::endl;
//                 } // end loop over global vars
                
//             } // "eos_global_vars"
//             else {
//                 std::cout << "ERROR: invalid input: " << a_word << std::endl;
//             }
            
            
//         } // end for words in material
        
//     } // end loop over materials

// } // end of function to parse material information

// wrap long lines