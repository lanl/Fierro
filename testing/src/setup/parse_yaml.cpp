
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <variant>
#include <algorithm>
#include <map>

#include "matar.h"
#include "parse_yaml.h"




#define PI 3.141592653589793

using namespace mtr;


bool VERBOSE = false;




//==============================================================================
//   Function Definitions
//==============================================================================


// modified code from stackover flow for string delimiter parsing
std::vector<std::string> exact_array_values (std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    // remove first and last char in the string, which are [ and ]
    s.erase(s.begin());
    s.erase(s.end()-1);

    // now parse the values in the array into a vector
    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
    
} // end of extract_array_values

// a function to print a yaml file to 6 levels
void print_yaml(Yaml::Node root){

    Yaml::Node & layer0_items = root;

    if (layer0_items.Size()!=0){

        std::cout << "\n";

        if(VERBOSE) std::cout << "Layer 0 member size = " << layer0_items.Size() << "\n";

        for(auto layer0_item = layer0_items.Begin(); layer0_item != layer0_items.End(); layer0_item++)
        {

            // print the outlayer variable
            std::cout << (*layer0_item).first << "\n";
        
            Yaml::Node & layer1_items = (*layer0_item).second;

            // layer 1
            if (layer1_items.Size()!=0){

                size_t count_layer_1 = 0;

                if(VERBOSE) std::cout << "\t num_items in this layer = " << layer1_items.Size() << "\n";

                for(auto layer1_item = layer1_items.Begin(); layer1_item != layer1_items.End(); layer1_item++)
                {

                    std::string text_here1 = (*layer1_item).first;
                    if(text_here1.size()>0) std::cout << "\t " << text_here1 << "\n";
                    else std::cout << "\t " << count_layer_1 << "\n";

                    // layer 2
                    Yaml::Node & layer2_items = (*layer1_item).second;
                    
                    if (layer2_items.Size()!=0){

                        size_t count_layer_2 = 0; 

                        if(VERBOSE) std::cout << "\t\t num_items in this layer = " << layer2_items.Size() << "\n";

                        for(auto layer2_item = layer2_items.Begin(); layer2_item !=  layer2_items.End(); layer2_item++)
                        {

                            std::string text_here2 = (*layer2_item).first;
                            if(text_here2.size()>0) std::cout << "\t\t " << text_here2 << std::endl;
                            else std::cout << "\t\t " << count_layer_2 << "\n";

                            // layer 3
                            Yaml::Node & layer3_items = (*layer2_item).second;

                            if (layer3_items.Size()!=0){    

                                size_t count_layer_3 = 0; 

                                if(VERBOSE) std::cout << "\t\t\t num_items in this layer = " << layer3_items.Size() << "\n";

                                for(auto layer3_item = layer3_items.Begin(); layer3_item !=  layer3_items.End(); layer3_item++)
                                {

                                    std::string text_here3 = (*layer3_item).first;
                                    if(text_here3.size()>0) std::cout << "\t\t\t " << text_here3 << std::endl;
                                    else std::cout << "\t\t\t " << count_layer_3 << "\n";

                                    // layer 4
                                    Yaml::Node & layer4_items = (*layer3_item).second;
                                    
                                    if (layer4_items.Size()!=0){ 

                                        size_t count_layer_4 = 0; 

                                        if(VERBOSE) std::cout << "\t\t\t\t num_items in layer 4 = " << layer4_items.Size() << "\n";

                                        for(auto layer4_item = layer4_items.Begin(); layer4_item !=  layer4_items.End(); layer4_item++)
                                        {
                                        
                                            std::string text_here4 = (*layer4_item).first;
                                            if(text_here4.size()>0) std::cout << "\t\t\t\t " << text_here4 <<std::endl; 
                                            else std::cout << "\t\t\t\t " << count_layer_4 << "\n";
                                            
                                            // layer 5
                                            Yaml::Node & layer5_items = (*layer4_item).second;

                                            if (layer5_items.Size()!=0){    

                                                size_t count_layer_5 = 0; 
                                                if(VERBOSE) std::cout << "\t\t\t\t\t num_items in layer 5 = " << layer5_items.Size() << "\n";

                                                for(auto layer5_item = layer5_items.Begin(); layer5_item !=  layer5_items.End(); layer5_item++)
                                                {
                                                
                                                    std::string text_here5 = (*layer5_item).first;
                                                    if(text_here5.size()>0) std::cout << "\t\t\t\t\t " << text_here5 << std::endl;
                                                    else std::cout << "\t\t\t\t\t " << count_layer_5 << "\n";

                                                    // layer 6
                                                    Yaml::Node & layer6_items = (*layer5_item).second;

                                                    if (layer6_items.Size()!=0){  

                                                        size_t count_layer_6 = 0;  
                                                        if(VERBOSE) std::cout << "\t\t\t\t\t\t num_items in layer 6 = " << layer6_items.Size() << "\n";

                                                        for(auto layer6_item = layer6_items.Begin(); layer6_item !=  layer6_items.End(); layer6_item++)
                                                        {
                                                        
                                                            std::string text_here6 = (*layer6_item).first;
                                                            if(text_here6.size()>0) std::cout << "\t\t\t\t\t\t layer 6 = " << text_here6 << std::endl;
                                                            else std::cout << "\t\t\t\t\t\t " << count_layer_6 << "\n";

                                                            count_layer_6++;
                                                        } // end loop over layer 6
                                                    } // end of layer6

                                                    count_layer_5++;
                                                } // end loop over layer 5
                                            } // end if layer5 exists

                                            count_layer_4++;
                                        } // end loop over layer4
                                    } // end if layer4 exists

                                    count_layer_3 ++;
                                } // end loop over layer 3 items
                            } // end if layer 3 exists

                            count_layer_2 ++;
                        } // end if loop over layer 2
                    } // end if layer 2 exists

                    count_layer_1 ++;
                } // end loop over layer 1
            } // end if layer 1 exists

        } // end loop over layer 0
    } // end if layer0 exists

} // end print yaml function


// =================================================================================
//    Parse Mesh options
// =================================================================================
void parse_mesh_input(Yaml::Node &root, mesh_input_t &mesh_input){

    Yaml::Node & mesh_yaml = root["mesh_options"];
    
    // size_t num_meshes = mesh_yaml.Size();

    // if(num_meshes != 1){

    //     std::cout << "Num meshes =  "<< num_meshes<< std::endl;
    //     std::cout << "Fierro only operates on a single mesh (it does not have to be connected) "<< std::endl;
    //     exit(0);
    // }
    
    // get the mesh variables names set by the user
    std::vector <std::string> user_mesh_inputs;
    
    
    // extract words from the input file and validate they are correct
    for(auto item = mesh_yaml.Begin(); item != mesh_yaml.End(); item++)
    {
        
        std::string var_name = (*item).first;
        
        // print the variable
        std::cout << "This is var name = "<< var_name << "\n";
        
        user_mesh_inputs.push_back(var_name);
        
        // validate input: user_mesh_inputs match words in the str_mesh_inps
        if (std::find(str_mesh_inps.begin(), str_mesh_inps.end(), var_name) == str_mesh_inps.end())
        {
            std::cout << "ERROR: invalid input: " << var_name << std::endl;
        } // end if variable exists
        
    } // end for item in this yaml input
    

    // loop over the words in the material input definition
    for(auto &a_word : user_mesh_inputs){
        
        std::cout << a_word << std::endl;
    
        Yaml::Node & material_inps_yaml = root["mesh_options"][a_word];
        
        // get mesh source [generate or from file]
        if(a_word.compare("source") == 0){
            
            std::string source = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_source_map;
            
            // set the mesh source
            if(map.find(source) != map.end()){
                
                mesh_input.source = map[source];
                std::cout << "\tsource = " << source << std::endl;

                if(mesh_input.source == mesh_input::generate && !mesh_input.file_path.empty()){
                    std::cout << "ERROR: When the mesh source is set to generate, a mesh file cannot be passed in" << std::endl;
                    exit(0);
                }
            }
            else{
                std::cout << "ERROR: invalid mesh option input in YAML file: " << source << std::endl;
                std::cout << "Valid options are: "<< std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }

            } // end if
        } // source 

        // get mesh type for generation
        if(a_word.compare("type") == 0){
            
            std::string type = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_type_map;
            
            // set the mesh type
            if(map.find(type) != map.end()){
                
                mesh_input.type = map[type];
                std::cout << "\ttype = " << type << std::endl;
            }
            else{
                std::cout << "ERROR: invalid mesh option input in YAML file: " << type << std::endl;
                std::cout << "Valid options are: "<< std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
            } // end if
        } // type 
        
        // Get mesh file path 
        else if(a_word.compare("file_path") == 0){
            std::string path = root["mesh_options"][a_word].As<std::string>();
            std::cout << "\tfile_path = " << path << std::endl;
            
            mesh_input.file_path = path;

            if(mesh_input.source == mesh_input::file && mesh_input.file_path.empty()){
                std::cout << "ERROR: When the mesh source is a file, a file_path must be set to point to the mesh file" << std::endl;
                std::cout << "A mesh can either be generated or read in from a file, but not both"<< std::endl;
            }

            if(mesh_input.source == mesh_input::generate ){
                std::cout << "ERROR: When the mesh source is set to generate, a mesh file cannot be passed in" << std::endl;
                exit(0);
            }

        } // file path

        // Origin for the mesh
        else if(a_word.compare("origin") == 0){

            std::string origin = root["mesh_options"][a_word].As<std::string>();
            std::cout << "\torigin = " << origin << std::endl;
            
            // get the origin numbers, values are words
            std::vector<std::string> numbers = exact_array_values(origin, ",");


            std::vector<double> val;
            for(int i = 0; i < 3; i++){
                val.push_back(std::stod(numbers[i]));
            }

            mesh_input.origin = val;
        }

        // Extents of the mesh
        else if(a_word.compare("length") == 0){

            std::string origin = root["mesh_options"][a_word].As<std::string>();
            std::cout << "\tlength = " << origin << std::endl;
            
            // get the origin numbers, values are words
            std::vector<std::string> numbers = exact_array_values(origin, ",");


            std::vector<double> val;
            for(int i = 0; i < 3; i++){
                val.push_back(std::stod(numbers[i]));
            }

            mesh_input.length = val;
        }

        // Number of elements per direction
        else if(a_word.compare("num_elems") == 0){

            std::string origin = root["mesh_options"][a_word].As<std::string>();
            std::cout << "\tnum_elems = " << origin << std::endl;
            
            // get the origin numbers, values are words
            std::vector<std::string> numbers = exact_array_values(origin, ",");


            std::vector<int> val;
            for(int i = 0; i < 3; i++){
                val.push_back(std::stoi(numbers[i]));
            }

            mesh_input.num_elems = val;
        }

        // Polynomial order for the mesh
        else if(a_word.compare("polynomial_order") == 0){
            
            double p_order = root["mesh_options"][a_word].As<int>();
            std::cout << "\tPoly order = " << p_order << std::endl;
            
            mesh_input.p_order = p_order;
        
        } // polynomial order
} // end of function to parse region


// =================================================================================
//    Parse Fill regions
// =================================================================================
void parse_regions(Yaml::Node &root,
                   std::vector <reg_fill_t> &region_fills){

    Yaml::Node & region_yaml = root["regions"];
    
    size_t num_regions = region_yaml.Size();
    
    region_fills = std::vector <reg_fill_t>(num_regions);
    
    // loop over the fill regions specified
    for(int reg_id=0; reg_id<num_regions; reg_id++){
        
        
        // read the variables names
        Yaml::Node & inps_yaml = root["regions"][reg_id]["fill_volume"];
        
        
        // get the material variables names set by the user
        std::vector <std::string> user_str_region_inps;
        
        
        // extract words from the input file and validate they are correct
        for(auto item = inps_yaml.Begin(); item != inps_yaml.End(); item++)
        {
            
            std::string var_name = (*item).first;
            
            // print the variable
            std::cout << "this is var name = "<< var_name << "\n";
            
            user_str_region_inps.push_back(var_name);
            
            // validate input: user_str_region_inps match words in the str_region_inps
            if (std::find(str_region_inps.begin(), str_region_inps.end(), var_name) == str_region_inps.end())
            {
                std::cout << "ERROR: invalid input: " << var_name << std::endl;
            } // end if variable exists
            
        } // end for item in this yaml input
        
        
        
        // loop over the words in the material input definition
        for(auto &a_word : user_str_region_inps){
            
            std::cout << a_word << std::endl;
        
            Yaml::Node & material_inps_yaml = root["regions"][reg_id]["fill_volume"][a_word];
            

            // set the values
            if(a_word.compare("material_id") == 0){
                region_fills[reg_id].material_id = 
                    root["regions"][reg_id]["fill_volume"][a_word].As<int>();
            
            } // mat_id
            else if(a_word.compare("den") == 0){

                double den = root["regions"][reg_id]["fill_volume"]["den"].As<double>();
                
                // check for a valid density else save it
                if(den < 0.0){
                    std::cout << "ERROR: density is negative: " << den << std::endl;
                } else {
                    region_fills[reg_id].den = den;  // NOTE: GPUs will require a RUN({})
                }
            
            } // den
            else if(a_word.compare("sie") == 0){
                
                // specific internal energy

                double sie = root["regions"][reg_id]["fill_volume"]["sie"].As<double>();
                std::cout << "\tsie = " << sie << std::endl;
                
                region_fills[reg_id].sie = sie;
            
            } // sie
            else if(a_word.compare("ie") == 0){
                
                // extensive internal energy

                double ie = root["regions"][reg_id]["fill_volume"]["ie"].As<double>();
                std::cout << "\tie = " << ie << std::endl;
                
                region_fills[reg_id].ie = ie;
            
            } // ie
            else if(a_word.compare("speed") == 0){

                double speed = root["regions"][reg_id]["fill_volume"]["speed"].As<double>();
                std::cout << "\tspeed = " << speed << std::endl;
                
                region_fills[reg_id].speed = speed;
            
            } // speed
            else if(a_word.compare("u") == 0){

                // x-component of velocity
                
                double u = root["regions"][reg_id]["fill_volume"]["u"].As<double>();
                std::cout << "\tu = " << u << std::endl;
                
                region_fills[reg_id].u = u;
            
            } // u
            else if(a_word.compare("v") == 0){
                
                // y-component of velocity

                double v = root["regions"][reg_id]["fill_volume"]["v"].As<double>();
                std::cout << "\tie = " << v << std::endl;
                
                region_fills[reg_id].v = v;
            
            } // v
            else if(a_word.compare("w") == 0){
                
                // z-component of velocity

                double w = root["regions"][reg_id]["fill_volume"]["w"].As<double>();
                std::cout << "\tw = " << w << std::endl;
                
                region_fills[reg_id].w = w;
            
            } // w
            else if(a_word.compare("radius1") == 0){
                
                // inner radius of sphere/cylinder

                double radius1 = root["regions"][reg_id]["fill_volume"]["radius1"].As<double>();
                std::cout << "\tradius1 = " << radius1 << std::endl;
                
                region_fills[reg_id].radius1 = radius1;
            
            } // radius1
            else if(a_word.compare("radius2") == 0){
                
                // outer radius of sphere/cylinder

                double radius2 = root["regions"][reg_id]["fill_volume"]["radius2"].As<double>();
                std::cout << "\tradius2 = " << radius2 << std::endl;
                
                region_fills[reg_id].radius2 = radius2;
            
            } // radius2
            else if(a_word.compare("x1") == 0){
                
                // inner plane

                double x1 = root["regions"][reg_id]["fill_volume"]["x1"].As<double>();
                std::cout << "\tx1 = " << x1 << std::endl;
                
                region_fills[reg_id].x1 = x1;
            
            } // x1
            else if(a_word.compare("x2") == 0){
                
                // outer plane

                double x2 = root["regions"][reg_id]["fill_volume"]["x2"].As<double>();
                std::cout << "\tx2 = " << x2 << std::endl;
                
                region_fills[reg_id].x2 = x2;
            
            } // x2
            else if(a_word.compare("y1") == 0){
                
                // inner plane

                double y1 = root["regions"][reg_id]["fill_volume"]["y1"].As<double>();
                std::cout << "\ty1 = " << y1 << std::endl;
                
                region_fills[reg_id].y1 = y1;
            
            } // y1
            else if(a_word.compare("y2") == 0){
                
                // outer plane

                double y2 = root["regions"][reg_id]["fill_volume"]["y2"].As<double>();
                std::cout << "\ty2 = " << y2 << std::endl;
                
                region_fills[reg_id].y2 = y2;
            
            } // y2
            else if(a_word.compare("z1") == 0){
                
                // inner plane

                double z1 = root["regions"][reg_id]["fill_volume"]["z1"].As<double>();
                std::cout << "\tz1 = " << z1 << std::endl;
                
                region_fills[reg_id].z1 = z1;
            
            } // z1
            else if(a_word.compare("z2") == 0){
                
                // outer plane

                double z2 = root["regions"][reg_id]["fill_volume"]["z2"].As<double>();
                std::cout << "\tz2 = " << z2 << std::endl;
                
                region_fills[reg_id].z2 = z2;
            
            } // z2
            else if(a_word.compare("type") == 0){

                std::string type = root["regions"][reg_id]["fill_volume"]["type"].As<std::string>();
                std::cout << "\ttype = " << type << std::endl;
                
                // set the volume tag type
                if(region_type_map.find(type) != region_type_map.end()){
                    
                    region_fills[reg_id].volume = region_type_map[type];
                    std::cout << "\tvolume_fill = " << type << std::endl;
                    std::cout << region_fills[reg_id].volume << std::endl;
                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                } // end if
                    
            } // end volume fill type
            else if(a_word.compare("velocity") == 0){

                std::string type = root["regions"][reg_id]["fill_volume"]["velocity"].As<std::string>();
                std::cout << "\tvelocity = " << type << std::endl;
                
                
                // set the volume tag type
                if(velocity_type_map.find(type) != velocity_type_map.end()){
                    
                    region_fills[reg_id].velocity = velocity_type_map[type];
                    std::cout << "\tvelocity_fill = " << type << std::endl;
                    std::cout << region_fills[reg_id].velocity << std::endl;
                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                } // end if
                    
            } // end velocity
            else if(a_word.compare("origin") == 0){

                std::string origin = root["regions"][reg_id]["fill_volume"]["origin"].As<std::string>();
                std::cout << "\torigin = " << origin << std::endl;
                
                // get the origin numbers, values are words
                std::vector<std::string> numbers = exact_array_values(origin, ",");
                
                double x1 = std::stod(numbers[0]);
                double y1 = std::stod(numbers[1]);
                double z1 = std::stod(numbers[2]);
                
                std::cout << "\tx1 = " << x1 << std::endl;
                std::cout << "\ty1 = " << y1 << std::endl;
                std::cout << "\tz1 = " << z1 << std::endl;
                
                // storing the origin values as (x1,y1,z1)
                region_fills[reg_id].x1 = x1;
                region_fills[reg_id].y1 = y1;
                region_fills[reg_id].z1 = z1;
                    
            } // origin
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
            }
            
            
        } // end for words in material
    
    } // end loop over regions
    
} // end of function to parse region


// =================================================================================
//    Parse Material Definitions
// =================================================================================
void parse_materials(Yaml::Node &root,
                     std::vector <material_t> &materials,
                     std::vector <std::vector <double>> &eos_global_vars){


    Yaml::Node & material_yaml = root["materials"];

    size_t num_materials = material_yaml.Size();
    
    materials = std::vector <material_t>(num_materials);
    
    // allocate room for each material to store eos_global_vars
    eos_global_vars = std::vector <std::vector <double>>(num_materials);

    // loop over the materials specified
    for(int mat_id=0; mat_id<num_materials; mat_id++){
        
        
        // read the variables names
        Yaml::Node & inps_yaml = root["materials"][mat_id]["material"];
        
        size_t num_vars_set = inps_yaml.Size();
        
        // get the material variables names set by the user
        std::vector <std::string> user_str_material_inps;
        
        
        // extract words from the input file and validate they are correct
        for(auto item = inps_yaml.Begin(); item != inps_yaml.End(); item++)
        {
            
            std::string var_name = (*item).first;
            
            // print the variable
            std::cout << "this is var name = "<< var_name << "\n";
            
            user_str_material_inps.push_back(var_name);
            
            // validate input: user_str_material_inps match words in the str_material_inps
            if (std::find(str_material_inps.begin(), str_material_inps.end(), var_name) == str_material_inps.end())
            {
                std::cout << "ERROR: invalid input: " << var_name << std::endl;
            }
            
        } // end for item in this yaml input
        
        
        // loop over the words in the material input definition
        for(auto &a_word : user_str_material_inps){
            
            std::cout << a_word << std::endl;
        
            Yaml::Node & material_inps_yaml = root["materials"][mat_id]["material"][a_word];
            
            
            // set the values in the input for this word
            
            
            // set the values
            if(a_word.compare("q1") == 0){
                
                // outer plane

                double q1 = root["materials"][mat_id]["material"]["q1"].As<double>();
                std::cout << "\tq1 = " << q1 << std::endl;
                
                materials[mat_id].q1 = q1;
            
            } // q1
            else if(a_word.compare("q1ex") == 0){
                
                // outer plane

                double q1ex = root["materials"][mat_id]["material"]["q1ex"].As<double>();
                std::cout << "\tq1ex = " << q1ex << std::endl;
                
                materials[mat_id].q1ex = q1ex;
            
            } // q1ex
            else if(a_word.compare("q2") == 0){
                
                // outer plane

                double q2 = root["materials"][mat_id]["material"]["q2"].As<double>();
                std::cout << "\tq2 = " << q2 << std::endl;
                
                materials[mat_id].q2 = q2;
            
            } // q1
            else if(a_word.compare("q2ex") == 0){
                
                // outer plane

                double q2ex = root["materials"][mat_id]["material"]["q2ex"].As<double>();
                std::cout << "\tq2ex = " << q2ex << std::endl;
                
                materials[mat_id].q2ex = q2ex;
            
            } // q1ex
            else if(a_word.compare("id") == 0){
                  
                int m_id = root["materials"][mat_id]["material"]["id"].As<int>();
                std::cout << "\tid = " << m_id << std::endl;
            
                materials[mat_id].id = m_id;
            
            } // id
            else if(a_word.compare("eos_model") == 0){
                  
                std::string eos = root["materials"][mat_id]["material"]["eos_model"].As<std::string>();
                
                // set the EOS
                if(eos_map.find(eos) != eos_map.end()){
                    materials[mat_id].eos_model = eos_map[eos];
                    std::cout << "\teos_model = " << eos << std::endl;
                    materials[mat_id].eos_model(0.,1.,2.);
                }
                else{
                    std::cout << "ERROR: invalid input: " << eos << std::endl;
                } // end if
            
            } // id
            // exact the eos_global_vars
            else if(a_word.compare("eos_global_vars") == 0){
                
                size_t num_global_vars = material_inps_yaml.Size();
                
                std::cout << "num global eos vars = " << num_global_vars << std::endl;
                
                for(int global_var_id=0; global_var_id<num_global_vars; global_var_id++){
                    
                    double eos_var = root["materials"][mat_id]["material"]["eos_global_vars"][global_var_id].As<double>();
                    
                    eos_global_vars[mat_id].push_back(eos_var);
                    
                    std::cout << "\t var = " << eos_var<< std::endl;
                } // end loop over global vars
                
            } // "eos_global_vars"
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
            }
            
            
        } // end for words in material
        
    } // end loop over materials

} // end of function to parse material information