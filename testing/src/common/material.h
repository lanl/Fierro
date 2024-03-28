

#ifndef FIERRO_MATERIAL_H
#define FIERRO_MATERIAL_H

#include <stdio.h>

#include "matar.h"



namespace model
{

    // strength model types
    enum strength_tag
    {
        none = 0,
        hypo = 1,     // hypoelastic plastic model
        hyper = 2,    // hyperelastic plastic model
    };

} // end of namespace


namespace model_init
{

    // strength model setup
    enum strength_setup_tag
    {
        input = 0,
        user_init = 1,
    };

} // end of namespace




// material model parameters
struct material_t {

    size_t id;
    
    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy
    
    // eos fcn pointer
    void (*eos_model)(double, double, double); // WARNING: a placeholder
    
    // strength fcn pointer
    void (*strength_model)(double, double, double); // WARNING: a placeholder
    
    // hypo or hyper elastic plastic model
    model::strength_tag strength_type;
    
    // setup the strength model via the input file for via a user_setup
    model_init::strength_setup_tag strength_setup=model_init::input;
    
    size_t num_eos_state_vars;
    size_t num_strength_state_vars;
    size_t num_eos_global_vars;
    size_t num_strength_global_vars;
    
    double q1;    // acoustic coefficient in Riemann solver for compresion
    double q1ex;  // acoustic coefficient in Riemann solver for expansion
    double q2;    // linear coefficient in Riemann solver for compression
    double q2ex;  // linear coefficient in Riemann solver for expansion
}; // end material_t



// ----------------------------------
// valid inputs for a material fill
// 
//   materials_text_inp["words"]
//
static std::vector <std::string> str_material_inps
{
    "id",
    "eos_model",
    "strength_model",
    "q1",
    "q2",
    "q1ex",
    "q2ex",
    "eos_global_vars"
};





//WARNING: placeholder
static void ideal_gas(double pres, double den, double sie){
    // do nothing
    std::cout << "hello from ideal_gas! Replace with actual EOS!" << std::endl;
};

//WARNING: placeholder
static void elastic_plastic(double stress, double strain){
    // do nothing
    std::cout << "hello from elastic_plastic! Replace with actual strength model!" << std::endl;
}

// add the eos models here
typedef void (*eos_type)(double, double, double);
static std::map <std::string, eos_type> eos_map
{
    {"ideal_gas", ideal_gas}
};

// add the strength models here
typedef void (*strength_type)(double, double);
static std::map <std::string, strength_type> strength_map
{
    {"elastic_plastic", elastic_plastic}
};

#endif // end Header Guard