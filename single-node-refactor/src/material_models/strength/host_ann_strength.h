/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/

#ifndef HOST_ANN_STRENGTH_H
#define HOST_ANN_STRENGTH_H


/////////////////////////////////////////////////////////////////////////////
// helper functions
/////////////////////////////////////////////////////////////////////////////

// the number of nodes in each layer of the ANN
std::vector <size_t> num_nodes_in_layer = {6, 3, 1}; //{64000, 30000, 8000, 4000, 100, 25, 6};

// array of ANN structs
struct ANNLayer_t{

    DCArrayKokkos <float> outputs;  // dims = [layer]
    DFArrayKokkos <float> weights;  // dims = [layer-1, layer]
    DCArrayKokkos <float> biases;   // dims = [layer]  

}; // end struct

void set_biases_host(const DCArrayKokkos <float> &biases){
    const size_t num_j = biases.size();

    FOR_ALL(j,0,num_j, {
		    biases(j) = 0.0;
	}); // end parallel for

}; // end function

void set_weights_host(const DFArrayKokkos <float> &weights){

    const size_t num_i = weights.dims(0);
    const size_t num_j = weights.dims(1);
    
	FOR_ALL(i,0,num_i,
	        j,0,num_j, {
		    
		    weights(i,j) = 1.0;
	}); // end parallel for

}; // end function

void forward_propagate_layer_host(const DCArrayKokkos <float> &inputs,
                                  const DCArrayKokkos <float> &outputs, 
                                  const DFArrayKokkos <float> &weights,
                                  const DCArrayKokkos <float> &biases){
    


    const size_t num_i = inputs.size();
    const size_t num_j = outputs.size();



    // For a GPU, use the nested parallelism below here
    using team_t = typename Kokkos::TeamPolicy<>::member_type;
    Kokkos::parallel_for ("MatVec", Kokkos::TeamPolicy<> (num_j, Kokkos::AUTO),
                 KOKKOS_LAMBDA (const team_t& team_h) {

        float sum = 0;
        int j = team_h.league_rank();
        Kokkos::parallel_reduce (Kokkos::TeamThreadRange (team_h, num_i),
                        [&] (int i, float& lsum) {

            lsum += inputs(i)*weights(i,j) + biases(j);

            
        }, sum); // end parallel reduce

        outputs(j) = 1.0/(1.0 + exp(-sum)); 

        

    }); // end parallel for
    


    return;

}; // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn ANNStrengthModelHost
///
/// \brief Large Artificial Neural Network strength model
///
///  This is the ANN model that returns the stress tensor, the parallelism is
///  inside this model, thus it is launched from the CPU.  The loop is serial
///  over the elements of the mesh, and  nested parallelism is used in the
///  material model (i.e., the ANN).  For performance, there must be more ANN
///  nodes than elements, otherwise, use parallelism over elements and 
///  launch the ANN material model in serial on the device (e.g., GPU).
///
/// \param Element pressure
/// \param Element stress
/// \param Global ID for the element
/// \param Material ID for the element
/// \param Element state variables
/// \param Element Sound speed
/// \param Material density
/// \param Material specific internal energy
/// \param Element velocity gradient
/// \param Element nodes IDs in the element
/// \param Node node coordinates
/// \param Noe velocity 
/// \param Element volume
/// \param Time time step
/// \param Time coefficient in the Runge Kutta time integration step
///
/////////////////////////////////////////////////////////////////////////////
namespace HostANNStrengthModel {

    CMatrix <ANNLayer_t> ANNLayers;
    DCArrayKokkos <float> inputs;
    size_t num_layers;

    static void init_strength_state_vars(
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {

        // =================================================================
        // allocate arrays
        // =================================================================

        // note: the num_nodes_in_layer has the inputs into the ANN, so subtract 1 for the layers
        num_layers = num_nodes_in_layer.size()-1;  

        ANNLayers = CMatrix <ANNLayer_t> (num_layers); // starts at 1 and goes to num_layers

        // input and ouput values to ANN
        inputs = DCArrayKokkos <float> (num_nodes_in_layer[0]);


        // set the strides
        // layer 0 are the inputs to the ANN
        // layer n-1 are the outputs from the ANN
        for (size_t layer=1; layer<=num_layers; layer++){

            // dimensions
            size_t num_i = num_nodes_in_layer[layer-1];
            size_t num_j = num_nodes_in_layer[layer];

            // allocate the weights in this layer
            ANNLayers(layer).weights = DFArrayKokkos <float> (num_i, num_j); 
            ANNLayers(layer).outputs = DCArrayKokkos <float> (num_j);
            ANNLayers(layer).biases = DCArrayKokkos <float> (num_j);

        } // end for


        // =================================================================
        // set weights, biases, and inputs
        // =================================================================
        
        // inputs to ANN
        for (size_t i=0; i<num_nodes_in_layer[0]; i++) {
            inputs.host(i) = 1.0;
        }
        inputs.update_device();  // copy inputs to device

        // weights of the ANN
        for (size_t layer=1; layer<=num_layers; layer++){

            // dimensions
            size_t num_i = num_nodes_in_layer[layer-1];
            size_t num_j = num_nodes_in_layer[layer];

            set_weights_host(ANNLayers(layer).weights);
            set_biases_host(ANNLayers(layer).biases);

        } // end for over layers


    }  // end of init_strength_state_vars


    // this model is launched from the CPU, coding inside is run on GPUS
    static void calc_stress(
        const DCArrayKokkos<double>  &GaussPoints_vel_grad,
        const DCArrayKokkos <double> &node_coords,
        const DCArrayKokkos <double> &node_vel,
        const DCArrayKokkos<size_t>  &nodes_in_elem,
        const DCArrayKokkos<double>  &MaterialPoints_pres,
        const DCArrayKokkos<double>  &MaterialPoints_stress,
        const DCArrayKokkos<double>  &MaterialPoints_sspd,
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const double vol,
        const double dt,
        const double rk_alpha,
        const double time,
        const size_t cycle,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {

        // layer 1, hidden layer 0, uses the inputs as the input values
        forward_propagate_layer_host(inputs,
                                     ANNLayers(1).outputs,
                                     ANNLayers(1).weights,
                                     ANNLayers(1).biases); 

        // layer 2 through n-1, layer n-1 goes to the output
        for (size_t layer=2; layer<=num_layers; layer++){

            // go through this layer, the fcn takes(inputs, outputs, weights, biases)
            forward_propagate_layer_host(ANNLayers(layer-1).outputs, 
                                         ANNLayers(layer).outputs,
                                         ANNLayers(layer).weights,
                                         ANNLayers(1).biases); 
        } // end for

        return;
    } // end of user mat

    
    void destroy(
        const DCArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DCArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const size_t num_material_points,
        const size_t mat_ids)
    {

    } // end destory

} // end namespace





#endif // end Header Guard