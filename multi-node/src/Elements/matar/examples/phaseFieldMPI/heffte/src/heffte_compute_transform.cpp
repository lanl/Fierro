/**
 * @class
 * heFFTe kernels for complex-to-complex transforms
 */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "heffte_compute_transform.h"

namespace heffte {

template<typename location_tag, typename index, typename scalar_type>
void compute_transform(typename backend::data_manipulator<location_tag>::stream_type stream,
                       int const batch_size,
                       scalar_type const input[], scalar_type output[], scalar_type workspace[],
                       size_t executor_buffer_offset, size_t size_comm_buffers,
                       std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                       std::array<executor_base*, 3> const &executor,
                       direction dir){

    /*
     * The logic is a bit messy, but the objective is:
     * - call all shaper and executor objects in the correct order
     * - assume that any or all of the shapers can be missing, i.e., null unique_ptr()
     * - do not allocate buffers if not needed
     * - never have more than 2 allocated buffers (input and output)
     */

    scalar_type *executor_workspace = (executor_buffer_offset == 0) ? nullptr : workspace + batch_size * executor_buffer_offset;

    auto apply_fft = [&](int i, scalar_type data[])
        ->void{
            add_trace name("fft-1d");
            if (dir == direction::forward){
                if (executor[i] != nullptr){
                    for(int j=0; j<batch_size; j++)
                        executor[i]->forward(data + j * executor[i]->box_size(), executor_workspace);
                }
            }else{
                if (executor[i] != nullptr){
                    for(int j=0; j<batch_size; j++)
                        executor[i]->backward(data + j * executor[i]->box_size(), executor_workspace);
                }
            }
        };

    int num_active = count_active(shaper);
    int last = get_last_active(shaper);

    if (last < 1){ // no extra buffer case
        add_trace name("reshape/copy");
        // move input -> output and apply all ffts
        // use either zeroth shaper or simple copy (or nothing in case of in-place transform)
        if (last == 0){
            shaper[0]->apply(batch_size, input, output, workspace);
        }else if (input != output){
            int valid_executor = (executor[0] != nullptr) ? 0 : ((executor[1] != nullptr) ? 1 : 2);
            backend::data_manipulator<location_tag>::copy_n(stream, input, batch_size * executor[valid_executor]->box_size(), output);
        }
        for(int i=0; i<3; i++)
            apply_fft(i, output);

        return;
    }

    // with only one reshape, the temp buffer would be used only if not doing in-place
    scalar_type *temp_buffer = workspace + batch_size * size_comm_buffers;
    if (num_active == 1){ // one active and not shaper 0
        scalar_type *effective_input = output;
        if (input != output){
            add_trace name("copy");
            if (executor[0] != nullptr)
                backend::data_manipulator<location_tag>::copy_n(stream, input, batch_size * executor[0]->box_size(), temp_buffer);
            effective_input = temp_buffer;
        }
        for(int i=0; i<last; i++)
            apply_fft(i, effective_input);
        { add_trace name("reshape");
        shaper[last]->apply(batch_size, effective_input, output, workspace);
        }
        for(int i=last; i<3; i++)
            apply_fft(i, output);

        return;
    }

    // with two or more reshapes, the first reshape must move to the temp_buffer and the last must move to output
    int active_shaper = 0;
    if (shaper[0] or input != output){
        if (shaper[0]){
            add_trace name("reshape");
            shaper[0]->apply(batch_size, input, temp_buffer, workspace);
        }else{
            add_trace name("copy");
            backend::data_manipulator<location_tag>::copy_n(stream, input, batch_size * executor[0]->box_size(), temp_buffer);
        }
        active_shaper = 1;
    }else{
        // in place transform and shaper[0] is not active
        while(not shaper[active_shaper]){
            // note, at least one shaper must be active, otherwise last will catch it
            apply_fft(active_shaper++, output);
        }
        { add_trace name("reshape");
        shaper[active_shaper]->apply(batch_size, output, temp_buffer, workspace);
        }
        active_shaper += 1;
    }
    apply_fft(active_shaper - 1, temp_buffer); // one reshape was applied above

    for(int i=active_shaper; i<last; i++){
        if (shaper[i]){
            add_trace name("reshape");
            shaper[i]->apply(batch_size, temp_buffer, temp_buffer, workspace);
        }
        apply_fft(i, temp_buffer);
    }
    { add_trace name("reshape");
    shaper[last]->apply(batch_size, temp_buffer, output, workspace);
    }

    for(int i=last; i<3; i++)
        apply_fft(i, output);
}

template<typename location_tag, typename index, typename scalar_type>
void compute_transform(typename backend::data_manipulator<location_tag>::stream_type,
                       int const batch_size,
                       scalar_type const input[], std::complex<scalar_type> output[],
                       std::complex<scalar_type> workspace[],
                       size_t executor_buffer_offset, size_t size_comm_buffers,
                       std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                       std::array<executor_base*, 3> const &executor, direction){
    /*
     * Follows logic similar to the complex-to-complex case but the first shaper and executor will be applied to real data.
     */
    int last = get_last_active(shaper);
    std::complex<scalar_type> *executor_workspace = (executor_buffer_offset == 0) ?
                                                    nullptr : workspace + batch_size * executor_buffer_offset;

    scalar_type* reshaped_input = reinterpret_cast<scalar_type*>(workspace);
    scalar_type const *effective_input = input; // either input or the result of reshape operation 0
    if (shaper[0]){
        add_trace name("reshape");
        shaper[0]->apply(batch_size, input, reshaped_input,
                         reinterpret_cast<scalar_type*>(workspace + batch_size * get_max_box_size(executor)));
        effective_input = reshaped_input;
    }

    if (last < 1){ // no reshapes after 0
        add_trace name("fft-1d x3");
        for(int j=0; j<batch_size; j++){
            if (executor[0] != nullptr) executor[0]->forward(effective_input + j * executor[0]->box_size(),
                                                             output + j * executor[0]->complex_size(), executor_workspace);
            if (executor[1] != nullptr) executor[1]->forward(output + j * executor[0]->box_size(), executor_workspace);
            if (executor[2] != nullptr) executor[2]->forward(output + j * executor[0]->box_size(), executor_workspace);
        }
        return;
    }

    // if there is messier combination of transforms, then we need internal buffers
    std::complex<scalar_type> *temp_buffer = workspace + batch_size * size_comm_buffers;
    { add_trace name("fft-1d");
        if (executor[0] != nullptr){
            for(int j=0; j<batch_size; j++)
                executor[0]->forward(effective_input + j * executor[0]->box_size(),
                                     temp_buffer + j * executor[0]->complex_size(), executor_workspace);
        }
    }

    for(int i=1; i<last; i++){
        if (shaper[i]){
            add_trace name("reshape");
            shaper[i]->apply(batch_size, temp_buffer, temp_buffer, workspace);
        }
        add_trace name("fft-1d");
        if (executor[i] != nullptr){
            for(int j=0; j<batch_size; j++)
                executor[i]->forward(temp_buffer + j * executor[i]->box_size(), executor_workspace);
        }
    }
    { add_trace name("reshape");
        shaper[last]->apply(batch_size, temp_buffer, output, workspace);
    }

    for(int i=last; i<3; i++){
        add_trace name("fft-1d");
        if (executor[i] != nullptr){
            for(int j=0; j<batch_size; j++)
                executor[i]->forward(output + j * executor[i]->box_size(), executor_workspace);
        }
    }
}
template<typename location_tag, typename index, typename scalar_type>
void compute_transform(typename backend::data_manipulator<location_tag>::stream_type stream,
                       int const batch_size,
                       std::complex<scalar_type> const input[], scalar_type output[],
                       std::complex<scalar_type> workspace[],
                       size_t executor_buffer_offset, size_t size_comm_buffers,
                       std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                       std::array<executor_base*, 3> const &executor, direction){
    /*
     * Follows logic similar to the complex-to-complex case but the last shaper and executor will be applied to real data.
     */
    std::complex<scalar_type> *temp_buffer = workspace + batch_size * size_comm_buffers;
    std::complex<scalar_type> *executor_workspace = (executor_buffer_offset == 0) ?
                                                     nullptr : workspace + batch_size * executor_buffer_offset;

    if (shaper[0]){
        add_trace name("reshape");
        shaper[0]->apply(batch_size, input, temp_buffer, workspace);
    }else{
        add_trace name("copy");
        int valid_executor = (executor[0] != nullptr) ? 0 : ((executor[1] != nullptr) ? 1 : 2);
        backend::data_manipulator<location_tag>::copy_n(stream, input,
                                                        batch_size * executor[valid_executor]->box_size(), temp_buffer);
    }

    for(int i=0; i<2; i++){ // apply the two complex-to-complex ffts
        { add_trace name("fft-1d x3");
            if (executor[i] != nullptr){
                for(int j=0; j<batch_size; j++)
                    executor[i]->backward(temp_buffer + j * executor[i]->box_size(), executor_workspace);
            }
        }
        if (shaper[i+1]){
            add_trace name("reshape");
            shaper[i+1]->apply(batch_size, temp_buffer, temp_buffer, workspace);
        }
    }

    // the result of the first two ffts and three reshapes is stored in temp_buffer
    // executor 2 must apply complex to real backward transform
    if (shaper[3]){
        // there is one more reshape left, transform into a real temporary buffer
        scalar_type* real_buffer = reinterpret_cast<scalar_type*>(workspace);
        { add_trace name("fft-1d");
            if (executor[2] != nullptr){
                for(int j=0; j<batch_size; j++)
                    executor[2]->backward(temp_buffer + j * executor[2]->complex_size(),
                                          real_buffer + j * executor[2]->box_size(), executor_workspace);
            }
        }
        add_trace name("reshape");
        shaper[3]->apply(batch_size, real_buffer, output, reinterpret_cast<scalar_type*>(workspace +
                         batch_size * ((executor[2] == nullptr) ? 0 : executor[2]->box_size()) ));
    }else{
        add_trace name("fft-1d");
        if (executor[2] != nullptr){
            for(int j=0; j<batch_size; j++)
                executor[2]->backward(temp_buffer + j * executor[2]->complex_size(),
                                      output + j * executor[2]->box_size(), executor_workspace);
        }
    }
}

#define heffte_instantiate_transform(location_tag, index) \
    template void compute_transform<location_tag, index, std::complex<float>>( \
                        typename backend::data_manipulator<location_tag>::stream_type, int const, \
                        std::complex<float> const input[], std::complex<float> output[], std::complex<float> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, \
                        direction dir); \
    template void compute_transform<location_tag, index, std::complex<double>>( \
                        typename backend::data_manipulator<location_tag>::stream_type, int const, \
                        std::complex<double> const input[], std::complex<double> output[], std::complex<double> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, \
                        direction dir); \
    template void compute_transform<location_tag, index, float>( \
                        typename backend::data_manipulator<location_tag>::stream_type, int const, \
                        float const input[], float output[], float workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, \
                        direction dir); \
    template void compute_transform<location_tag, index, double>( \
                        typename backend::data_manipulator<location_tag>::stream_type, int const, \
                        double const input[], double output[], double workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, \
                        direction dir); \
    template void compute_transform<location_tag, index, float>( \
                        typename backend::data_manipulator<location_tag>::stream_type stream, int const, \
                        float const input[], std::complex<float> output[], std::complex<float> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, direction); \
    template void compute_transform<location_tag, index, double>( \
                        typename backend::data_manipulator<location_tag>::stream_type stream, int const, \
                        double const input[], std::complex<double> output[], std::complex<double> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, direction); \
    template void compute_transform<location_tag, index, float>( \
                        typename backend::data_manipulator<location_tag>::stream_type stream, int const, \
                        std::complex<float> const input[], float output[], std::complex<float> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, direction); \
    template void compute_transform<location_tag, index, double>( \
                        typename backend::data_manipulator<location_tag>::stream_type stream, int const, \
                        std::complex<double> const input[], double output[], std::complex<double> workspace[], \
                        size_t executor_buffer_offset, size_t size_comm_buffers, \
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper, \
                        std::array<executor_base*, 3> const &executor, direction); \

heffte_instantiate_transform(tag::cpu, int)
heffte_instantiate_transform(tag::cpu, long long)

#ifdef Heffte_ENABLE_GPU
heffte_instantiate_transform(tag::gpu, int)
heffte_instantiate_transform(tag::gpu, long long)
#endif

}
