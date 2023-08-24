/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_PACK3D_H
#define HEFFTE_PACK3D_H

#include "heffte_common.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftepacking Packing/Unpacking operations
 *
 * MPI communications assume that the data is located in contiguous arrays;
 * however, the blocks that need to be transmitted in an FFT algorithm
 * correspond to sub-boxes of a three dimensional array, which is never contiguous.
 * Thus, packing and unpacking operations are needed to copy the sub-box into contiguous arrays.
 * Furthermore, some backends (e.g., fftw3) work much faster with contiguous FFT
 * transforms, thus it is beneficial to transpose the data between backend calls.
 * Combining unpack and transpose operations reduces data movement.
 */

namespace heffte {

/*!
 * \ingroup hefftepacking
 * \brief Holds the plan for a pack/unpack operation.
 */
template<typename index>
struct pack_plan_3d{
    //! \brief Number of elements in the three directions.
    std::array<index, 3> size;
    //! \brief Stride of the lines.
    index line_stride;
    //! \brief Stride of the planes.
    index plane_stride;
    //! \brief Stride of the lines in the received buffer (transpose packing only).
    index buff_line_stride;
    //! \brief Stride of the planes in the received buffer (transpose packing only).
    index buff_plane_stride;
    //! \brief Maps the i,j,k indexes from input to the output (transpose packing only).
    std::array<int, 3> map;
};

/*!
 * \ingroup hefftepacking
 * \brief Writes a plan to the stream, useful for debugging.
 */
template<typename index>
inline std::ostream & operator << (std::ostream &os, pack_plan_3d<index> const &plan){
    os << "nfast = " << plan.size[0] << "\n";
    os << "nmid  = " << plan.size[1] << "\n";
    os << "nslow = " << plan.size[2] << "\n";
    os << "line_stride = "  << plan.line_stride << "\n";
    os << "plane_stride = " << plan.plane_stride << "\n";
    if (plan.buff_line_stride > 0){
        os << "buff_line_stride = " << plan.buff_line_stride << "\n";
        os << "buff_plane_stride = " << plan.buff_plane_stride << "\n";
        os << "map = (" << plan.map[0] << ", " << plan.map[1] << ", " << plan.map[2] << ")\n";
    }
    os << "\n";
    return os;
}

/*!
 * \ingroup hefftepacking
 * \brief The packer needs to know whether the data will be on the CPU or GPU devices.
 *
 * Specializations of this template will define the type alias \b mode
 * that will be set to either the tag::cpu or tag::gpu.
 */
template<typename backend>
struct packer_backend{};

// typename struct packer_backend<cuda>{ using mode = tag::gpu; } // specialization can differentiate between gpu and cpu backends

/*!
 * \ingroup hefftepacking
 * \brief Defines the direct packer without implementation, use the specializations to get the CPU or GPU implementation.
 */
template<typename mode> struct direct_packer{};

/*!
 * \ingroup hefftepacking
 * \brief Simple packer that copies sub-boxes without transposing the order of the indexes.
 */
template<> struct direct_packer<tag::cpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(void*, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        scalar_type* buffer_iterator = buffer;
        for(index slow = 0; slow < plan.size[2]; slow++){
            for(index mid = 0; mid < plan.size[1]; mid++){
                buffer_iterator = std::copy_n(&data[slow * plan.plane_stride + mid * plan.line_stride], plan.size[0], buffer_iterator);
            }
        }
    }
    //! \brief Execute the planned unpack operation.
    template<typename scalar_type, typename index>
    void unpack(void*, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        for(index slow = 0; slow < plan.size[2]; slow++){
            for(index mid = 0; mid < plan.size[1]; mid++){
                std::copy_n(&buffer[(slow * plan.size[1] + mid) * plan.size[0]],
                            plan.size[0], &data[slow * plan.plane_stride + mid * plan.line_stride]);
            }
        }
    }
};

/*!
 * \ingroup hefftepacking
 * \brief Defines the transpose packer without implementation, use the specializations to get the CPU implementation.
 */
template<typename mode> struct transpose_packer{};

/*!
 * \ingroup hefftepacking
 * \brief Transpose packer that packs sub-boxes without transposing, but unpacks applying a transpose operation.
 */
template<> struct transpose_packer<tag::cpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(void *q, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        direct_packer<tag::cpu>().pack(q, plan, data, buffer); // packing is done the same way as the direct_packer
    }
    /*!
     * \brief Execute the planned unpack operation.
     *
     * Note that this will transpose the data in the process.
     * The transpose is done in blocks to maximize cache reuse.
     */
    template<typename scalar_type, typename index>
    void unpack(void*, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        constexpr index stride = 256 / sizeof(scalar_type);
        if (plan.map[0] == 0 and plan.map[1] == 1){
            for(index i=0; i<plan.size[2]; i++)
                for(index j=0; j<plan.size[1]; j++)
                    for(index k=0; k<plan.size[0]; k++)
                        data[i * plan.plane_stride + j * plan.line_stride + k]
                            = buffer[ i * plan.buff_plane_stride + j * plan.buff_line_stride + k ];

        }else if (plan.map[0] == 0 and plan.map[1] == 2){
            for(index bi=0; bi<plan.size[2]; bi+=stride)
                for(index bj=0; bj<plan.size[1]; bj+=stride)
                    for(index bk=0; bk<plan.size[0]; bk+=stride)
                        for(index i=bi; i<std::min(bi + stride, plan.size[2]); i++)
                            for(index j=bj; j<std::min(bj + stride, plan.size[1]); j++)
                                for(index k=bk; k<std::min(bk + stride, plan.size[0]); k++)
                                    data[i * plan.plane_stride + j * plan.line_stride + k]
                                        = buffer[ j * plan.buff_plane_stride + i * plan.buff_line_stride + k ];

        }else if (plan.map[0] == 1 and plan.map[1] == 0){
            for(index bi=0; bi<plan.size[2]; bi+=stride)
                for(index bj=0; bj<plan.size[1]; bj+=stride)
                    for(index bk=0; bk<plan.size[0]; bk+=stride)
                        for(index i=bi; i<std::min(bi + stride, plan.size[2]); i++)
                            for(index j=bj; j<std::min(bj + stride, plan.size[1]); j++)
                                for(index k=bk; k<std::min(bk + stride, plan.size[0]); k++)
                                    data[i * plan.plane_stride + j * plan.line_stride + k]
                                        = buffer[ i * plan.buff_plane_stride + k * plan.buff_line_stride + j ];

        }else if (plan.map[0] == 1 and plan.map[1] == 2){
            for(index bi=0; bi<plan.size[2]; bi+=stride)
                for(index bj=0; bj<plan.size[1]; bj+=stride)
                    for(index bk=0; bk<plan.size[0]; bk+=stride)
                        for(index i=bi; i<std::min(bi + stride, plan.size[2]); i++)
                            for(index j=bj; j<std::min(bj + stride, plan.size[1]); j++)
                                for(index k=bk; k<std::min(bk + stride, plan.size[0]); k++)
                                    data[i * plan.plane_stride + j * plan.line_stride + k]
                                        = buffer[ k * plan.buff_plane_stride + i * plan.buff_line_stride + j ];

        }else if (plan.map[0] == 2 and plan.map[1] == 0){
            for(index bi=0; bi<plan.size[2]; bi+=stride)
                for(index bj=0; bj<plan.size[1]; bj+=stride)
                    for(index bk=0; bk<plan.size[0]; bk+=stride)
                        for(index i=bi; i<std::min(bi + stride, plan.size[2]); i++)
                            for(index j=bj; j<std::min(bj + stride, plan.size[1]); j++)
                                for(index k=bk; k<std::min(bk + stride, plan.size[0]); k++)
                                    data[i * plan.plane_stride + j * plan.line_stride + k]
                                        = buffer[ j * plan.buff_plane_stride + k * plan.buff_line_stride + i ];

        }else{ // if (plan.map[0] == 2 and plan.map[1] == 1){
            for(index bi=0; bi<plan.size[2]; bi+=stride)
                for(index bj=0; bj<plan.size[1]; bj+=stride)
                    for(index bk=0; bk<plan.size[0]; bk+=stride)
                        for(index i=bi; i<std::min(bi + stride, plan.size[2]); i++)
                            for(index j=bj; j<std::min(bj + stride, plan.size[1]); j++)
                                for(index k=bk; k<std::min(bk + stride, plan.size[0]); k++)
                                    data[i * plan.plane_stride + j * plan.line_stride + k]
                                        = buffer[ k * plan.buff_plane_stride + j * plan.buff_line_stride + i ];

        }

    }
};

/*!
 * \ingroup hefftepacking
 * \brief Specialization for the CPU case.
 */
namespace data_scaling {
    /*!
     * \ingroup hefftepacking
     * \brief Simply multiply the \b num_entries in the \b data by the \b scale_factor.
     */
    template<typename scalar_type, typename index>
    void apply(void*, index num_entries, scalar_type *data, double scale_factor){;
        for(index i=0; i<num_entries; i++) data[i] *= scale_factor;
    }
    /*!
     * \ingroup hefftepacking
     * \brief Complex by real scaling.
     *
     * Depending on the compiler and type of operation, C++ complex numbers can have bad
     * performance compared to float and double operations.
     * Since the scaling factor is always real, scaling can be performed
     * with real arithmetic which is easier to vectorize.
     */
    template<typename precision_type, typename index>
    void apply(void *stream, index num_entries, std::complex<precision_type> *data, double scale_factor){
        apply<precision_type>(stream, 2*num_entries, reinterpret_cast<precision_type*>(data), scale_factor);
    }
    /*!
     * \ingroup hefftepacking
     * \brief Helper method that omits the stream for the CPU case.
     */
    template<typename scalar_type, typename index>
    void apply(index num_entries, scalar_type *data, double scale_factor){
        apply(nullptr, num_entries, data, scale_factor);
    }
}

}

#endif
