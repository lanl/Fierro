/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_RESHAPE3D_H
#define HEFFTE_RESHAPE3D_H

#include "heffte_plan_logic.h"
#include "heffte_backends.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftereshape Reshape operations
 *
 * A reshape operation is one that modifies the distribution of the indexes
 * across an MPI communicator. In a special case, the reshape can correspond
 * to a simple in-node data transpose (i.e., no communication).
 *
 * The reshape operations inherit from a common heffte::reshape3d_base class
 * that defines the apply method for different data-types and the sizes
 * of the input, output, and scratch workspace.
 * Reshape objects are usually wrapped in std::unique_ptr containers,
 * which handles the polymorphic calls at runtime and also indicates
 * the special case of no-reshape when the container is empty.
 */

namespace heffte {

#ifdef Heffte_ENABLE_CUDA
namespace gpu { using namespace cuda; }
#else
#ifdef Heffte_ENABLE_ROCM
namespace gpu { using namespace rocm; }
#endif
#ifdef Heffte_ENABLE_ONEAPI
namespace gpu { using namespace oneapi; }
#endif
#endif

/*!
 * \ingroup hefftereshape
 * \brief Generates an unpack plan where the boxes and the destination do not have the same order.
 *
 * This method does not make any MPI calls, but it uses the set of boxes the define the current distribution of the indexes
 * and computes the overlap and the proc, offset, and sizes vectors for the receive stage of an all-to-all-v communication patterns.
 * In addition, a set of unpack plans is created where the order of the boxes and the destination are different,
 * which will transpose the data. The plan has to be used in conjunction with the transpose packer.
 */
template<typename index>
void compute_overlap_map_transpose_pack(int me, int nprocs, box3d<index> const destination, std::vector<box3d<index>> const &boxes,
                                        std::vector<int> &proc, std::vector<int> &offset, std::vector<int> &sizes, std::vector<pack_plan_3d<index>> &plans);

/*!
 * \ingroup hefftereshape
 * \brief Base reshape interface.
 */
template<typename index>
class reshape3d_base{
public:
    //! \brief Constructor that sets the input and output sizes.
    reshape3d_base(index cinput_size, index coutput_size) : input_size(cinput_size), output_size(coutput_size){};
    //! \brief Default virtual destructor.
    virtual ~reshape3d_base() = default;
    //! \brief Apply the reshape, single precision.
    virtual void apply(int batch_size, float const source[], float destination[], float workspace[]) const = 0;
    //! \brief Apply the reshape, double precision.
    virtual void apply(int batch_size, double const source[], double destination[], double workspace[]) const = 0;
    //! \brief Apply the reshape, single precision complex.
    virtual void apply(int batch_size, std::complex<float> const source[], std::complex<float> destination[], std::complex<float> workspace[]) const = 0;
    //! \brief Apply the reshape, double precision complex.
    virtual void apply(int batch_size, std::complex<double> const source[], std::complex<double> destination[], std::complex<double> workspace[]) const = 0;

    //! \brief Returns the input size.
    index size_intput() const{ return input_size; }
    //! \brief Returns the output size.
    index size_output() const{ return output_size; }
    //! \brief Returns the workspace size.
    virtual size_t size_workspace() const{ return input_size + output_size; }

protected:
    //! \brief Stores the size of the input.
    index const input_size;
    //! \brief Stores the size of the output.
    index const output_size;

    // buffers to be used in the no-gpu-aware algorithm for the temporary cpu storage
    // the no-gpu-aware version alleviate the latency when working with small FFTs
    // hence the cpu buffers will be small and will not cause issues
    // note that the main API accepts a GPU buffer for scratch work and cannot be used here
    //! \brief Allocates and returns a CPU buffer when GPU-Aware communication has been disabled.
    template<typename scalar_type> scalar_type* cpu_send_buffer(size_t num_entries) const{
        size_t float_entries = num_entries * sizeof(scalar_type) / sizeof(float);
        send_unaware.resize(float_entries);
        return reinterpret_cast<scalar_type*>(send_unaware.data());
    }
    //! \brief Allocates and returns a CPU buffer when GPU-Aware communication has been disabled.
    template<typename scalar_type> scalar_type* cpu_recv_buffer(size_t num_entries) const{
        size_t float_entries = num_entries * sizeof(scalar_type) / sizeof(float);
        recv_unaware.resize(float_entries);
        return reinterpret_cast<scalar_type*>(recv_unaware.data());
    }
    //! \brief Temp buffers for the gpu-unaware algorithms.
    mutable std::vector<float> send_unaware;
    //! \brief Temp buffers for the gpu-unaware algorithms.
    mutable std::vector<float> recv_unaware;
};

/*!
 * \ingroup hefftereshape
 * \brief Returns the maximum workspace size used by the shapers.
 */
template<typename index>
inline size_t get_workspace_size(std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shapers){
    size_t max_size = 0;
    for(auto const &s : shapers) if (s) max_size = std::max(max_size, s->size_workspace());
    return max_size;
}

/*!
 * \ingroup hefftereshape
 * \brief Reshape algorithm based on the MPI_Alltoall() method.
 *
 * The communication plan for the reshape requires complex initialization,
 * which is put outside of the class into a factory method.
 *
 * \tparam location_tag is tag::cpu or tag::gpu indicating the location of the input/output buffers
 * \tparam packer the packer algorithms to use in arranging the sub-boxes into the global send/recv buffer,
 *         will work with either heffte::direct_packer or heffte::transpose_packer
 */
template<typename location_tag, template<typename device> class packer, typename index>
class reshape3d_alltoall : public reshape3d_base<index>, public backend::device_instance<location_tag>{
public:
    //! \brief Destructor, frees the comm generated by the constructor.
    ~reshape3d_alltoall(){ mpi::comm_free(comm); }
    //! \brief Factory method, use to construct instances of the class.
    template<typename b, template<typename d> class p, typename i> friend std::unique_ptr<reshape3d_alltoall<b, p, i>>
    make_reshape3d_alltoall(typename backend::device_instance<b>::stream_type, std::vector<box3d<i>> const&, std::vector<box3d<i>> const&, bool, MPI_Comm const);

    //! \brief Apply the reshape operations, single precision overload.
    void apply(int batch_size, float const source[], float destination[], float workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision overload.
    void apply(int batch_size, double const source[], double destination[], double workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, single precision complex overload.
    void apply(int batch_size, std::complex<float> const source[], std::complex<float> destination[], std::complex<float> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision complex overload.
    void apply(int batch_size, std::complex<double> const source[], std::complex<double> destination[], std::complex<double> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }

    //! \brief Templated reshape3d_alltoallv::apply() algorithm for all scalar types.
    template<typename scalar_type>
    void apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const;

    //! \brief The size of the workspace must include padding.
    size_t size_workspace() const override { return 2 * num_entries * packplan.size(); }

private:
    /*!
     * \brief Private constructor that accepts a set of arrays that have been pre-computed by the factory.
     */
    reshape3d_alltoall(typename backend::device_instance<location_tag>::stream_type q,
                       int input_size, int output_size, bool gpu_aware, MPI_Comm ccomm,
                       std::vector<pack_plan_3d<index>>&&, std::vector<pack_plan_3d<index>>&&,
                       std::vector<int>&&, std::vector<int>&&, int);

    MPI_Comm const comm;
    int const me, nprocs;
    bool const use_gpu_aware;

    std::vector<pack_plan_3d<index>> packplan, unpackplan;
    std::vector<int> send_offset, recv_offset;
    int const num_entries;
};

/*!
 * \ingroup hefftereshape
 * \brief Factory method that all the necessary work to establish the communication patterns.
 *
 * The purpose of the factory method is to isolate the initialization code and ensure that the internal
 * state of the class is minimal and const-correct, i.e., objects do not hold onto data that will not be
 * used in a reshape apply and the data is labeled const to prevent accidental corruption.
 *
 * \tparam location_tag the location for the input/output buffers for the reshape operation, tag::cpu or tag::gpu
 * \tparam packer is the packer to use to parts of boxes into global send/recv buffer
 *
 * \param q device stream
 * \param input_boxes list of all input boxes across all ranks in the comm
 * \param output_boxes list of all output boxes across all ranks in the comm
 * \param uses_gpu_aware use MPI calls directly from the GPU (GPU backends only)
 * \param comm the communicator associated with all the boxes
 *
 * \returns unique_ptr containing an instance of the heffte::reshape3d_alltoall
 *
 * Note: the input and output boxes associated with this rank are located at position
 * mpi::comm_rank() in the respective lists.
 */
template<typename location_tag, template<typename device> class packer = direct_packer, typename index>
std::unique_ptr<reshape3d_alltoall<location_tag, packer, index>>
make_reshape3d_alltoall(typename backend::device_instance<location_tag>::stream_type q,
                        std::vector<box3d<index>> const &input_boxes, std::vector<box3d<index>> const &output_boxes,
                        bool uses_gpu_aware, MPI_Comm const comm);

/*!
 * \ingroup hefftereshape
 * \brief Reshape algorithm based on the MPI_Alltoallv() method.
 *
 * The communication plan for the reshape requires complex initialization,
 * which is put outside of the class into a factory method.
 * An instance of the class can be created only via the factory method
 * heffte::make_reshape3d_alltoallv()
 * which allows for stronger const correctness and reduces memory footprint.
 *
 * \tparam location_tag is the location of the input/output buffers, tag::cpu or tag::gpu
 * \tparam packer the packer algorithms to use in arranging the sub-boxes into the global send/recv buffer,
 *         will work with either heffte::direct_packer or heffte::transpose_packer
 */
template<typename location_tag, template<typename device> class packer, typename index>
class reshape3d_alltoallv : public reshape3d_base<index>, public backend::device_instance<location_tag>{
public:
    //! \brief Destructor, frees the comm generated by the constructor.
    ~reshape3d_alltoallv(){ mpi::comm_free(comm); }
    //! \brief Factory method, use to construct instances of the class.
    template<typename b, template<typename d> class p, typename i> friend std::unique_ptr<reshape3d_alltoallv<b, p, i>>
    make_reshape3d_alltoallv(typename backend::device_instance<b>::stream_type, std::vector<box3d<i>> const&, std::vector<box3d<i>> const&, bool, MPI_Comm const);

    //! \brief Apply the reshape operations, single precision overload.
    void apply(int batch_size, float const source[], float destination[], float workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision overload.
    void apply(int batch_size, double const source[], double destination[], double workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, single precision complex overload.
    void apply(int batch_size, std::complex<float> const source[], std::complex<float> destination[], std::complex<float> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision complex overload.
    void apply(int batch_size, std::complex<double> const source[], std::complex<double> destination[], std::complex<double> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }

    //! \brief Templated reshape3d_alltoallv::apply() algorithm for all scalar types.
    template<typename scalar_type>
    void apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const;

private:
    /*!
     * \brief Private constructor that accepts a set of arrays that have been pre-computed by the factory.
     */
    reshape3d_alltoallv(typename backend::device_instance<location_tag>::stream_type q,
                        int input_size, int output_size,
                        bool gpu_aware, MPI_Comm new_comm, std::vector<int> const &pgroup,
                        std::vector<int> &&send_offset, std::vector<int> &&send_size, std::vector<int> const &send_proc,
                        std::vector<int> &&recv_offset, std::vector<int> &&recv_size, std::vector<int> const &recv_proc,
                        std::vector<pack_plan_3d<index>> &&packplan, std::vector<pack_plan_3d<index>> &&unpackplan);

    MPI_Comm const comm;
    int const me, nprocs;
    bool const use_gpu_aware;

    std::vector<int> const send_offset;   // extraction loc for each send
    std::vector<int> const send_size;     // size of each send message
    std::vector<int> const recv_offset;   // insertion loc for each recv
    std::vector<int> const recv_size;     // size of each recv message
    int const send_total, recv_total;

    std::vector<pack_plan_3d<index>> const packplan, unpackplan;

    struct iotripple{
        std::vector<int> counts, displacements, map;
        iotripple(std::vector<int> const &pgroup, std::vector<int> const &proc, std::vector<int> const &sizes) :
            counts(pgroup.size(), 0), displacements(pgroup.size(), 0), map(pgroup.size(), -1)
        {
            int offset = 0;
            for(size_t src = 0; src < pgroup.size(); src++){
                for(size_t i=0; i<proc.size(); i++){
                    if (proc[i] != pgroup[src]) continue;
                    counts[src] = sizes[i];
                    displacements[src] = offset;
                    offset += sizes[i];
                    map[src] = i;
                }
            }
        }

    };

    iotripple const send, recv;
};

/*!
 * \ingroup hefftereshape
 * \brief Factory method that all the necessary work to establish the communication patterns.
 *
 * The purpose of the factory method is to isolate the initialization code and ensure that the internal
 * state of the class is minimal and const-correct, i.e., objects do not hold onto data that will not be
 * used in a reshape apply and the data is labeled const to prevent accidental corruption.
 *
 * \tparam location_tag the location of the input/output buffers, tag::cpu or tag::gpu
 * \tparam packer is the packer to use to parts of boxes into global send/recv buffer
 *
 * \param q device stream
 * \param input_boxes list of all input boxes across all ranks in the comm
 * \param output_boxes list of all output boxes across all ranks in the comm
 * \param use_gpu_aware use MPI calls directly from the GPU (GPU backends only)
 * \param comm the communicator associated with all the boxes
 *
 * \returns unique_ptr containing an instance of the heffte::reshape3d_alltoallv
 *
 * Note: the input and output boxes associated with this rank are located at position
 * mpi::comm_rank() in the respective lists.
 */
template<typename location_tag, template<typename device> class packer = direct_packer, typename index>
std::unique_ptr<reshape3d_alltoallv<location_tag, packer, index>>
make_reshape3d_alltoallv(typename backend::device_instance<location_tag>::stream_type q,
                         std::vector<box3d<index>> const &input_boxes,
                         std::vector<box3d<index>> const &output_boxes,
                         bool use_gpu_aware,
                         MPI_Comm const comm);

/*!
 * \ingroup hefftereshape
 * \brief Reshape algorithm based on the MPI_Send() and MPI_Irecv() methods.
 *
 * Similar to heffte::reshape3d_alltoallv, this class handles a point-to-point reshape
 * and the initialization can be done only with the heffte::make_reshape3d_pointtopoint() factory.
 *
 * \tparam location_tag is tag::cpu or tag::gpu, see the alltoall for details
 * \tparam packer the packer algorithms to use in arranging the sub-boxes into the global send/recv buffer
 */
template<typename location_tag, template<typename device> class packer, typename index>
class reshape3d_pointtopoint : public reshape3d_base<index>, public backend::device_instance<location_tag>{
public:
    //! \brief Destructor, frees the comm generated by the constructor.
    ~reshape3d_pointtopoint() = default;
    //! \brief Factory method, use to construct instances of the class.
    template<typename b, template<typename d> class p, typename i> friend std::unique_ptr<reshape3d_pointtopoint<b, p, i>>
    make_reshape3d_pointtopoint(typename backend::device_instance<b>::stream_type, std::vector<box3d<i>> const&, std::vector<box3d<i>> const&, reshape_algorithm, bool, MPI_Comm const);

    //! \brief Apply the reshape operations, single precision overload.
    void apply(int batch_size, float const source[], float destination[], float workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision overload.
    void apply(int batch_size, double const source[], double destination[], double workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, single precision complex overload.
    void apply(int batch_size, std::complex<float> const source[], std::complex<float> destination[], std::complex<float> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision complex overload.
    void apply(int batch_size, std::complex<double> const source[], std::complex<double> destination[], std::complex<double> workspace[]) const override final{
        apply_base(batch_size, source, destination, workspace);
    }

    //! \brief Templated reshape3d_pointtopoint::apply() algorithm for all scalar types.
    template<typename scalar_type>
    void apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const;

    //! \brief Templated reshape3d_pointtopoint::apply() algorithm that does not use GPU-Aware MPI.
    template<typename scalar_type>
    void no_gpuaware_send_recv(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const;

private:
    /*!
     * \brief Private constructor that accepts a set of arrays that have been pre-computed by the factory.
     */
    reshape3d_pointtopoint(typename backend::device_instance<location_tag>::stream_type stream,
                           int input_size, int output_size, reshape_algorithm alg, bool gpu_aware,  MPI_Comm ccomm,
                           std::vector<int> &&send_offset, std::vector<int> &&send_size, std::vector<int> &&send_proc,
                           std::vector<int> &&recv_offset, std::vector<int> &&recv_size, std::vector<int> &&recv_proc,
                           std::vector<int> &&recv_loc,
                           std::vector<pack_plan_3d<index>> &&packplan, std::vector<pack_plan_3d<index>> &&unpackplan);

    MPI_Comm const comm;
    int const me, nprocs;
    bool const self_to_self;
    reshape_algorithm const algorithm;
    bool const use_gpu_aware;
    mutable std::vector<MPI_Request> requests; // recv_proc.size() requests, but remove one if using self_to_self communication
    mutable std::vector<MPI_Request> isends;

    std::vector<int> const send_proc;     // processor to send towards
    std::vector<int> const send_offset;   // extraction loc for each send
    std::vector<int> const send_size;     // size of each send message
    std::vector<int> const recv_proc;     // processor to receive from
    std::vector<int> const recv_offset;   // insertion loc for each recv
    std::vector<int> const recv_size;     // size of each recv message
    std::vector<int> const recv_loc;      // offset in the receive buffer (recv_offset refers to the the destination buffer)
    int const send_total, recv_total;

    std::vector<pack_plan_3d<index>> const packplan, unpackplan;
    int max_send_size;
};

/*!
 * \ingroup hefftereshape
 * \brief Factory method that all the necessary work to establish the communication patterns.
 *
 * The purpose of the factory method is to isolate the initialization code and ensure that the internal
 * state of the class is minimal and const-correct, i.e., objects do not hold onto data that will not be
 * used in a reshape apply and the data is labeled const to prevent accidental corruption.
 *
 * \tparam location_tag the tag for the input/output buffers, tag::cpu or tag::gpu
 * \tparam packer is the packer to use to parts of boxes into global send/recv buffer
 *
 * \param q device stream
 * \param input_boxes list of all input boxes across all ranks in the comm
 * \param output_boxes list of all output boxes across all ranks in the comm
 * \param algorithm must be either reshape_algorithm::p2p or reshape_algorithm::p2p_plined
 * \param use_gpu_aware use MPI calls directly from the GPU (GPU backends only)
 * \param comm the communicator associated with all the boxes
 *
 * \returns unique_ptr containing an instance of the heffte::reshape3d_pointtopoint
 *
 * Note: the input and output boxes associated with this rank are located at position
 * mpi::comm_rank() in the respective lists.
 */
template<typename location_tag, template<typename device> class packer = direct_packer, typename index>
std::unique_ptr<reshape3d_pointtopoint<location_tag, packer, index>>
make_reshape3d_pointtopoint(typename backend::device_instance<location_tag>::stream_type q,
                            std::vector<box3d<index>> const &input_boxes,
                            std::vector<box3d<index>> const &output_boxes,
                            reshape_algorithm algorithm, bool use_gpu_aware,
                            MPI_Comm const comm);

/*!
 * \ingroup hefftereshape
 * \brief Special case of the reshape that does not involve MPI communication but applies a transpose instead.
 *
 * The operations is implemented as a single unpack operation using the transpose_packer with the same location tag.
 */
template<typename location_tag, typename index>
class reshape3d_transpose : public reshape3d_base<index>, public backend::device_instance<location_tag>{
public:
    //! \brief Constructor using the provided unpack plan.
    reshape3d_transpose(typename backend::device_instance<location_tag>::stream_type q,
                        pack_plan_3d<index> const cplan) :
        reshape3d_base<index>(cplan.size[0] * cplan.size[1] * cplan.size[2], cplan.size[0] * cplan.size[1] * cplan.size[2]),
        backend::device_instance<location_tag>(q),
        plan(cplan)
        {}

    //! \brief Apply the reshape operations, single precision overload.
    void apply(int batch_size, float const source[], float destination[], float workspace[]) const override final{
        transpose(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision overload.
    void apply(int batch_size, double const source[], double destination[], double workspace[]) const override final{
        transpose(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, single precision complex overload.
    void apply(int batch_size, std::complex<float> const source[], std::complex<float> destination[], std::complex<float> workspace[]) const override final{
        transpose(batch_size, source, destination, workspace);
    }
    //! \brief Apply the reshape operations, double precision complex overload.
    void apply(int batch_size, std::complex<double> const source[], std::complex<double> destination[], std::complex<double> workspace[]) const override final{
        transpose(batch_size, source, destination, workspace);
    }

private:
    template<typename scalar_type>
    void transpose(int batch_size, scalar_type const *source, scalar_type *destination, scalar_type *workspace) const{
        if (source == destination){ // in-place transpose will need workspace
            backend::data_manipulator<location_tag>::copy_n(this->stream(), source, batch_size * this->input_size, workspace);
            for(int j=0; j<batch_size; j++)
                transpose_packer<location_tag>().unpack(this->stream(), plan, workspace + j * this->input_size,
                                                        destination + j * this->input_size);
        }else{
            for(int j=0; j<batch_size; j++)
                transpose_packer<location_tag>().unpack(this->stream(), plan, source + j * this->input_size,
                                                        destination + j * this->input_size);
        }
    }

    pack_plan_3d<index> const plan;
};

/*!
 * \ingroup hefftereshape
 * \brief Factory method to create a reshape3d instance.
 *
 * Creates a reshape operation from the geometry defined by the input boxes to the geometry defined but the output boxes.
 * The boxes are spread across the given MPI communicator where the boxes associated with the current MPI rank is located
 * at input_boxes[mpi::comm_rank(comm)] and output_boxes[mpi::comm_rank(comm)].
 *
 * - If the input and output are the same, then an empty unique_ptr is created.
 * - If the geometries differ only in the order, then a reshape3d_transpose instance is created.
 * - In all other cases, a reshape3d_alltoallv instance is created using either direct_packer or transpose_packer.
 *
 * Assumes that the order of the input and output geometries are consistent, i.e.,
 * input_boxes[i].order == input_boxes[j].order for all i, j.
 */
template<typename backend_tag, typename index>
std::unique_ptr<reshape3d_base<index>> make_reshape3d(typename backend::device_instance<typename backend::buffer_traits<backend_tag>::location>::stream_type stream,
                                               std::vector<box3d<index>> const &input_boxes,
                                               std::vector<box3d<index>> const &output_boxes,
                                               MPI_Comm const comm,
                                               plan_options const options){
    using location_tag = typename backend::buffer_traits<backend_tag>::location;

    if (match(input_boxes, output_boxes)){
        if (input_boxes[0].ordered_same_as(output_boxes[0])){
            return std::unique_ptr<reshape3d_base<index>>();
        }else{
            int const me = mpi::comm_rank(comm);
            std::vector<int> proc, offset, sizes;
            std::vector<pack_plan_3d<index>> plans;

            compute_overlap_map_transpose_pack(0, 1, output_boxes[me], {input_boxes[me]}, proc, offset, sizes, plans);

            if (not plans.empty()){
                return std::unique_ptr<reshape3d_base<index>>(new reshape3d_transpose<location_tag, index >(stream, plans[0]));
            }else{
                // when the number of indexes is very small, the current box can be empty
                return std::unique_ptr<reshape3d_base<index>>();
            }
        }
    }else{
        if (options.algorithm == reshape_algorithm::alltoallv){
            if (input_boxes[0].ordered_same_as(output_boxes[0])){
                return make_reshape3d_alltoallv<location_tag, direct_packer, index>(stream, input_boxes, output_boxes,
                                                                                    options.use_gpu_aware, comm);
            }else{
                return make_reshape3d_alltoallv<location_tag, transpose_packer, index>(stream, input_boxes, output_boxes,
                                                                                       options.use_gpu_aware, comm);
            }
        }else if (options.algorithm == reshape_algorithm::alltoall){
            if (input_boxes[0].ordered_same_as(output_boxes[0])){
                return make_reshape3d_alltoall<location_tag, direct_packer, index>(stream, input_boxes, output_boxes,
                                                                                   options.use_gpu_aware, comm);
            }else{
                return make_reshape3d_alltoall<location_tag, transpose_packer, index>(stream, input_boxes, output_boxes,
                                                                                      options.use_gpu_aware, comm);
            }
        }else{
            if (input_boxes[0].ordered_same_as(output_boxes[0])){
                return make_reshape3d_pointtopoint<location_tag, direct_packer, index>(stream, input_boxes, output_boxes,
                                                                                       options.algorithm, options.use_gpu_aware, comm);
            }else{
                return make_reshape3d_pointtopoint<location_tag, transpose_packer, index>(stream, input_boxes, output_boxes,
                                                                                          options.algorithm, options.use_gpu_aware, comm);
            }
        }
    }
}

}

#endif
