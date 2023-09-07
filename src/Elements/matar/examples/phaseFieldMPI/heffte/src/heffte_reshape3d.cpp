/**
 * @class
 * heFFTe kernels for data reshape - communication frameworks
 */
 /*
    -- heFFTe --
        Univ. of Tennessee, Knoxville
        @date
 */

#include "heffte_reshape3d.h"

namespace heffte {

#ifdef Heffte_ENABLE_TRACING

    std::deque<event> event_log;
    std::string log_filename;

#endif

    /*!
     * \internal
     * \ingroup hefftereshape
     * \brief Struct type to convert pre-processor define Heffte_DISABLE_GPU_AWARE_MPI to boolean.
     *
     * \endinternal
     */
#ifdef Heffte_DISABLE_GPU_AWARE_MPI
    struct disable_gpu_aware : std::true_type {};
#else
    struct disable_gpu_aware : std::false_type {};
#endif


/*!
 * \brief Counts how many boxes from the list have a non-empty intersection with the reference box.
 */
template<typename index>
int count_collisions(std::vector<box3d<index>> const &boxes, box3d<index> const reference){
    return std::count_if(boxes.begin(), boxes.end(), [&](box3d<index> const b)->bool{ return not reference.collide(b).empty(); });
}

/*!
 * \brief Returns the ranks that will participate in an all-to-all communication.
 *
 * In a reshape algorithm, consider all ranks and connected them into a graph, where each edge
 * corresponds to a piece of data that must be communicated (send or receive).
 * Then take this rank (defined by the list of send and recv procs) and find the larges connected sub-graph.
 * That corresponds to all the processes that need to participate in an all-to-all communication pattern.
 *
 * \param send_proc is the list of ranks that need data from this rank
 * \param recv_proc is the list of ranks that need to send data to this rank
 * \param input_boxes is the list of all boxes held currently across the comm
 * \param output_boxes is the list of all boxes at the end of the communication
 *
 * \returns a list of ranks that must participate in an all-to-all communication
 *          the list is sorted by the order of the rank
 */
template<typename index>
std::vector<int> a2a_group(std::vector<int> const &send_proc, std::vector<int> const &recv_proc,
                           std::vector<box3d<index>> const &input_boxes, std::vector<box3d<index>> const &output_boxes){
    assert(input_boxes.size() == output_boxes.size());
    std::vector<int> result;
    std::vector<bool> marked(input_boxes.size(), false);

    // start with the processes that are connected to this rank
    for(auto p : send_proc){
        if (marked[p]) continue;
        marked[p] = true;
        result.push_back(p);
    }
    for(auto p : recv_proc){
        if (marked[p]) continue;
        marked[p] = true;
        result.push_back(p);
    }

    // loop over procs in result
    // collide each input_boxes extent with all Nprocs output extents
    // collide each output_boxes extent with all Nprocs input extents
    // add any new collision to result
    // keep iterating until nothing is added to result
    bool adding = true;
    while(adding){
        size_t num_current = result.size();
        for(size_t i=0; i<num_current; i++){
            int iproc = result[i];
            // note the O(n^2) graph search, but should be OK for now
            for(size_t j=0; j<input_boxes.size(); j++){
                if (not marked[j] and not input_boxes[iproc].collide(output_boxes[j]).empty()){
                    result.push_back(j);
                    marked[j] = true;
                }
                if (not marked[j] and not output_boxes[iproc].collide(input_boxes[j]).empty()){
                    result.push_back(j);
                    marked[j] = true;
                }
            }
        }
        adding = (num_current != result.size()); // if nothing got added
    }

    // sort based on the flag
    result.resize(0);
    for(size_t i=0; i<input_boxes.size(); i++)
        if (marked[i]) result.push_back(i);

    return result;
}

template<typename index>
int get_max_size(std::vector<int> const &group, std::vector<box3d<index>> const &input_boxes, std::vector<box3d<index>> const &output_boxes){
    long long max_overlap = 0;
    for(auto isend : group)
        for(auto irecv : group)
            max_overlap = std::max(max_overlap, input_boxes[isend].collide(output_boxes[irecv]).count());
    return static_cast<int>(max_overlap);
}

/*
 * Assumes that all boxes have the same order which may be different from (0, 1, 2).
 * The data-movement will be done from a contiguous buffer into the lines of a box.
 */
template<typename index>
void compute_overlap_map_direct_pack(int me, int nprocs, box3d<index> const source, std::vector<box3d<index>> const &boxes,
                                     std::vector<int> &proc, std::vector<int> &offset, std::vector<int> &sizes,
                                     std::vector<pack_plan_3d<index>> &plans){
    for(int i=0; i<nprocs; i++){
        int iproc = (i + me + 1) % nprocs;
        box3d<index> overlap = source.collide(boxes[iproc]);
        if (not overlap.empty()){
            proc.push_back(iproc);
            offset.push_back((overlap.low[source.order[2]] - source.low[source.order[2]]) * source.osize(0) * source.osize(1)
                              + (overlap.low[source.order[1]] - source.low[source.order[1]]) * source.osize(0)
                              + (overlap.low[source.order[0]] - source.low[source.order[0]]));

            plans.push_back({{overlap.osize(0), overlap.osize(1), overlap.osize(2)}, // fast, mid, and slow sizes
                             source.osize(0), source.osize(1) * source.osize(0), // line and plane strides
                             0, 0, {0, 0, 0}});  // ignore the transpose parameters
            sizes.push_back(overlap.count());
        }
    }
}
template<typename index, bool transpose>
void compute_overlap_map_all2all_pack(std::vector<int> const &group, box3d<index> const target, std::vector<box3d<index>> const &boxes,
                                      std::vector<int> &offset, std::vector<pack_plan_3d<index>> &plans){
    for(auto iproc : group){
        box3d<index> overlap = target.collide(boxes[iproc]);
        if (not overlap.empty()){
            if (transpose){
                // figure out the map between the fast-mid-slow directions of the destination and the fast-mid-slow directions of the source
                std::array<int, 3> map = {-1, -1, -1};
                for(int j=0; j<3; j++)
                    for(int k=0; k<3; k++)
                        if (target.order[k] == boxes[iproc].order[j])
                            map[j] = k;

                plans.push_back({{overlap.osize(0), overlap.osize(1), overlap.osize(2)}, // fast, mid, and slow sizes
                                target.osize(0), target.osize(1) * target.osize(0), // target line and plane strides
                                overlap.size[boxes[iproc].order[0]], // strides for the buffer from the received data
                                overlap.size[boxes[iproc].order[0]] * overlap.size[boxes[iproc].order[1]],
                                map});  // map of the sizes of the overlap to the fast, mid and slow directions of the input
            }else{
                plans.push_back({{overlap.osize(0), overlap.osize(1), overlap.osize(2)}, // fast, mid, and slow sizes
                                target.osize(0), target.osize(1) * target.osize(0), // line and plane strides
                                0, 0, {0, 0, 0}});  // ignore the transpose parameters
            }
            offset.push_back((overlap.low[target.order[2]] - target.low[target.order[2]]) * target.osize(0) * target.osize(1)
                              + (overlap.low[target.order[1]] - target.low[target.order[1]]) * target.osize(0)
                              + (overlap.low[target.order[0]] - target.low[target.order[0]]));
        }else{
            plans.push_back({{0, 0, 0}, 0, 0, 0, 0, {0, 0, 0}});
            offset.push_back(0);
        }
    }
}

template<typename index>
void compute_overlap_map_transpose_pack(int me, int nprocs, box3d<index> const destination, std::vector<box3d<index>> const &boxes,
                                        std::vector<int> &proc, std::vector<int> &offset, std::vector<int> &sizes, std::vector<pack_plan_3d<index>> &plans){
    for(int i=0; i<nprocs; i++){
        int iproc = (i + me + 1) % nprocs;
        box3d<index> overlap = destination.collide(boxes[iproc]);
        if (not overlap.empty()){
            proc.push_back(iproc);
            offset.push_back((overlap.low[destination.order[2]] - destination.low[destination.order[2]]) * destination.osize(0) * destination.osize(1)
                              + (overlap.low[destination.order[1]] - destination.low[destination.order[1]]) * destination.osize(0)
                              + (overlap.low[destination.order[0]] - destination.low[destination.order[0]]));

            // figure out the map between the fast-mid-slow directions of the destination and the fast-mid-slow directions of the source
            std::array<int, 3> map = {-1, -1, -1};
            for(int j=0; j<3; j++)
                for(int k=0; k<3; k++)
                    if (destination.order[k] == boxes[iproc].order[j])
                        map[j] = k;

            plans.push_back({{overlap.osize(0), overlap.osize(1), overlap.osize(2)}, // fast, mid, and slow sizes
                             destination.osize(0), destination.osize(1) * destination.osize(0), // destination line and plane strides
                             overlap.size[boxes[iproc].order[0]], // strides for the buffer from the received data
                             overlap.size[boxes[iproc].order[0]] * overlap.size[boxes[iproc].order[1]],
                             map});  // map of the sizes of the overlap to the fast, mid and slow directions of the input
            sizes.push_back(overlap.count());
        }
    }
}

template
void compute_overlap_map_transpose_pack<int>(int me, int nprocs, box3d<int> const destination, std::vector<box3d<int>> const &boxes,
                                         std::vector<int> &proc, std::vector<int> &offset, std::vector<int> &sizes, std::vector<pack_plan_3d<int>> &plans);
template
void compute_overlap_map_transpose_pack<long long>(int me, int nprocs, box3d<long long> const destination,
                                                   std::vector<box3d<long long>> const &boxes,
                                                   std::vector<int> &proc, std::vector<int> &offset, std::vector<int> &sizes, std::vector<pack_plan_3d<long long>> &plans);

template<typename device_type> void sync_if_not_default(device_type const *device){
    if (device->stream() != nullptr) device->synchronize_device();
}

template<typename location_tag, template<typename device> class packer, typename index>
reshape3d_alltoall<location_tag, packer, index>::reshape3d_alltoall(
                   typename backend::device_instance<location_tag>::stream_type q,
                   int cinput_size, int coutput_size, bool gpu_aware, MPI_Comm ccomm,
                   std::vector<pack_plan_3d<index>> &&cpackplan, std::vector<pack_plan_3d<index>> &&cunpackplan,
                   std::vector<int> &&csend_offset, std::vector<int> &&crecv_offset,
                   int cnum_entries
                                                                  ) :
                   reshape3d_base<index>(cinput_size, coutput_size),
                   backend::device_instance<location_tag>(q),
                   comm(ccomm), me(mpi::comm_rank(comm)), nprocs(mpi::comm_size(comm)),
                   use_gpu_aware( (disable_gpu_aware::value) ? false : gpu_aware ),
                   packplan(std::move(cpackplan)), unpackplan(std::move(cunpackplan)),
                   send_offset(std::move(csend_offset)), recv_offset(std::move(crecv_offset)),
                   num_entries(cnum_entries)
{}

template<typename location_tag, template<typename device> class packer, typename index>
template<typename scalar_type>
void reshape3d_alltoall<location_tag, packer, index>::apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const{

    scalar_type *send_buffer = workspace;
    scalar_type *recv_buffer = workspace + batch_size * num_entries * packplan.size();

    packer<location_tag> packit;

    int offset = 0;

    { add_trace name("packing");
        for(size_t i=0; i<packplan.size(); i++){
            if (packplan[i].size[0] > 0){
                for(int j=0; j<batch_size; j++){
                    packit.pack(this->stream(), packplan[i], source + send_offset[i] + j * this->input_size, send_buffer + offset);
                    offset += num_entries;
                }
            }else{
                offset += batch_size * num_entries;
            }
        }
    }

    this->synchronize_device();

    #ifdef Heffte_ENABLE_GPU
    if (std::is_same<location_tag, tag::gpu>::value and not use_gpu_aware){
        scalar_type *temp = this->template cpu_send_buffer<scalar_type>(batch_size * num_entries * packplan.size());
        gpu::transfer::unload(this->stream(), send_buffer, batch_size * num_entries * packplan.size(), temp);
        send_buffer = temp;
        recv_buffer = this->template cpu_recv_buffer<scalar_type>(batch_size * num_entries * packplan.size());
    }
    #endif

    { add_trace name("all2all");
        MPI_Alltoall(send_buffer, batch_size * num_entries, mpi::type_from<scalar_type>(),
                     recv_buffer, batch_size * num_entries, mpi::type_from<scalar_type>(),
                     comm);
    }

    #ifdef Heffte_ENABLE_GPU
    if (std::is_same<location_tag, tag::gpu>::value and not use_gpu_aware){
        scalar_type* temp = workspace + batch_size * num_entries * packplan.size();
        gpu::transfer::load(this->stream(), recv_buffer, batch_size * num_entries * packplan.size(), temp);
        recv_buffer = temp;
    }
    #endif

    offset = 0;
    { add_trace name("unpacking");
        for(size_t i=0; i<unpackplan.size(); i++){
            if (unpackplan[i].size[0] > 0){
                for(int j=0; j<batch_size; j++){
                    packit.unpack(this->stream(), unpackplan[i],
                                  recv_buffer + offset,
                                  destination + recv_offset[i] + j * this->output_size);
                    offset += num_entries;
                }
            }else{
                offset += batch_size * num_entries;
            }
        }
    }
}

template<typename location_tag, template<typename device> class packer, typename index> std::unique_ptr<reshape3d_alltoall<location_tag, packer, index>>
make_reshape3d_alltoall(typename backend::device_instance<location_tag>::stream_type q,
                        std::vector<box3d<index>> const &input_boxes, std::vector<box3d<index>> const &output_boxes,
                        bool uses_gpu_aware, MPI_Comm const comm){
    int const me = mpi::comm_rank(comm);
    std::vector<int> group = a2a_group({me}, {me}, input_boxes, output_boxes);

    constexpr bool transpose = true;
    constexpr bool non_transpose = false;

    std::vector<pack_plan_3d<index>> packplans, unpackplans;
    std::vector<int> send_offset, recv_offset;

    compute_overlap_map_all2all_pack<index, non_transpose>(group, input_boxes[me], output_boxes, send_offset, packplans);
    if (std::is_same<packer<location_tag>, direct_packer<location_tag>>::value){
        compute_overlap_map_all2all_pack<index, non_transpose>(group, output_boxes[me], input_boxes, recv_offset, unpackplans);
    }else{
        compute_overlap_map_all2all_pack<index, transpose>(group, output_boxes[me], input_boxes, recv_offset, unpackplans);
    }

    return std::unique_ptr<reshape3d_alltoall<location_tag, packer, index>>(new reshape3d_alltoall<location_tag, packer, index>(
        q, input_boxes[me].count(), output_boxes[me].count(),
        uses_gpu_aware, mpi::new_comm_from_group(group, comm),
        std::move(packplans), std::move(unpackplans), std::move(send_offset), std::move(recv_offset),
        get_max_size(group, input_boxes, output_boxes)
                                                       ));

}

template<typename location_tag, template<typename device> class packer, typename index>
reshape3d_alltoallv<location_tag, packer, index>::reshape3d_alltoallv(
                        typename backend::device_instance<location_tag>::stream_type q,
                        int cinput_size, int coutput_size,
                        bool gpu_aware, MPI_Comm new_comm, std::vector<int> const &pgroup,
                        std::vector<int> &&csend_offset, std::vector<int> &&csend_size, std::vector<int> const &send_proc,
                        std::vector<int> &&crecv_offset, std::vector<int> &&crecv_size, std::vector<int> const &recv_proc,
                        std::vector<pack_plan_3d<index>> &&cpackplan, std::vector<pack_plan_3d<index>> &&cunpackplan
                                                                ) :
    reshape3d_base<index>(cinput_size, coutput_size),
    backend::device_instance<location_tag>(q),
    comm(new_comm), me(mpi::comm_rank(comm)), nprocs(mpi::comm_size(comm)),
    use_gpu_aware( (disable_gpu_aware::value) ? false : gpu_aware ),
    send_offset(std::move(csend_offset)), send_size(std::move(csend_size)),
    recv_offset(std::move(crecv_offset)), recv_size(std::move(crecv_size)),
    send_total(std::accumulate(send_size.begin(), send_size.end(), 0)),
    recv_total(std::accumulate(recv_size.begin(), recv_size.end(), 0)),
    packplan(std::move(cpackplan)), unpackplan(std::move(cunpackplan)),
    send(pgroup, send_proc, send_size),
    recv(pgroup, recv_proc, recv_size)
{}

template<typename location_tag, template<typename device> class packer, typename index>
template<typename scalar_type>
void reshape3d_alltoallv<location_tag, packer, index>::apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const{

    scalar_type *send_buffer = workspace;
    scalar_type *recv_buffer = workspace + batch_size * this->input_size;

    packer<location_tag> packit;

    int offset = 0;

    { add_trace name("packing");
    for(auto isend : send.map){
        if (isend >= 0){ // something to send
            for(int j=0; j<batch_size; j++){
                packit.pack(this->stream(), packplan[isend],
                            source + send_offset[isend] + j * this->input_size,
                            send_buffer + offset);
                offset += send_size[isend];
            }
        }
    }
    }

    this->synchronize_device();

    #ifdef Heffte_ENABLE_GPU
    // the device_synchronize() is needed to flush the kernels of the asynchronous packing
    if (std::is_same<location_tag, tag::gpu>::value and not use_gpu_aware){
        scalar_type *temp = this->template cpu_send_buffer<scalar_type>(batch_size * this->input_size);
        gpu::transfer::unload(this->stream(), send_buffer, batch_size * this->input_size, temp);
        send_buffer = temp;
        recv_buffer = this->template cpu_recv_buffer<scalar_type>(batch_size * this->output_size);
    }
    #endif

    { add_trace name("all2allv");
        if (batch_size == 1){
            MPI_Alltoallv(send_buffer, send.counts.data(), send.displacements.data(), mpi::type_from<scalar_type>(),
                          recv_buffer, recv.counts.data(), recv.displacements.data(), mpi::type_from<scalar_type>(),
                          comm);
        }else{
            std::vector<int> send_counts = send.counts;
            std::vector<int> send_displacements = send.displacements;
            std::vector<int> recv_counts = recv.counts;
            std::vector<int> recv_displacements = recv.displacements;
            for(size_t i=0; i<send_counts.size(); i++){
                send_counts[i] *= batch_size;
                send_displacements[i] *= batch_size;
                recv_counts[i] *= batch_size;
                recv_displacements[i] *= batch_size;
            }
            MPI_Alltoallv(send_buffer, send_counts.data(), send_displacements.data(), mpi::type_from<scalar_type>(),
                          recv_buffer, recv_counts.data(), recv_displacements.data(), mpi::type_from<scalar_type>(),
                          comm);
        }
    }

    #ifdef Heffte_ENABLE_GPU
    if (std::is_same<location_tag, tag::gpu>::value and not use_gpu_aware){
        scalar_type* temp = workspace + batch_size * this->input_size;
        gpu::transfer::load(this->stream(), recv_buffer, batch_size * this->output_size, temp);
        recv_buffer = temp;
    }
    #endif

    offset = 0;
    { add_trace name("unpacking");
    for(auto irecv : recv.map){
        if (irecv >= 0){ // something received
            for(int j=0; j<batch_size; j++){
                packit.unpack(this->stream(), unpackplan[irecv],
                              recv_buffer + offset,
                              destination + recv_offset[irecv] + j * this->output_size);
                offset += recv_size[irecv];
            }
        }
    }
    }
}

template<typename location_tag, template<typename device> class packer, typename index>
std::unique_ptr<reshape3d_alltoallv<location_tag, packer, index>>
make_reshape3d_alltoallv(typename backend::device_instance<location_tag>::stream_type q,
                         std::vector<box3d<index>> const &input_boxes,
                         std::vector<box3d<index>> const &output_boxes,
                         bool uses_gpu_aware,
                         MPI_Comm const comm){

    int const me = mpi::comm_rank(comm);
    int const nprocs = mpi::comm_size(comm);

    std::vector<pack_plan_3d<index>> packplan, unpackplan; // will be moved into the class
    std::vector<int> send_offset;
    std::vector<int> send_size;
    std::vector<int> send_proc;
    std::vector<int> recv_offset;
    std::vector<int> recv_size;
    std::vector<int> recv_proc;

    box3d<index> outbox = output_boxes[me];
    box3d<index> inbox  = input_boxes[me];

    // number of ranks that need data from me
    int nsend = count_collisions(output_boxes, inbox);

    if (nsend > 0) // if others need something from me, prepare the corresponding sizes and plans
        compute_overlap_map_direct_pack(me, nprocs, input_boxes[me], output_boxes, send_proc, send_offset, send_size, packplan);

    // number of ranks that I need data from
    int nrecv = count_collisions(input_boxes, outbox);

    if (nrecv > 0){ // if I need something from others, prepare the corresponding sizes and plans
        // the transpose logic is included in the unpack procedure, direct_packer does not transpose
        if (std::is_same<packer<location_tag>, direct_packer<location_tag>>::value){
            compute_overlap_map_direct_pack(me, nprocs, output_boxes[me], input_boxes, recv_proc, recv_offset, recv_size, unpackplan);
        }else{
            compute_overlap_map_transpose_pack(me, nprocs, output_boxes[me], input_boxes, recv_proc, recv_offset, recv_size, unpackplan);
        }
    }

    std::vector<int> pgroup = a2a_group(send_proc, recv_proc, input_boxes, output_boxes);

    MPI_Comm new_comm = mpi::new_comm_from_group(pgroup, comm);

    if (pgroup.empty())
        return std::unique_ptr<reshape3d_alltoallv<location_tag, packer, index>>();
    else
        return std::unique_ptr<reshape3d_alltoallv<location_tag, packer, index>>(new reshape3d_alltoallv<location_tag, packer, index>(
        q, inbox.count(), outbox.count(),
        uses_gpu_aware, new_comm, pgroup,
        std::move(send_offset), std::move(send_size), send_proc,
        std::move(recv_offset), std::move(recv_size), recv_proc,
        std::move(packplan), std::move(unpackplan)
                                                       ));
}

template<typename location_tag, template<typename device> class packer, typename index>
reshape3d_pointtopoint<location_tag, packer, index>::reshape3d_pointtopoint(
                        typename backend::device_instance<location_tag>::stream_type q,
                        int cinput_size, int coutput_size, reshape_algorithm alg, bool gpu_aware, MPI_Comm ccomm,
                        std::vector<int> &&csend_offset, std::vector<int> &&csend_size, std::vector<int> &&csend_proc,
                        std::vector<int> &&crecv_offset, std::vector<int> &&crecv_size, std::vector<int> &&crecv_proc,
                        std::vector<int> &&crecv_loc,
                        std::vector<pack_plan_3d<index>> &&cpackplan, std::vector<pack_plan_3d<index>> &&cunpackplan
                                                                ) :
    reshape3d_base<index>(cinput_size, coutput_size),
    backend::device_instance<location_tag>(q),
    comm(ccomm), me(mpi::comm_rank(comm)), nprocs(mpi::comm_size(comm)),
    self_to_self(not crecv_proc.empty() and (crecv_proc.back() == me)), // check whether we should include "me" in the communication scheme
    algorithm(alg),
    use_gpu_aware( (disable_gpu_aware::value) ? false : gpu_aware ),
    requests(crecv_proc.size() + ((self_to_self) ? -1 : 0)), // remove 1 if using self-to-self
    isends(csend_proc.size() + ((self_to_self) ? -1 : 0)), // remove 1 if using self-to-self
    send_proc(std::move(csend_proc)), send_offset(std::move(csend_offset)), send_size(std::move(csend_size)),
    recv_proc(std::move(crecv_proc)), recv_offset(std::move(crecv_offset)), recv_size(std::move(crecv_size)),
    recv_loc(std::move(crecv_loc)),
    send_total(std::accumulate(send_size.begin(), send_size.end(), 0)),
    recv_total(std::accumulate(recv_size.begin(), recv_size.end(), 0)),
    packplan(std::move(cpackplan)), unpackplan(std::move(cunpackplan))
{
    if (algorithm == reshape_algorithm::p2p_plined){
        max_send_size = this->input_size;
    }else{
        max_send_size = 0;
        for(auto s : send_size) if (max_send_size < s) max_send_size = s;
    }
}

#ifdef Heffte_ENABLE_GPU
template<typename location_tag, template<typename device> class packer, typename index>
template<typename scalar_type>
void reshape3d_pointtopoint<location_tag, packer, index>::no_gpuaware_send_recv(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const{
    scalar_type *send_buffer = workspace;
    scalar_type *recv_buffer = workspace + batch_size * this->input_size;

    scalar_type *cpu_send = this->template cpu_send_buffer<scalar_type>(batch_size * max_send_size);
    scalar_type *cpu_recv = this->template cpu_recv_buffer<scalar_type>(batch_size * this->output_size);

    packer<location_tag> packit;

    // synchronize before starting the receives, because kernels from an other reshape might
    // still be running, using the workspace
    this->synchronize_device();
    // queue the receive messages, using asynchronous receive
    for(size_t i=0; i<requests.size(); i++){
        heffte::add_trace name("irecv " + std::to_string(batch_size * recv_size[i]) + " from " + std::to_string(recv_proc[i]));
        MPI_Irecv(cpu_recv + batch_size * recv_loc[i], batch_size * recv_size[i], mpi::type_from<scalar_type>(),
                  recv_proc[i], 0, comm, &requests[i]);
    }

    // perform the send commands, using blocking send
    if (algorithm == reshape_algorithm::p2p_plined){
        size_t offset = 0;
        for(size_t i=0; i<send_proc.size() + ((self_to_self) ? -1 : 0); i++){
            { heffte::add_trace name("packing");
                for(int j=0; j<batch_size; j++){
                    packit.pack(this->stream(), packplan[i], source + j * this->input_size + send_offset[i],
                                send_buffer + offset + j * send_size[i]);
                }
            }

            gpu::transfer::unload(this->stream(), send_buffer + offset, batch_size * send_size[i], cpu_send + offset);

            { heffte::add_trace name("isend " + std::to_string(send_size[i]) + " for " + std::to_string(send_proc[i]));
            MPI_Isend(cpu_send + offset, batch_size * send_size[i], mpi::type_from<scalar_type>(),
                      send_proc[i], 0, comm, &isends[i]);
            }
            offset += batch_size * send_size[i];
        }
    }else{
        for(size_t i=0; i<send_proc.size() + ((self_to_self) ? -1 : 0); i++){
            { heffte::add_trace name("packing");
                for(int j=0; j<batch_size; j++){
                    packit.pack(this->stream(), packplan[i], source + j * this->input_size + send_offset[i],
                                send_buffer + j * send_size[i]);
                }
            }

            gpu::transfer::unload(this->stream(), send_buffer, batch_size * send_size[i], cpu_send);

            { heffte::add_trace name("send " + std::to_string(batch_size * send_size[i]) + " for " + std::to_string(send_proc[i]));
            MPI_Send(cpu_send, batch_size * send_size[i], mpi::type_from<scalar_type>(), send_proc[i], 0, comm);
            }
        }
    }

    if (self_to_self){ // if using self-to-self, do not invoke an MPI command
        { heffte::add_trace name("self packing");
            for(int j=0; j<batch_size; j++){
                packit.pack(this->stream(), packplan.back(), source + send_offset.back() + j * this->input_size,
                            recv_buffer + batch_size * recv_loc.back() + j * send_size.back());
            }
        }

        { heffte::add_trace name("self unpacking");
            for(int j=0; j<batch_size; j++){
                packit.unpack(this->stream(), unpackplan.back(), recv_buffer + batch_size * recv_loc.back() + j * send_size.back(),
                              destination + j * this->output_size + recv_offset.back());
            }
        }
    }

    for(size_t i=0; i<requests.size(); i++){
        int irecv;
        { heffte::add_trace name("waitany");
        MPI_Waitany(requests.size(), requests.data(), &irecv, MPI_STATUS_IGNORE);
        }
        gpu::transfer::load(this->stream(), cpu_recv + batch_size * recv_loc[irecv], batch_size * recv_size[irecv],
                            recv_buffer + batch_size * recv_loc[irecv]);
        { heffte::add_trace name("unpacking from " + std::to_string(recv_proc[irecv]));
            for(int j=0; j<batch_size; j++){
                packit.unpack(this->stream(), unpackplan[irecv], recv_buffer + batch_size * recv_loc[irecv] + j * recv_size[irecv],
                              destination + j * this->output_size + recv_offset[irecv]);
            }
        }
    }

    if (algorithm == reshape_algorithm::p2p_plined)
        MPI_Waitall(isends.size(), isends.data(), MPI_STATUS_IGNORE);
}
#endif

template<typename location_tag, template<typename device> class packer, typename index>
template<typename scalar_type>
void reshape3d_pointtopoint<location_tag, packer, index>::apply_base(int batch_size, scalar_type const source[], scalar_type destination[], scalar_type workspace[]) const{

    #ifdef Heffte_ENABLE_GPU
    if (std::is_same<location_tag, tag::gpu>::value and not use_gpu_aware){
        no_gpuaware_send_recv(batch_size, source, destination, workspace);
        return;
    }
    #endif

    scalar_type *send_buffer = workspace;
    scalar_type *recv_buffer = workspace + batch_size * this->input_size;

    packer<location_tag> packit;

    // synchronize before starting the receives, because otherwise kernels could be still using
    // the workspace
    this->synchronize_device();
    // queue the receive messages, using asynchronous receive
    for(size_t i=0; i<requests.size(); i++){
        heffte::add_trace name("irecv " + std::to_string(batch_size * recv_size[i]) + " from " + std::to_string(recv_proc[i]));
        MPI_Irecv(recv_buffer + batch_size * recv_loc[i], batch_size * recv_size[i], mpi::type_from<scalar_type>(),
                  recv_proc[i], 0, comm, &requests[i]);
    }

    // perform the send commands, using blocking send
    size_t offset = 0;
    for(size_t i=0; i<send_proc.size() + ((self_to_self) ? -1 : 0); i++){
        { heffte::add_trace name("packing");
            for(int j=0; j<batch_size; j++){
                packit.pack(this->stream(), packplan[i], source + j * this->input_size + send_offset[i],
                            send_buffer + offset + j * send_size[i]);
            }
        }

        this->synchronize_device();

        if (algorithm == reshape_algorithm::p2p_plined){
            heffte::add_trace name("isend " + std::to_string(batch_size * send_size[i]) + " for " + std::to_string(send_proc[i]));
            MPI_Isend(send_buffer + offset, batch_size * send_size[i], mpi::type_from<scalar_type>(),
                      send_proc[i], 0, comm, &isends[i]);
        }else{
            heffte::add_trace name("send " + std::to_string(batch_size* send_size[i]) + " for " + std::to_string(send_proc[i]));
            MPI_Send(send_buffer + offset, batch_size * send_size[i], mpi::type_from<scalar_type>(),
                     send_proc[i], 0, comm);
        }
        offset += batch_size * send_size[i];
    }

    if (self_to_self){ // if using self-to-self, do not invoke an MPI command
        { heffte::add_trace name("self packing");
            for(int j=0; j<batch_size; j++){
                packit.pack(this->stream(), packplan.back(), source + send_offset.back() + j * this->input_size,
                            recv_buffer + batch_size * recv_loc.back() + j * send_size.back());
            }
        }

        { heffte::add_trace name("self unpacking");
            for(int j=0; j<batch_size; j++){
                packit.unpack(this->stream(), unpackplan.back(), recv_buffer + batch_size * recv_loc.back() + j * send_size.back(),
                              destination + j * this->output_size + recv_offset.back());
            }
        }
    }

    for(size_t i=0; i<requests.size(); i++){
        int irecv;
        { heffte::add_trace name("waitany");
        MPI_Waitany(requests.size(), requests.data(), &irecv, MPI_STATUS_IGNORE);
        }

        #ifdef Heffte_ENABLE_ROCM // this synch is not needed under CUDA
        if (std::is_same<location_tag, tag::gpu>::value)
            gpu::synchronize_default_stream();
        #endif

        { heffte::add_trace name("unpacking from " + std::to_string(irecv));
            for(int j=0; j<batch_size; j++){
                packit.unpack(this->stream(), unpackplan[irecv], recv_buffer + batch_size * recv_loc[irecv] + j * recv_size[irecv],
                              destination + j * this->output_size + recv_offset[irecv]);
            }
        }
    }

    if (algorithm == reshape_algorithm::p2p_plined)
        MPI_Waitall(isends.size(), isends.data(), MPI_STATUS_IGNORE);
}

template<typename location_tag, template<typename device> class packer, typename index>
std::unique_ptr<reshape3d_pointtopoint<location_tag, packer, index>>
make_reshape3d_pointtopoint(typename backend::device_instance<location_tag>::stream_type stream,
                         std::vector<box3d<index>> const &input_boxes,
                         std::vector<box3d<index>> const &output_boxes,
                         reshape_algorithm algorithm, bool uses_gpu_aware,
                         MPI_Comm const comm){

    int const me = mpi::comm_rank(comm);
    int const nprocs = mpi::comm_size(comm);

    std::vector<pack_plan_3d<index>> packplan, unpackplan; // will be moved into the class
    std::vector<int> send_offset;
    std::vector<int> send_size;
    std::vector<int> send_proc;
    std::vector<int> recv_offset;
    std::vector<int> recv_size;
    std::vector<int> recv_proc;
    std::vector<int> recv_loc;

    box3d<index> outbox = output_boxes[me];
    box3d<index> inbox  = input_boxes[me];

    // number of ranks that need data from me
    int nsend = count_collisions(output_boxes, inbox);

    if (nsend > 0) // if others need something from me, prepare the corresponding sizes and plans
        compute_overlap_map_direct_pack(me, nprocs, input_boxes[me], output_boxes, send_proc, send_offset, send_size, packplan);

    // number of ranks that I need data from
    int nrecv = count_collisions(input_boxes, outbox);

    if (nrecv > 0){ // if I need something from others, prepare the corresponding sizes and plans
        // the transpose logic is included in the unpack procedure, direct_packer does not transpose
        if (std::is_same<packer<location_tag>, direct_packer<location_tag>>::value){
            compute_overlap_map_direct_pack(me, nprocs, output_boxes[me], input_boxes, recv_proc, recv_offset, recv_size, unpackplan);
        }else{
            compute_overlap_map_transpose_pack(me, nprocs, output_boxes[me], input_boxes, recv_proc, recv_offset, recv_size, unpackplan);
        }
    }

    recv_loc.push_back(0);
    if (not recv_size.empty())
        for(size_t i=0; i<recv_size.size() - 1; i++)
            recv_loc.push_back(recv_loc.back() + recv_size[i]);

    return std::unique_ptr<reshape3d_pointtopoint<location_tag, packer, index>>(new reshape3d_pointtopoint<location_tag, packer, index>(
        stream, inbox.count(), outbox.count(), algorithm, uses_gpu_aware, comm,
        std::move(send_offset), std::move(send_size), std::move(send_proc),
        std::move(recv_offset), std::move(recv_size), std::move(recv_proc),
        std::move(recv_loc),
        std::move(packplan), std::move(unpackplan)
                                                       ));
}

#define heffte_instantiate_reshape3d_algorithm(alg, make_alg, some_backend, index) \
template void alg<some_backend, direct_packer, index>::apply_base<float>(int, float const[], float[], float[]) const; \
template void alg<some_backend, direct_packer, index>::apply_base<double>(int, double const[], double[], double[]) const; \
template void alg<some_backend, direct_packer, index>::apply_base<std::complex<float>>(int, std::complex<float> const[], std::complex<float>[], std::complex<float>[]) const; \
template void alg<some_backend, direct_packer, index>::apply_base<std::complex<double>>(int, std::complex<double> const[], std::complex<double> [], std::complex<double> []) const; \
template void alg<some_backend, transpose_packer, index>::apply_base<float>(int, float const[], float[], float[]) const; \
template void alg<some_backend, transpose_packer, index>::apply_base<double>(int, double const[], double[], double[]) const; \
template void alg<some_backend, transpose_packer, index>::apply_base<std::complex<float>>(int, std::complex<float> const[], std::complex<float>[], std::complex<float>[]) const; \
template void alg<some_backend, transpose_packer, index>::apply_base<std::complex<double>>(int, std::complex<double> const[], std::complex<double> [], std::complex<double> []) const; \
 \
template std::unique_ptr<alg<some_backend, direct_packer, index>> \
make_alg<some_backend, direct_packer, index>(typename backend::device_instance<some_backend>::stream_type, std::vector<box3d<index>> const&, \
                                                           std::vector<box3d<index>> const&, bool, MPI_Comm const); \
template std::unique_ptr<alg<some_backend, transpose_packer, index>> \
make_alg<some_backend, transpose_packer, index>(typename backend::device_instance<some_backend>::stream_type, std::vector<box3d<index>> const&, \
                                                              std::vector<box3d<index>> const&, bool, MPI_Comm const); \


#define heffte_instantiate_reshape3d(some_backend, index) \
heffte_instantiate_reshape3d_algorithm(reshape3d_alltoall, make_reshape3d_alltoall, some_backend, index) \
heffte_instantiate_reshape3d_algorithm(reshape3d_alltoallv, make_reshape3d_alltoallv, some_backend, index) \
 \
template void reshape3d_pointtopoint<some_backend, direct_packer, index>::apply_base<float>(int, float const[], float[], float[]) const; \
template void reshape3d_pointtopoint<some_backend, direct_packer, index>::apply_base<double>(int, double const[], double[], double[]) const; \
template void reshape3d_pointtopoint<some_backend, direct_packer, index>::apply_base<std::complex<float>>(int, std::complex<float> const[], std::complex<float>[], std::complex<float>[]) const; \
template void reshape3d_pointtopoint<some_backend, direct_packer, index>::apply_base<std::complex<double>>(int, std::complex<double> const[], std::complex<double> [], std::complex<double> []) const; \
template void reshape3d_pointtopoint<some_backend, transpose_packer, index>::apply_base<float>(int, float const[], float[], float[]) const; \
template void reshape3d_pointtopoint<some_backend, transpose_packer, index>::apply_base<double>(int, double const[], double[], double[]) const; \
template void reshape3d_pointtopoint<some_backend, transpose_packer, index>::apply_base<std::complex<float>>(int, std::complex<float> const[], std::complex<float>[], std::complex<float>[]) const; \
template void reshape3d_pointtopoint<some_backend, transpose_packer, index>::apply_base<std::complex<double>>(int, std::complex<double> const[], std::complex<double> [], std::complex<double> []) const; \
 \
template std::unique_ptr<reshape3d_pointtopoint<some_backend, direct_packer, index>> \
make_reshape3d_pointtopoint<some_backend, direct_packer, index>(typename backend::device_instance<some_backend>::stream_type, \
                                                              std::vector<box3d<index>> const&, \
                                                              std::vector<box3d<index>> const&, reshape_algorithm, bool, MPI_Comm const); \
template std::unique_ptr<reshape3d_pointtopoint<some_backend, transpose_packer, index>> \
make_reshape3d_pointtopoint<some_backend, transpose_packer, index>(typename backend::device_instance<some_backend>::stream_type, \
                                                                 std::vector<box3d<index>> const&, \
                                                                 std::vector<box3d<index>> const&, reshape_algorithm, bool, MPI_Comm const); \

heffte_instantiate_reshape3d(tag::cpu, int)
heffte_instantiate_reshape3d(tag::cpu, long long)

#ifdef Heffte_ENABLE_GPU
heffte_instantiate_reshape3d(tag::gpu, int)
heffte_instantiate_reshape3d(tag::gpu, long long)
#endif


}
