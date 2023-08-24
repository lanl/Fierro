/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_GEOMETRY_H
#define HEFFTE_GEOMETRY_H

#include <array>
#include <valarray>
#include "heffte_utils.h"

namespace heffte {

/*!
 * \ingroup fft3d
 * \addtogroup fft3dgeometry Box-geometry operations
 *
 * HeFFTe operates with indexes that are distributed in boxes across the mpi ranks.
 * Several methods help manipulate such vectors of boxes,
 * note that in each instance the order of the boxes in a single vector should always match.
 */

/*!
 * \ingroup fft3d
 * \brief A generic container that describes a 3d box of indexes.
 *
 * The box is defined by three low and three high indexes and holds all indexes from low to high
 * including the low and high. For example,
 * \code
 *  box3d box({0, 0, 0}, {0, 1, 2});
 *  // box will hold the 6 indexes
 *  // (0, 0, 0), (0, 0, 1), (0, 0, 2)
 *  // (0, 1, 0), (0, 1, 1), (0, 1, 2)
 * \endcode
 *
 * The box3d also defines a size field that holds the number of indexes in each direction,
 * for the example above size will be {1, 2, 3}.
 *
 * In addition to the low and high indexes, the box has orientation with respect to the i j k indexing,
 * indicating which are the fast, middle and slow growing dimensions.
 * For example, suppose we want to access the entries of a \b box stored in a \b data array
 * \code
 *  // if box.order is {0, 1, 2}
 *  int const plane = box.size[0] * box.size[1];
 *  int const lane  = box.size[0];
 *  for(int i = box.low[0]; i <= box.high[0]; i++)
 *      for(int j = box.low[1]; j <= box.high[1]; j++)
 *          for(int k = box.low[2]; j <= box.high[2]; j++)
 *              std::cout << data[k * plane + j * line + i] << "\n";
 * \endcode
 * However,
 * \code
 *  // if box.order is {2, 0, 1}
 *  int const plane = box.size[box.order[0]] * box.size[box.order[1]];
 *  int const lane  = box.size[box.order[0]];
 *  for(int i = box.low[0]; i <= box.high[0]; i++)
 *      for(int j = box.low[1]; j <= box.high[1]; j++)
 *          for(int k = box.low[2]; j <= box.high[2]; j++)
 *              std::cout << data[j * plane + i * line + k] << "\n";
 * \endcode
 *
 * The order must always hold a permutation of the entries 0, 1, and 2.
 */
template<typename index = int>
struct box3d{
    //! \brief Constructs a box from the low and high indexes, the span in each direction includes the low and high (uses default order).
    box3d(std::array<index, 3> clow, std::array<index, 3> chigh) :
        low(clow), high(chigh), size({high[0] - low[0] + 1, high[1] - low[1] + 1, high[2] - low[2] + 1}), order({0, 1, 2})
    {}
    //! \brief Constructs a box from the low and high indexes, also sets an order.
    box3d(std::array<index, 3> clow, std::array<index, 3> chigh, std::array<int, 3> corder) :
        low(clow), high(chigh), size({high[0] - low[0] + 1, high[1] - low[1] + 1, high[2] - low[2] + 1}), order(corder)
    {}
    //! \brief Constructor for the two dimensional case.
    box3d(std::array<index, 2> clow, std::array<index, 2> chigh) :
        box3d(std::array<index, 3>{clow[0], clow[1], 0}, std::array<index, 3>{chigh[0], chigh[1], 0}){}
    //! \brief Constructor for the two dimensional case, order is either (0, 1) or (1, 0).
    box3d(std::array<index, 2> clow, std::array<index, 2> chigh, std::array<int, 2> corder) :
        box3d(std::array<index, 3>{clow[0], clow[1], 0}, std::array<index, 3>{chigh[0], chigh[1], 0},
              std::array<index, 3>{corder[0], corder[1], 2}){}
    //! \brief Returns true if the box contains no indexes.
    bool empty() const{ return (size[0] <= 0 or size[1] <= 0 or size[2] <= 0); }
    //! \brief Counts all indexes in the box, i.e., the volume.
    long long count() const{ return (empty()) ? 0 :
        (static_cast<long long>(size[0]) * static_cast<long long>(size[1]) * static_cast<long long>(size[2])); }
    //! \brief Creates a box that holds the intersection of this box and the \b other.
    box3d collide(box3d const other) const{
        if (empty() or other.empty()) return box3d({0, 0, 0}, {-1, -1, -1});
        return box3d({std::max(low[0], other.low[0]), std::max(low[1], other.low[1]), std::max(low[2], other.low[2])},
                     {std::min(high[0], other.high[0]), std::min(high[1], other.high[1]), std::min(high[2], other.high[2])}, order);
    }
    //! \brief Returns the box that is reduced in the given dimension according to the real-to-complex symmetry.
    box3d r2c(int dimension) const{
        if (empty()) return box3d({0, 0, 0}, {-1, -1, -1});
        switch(dimension){
            case 0: return box3d(low, {low[0] + size[0] / 2, high[1], high[2]}, order);
            case 1: return box3d(low, {high[0], low[1] + size[1] / 2, high[2]}, order);
            default: // dimension == 2
                return box3d(low, {high[0], high[1], low[2] + size[2] / 2}, order);
        }
    }
    //! \brief Returns \b true if either dimension contains only a single index.
    bool is2d() const{ return (size[0] == 1 or size[1] == 1 or size[2] == 1); }
    //! \brief Compares two boxes, ignoring the order, returns \b true if all sizes and boundaries match.
    bool operator == (box3d const &other) const{
        return not (*this != other);
    }
    //! \brief Compares two boxes, ignoring the order, returns \b true if either of the box boundaries do not match.
    bool operator != (box3d const &other) const{
        for(int i=0; i<3; i++)
            if (low[i] != other.low[i] or high[i] != other.high[i]) return true;
        return false;
    }
    //! \brief Compares the order of two boxes, ignores the dimensions, returns \b true if the order is the same.
    bool ordered_same_as(box3d const &other) const{
        return (order[0] == other.order[0]) and (order[1] == other.order[1]) and (order[2] == other.order[2]);
    }
    //! \brief Returns the effective order of the direction (dir), 0 -> fast, 1 -> mid, 2 -> slow.
    int find_order(int dir) const{ return ((dir == order[0]) ? 0 : ((dir == order[1]) ? 1 : 2)); }
    //! \brief Get the ordered size of the dimension, i.e., size[order[dimension]].
    index osize(int dimension) const{ return size[order[dimension]]; }
    //! \brief The three lowest indexes.
    std::array<index, 3> const low;
    //! \brief The three highest indexes.
    std::array<index, 3> const high;
    //! \brief The number of indexes in each direction.
    std::array<index, 3> const size;
    //! \brief The order of the dimensions in the k * plane_stride + j * line_stride + i indexing.
    std::array<int, 3> const order;
};

/*!
 * \ingroup fft3d
 * \brief Alias for expressive calls to heffte::fft2d and heffte::fft2d_r2c.
 */
template<typename index = int>
using box2d = box3d<index>;

/*!
 * \ingroup fft3dmisc
 * \brief Debugging info, writes out the box to a stream.
 */
template<typename index>
inline std::ostream & operator << (std::ostream &os, box3d<index> const box){
    for(int i=0; i<3; i++)
        os << box.low[i] << "  " << box.high[i] << "  (" << box.size[i] << ")\n";
    os << "(" << box.order[0] << "," << box.order[1] << "," << box.order[2] << ")\n";
    os << "\n";
    return os;
}

/*!
 * \ingroup fft3dbackend
 * \brief Return the number of 1-D ffts contained in the box in the given dimension.
 */
template<typename index>
inline int fft1d_get_howmany(box3d<index> const box, int const dimension){
    if (dimension == box.order[0]) return box.osize(1) * box.osize(2);
    if (dimension == box.order[1]) return box.osize(0);
    return box.osize(0) * box.osize(1);
}
/*!
 * \ingroup fft3dbackend
 * \brief Return the stride of the 1-D ffts contained in the box in the given dimension.
 */
template<typename index>
inline int fft1d_get_stride(box3d<index> const box, int const dimension){
    if (dimension == box.order[0]) return 1;
    if (dimension == box.order[1]) return box.osize(0);
    return box.osize(0) * box.osize(1);
}

/*!
 * \ingroup fft3dgeometry
 * \brief Keeps the local rank and the map from the global rank to the sub-ranks used in the work.
 */
struct rank_remap{
    //! \brief Default constructor.
    rank_remap() : mpi_rank(-1), size_subcomm(0){}
    //! \brief Creates a new remap and initializes only the mpi_rank.
    rank_remap(int rank) : mpi_rank(rank), size_subcomm(0){}
    //! \brief The local rank.
    int const mpi_rank;
    //! \brief Global to local map.
    std::vector<int> map;
    //! \brief Build the map when using fixed number of sub-ranks.
    void set_subranks(size_t all_ranks, int sub_ranks){
        size_subcomm = sub_ranks;
        map = std::vector<int>(all_ranks, -1);
        for(int i=0; i<sub_ranks; i++) map[i] = i;
    }
    //! \brief Build the map when using user provided sub-communicator.
    void set_subranks(MPI_Comm global, MPI_Comm subcomm){
        int subcomm_rank = (subcomm == MPI_COMM_NULL) ? -1 : mpi::comm_rank(subcomm);
        map = std::vector<int>(mpi::comm_size(global));
        MPI_Allgather(&subcomm_rank, 1, MPI_INT, map.data(), 1, MPI_INT, global);
    }
    //! \brief The number of ranks in the sub-comm.
    int size_subcomm;
    //! \brief Indicates whether the remap is empty.
    bool empty() const{ return map.empty(); }
};

/*!
 * \ingroup fft3dgeometry
 * \brief Pair of lists of input-output boxes as used by the heffte::fft3d.
 *
 * The strict contains all inboxes and outboxes for the ranks associated with the geometry of a heffte::fft3d transformation.
 */
template<typename index = int>
struct ioboxes{
    //! \brief Inboxes for all ranks across the comm.
    std::vector<box3d<index>> in;
    //! \brief Outboxes for all ranks across the comm.
    std::vector<box3d<index>> out;
};

/*!
 * \ingroup fft3dgeometry
 * \brief Returns the box that encapsulates all other boxes.
 *
 * Searches through the world.in boxes and computes the highest and lowest of all entries.
 *
 * \param boxes the collection of all input and output boxes.
 */
template<typename index>
inline box3d<index> find_world(std::vector<box3d<index>> const &boxes){
    std::array<index, 3> low  = boxes[0].low;
    std::array<index, 3> high = boxes[0].high;
    for(auto b : boxes){
        for(index i=0; i<3; i++)
            low[i] = std::min(low[i], b.low[i]);
        for(index i=0; i<3; i++)
            high[i] = std::max(high[i], b.high[i]);
    }
    return {low, high};
}

/*!
 * \ingroup fft3dgeometry
 * \brief Compares two vectors of boxes, returns true if all boxes match.
 */
template<typename index>
inline bool match(std::vector<box3d<index>> const &shape0, std::vector<box3d<index>> const &shape1){
    if (shape0.size() != shape1.size()) return false;
    for(size_t i=0; i<shape0.size(); i++)
        if (shape0[i] != shape1[i])
            return false;
    return true;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Returns true if the geometry of the world is as expected.
 *
 * Runs simple checks to ensure that the inboxes will fill the world.
 * \param boxes is the collection of all world boxes
 * \param world the box that incorporates all other boxes
 *
 * The check is not very rigorous at the moment, a true rigorous test
 * will probably be too expensive unless lots of thought is put into it.
 */
template<typename index>
inline bool world_complete(std::vector<box3d<index>> const &boxes, box3d<index> const world){
    long long wsize = 0;
    for(auto b : boxes) wsize += b.count();
    if (wsize < world.count())
        throw std::invalid_argument("The provided input boxes do not fill the world box!");

    for(size_t i=0; i<3; i++)
        if (world.low[i] != 0)
            throw std::invalid_argument("Global box indexing must start from 0!");

    for(size_t i=0; i<boxes.size(); i++)
        for(size_t j=0; j<boxes.size(); j++)
            if (i != j and not boxes[i].collide(boxes[j]).empty())
                throw std::invalid_argument("Input boxes cannot overlap!");

    return true;
}

/*!
 * \brief fft3dmisc
 * \brief Generate all possible factors of a number.
 *
 * Taking an integer \b n, generates a list of all possible way that
 * \b n can be split into a product of two numbers.
 * This is using a brute force method.
 *
 * \param n is the number to factor
 * \returns list of integer factors that multiply to \b n, i.e.,
 * \code
 *  std::vector<std::array<int, 2>> factors = get_factors(n);
 *  for(auto f : factors) assert( f[0] * f[1] == n );
 * \endcode
 *
 * Note: the result is never empty as it always contains (1, n) and (n, 1).
 */
inline std::vector<std::array<int, 2>> get_factors(int const n){
    std::vector<std::array<int, 2>> result;
    for(int i=1; i<=n; i++){
        if (n % i == 0)
            result.push_back({i, n / i});
    }
    if (n == 1) {
        result.push_back({1, 1});
    }
    return result;
}

/*!
 * \ingroup fft3dmisc
 * \brief Get the surface area of a processor grid.
 *
 * For a three dimensional grid with size dims[0] by dims[1] by 1,
 * returns the surface area.
 * Useful for optimizing average communication cost.
 */
inline int get_area(std::array<int, 2> const &dims){ return dims[0] * dims[1] + dims[0] + dims[1]; }


/*!
 * \ingroup fft3dmisc
 * \brief Factorize the MPI ranks into a 2D grid.
 *
 * Considers all possible factorizations of the total number of processors
 * and select the one with the lowest area which heuristically reduces the
 * number of ranks that need to communicate in each of the all-to-all operations.
 *
 * \param num_procs is the total number of processors to factorize
 *
 * \returns the two dimensions of the grid
 */
inline std::array<int, 2> make_procgrid(int const num_procs){
    auto factors = get_factors(num_procs);
    std::array<int, 2> min_array = factors.front();
    int min_area = get_area(min_array);
    for(auto f : factors){
        int farea = get_area(f);
        if (farea < min_area){
            min_array = f;
            min_area = farea;
        }
    }
    return min_array;
}

/*!
 * \ingroup fft3dmisc
 * \brief Factorize the MPI ranks into a 2D grid with specific constraints.
 *
 * The constraints satisfied by the grid will be as follow:
 * - result[i] <= world.size[i], i.e., we don't use more processors than there are indexes
 * - result[direction_1d] == 1, i.e., the grid is one-dimensional in the given \b direction_1d
 * - the product of result[i] will be equal to the product of candidate_grid[0] * candidate_grid[1]
 * - if possible, the \b candidate_grid factorization will be used with the implicit assumption
 *   that it will be best or optimal except when the fist constraint is violated
 */
template<typename index>
inline std::array<int, 3> make_procgrid2d(box3d<index> const world, int direction_1d, std::array<int, 2> const candidate_grid){
    auto make_grid = [&](std::array<int, 2> const &grid)
                     ->std::array<int, 3>{
                         return  (direction_1d == 0) ?  std::array<int, 3>{1, grid[0], grid[1]} :
                                ((direction_1d == 1) ? std::array<int, 3>{grid[0], 1, grid[1]} :
                                                       std::array<int, 3>{grid[0], grid[1], 1});
                     };
    std::array<int, 3> result = make_grid(candidate_grid);

    auto valid = [&](std::array<int, 3> const &grid)
                 ->bool{
                     for(int i=0; i<3; i++) if (grid[i] > world.size[i]) return false;
                     return true;
                 };
    if (valid(result)) return result; // if the first constraint is OK

    // otherwise seek a new factorization
    auto factors = get_factors(candidate_grid[0] * candidate_grid[1]);

    int min_area = get_area(factors.front());
    // a bit misleading, but here we initialize min_area with the largest possible area
    for(auto const &g : factors) min_area = std::max(min_area, get_area(g));

    for(auto const &g : factors){
        if (valid(make_grid(g)) and get_area(g) <= min_area){
            result = make_grid(g);
            min_area = get_area(g);
        }
    }

    if (not valid(result)) throw std::runtime_error("Cannot split the given number of indexes into the given set of mpi-ranks. Most liklely, the number of indexes is too small compared to the number of mpi-ranks.");

    return result;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Splits the world box into a set of boxes that will be assigned to a process in the process grid.
 *
 * \param world is a box describing all indexes of consideration,
 *              there is no assumption on the lower or upper bound of the world
 * \param proc_grid describes the number of boxes in each dimension
 *
 * \returns a list of non-overlapping boxes with union that fills the world where each box contains
 *          approximately the same number of indexes
 */
template<typename index>
inline std::vector<box3d<index>> split_world(box3d<index> const world, std::array<int, 3> const proc_grid, rank_remap const &remap = rank_remap()){

    auto fast = [=](index i)->index{ return world.low[0] + i * (world.size[0] / proc_grid[0]) + std::min(i, (world.size[0] % proc_grid[0])); };
    auto mid  = [=](index i)->index{ return world.low[1] + i * (world.size[1] / proc_grid[1]) + std::min(i, (world.size[1] % proc_grid[1])); };
    auto slow = [=](index i)->index{ return world.low[2] + i * (world.size[2] / proc_grid[2]) + std::min(i, (world.size[2] % proc_grid[2])); };

    std::vector<box3d<index>> result;
    result.reserve(proc_grid[0] * proc_grid[1] * proc_grid[2]);
    for(index k = 0; k < proc_grid[2]; k++){
        for(index j = 0; j < proc_grid[1]; j++){
            for(index i = 0; i < proc_grid[0]; i++){
                result.push_back(box3d<index>({fast(i), mid(j), slow(k)}, {fast(i+1)-1, mid(j+1)-1, slow(k+1)-1}, world.order));
            }
        }
    }
    if (remap.empty()){
        return result;
    }else{
        std::vector<box3d<index>> remapped_result;
        for(size_t i=0; i<remap.map.size(); i++)
            remapped_result.push_back( (remap.map[i] == -1) ?
                                       box3d<index>(std::array<index, 3>{0, 0, 0}, std::array<index, 3>{-1, -1, -1}, world.order) :
                                       result[remap.map[i]]
                            );
        return remapped_result;
    }
}

/*!
 * \ingroup fft3dgeometry
 * \brief Returns true if the shape forms pencils in the given direction.
 */
template<typename index>
inline bool is_pencils(box3d<index> const world, std::vector<box3d<index>> const &shape, int direction){
    for(auto s : shape)
        if (s.size[direction] != world.size[direction])
            return false;
    return true;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Returns true if the shape forms slabs in the given directions.
 */
template<typename index>
inline bool is_slab(box3d<index> const world, std::vector<box3d<index>> const &shape, int direction1, int direction2){
    for(auto s : shape)
        if (s.size[direction1] != world.size[direction1] or s.size[direction2] != world.size[direction2])
            return false;
    return true;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Returns the same shape, but sets a different order for each box.
 */
template<typename index>
inline std::vector<box3d<index>> reorder(std::vector<box3d<index>> const &shape, std::array<int, 3> order){
    std::vector<box3d<index>> result;
    result.reserve(shape.size());
    for(auto const &b : shape)
        result.push_back(box3d<index>(b.low, b.high, order));
    return result;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Shuffle the new boxes to maximize the overlap with the old boxes
 *
 * A reshape operation from an old to a new configuration will require as much
 * MPI communication as the lack of overlap between the two box sets,
 * hence using a heuristic algorithms, in an attempt to find a reordering
 * of the boxes to increase the overlap with the old and new.
 * Also assigns the given order to the result.
 *
 * \param new_boxes is a set of new boxes to be reordered
 * \param old_boxes is the current box configuration
 * \param order is the new order to be assigned to the result boxes
 */
template<typename index>
inline std::vector<box3d<index>> maximize_overlap(std::vector<box3d<index>> const &new_boxes,
                                                  std::vector<box3d<index>> const &old_boxes,
                                                  std::array<int, 3> const order,
                                                  rank_remap const &remap){
    std::vector<box3d<index>> result;
    result.reserve(new_boxes.size());
    std::vector<bool> taken(new_boxes.size(), false);
    if (not remap.empty()){
        for(size_t i=0; i<remap.map.size(); i++)
            if (remap.map[i] == -1) taken[i] = true;
    }

    for(size_t i=0; i<new_boxes.size(); i++){
        if (not remap.empty() and remap.map[i] == -1){
            result.push_back(box3d<index>(new_boxes[i].low, new_boxes[i].high, order));
            continue;
        }
        // for each box in the result, find the box among the new_boxes
        // that has not been taken and has the largest overlap with the corresponding old box
        int max_overlap = -1;
        size_t max_index = new_boxes.size();
        for(size_t j=0; j<new_boxes.size(); j++){
            int overlap = old_boxes[i].collide(new_boxes[j]).count();
            if (not taken[j] and overlap > max_overlap){
                max_overlap = overlap;
                max_index = j;
            }
        }
        assert( max_index < new_boxes.size() ); // if we found a box
        taken[max_index] = true;
        result.push_back(box3d<index>(new_boxes[max_index].low, new_boxes[max_index].high, order));
    }

    return result;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Counts the number of point-to-point connections between the old and new box geometries.
 *
 * Given a grid factorization of (a, b) the split into pencils can be done as either (1, a, b) or (1, b, a),
 * the heffte::make_pencils() method computes both factorizations and selects the one that will lead to fewer
 * point-to-point communications. This allows the pencil rotations to be done within the rows/columns of the 2D grid
 * as opposed to communicating between all available ranks.
 */
template<typename index>
inline long long count_connections(std::vector<box3d<index>> const &new_boxes, std::vector<box3d<index>> const &old_boxes){
    long long count = 0;

    for(auto &nbox : new_boxes)
        for(auto &obox : old_boxes)
            if (not nbox.collide(obox).empty())
                count++;

    return count;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Breaks the world into a grid of pencils and orders the pencils to the ranks that will minimize communication
 *
 * A pencil is a box with one dimension that matches the entire world,
 * a pencils grid is a two-dimensional grid of boxes that captures a three dimensional world box.
 *
 * This calls heffte::split_world() and then rearranges the list so that performing a reshape operation
 * from the source to the resulting list will minimize communication.
 * \param world is a box describing all indexes of consideration,
 *              it is assumed that the world is the union of the \b source boxes
 * \param proc_grid gives the number of boxes to use for the two-by-two grid
 * \param dimension is 0, 1, or 2, indicating the direction of orientation of the pencils,
 *                  e.g., dimension 1 means that pencil.size[1] == world.size[1]
 *                  for each pencil in the output list
 * \param source is the current distribution of boxes across MPI ranks,
 *               and will be used as a reference when remapping boxes to ranks
 * \param order is the box index order (fast, mid, slow) that will be assigned to the result.
 *
 * \returns a sorted list of non-overlapping pencils which union is the \b world box
 */
template<typename index>
inline std::vector<box3d<index>> make_pencils(box3d<index> const world,
                                              std::array<int, 2> const proc_grid,
                                              int const dimension,
                                              std::vector<box3d<index>> const &source,
                                              std::array<int, 3> const order,
                                              rank_remap const &remap = rank_remap()
                                             ){
    // trivial case, the grid is already in a pencil format
    if (is_pencils(world, source, dimension))
        return reorder(source, order);

    // create a list of boxes ordered in column major format (following the proc_grid box)
    // using two boxes corresponding to the two dimensions of proc-grid in both variants
    std::vector<box3d<index>> pencilsA =
        maximize_overlap(
            split_world(world, make_procgrid2d(world, dimension, proc_grid), remap), source, order, remap);

    std::vector<box3d<index>> pencilsB =
        maximize_overlap(
            split_world(world, make_procgrid2d(world, dimension, std::array<int, 2>{proc_grid[1], proc_grid[0]}), remap), source, order, remap);

    return (count_connections(pencilsB, source) < count_connections(pencilsA, source)) ? pencilsB : pencilsA;
}

/*!
 * \ingroup fft3dgeometry
 * \brief Breaks the world into a set of slabs that span the given dimensions
 *
 * The method is near identical to make_pencils, but the slabs span two dimensions.
 */
template<typename index>
inline std::vector<box3d<index>> make_slabs(box3d<index> const world, int num_slabs,
                                            int const dimension1, int const dimension2,
                                            std::vector<box3d<index>> const &source,
                                            std::array<int, 3> const order,
                                            rank_remap const &remap
                                           ){
    assert( dimension1 != dimension2 );
    std::vector<box3d<index>> slabs;
    if (dimension1 == 0){
        if (dimension2 == 1){
            slabs = split_world(world, {1, 1, num_slabs}, remap);
        }else{
            slabs = split_world(world, {1, num_slabs, 1}, remap);
        }
    }else if (dimension1 == 1){
        if (dimension2 == 0){
            slabs = split_world(world, {1, 1, num_slabs}, remap);
        }else{
            slabs = split_world(world, {num_slabs, 1, 1}, remap);
        }
    }else{ // third dimension
        if (dimension2 == 0){
            slabs = split_world(world, {1, num_slabs, 1}, remap);
        }else{
            slabs = split_world(world, {num_slabs, 1, 1}, remap);
        }
    }

    return maximize_overlap(slabs, source, order, remap);
}

/*!
 * \ingroup fft3dgeometry
 * \brief Creates a grid of mpi-ranks that will minimize the area of each of the boxes.
 *
 * Given the world box of indexes, generate the dimensions of a 3d grid of mpi-ranks,
 * where the grid is chosen to minimize the total surface area of each of the boxes.
 *
 * \param world is the box of all indexes starting from 0.
 * \param num_procs is the total number of mpi-ranks to use for the process grid.
 *
 * \returns the dimensions of the 3d grid that will minimize the size of each box.
 */
template<typename index>
inline std::array<int, 3> proc_setup_min_surface(box3d<index> const world, int num_procs){
    assert(world.count() > 0); // make sure the world is not empty

    // using valarrays that work much like vectors, but can perform basic
    // point-wise operations such as addition, multiply, and division
    std::valarray<index> all_indexes = {world.size[0], world.size[1], world.size[2]};
    // set initial guess, probably the worst grid but a valid one
    std::valarray<index> best_grid = {1, 1, num_procs};

    // internal helper method to compute the surface
    auto surface = [&](std::valarray<index> const &proc_grid)->
        index{
            auto box_size = all_indexes / proc_grid;
            return ( box_size * box_size.cshift(1) ).sum();
        };

    index best_surface = surface({1, 1, num_procs});

    for(int i=1; i<=num_procs; i++){
        if (num_procs % i == 0){
            int const remainder = num_procs / i;
            for(int j=1; j<=remainder; j++){
                if (remainder % j == 0){
                    std::valarray<index> candidate_grid = {i, j, remainder / j};
                    index const candidate_surface = surface(candidate_grid);
                    if (candidate_surface < best_surface){
                        best_surface = candidate_surface;
                        best_grid    = candidate_grid;
                    }
                }
            }
        }
    }

    assert(best_grid[0] * best_grid[1] * best_grid[2] == num_procs);

    return {static_cast<int>(best_grid[0]), static_cast<int>(best_grid[1]), static_cast<int>(best_grid[2])};
}

namespace mpi {
/*!
 * \ingroup hefftempi
 * \brief Gather all boxes across all ranks in the comm.
 *
 * Constructs an \b ioboxes struct with input and output boxes collected from all ranks.
 * \param my_inbox is the input box on this rank
 * \param my_outbox is the output box on this rank
 * \param comm is the communicator with all ranks
 *
 * \returns an \b ioboxes struct that holds all boxes across all ranks in the comm
 *
 * Uses MPI_Allgather().
 */
template<typename index>
inline ioboxes<index> gather_boxes(box3d<index> const my_inbox, box3d<index> const my_outbox, MPI_Comm const comm){
    std::array<box3d<index>, 2> my_data = {my_inbox, my_outbox};
    std::vector<box3d<index>> all_boxes(2 * mpi::comm_size(comm), box3d<index>({0, 0, 0}, {0, 0, 0}));
    MPI_Allgather(&my_data, 2 * sizeof(box3d<index>), MPI_BYTE, all_boxes.data(), 2 * sizeof(box3d<index>), MPI_BYTE, comm);
    ioboxes<index> result;
    for(auto i = all_boxes.begin(); i < all_boxes.end(); i += 2){
        result.in.push_back(*i);
        result.out.push_back(*(i+1));
    }
    return result;
}

}

}


#endif /* HEFFTE_GEOMETRY_H */
