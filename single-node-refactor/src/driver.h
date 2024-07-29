/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
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

#include "io_utils.h"
#include "parse_yaml.h"
#include "solver.h"
#include "simulation_parameters.h"

// Headers for solver classes
#include "sgh_solver.h"

// Physical state data
#include "state.h"


void fill_regions(DCArrayKokkos<reg_fill_t>&, 
                  Material_t&, 
                  mesh_t&, 
                  State_t&, 
                  size_t);

// ==============================================================================
//   Function that returns 1 or 0 if the mesh location is inside an object
// ==============================================================================
size_t fill_geometric_region(const mesh_t& mesh,
                             const DCArrayKokkos<size_t>& voxel_elem_mat_id,
                             const DCArrayKokkos<reg_fill_t>& region_fills,
                             const ViewCArrayKokkos <double>& mesh_coords,
                             const double voxel_dx, 
                             const double voxel_dy, 
                             const double voxel_dz,
                             const double orig_x, 
                             const double orig_y, 
                             const double orig_z,
                             const size_t voxel_num_i, 
                             const size_t voxel_num_j, 
                             const size_t voxel_num_k,
                             const size_t f_id);


// ==============================================================================
//   SGH related fill functions
// ==============================================================================
void fill_regions_sgh(const Material_t& Materials,
                      const mesh_t& mesh,
                      const DCArrayKokkos <double>& node_coords,
                      const DCArrayKokkos <double>& node_vel,
                      DCArrayKokkos <double>& GaussPoint_den,
                      DCArrayKokkos <double>& GaussPoint_sie,
                      DCArrayKokkos <size_t>& elem_mat_id,
                      DCArrayKokkos <reg_fill_t>& region_fills,
                      DCArrayKokkos <size_t>& voxel_elem_mat_id,
                      const DCArrayKokkos <size_t>& read_voxel_file,
                      const size_t num_fills,
                      const size_t num_elems,
                      const size_t num_nodes,
                      const size_t rk_num_bins);

void init_press_sspd_stress(const Material_t& Materials,
                            const mesh_t& mesh,
                            const DCArrayKokkos<double>& MaterialPoints_den,
                            const DCArrayKokkos<double>& MaterialPoints_pres,
                            const DCArrayKokkos<double>& MaterialPoints_stress,
                            const DCArrayKokkos<double>& MaterialPoints_sspd,
                            const DCArrayKokkos<double>& MaterialPoints_sie,
                            const DCArrayKokkos<double>& MaterialPoints_statev,
                            const size_t rk_num_bins,
                            const size_t num_mat_pts,
                            const size_t mat_id);

void init_corner_node_masses_zero(const mesh_t& mesh,
                                  const DCArrayKokkos<double>& node_mass,
                                  const DCArrayKokkos<double>& corner_mass);

void calc_corner_node_masses(const Material_t& Materials,
                      const mesh_t& mesh,
                      const DCArrayKokkos<double>& node_coords,
                      const DCArrayKokkos<double>& node_mass,
                      const DCArrayKokkos<double>& corner_mass,
                      const DCArrayKokkos<double>& MaterialPoints_mass,
                      const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                      const size_t num_mat_elems);        


// ==============================================================================
//   Functions to paint nodal fields onto the mesh
// ==============================================================================
KOKKOS_FUNCTION
void paint_node_vel(const DCArrayKokkos<reg_fill_t>& region_fills,
                    const DCArrayKokkos<double>& node_vel,
                    const DCArrayKokkos<double>& node_coords,
                    const double node_gid,
                    const double num_dims,
                    const size_t f_id,
                    const size_t rk_num_bins);             



// ==============================================================================
//   Functions to fields on the gauss points of the mesh
// ==============================================================================
KOKKOS_FUNCTION
void paint_gauss_den_sie(const Material_t& Materials,
                         const mesh_t& mesh,
                         const DCArrayKokkos <double>& node_coords,
                         const DCArrayKokkos <double>& GaussPoint_den,
                         const DCArrayKokkos <double>& GaussPoint_sie,
                         const DCArrayKokkos <size_t>& elem_mat_id,
                         const DCArrayKokkos<reg_fill_t>& region_fills,
                         const ViewCArrayKokkos <double> elem_coords,
                         const double elem_gid,
                         const size_t f_id);


// ==============================================================================
//   Functions to read voxel mesh
// ==============================================================================
void user_voxel_init(DCArrayKokkos<size_t>& elem_values,
                     double& dx,
                     double& dy,
                     double& dz,
                     double& orig_x,
                     double& orig_y,
                     double& orig_z,
                     size_t& voxel_num_i,
                     size_t& voxel_num_j,
                     size_t& voxel_num_k,
                     double scale_x,
                     double scale_y,
                     double scale_z,
                     std::string mesh_file);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_id_device
///
/// \brief This gives the index value of the point or the elem
///
/// Assumes that the grid has an i,j,k structure
/// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
/// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
///
/// \param i index
/// \param j index
/// \param k index
/// \param Number of i indices
/// \param Number of j indices
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
int get_id_device(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
}

// for string delimiter parsing
std::vector<std::string> split(std::string s, std::string delimiter);

// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str);

const std::string WHITESPACE = " ";

std::string ltrim(const std::string& s);

std::string rtrim(const std::string& s);

std::string trim(const std::string& s);



class Driver
{
public:

    char* mesh_file;
    char* yaml_file;

    // ---------------------------------------------------------------------
    //    input type declarations
    // ---------------------------------------------------------------------

    MeshReader  mesh_reader;
    MeshBuilder mesh_builder;

    SimulationParameters_t SimulationParamaters; ///< the input simulation parameters

    // ---------------------------------------------------------------------
    //    Material and Boundary declarations
    // ---------------------------------------------------------------------

    Material_t Materials;                   ///< Material data for simulation
    BoundaryCondition_t BoundaryConditions; ///< Simulation boundary conditions

    int num_dims = 3;

    // ---------------------------------------------------------------------
    //    mesh data type declarations
    // ---------------------------------------------------------------------
    mesh_t mesh;

    // ---------------------------------------------------------------------
    //    state data type declaration
    // ---------------------------------------------------------------------
    State_t  State;

    int num_solvers = 0;

    // set of enabled solvers
    std::vector<Solver*> solvers;

    Driver(char* YAML)
    {
        yaml_file = YAML;
    };
    ~Driver() {};

    // Initialize driver data.  Solver type, number of solvers
    // Will be parsed from YAML input
    void initialize()
    {
        std::cout << "Initializing Driver" << std::endl;
        Yaml::Node root;
        try
        {
            Yaml::Parse(root, yaml_file);
        }
        catch (const Yaml::Exception e)
        {
            std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
            exit(0);
        }

        parse_yaml(root, SimulationParamaters, Materials, BoundaryConditions);
        std::cout << "Finished  parsing YAML file" << std::endl;

        if (SimulationParamaters.mesh_input.source == mesh_input::file) {
            // Create and/or read mesh
            std::cout << "Mesh file path: " << SimulationParamaters.mesh_input.file_path << std::endl;
            mesh_reader.set_mesh_file(SimulationParamaters.mesh_input.file_path.data());
            mesh_reader.read_mesh(mesh, 
                                  State.GaussPoints, 
                                  State.node, 
                                  State.corner, 
                                  num_dims, 
                                  SimulationParamaters.dynamic_options.rk_num_bins);
        }
        else if (SimulationParamaters.mesh_input.source == mesh_input::generate) {
            mesh_builder.build_mesh(mesh, 
                                    State.GaussPoints, 
                                    State.node, 
                                    State.corner, 
                                    SimulationParamaters);
        }
        else{
            throw std::runtime_error("**** NO MESH INPUT OPTIONS PROVIDED IN YAML ****");
            return;
        }

        // Build boundary conditions
        int num_bcs = BoundaryConditions.num_bcs;
        printf("Num BC's = %d\n", num_bcs);

        // --- calculate bdy sets ---//
        mesh.init_bdy_sets(num_bcs);
        tag_bdys(BoundaryConditions, mesh, State.node.coords);
        mesh.build_boundry_node_sets(mesh);

        // Calculate element volume
        geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);

        
        //fill_regions();
        fill_regions(SimulationParamaters.region_fills, 
                     Materials, 
                     mesh, 
                     State,
                     SimulationParamaters.dynamic_options.rk_num_bins);


        // --- Move the following sovler setup to yaml parsing routine
        // Create solvers
        for (int solver_id = 0; solver_id < SimulationParamaters.solver_inputs.size(); solver_id++) {
            if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGH) {
                SGH* sgh_solver = new SGH(); // , mesh, node, MaterialPoints, corner
                sgh_solver->initialize(SimulationParamaters, Materials, BoundaryConditions, State);
                solvers.push_back(sgh_solver);
            }
        }

    } // end initialize

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls the setup function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup()
    {
        std::cout << "Inside driver setup" << std::endl;
        for (auto& solver : solvers) {
            solver->setup(SimulationParamaters, 
                          Materials, 
                          BoundaryConditions, 
                          mesh, 
                          State);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn run
    ///
    /// \brief Calls the exectue function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void run()
    {
        std::cout << "Inside driver run" << std::endl;
        for (auto& solver : solvers) {
            solver->execute(SimulationParamaters, 
                            Materials, 
                            BoundaryConditions, 
                            mesh, 
                            State);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn finalize
    ///
    /// \brief Calls the finalize function of each of the solvers assuming the 
    ///        finalize function exists and deletes the solver
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void finalize()
    {
        std::cout << "Inside driver finalize" << std::endl;
        for (auto& solver : solvers) {
            if (solver->finalize_flag) {
                solver->finalize(SimulationParamaters, 
                                 Materials, 
                                 BoundaryConditions);
            }
        }
        // destroy FEA modules
        for (auto& solver : solvers) {
            std::cout << "Deleting solver" << std::endl;
            delete solver;
        }
    }

    
}; // end driver class

/////////////////////////////////////////////////////////////////////////////
///
/// \fn fill_regions
///
/// \brief Fills mesh regions based on YAML input
///
/////////////////////////////////////////////////////////////////////////////
void fill_regions(DCArrayKokkos<reg_fill_t>& region_fills, 
                  Material_t& Materials, 
                  mesh_t& mesh, 
                  State_t& State,
                  size_t rk_num_bins)
{


    size_t num_fills = region_fills.size();
    printf("Num Fills's = %zu\n", num_fills);

    // the number of elems and nodes in the mesh
    const size_t num_elems = mesh.num_elems;
    const size_t num_nodes = mesh.num_nodes;



    // create temporary state fields
    // Painting routine requires only 1 material per GaussPoint
    DCArrayKokkos <double> GaussPoint_den(num_elems);
    DCArrayKokkos <double> GaussPoint_sie(num_elems);
    DCArrayKokkos <size_t> elem_mat_id(num_elems); // the mat_id in the elem


    // ---------------------------------------------
    // variables from a voxel file
    // ---------------------------------------------
    DCArrayKokkos<size_t> voxel_elem_mat_id;      // 1 or 0 if material exist, or it is the material_id
    
    DCArrayKokkos<size_t> read_voxel_file(num_fills); // check to see if readVoxelFile
    FOR_ALL(f_id, 0, num_fills, {
        if (region_fills(f_id).volume == region::readVoxelFile)
        {
            read_voxel_file(f_id) = region::readVoxelFile;  // read the  voxel file
        }
        // add other mesh voxel files
        else
        {
            read_voxel_file(f_id) = 0;
        }
    }); // end parallel for
    read_voxel_file.update_host(); // copy to CPU if code is to read a file
    Kokkos::fence();
    // ---------------------------------------------



    // ---------------------------------------------
    // fill den, sie, and velocity on the mesh
    // ---------------------------------------------
    fill_regions_sgh(Materials,
                     mesh,
                     State.node.coords,
                     State.node.vel,
                     GaussPoint_den,
                     GaussPoint_sie,
                     elem_mat_id,
                     region_fills,
                     voxel_elem_mat_id,
                     read_voxel_file,
                     num_fills,
                     num_elems,
                     num_nodes,
                     rk_num_bins);


    // note device and host are updated in the above function
    // ---------------------------------------------

    
    // ----------------------------------------------------------------
    //  Walk over the mesh and find dimensions of material arrays
    // ----------------------------------------------------------------
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // a counter for the Material index spaces
    DCArrayKokkos <size_t> num_elems_saved_for_mat(num_mats);  

    for(int mat_id=0; mat_id<num_mats; mat_id++){

        size_t sum_local;
        size_t sum_total;

        REDUCE_SUM(elem_gid, 0, num_elems, sum_local,{

            if(elem_mat_id(elem_gid) == mat_id){
                // increment the number of elements the materials live in
                sum_local++;
            } // end if    

        }, sum_total);

        // material index space size
        num_elems_saved_for_mat.host(mat_id) = sum_total;

    } // end for



    // ---------------------------------------
    //  SGH allocation of maps and state
    // ---------------------------------------
    State.MaterialToMeshMaps = CArray<MaterialToMeshMap_t> (num_mats);

    State.MaterialPoints  = CArray<MaterialPoint_t> (num_mats);
    State.MaterialCorners = CArray<MaterialCorner_t> (num_mats);
    // zones not needed with SGH
    
    
    // for ALE SGH, add a buffer to num_elems_for_mat, like 10% of num_elems up to num_elems.
    for(int mat_id=0; mat_id<num_mats; mat_id++){

        const size_t num_mat_pts_in_elem = 1; // 1 mat_point per elem with SGH

        // the following always have the exact memory needed, they omit the buffer
        State.MaterialToMeshMaps(mat_id).num_material_elems = num_elems_saved_for_mat.host(mat_id); 
        State.MaterialPoints(mat_id).num_material_points = num_elems_saved_for_mat.host(mat_id)*num_mat_pts_in_elem;
        State.MaterialCorners(mat_id).num_material_corners = num_elems_saved_for_mat.host(mat_id)*mesh.num_nodes_in_elem;

        // -----
        //  Allocation after here will include a buffer
        // -----
        size_t buffer = 0; // memory buffer to push back into
        size_t num_elems_for_mat = num_elems_saved_for_mat.host(mat_id)+buffer; // has a memory buffer for ALE

        size_t num_points_for_mat = num_elems_for_mat*num_mat_pts_in_elem;  
        size_t num_corners_for_mat = num_elems_for_mat*mesh.num_nodes_in_elem;

        State.MaterialToMeshMaps(mat_id).initialize(num_elems_for_mat); 
        State.MaterialPoints(mat_id).initialize(rk_num_bins, num_points_for_mat, 3); // aways 3D, even for 2D-RZ calcs
        State.MaterialCorners(mat_id).initialize(num_corners_for_mat, mesh.num_dims); 
        // zones are not used
    
    } // end for mat_id
    
    // data structures to access indices in other material index spaces
    State.corners_in_mat_elem = corners_in_mat_t(mesh.num_nodes_in_elem); 
    State.points_in_mat_elem  = points_in_mat_t(1);  // 1 material point per element
    // zones_in_mat_elem is not used with SGH



    // now a counter for how many elems have been saved
    for(int mat_id=0; mat_id<num_mats; mat_id++){
        num_elems_saved_for_mat.host(mat_id) = 0; // initializing to zero
    }



    // ---------------------------------------
    //  SGH save data, maps, and state
    // ---------------------------------------
    State.GaussPoints.vol.update_host(); 
    Kokkos::fence();

    // the following loop is not thread safe
    for(size_t elem_gid=0; elem_gid<num_elems; elem_gid++){

        // get the material_id in this element
        size_t mat_id = elem_mat_id.host(elem_gid);

        // mat elem lid (compressed storage) to save the data to, for this material mat_id
        size_t mat_elem_lid = num_elems_saved_for_mat.host(mat_id); 

        // --- mapping from material elem lid to elem ---
        State.MaterialToMeshMaps(mat_id).elem.host(mat_elem_lid) = elem_gid;

        // -----------------------
        // Save MaterialPoints
        // -----------------------

        // LOOP OVER Guass points in the element
        {
            size_t gauss_gid = elem_gid;  // 1 gauss point per element

            size_t mat_point_lid = mat_elem_lid; // for more than 1 gauss point, this must increment

            // --- density and mass ---
            State.MaterialPoints(mat_id).den.host(mat_point_lid)  = GaussPoint_den.host(gauss_gid); 
            State.MaterialPoints(mat_id).mass.host(mat_point_lid) = GaussPoint_den.host(gauss_gid) * State.GaussPoints.vol.host(gauss_gid);

            // --- set eroded flag to false ---
            State.MaterialPoints(mat_id).eroded.host(mat_point_lid) = false;

            // --- specific internal energy ---
            // save state, that is integrated in time, at the RK levels
            for(size_t rk_level=0; rk_level<rk_num_bins; rk_level++){
                State.MaterialPoints(mat_id).sie.host(rk_level, mat_point_lid) = GaussPoint_sie.host(gauss_gid);
            }
        } // end loop over gauss points in element
    

        // -----------------------
        // Save MaterialZones
        // -----------------------
        // For higher-order FE, least squares fit the sie at gauss points to get zone values

        
        // update counter for how many mat_elem_lid values have been saved
        num_elems_saved_for_mat.host(mat_id)++;

    } // end serial for loop over all elements

    // copy the state to the device
    for(int mat_id=0; mat_id<num_mats; mat_id++){
        State.MaterialPoints(mat_id).den.update_device();
        State.MaterialPoints(mat_id).mass.update_device();
        State.MaterialPoints(mat_id).sie.update_device();
        State.MaterialPoints(mat_id).eroded.update_device();

        State.MaterialToMeshMaps(mat_id).elem.update_device();
    } // end for
    Kokkos::fence();


    // calculate pressure, sound speed, and stress for each material
    for(int mat_id=0; mat_id<num_mats; mat_id++){

        size_t num_mat_points = State.MaterialPoints(mat_id).num_material_points;

        init_press_sspd_stress(Materials,
                            mesh,
                            State.MaterialPoints(mat_id).den,
                            State.MaterialPoints(mat_id).pres,
                            State.MaterialPoints(mat_id).stress,
                            State.MaterialPoints(mat_id).sspd,
                            State.MaterialPoints(mat_id).sie,
                            State.MaterialPoints(mat_id).statev,
                            rk_num_bins,
                            num_mat_points,
                            mat_id);

    } // for loop over mat_id


    // set corner and node masses to zero
    init_corner_node_masses_zero(mesh, State.node.mass, State.corner.mass);


    // calculate corner and node masses on the mesh
    if (mesh.num_dims == 3) {
        for(int mat_id=0; mat_id<num_mats; mat_id++){

            size_t num_mat_elems = State.MaterialToMeshMaps(mat_id).num_material_elems;
            
            calc_corner_node_masses(Materials,
                                    mesh,
                                    State.node.coords,
                                    State.node.mass,
                                    State.corner.mass,
                                    State.MaterialPoints(mat_id).mass,
                                    State.MaterialToMeshMaps(mat_id).elem,
                                    num_mat_elems);
        } // end for mat_id
    }
    else{
        // 2D RZ
        // // calculate the corner massess if 2D
        // if (mesh.num_dims == 2) {
        //     FOR_ALL(elem_gid, 0, mesh.num_elems, {
        //         // facial area of the corners
        //         double corner_areas_array[4];

        //         ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);
        //         ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);

        //         geometry::get_area_weights2D(corner_areas, elem_gid, node_coords, elem_node_gids);

        //         // loop over the corners of the element and calculate the mass
        //         for (size_t corner_lid = 0; corner_lid < 4; corner_lid++) {
        //             size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
        //             corner_mass(corner_gid) = corner_areas(corner_lid) * MaterialPoints.den(elem_gid); // node radius is added later
        //         } // end for over corners
        //     });
        //        
        // 
        //    FOR_ALL(nodes_gid=0; nodes_gid<mesh.num_nodes; nodes_gid++){
        //        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
        //            size_t corner_gid    = mesh.corners_in_node(node_gid, corner_lid);
        //            State.node.mass(node_gid) += corner.mass(corner_gid);  // sans the radius so it is areal node mass
        //
        //            corner.mass(corner_gid) *= State.node.coords(1, node_gid, 1); // true corner mass now
        //        } // end for elem_lid
        //    });
        // } // end of
    } // end if 2D



} // end fill regions


// -----------------------------------------------------------------------------
// The function to read a voxel vtk file from Dream3d and intialize the mesh
// ------------------------------------------------------------------------------
void user_voxel_init(DCArrayKokkos<size_t>& elem_values,
                     double& dx,
                     double& dy,
                     double& dz,
                     double& orig_x,
                     double& orig_y,
                     double& orig_z,
                     size_t& num_elems_i,
                     size_t& num_elems_j,
                     size_t& num_elems_k,
                     double scale_x,
                     double scale_y,
                     double scale_z,
                     std::string mesh_file)
{
    std::string MESH = mesh_file; // user specified

    std::ifstream in;  // FILE *in;
    in.open(MESH);

    // check to see of a mesh was supplied when running the code
    if (in)
    {
        printf("\nReading the 3D voxel mesh: ");
        std::cout << MESH << std::endl;
    }
    else
    {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Voxel vtk input does not exist \n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if

    size_t i;           // used for writing information to file
    size_t point_id;    // the global id for the point
    size_t elem_id;     // the global id for the elem
    size_t this_point;   // a local id for a point in a elem (0:7 for a Hexahedral elem)

    size_t num_points_i;
    size_t num_points_j;
    size_t num_points_k;

    size_t num_dims = 3;

    std::string token;

    bool found = false;

    // look for POINTS
    i = 0;
    while (found == false) {
        std::string str;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        //      POINTS %d float
        if (v[0] == "DIMENSIONS")
        {
            num_points_i = std::stoi(v[1]);
            num_points_j = std::stoi(v[2]);
            num_points_k = std::stoi(v[3]);
            printf("Num voxel nodes read in = %zu, %zu, %zu\n", num_points_i, num_points_j, num_points_k);

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find POINTS \n");
            break;
        } // end if

        i++;
    } // end while

    found = false;

    int            num_points = num_points_i * num_points_j * num_points_k;
    CArray<double> pt_coords_x(num_points_i);
    CArray<double> pt_coords_y(num_points_j);
    CArray<double> pt_coords_z(num_points_k);

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "X_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_i - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_x(num_saved) = scale_x*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find X_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "Y_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_j - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_y(num_saved) = scale_y*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find Y_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    while (found == false) {
        std::string str;
        std::string str0;
        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "Z_COORDINATES")
        {
            size_t num_saved = 0;

            while (num_saved < num_points_k - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_coords = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_point = 0; this_point < v_coords.size(); this_point++)
                {
                    pt_coords_z(num_saved) = scale_z*std::stod(v_coords[this_point]);
                    num_saved++;
                } // end for
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find Z_COORDINATES \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    size_t num_elems;
    num_elems_i = num_points_i - 1;
    num_elems_j = num_points_j - 1;
    num_elems_k = num_points_k - 1;

    // center to center distance between first and last elem along each edge
    double Lx = (pt_coords_x(num_points_i - 2) - pt_coords_x(0));
    double Ly = (pt_coords_y(num_points_j - 2) - pt_coords_y(0));
    double Lz = (pt_coords_z(num_points_k - 2) - pt_coords_z(0));

    // spacing between elems
    dx = Lx / ((double) num_elems_i);
    dy = Ly / ((double) num_elems_j);
    dz = Lz / ((double) num_elems_k);

    // element mesh origin
    orig_x = 0.5 * (pt_coords_x(0) + pt_coords_x(1)),
    orig_y = 0.5 * (pt_coords_y(0) + pt_coords_y(1)),
    orig_z = 0.5 * (pt_coords_z(0) + pt_coords_z(1)),

    // look for CELLS
    i = 0;
    while (found == false) {
        std::string str;
        std::getline(in, str);

        std::string              delimiter = " ";
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        //      CELLS num_elems size
        if (v[0] == "CELL_DATA")
        {
            num_elems = std::stoi(v[1]);
            printf("Num voxel elements read in %zu\n", num_elems);

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find CELL_DATA \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    // allocate memory for element voxel values
    elem_values = DCArrayKokkos<size_t>(num_elems);

    // reading the cell data
    while (found == false) {
        std::string str;
        std::string str0;

        std::string delimiter = " ";
        std::getline(in, str);
        std::vector<std::string> v = split(str, delimiter);

        // looking for the following text:
        if (v[0] == "LOOKUP_TABLE")
        {
            size_t num_saved = 0;

            while (num_saved < num_elems - 1) {
                // get next line
                std::getline(in, str0);

                // remove starting and trailing spaces
                str = trim(str0);
                std::vector<std::string> v_values = split(str, delimiter);

                // loop over the contents of the vector v_coords
                for (size_t this_elem = 0; this_elem < v_values.size(); this_elem++)
                {
                    // save integers (0 or 1) to host side
                    elem_values.host(num_saved) = std::stoi(v_values[this_elem]);
                    num_saved++;
                } // end for

                // printf(" done with one row of data \n");
            } // end while

            found = true;
        } // end if

        if (i > 1000)
        {
            printf("ERROR: Failed to find LOOKUP_TABLE data \n");
            break;
        } // end if

        i++;
    } // end while
    found = false;

    printf("\n");

    in.close();
} // end routine


// Code from stackover flow for string delimiter parsing
std::vector<std::string> split(std::string s, std::string delimiter)
{
    size_t                   pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string              token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token     = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
} // end of split


// retrieves multiple values between [ ]
std::vector<double> extract_list(std::string str)
{
    // replace '[' with a space and ']' with a space
    std::replace(str.begin(), str.end(), '[', ' ');
    std::replace(str.begin(), str.end(), ']', ' ');

    std::vector<std::string> str_values;
    std::vector<double>      values;

    // exact the str values into a vector
    str_values = split(str, ",");

    // convert the text values into double values
    for (auto& word : str_values)
    {
        values.push_back(atof(word.c_str()) );
    } // end for

    return values;
}  // end of extract_list


std::string ltrim(const std::string& s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}


std::string rtrim(const std::string& s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string& s)
{
    return rtrim(ltrim(s));
}



/////////////////////////////////////////////////////////////////////////////
///
/// \fn fill_geometric_region
///
/// \brief a function to calculate whether to fill this element based on the 
/// input instructions.  The output is
///  = 0 then no, do not fill this element
///  = 1 then yes, fill this element
///
/// \param mesh is the simulation mesh
/// \param node_coords is the nodal position array
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh 
/// \param region_fills are the instructures to paint state on the mesh
/// \param mesh_coords is the geometric center of the element or a node coordinates
///
/////////////////////////////////////////////////////////////////////////////
size_t fill_geometric_region(const mesh_t& mesh,
                             const DCArrayKokkos<size_t>& voxel_elem_mat_id,
                             const DCArrayKokkos<reg_fill_t>& region_fills,
                             const ViewCArrayKokkos <double>& mesh_coords,
                             const double voxel_dx, 
                             const double voxel_dy, 
                             const double voxel_dz,
                             const double orig_x, 
                             const double orig_y, 
                             const double orig_z,
                             const size_t voxel_num_i, 
                             const size_t voxel_num_j, 
                             const size_t voxel_num_k,
                             const size_t f_id){

    // default is not to fill the element
    size_t fill_this = 0;


    // for shapes with an origin (e.g., sphere and circle), accounting for the origin
    double dist_x = mesh_coords(0) - region_fills(f_id).origin[0];
    double dist_y = mesh_coords(1) - region_fills(f_id).origin[1];
    double dist_z = mesh_coords(2) - region_fills(f_id).origin[2];

    // spherical radius 
    double radius = sqrt(dist_x * dist_x +
                         dist_y * dist_y +
                         dist_z * dist_z);

    // cylindrical radius
    double radius_cyl = sqrt(dist_x * dist_x +
                             dist_y * dist_y);


    // check to see if this element should be filled
    switch (region_fills(f_id).volume) {
        case region::global:
            {
                fill_this = 1;
                break;
            }
        case region::box:
            {

                double x_lower_bound = region_fills(f_id).x1;
                double x_upper_bound = region_fills(f_id).x2;

                double y_lower_bound = region_fills(f_id).y1;
                double y_upper_bound = region_fills(f_id).y2;

                double z_lower_bound = region_fills(f_id).z1;
                double z_upper_bound = region_fills(f_id).z2;


                if (mesh_coords(0) >= x_lower_bound && mesh_coords(0) <= x_upper_bound &&
                    mesh_coords(1) >= y_lower_bound && mesh_coords(1) <= y_upper_bound &&
                    mesh_coords(2) >= z_lower_bound && mesh_coords(2) <= z_upper_bound) {
                    fill_this = 1;
                }
                break;
            }
        case region::cylinder:
            {
                if (radius_cyl >= region_fills(f_id).radius1
                    && radius_cyl <= region_fills(f_id).radius2) {
                    fill_this = 1;
                }
                break;
            }
        case region::sphere:
            {
                if (radius >= region_fills(f_id).radius1
                    && radius <= region_fills(f_id).radius2) {
                    fill_this = 1;
                }
                break;
            }

        case region::readVoxelFile:
            {

                fill_this = 0; // default is no, don't fill it

                // find the closest element in the voxel mesh to this element
                double i0_real = (mesh_coords(0) - orig_x - region_fills(f_id).origin[0]) / (voxel_dx);
                double j0_real = (mesh_coords(1) - orig_y - region_fills(f_id).origin[1]) / (voxel_dy);
                double k0_real = (mesh_coords(2) - orig_z - region_fills(f_id).origin[2]) / (voxel_dz);

                int i0 = (int)i0_real;
                int j0 = (int)j0_real;
                int k0 = (int)k0_real;

                // look for the closest element in the voxel mesh
                int elem_id0 = get_id_device(i0, j0, k0, voxel_num_i, voxel_num_j);

                // if voxel mesh overlaps this mesh, then fill it if =1
                if (elem_id0 < voxel_elem_mat_id.size() && elem_id0 >= 0 &&
                    i0 >= 0 && j0 >= 0 && k0 >= 0 &&
                    i0 < voxel_num_i && j0 < voxel_num_j && k0 < voxel_num_k) {
                    // voxel mesh elem values = 0 or 1
                    fill_this = voxel_elem_mat_id(elem_id0); // values from file

                } // end if

                break;

            } // end case
        case region::no_volume:
            {
                fill_this = 0; // default is no, don't fill it

                break;
            }
        default:
            {
                fill_this = 0; // default is no, don't fill it

                break;
            }

    } // end of switch


    return fill_this;

} // end function



/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_gauss_den_sie
///
/// \brief a function to paint den and sie on the Gauss points of the mesh 
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the node coordinates of the element
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param region_fills are the instructures to paint state on the mesh
/// \param elem_coords is the geometric center of the element
/// \param elem_gid is the element global mesh index
/// \param f_id is fill instruction
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_gauss_den_sie(const Material_t& Materials,
                         const mesh_t& mesh,
                         const DCArrayKokkos <double>& node_coords,
                         const DCArrayKokkos <double>& GaussPoint_den,
                         const DCArrayKokkos <double>& GaussPoint_sie,
                         const DCArrayKokkos <size_t>& elem_mat_id,
                         const DCArrayKokkos<reg_fill_t>& region_fills,
                         const ViewCArrayKokkos <double> elem_coords,
                         const double elem_gid,
                         const size_t f_id){

    // the material id
    size_t mat_id = region_fills(f_id).material_id;

    // --- material_id in elem ---
    elem_mat_id(elem_gid) = mat_id;

    // loop over the Gauss points in the element
    {
        
        const size_t gauss_gid = elem_gid;  // 1 gauss point per element

        // add test problem state setups here
        if (region_fills(f_id).velocity == init_conds::tg_vortex) {

            GaussPoint_den(gauss_gid) = 1.0;    

            // note: elem_coords are the gauss_coords, higher quadrature requires ref elem data
            double pres = 0.25 * (cos(2.0 * PI * elem_coords(0)) + 
                                  cos(2.0 * PI * elem_coords(1)) ) + 1.0;

            // p = rho*ie*(gamma - 1)
            // makes sure index 0 matches the gamma in the gamma law function 
            double gamma  = Materials.eos_global_vars(mat_id,0); 
            GaussPoint_sie(gauss_gid) =
                pres / (GaussPoint_den(gauss_gid) * (gamma - 1.0));
        } // end
        // add user initialization here
        else{
            
            // --- density ---
            GaussPoint_den(gauss_gid) = region_fills(f_id).den;

            // --- specific internal energy ---
            GaussPoint_sie(gauss_gid) = region_fills(f_id).sie;

        }  // end if 
        
    } // end loop over gauss points in element'

    // done setting the element state

} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn paint_node_vel
///
/// \brief a function to paint a velocity field on the nodes of the mesh 
///
/// \param mesh is the simulation mesh
/// \param node_vel is the nodal velocity array
/// \param node_coords are the coordinates of the nodes
/// \param elem_gid is the element global mesh index
/// \param f_id is fill instruction
/// \param rk_num_bins is time integration storage level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void paint_node_vel(const DCArrayKokkos<reg_fill_t>& region_fills,
                    const DCArrayKokkos<double>& node_vel,
                    const DCArrayKokkos<double>& node_coords,
                    const double node_gid,
                    const double num_dims,
                    const size_t f_id,
                    const size_t rk_num_bins){

    // save velocity at all rk_levels
    for(size_t rk_level=0; rk_level<rk_num_bins; rk_level++){

        // --- Velocity ---
        switch (region_fills(f_id).velocity) {
            case init_conds::cartesian:
                {
                    node_vel(rk_level, node_gid, 0) = region_fills(f_id).u;
                    node_vel(rk_level, node_gid, 1) = region_fills(f_id).v;
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = region_fills(f_id).w;
                    }

                    break;
                }
            // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
            case init_conds::radial:
                {
                    // Setting up radial
                    double dir[2];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 2; dim++) {
                        dir[dim]    = node_coords(rk_level, node_gid, dim);
                        radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 2; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    node_vel(rk_level, node_gid, 0) = region_fills(f_id).speed * dir[0];
                    node_vel(rk_level, node_gid, 1) = region_fills(f_id).speed * dir[1];
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = 0.0;
                    }

                    break;
                }
            case init_conds::spherical:
                {
                    // Setting up spherical
                    double dir[3];
                    dir[0] = 0.0;
                    dir[1] = 0.0;
                    dir[2] = 0.0;
                    double radius_val = 0.0;

                    for (int dim = 0; dim < 3; dim++) {
                        dir[dim]    = node_coords(rk_level, node_gid, dim);
                        radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                    } // end for
                    radius_val = sqrt(radius_val);

                    for (int dim = 0; dim < 3; dim++) {
                        if (radius_val > 1.0e-14) {
                            dir[dim] /= (radius_val);
                        }
                        else{
                            dir[dim] = 0.0;
                        }
                    } // end for

                    node_vel(rk_level, node_gid, 0) = region_fills(f_id).speed * dir[0];
                    node_vel(rk_level, node_gid, 1) = region_fills(f_id).speed * dir[1];
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = region_fills(f_id).speed * dir[2];
                    }

                    break;
                }
            case init_conds::radial_linear:
                {
                    printf("**** Radial_linear initial conditions not yet supported ****\n");
                    break;
                }
            case init_conds::spherical_linear:
                {
                    printf("**** spherical_linear initial conditions not yet supported ****\n");
                    break;
                }
            case init_conds::tg_vortex:
                {
                    node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level, node_gid, 0)) * 
                                                        cos(PI * node_coords(rk_level, node_gid, 1));
                    node_vel(rk_level, node_gid, 1) = -1.0 * cos(PI * node_coords(rk_level, node_gid, 0)) * 
                                                        sin(PI * node_coords(rk_level, node_gid, 1));
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = 0.0;
                    }

                    break;
                }

            case init_conds::no_ic_vel:
                {
                    // no velocity
                    node_vel(rk_level, node_gid, 0) = 0.0;
                    node_vel(rk_level, node_gid, 1) = 0.0;
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = 0.0;
                    }

                    break;
                }
            default:
                {
                    // no velocity
                    node_vel(rk_level, node_gid, 0) = 0.0;
                    node_vel(rk_level, node_gid, 1) = 0.0;
                    if (num_dims == 3) {
                        node_vel(rk_level, node_gid, 2) = 0.0;
                    }

                    break;
                }
        } // end of switch

    } // end loop over rk_num_bins


    // done setting the velocity
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn fill_regions_sgh
///
/// \brief a function to paint den, sie, vel, and mat_ids on the mesh 
/// The arrays populated (on host and device) are:
///       elem_mat_id
///       GaussPoint_den
///       GaussPoint_sie
///       node_vel
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the coordinates of the nodes
/// \param node_vel is the nodal velocity array
/// \param region_fills are the instructures to paint state on the mesh
/// \param voxel_elem_mat_id are the voxel values on a structured i,j,k mesh 
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param elem_mat_id is the material id in an element
/// \param num_fills is number of fill instruction
/// \param num_elems is number of elements on the mesh
/// \param num_nodes is number of nodes on the mesh
/// \param rk_num_bins is number of time integration storage bins
///
/////////////////////////////////////////////////////////////////////////////
void fill_regions_sgh(const Material_t& Materials,
                      const mesh_t& mesh,
                      const DCArrayKokkos <double>& node_coords,
                      const DCArrayKokkos <double>& node_vel,
                      DCArrayKokkos <double>& GaussPoint_den,
                      DCArrayKokkos <double>& GaussPoint_sie,
                      DCArrayKokkos <size_t>& elem_mat_id,
                      DCArrayKokkos <reg_fill_t>& region_fills,
                      DCArrayKokkos <size_t>& voxel_elem_mat_id,
                      const DCArrayKokkos <size_t>& read_voxel_file,
                      const size_t num_fills,
                      const size_t num_elems,
                      const size_t num_nodes,
                      const size_t rk_num_bins){


    double voxel_dx, voxel_dy, voxel_dz;          // voxel mesh resolution, set by input file
    double orig_x, orig_y, orig_z;                // origin of voxel elem center mesh, set by input file
    size_t voxel_num_i, voxel_num_j, voxel_num_k; // num voxel elements in each direction, set by input file


    // loop over the fill instructions
    for (size_t f_id = 0; f_id < num_fills; f_id++) {

        // ----
        // voxel mesh setup
        if (read_voxel_file.host(f_id) == region::readVoxelFile)
        {
            // read voxel mesh to get the values in the fcn interface
            user_voxel_init(voxel_elem_mat_id,
                            voxel_dx, 
                            voxel_dy, 
                            voxel_dz,
                            orig_x, 
                            orig_y, 
                            orig_z,
                            voxel_num_i, 
                            voxel_num_j, 
                            voxel_num_k,
                            region_fills(f_id).scale_x,
                            region_fills(f_id).scale_y,
                            region_fills(f_id).scale_z,
                            region_fills(f_id).file_path);

            // copy values read from file to device
            voxel_elem_mat_id.update_device();
        } // endif
        // add else if for other mesh reads including STL-2-voxel


        // parallel loop over elements in mesh
        FOR_ALL(elem_gid, 0, num_elems, {

            // calculate the coordinates and radius of the element
            double elem_coords_1D[3]; // note:initialization with a list won't work
            ViewCArrayKokkos <double> elem_coords(&elem_coords_1D[0], 3);
            elem_coords(0) = 0.0;
            elem_coords(1) = 0.0;
            elem_coords(2) = 0.0;

            // get the coordinates of the element center (using rk_level=1 or node coords)
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_coords(0) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords(1) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_coords(2) += node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_coords(2) = 0.0;
                }
            } // end loop over nodes in element 
            elem_coords(0) = (elem_coords(0) / mesh.num_nodes_in_elem);
            elem_coords(1) = (elem_coords(1) / mesh.num_nodes_in_elem);
            elem_coords(2) = (elem_coords(2) / mesh.num_nodes_in_elem);

            
            // calc if we are to fill this element
            size_t fill_this = fill_geometric_region(mesh, 
                                                     voxel_elem_mat_id, 
                                                     region_fills, 
                                                     elem_coords, 
                                                     voxel_dx, 
                                                     voxel_dy, 
                                                     voxel_dz,
                                                     orig_x, 
                                                     orig_y, 
                                                     orig_z,
                                                     voxel_num_i, 
                                                     voxel_num_j, 
                                                     voxel_num_k,
                                                     f_id);


            // paint the material state on the element if fill_this=1
            if (fill_this == 1) {

                // default sgh paint
                paint_gauss_den_sie(Materials,
                                    mesh,
                                    node_coords,
                                    GaussPoint_den,
                                    GaussPoint_sie,
                                    elem_mat_id,
                                    region_fills,
                                    elem_coords,
                                    elem_gid,
                                    f_id);

                // add user defined paint here
                // user_defined_sgh_state();

            } // end if fill

        }); // end FOR_ALL element loop
        Kokkos::fence();


        // parallel loop over nodes in mesh
        FOR_ALL(node_gid, 0, num_nodes, {

            // make a view to pass to fill and paint functions (using rk_level 1 for node coords)
            ViewCArrayKokkos <double> these_coords(&node_coords(1,node_gid,0), 3);

            
            // calc if we are to fill this element
            size_t fill_this = fill_geometric_region(mesh, 
                                                     voxel_elem_mat_id, 
                                                     region_fills, 
                                                     these_coords, 
                                                     voxel_dx, 
                                                     voxel_dy, 
                                                     voxel_dz,
                                                     orig_x, 
                                                     orig_y, 
                                                     orig_z,
                                                     voxel_num_i, 
                                                     voxel_num_j, 
                                                     voxel_num_k,
                                                     f_id);

            // paint the material state on the node if fill_this=1
            if (fill_this == 1) {

                // default sgh paint
                paint_node_vel(region_fills,
                               node_vel,
                               node_coords,
                               node_gid,
                               mesh.num_dims,
                               f_id,
                               rk_num_bins);

                // add user defined paint here
                // user_defined_vel_state();

            } // end if fill

        }); // end FOR_ALL node loop
        Kokkos::fence();

    } // end for loop over fills


    elem_mat_id.update_host();
    GaussPoint_den.update_host();
    GaussPoint_sie.update_host();
    Kokkos::fence();

} // end SGH fill regions


/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_press_sspd_stress
///
/// \brief a function to initialize pressure, sound speed and stress
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param GaussPoint_den is density at the GaussPoints on the mesh
/// \param GaussPoint_pres is pressure at the GaussPoints on the mesh
/// \param GaussPoint_stress is stress at the GaussPoints on the mesh
/// \param GaussPoint_sspd is sound speed at the GaussPoints on the mesh
/// \param GaussPoint_sie is specific internal energy at the GaussPoints on the mesh
/// \param GaussPoint_statev are the state variables at the GaussPoints on the mesh
/// \param num_mat_pts is the number of material points for mat_id
/// \param mat_id is material id
/// \param rk_num_bins is number of time integration storage bins
///
/////////////////////////////////////////////////////////////////////////////
void init_press_sspd_stress(const Material_t& Materials,
                            const mesh_t& mesh,
                            const DCArrayKokkos<double>& MaterialPoints_den,
                            const DCArrayKokkos<double>& MaterialPoints_pres,
                            const DCArrayKokkos<double>& MaterialPoints_stress,
                            const DCArrayKokkos<double>& MaterialPoints_sspd,
                            const DCArrayKokkos<double>& MaterialPoints_sie,
                            const DCArrayKokkos<double>& MaterialPoints_statev,
                            const size_t rk_num_bins,
                            const size_t num_mat_pts,
                            const size_t mat_id){


    // -------
    // the call to the model initialization goes here
    // -------

    // --- pressure and sound speed ---
    // loop over the material points
    FOR_ALL(mat_point_lid, 0, num_mat_pts, {

        // --- Pressure ---
        Materials.MaterialFunctions(mat_id).calc_pressure(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_lid,
                                        mat_id,
                                        MaterialPoints_statev,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(0, mat_point_lid),
                                        Materials.eos_global_vars);   

        // --- Sound Speed ---                               
        Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_lid,
                                        mat_id,
                                        MaterialPoints_statev,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(0, mat_point_lid),
                                        Materials.eos_global_vars);
    }); // end pressure and sound speed


    // --- stress tensor ---
    for(size_t rk_level=0; rk_level<rk_num_bins; rk_level++){                

        FOR_ALL(mat_point_lid, 0, num_mat_pts, {

            // always 3D even for 2D-RZ
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    MaterialPoints_stress(rk_level, mat_point_lid, i, j) = 0.0;
                }
            }  // end for i,j
                             
        }); // end parallel for over matpt storage

    }// end for rk_level

} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_corner_node_masses_zero
///
/// \brief a function to initialize corner and node masses to zero
///
/// \param mesh is the simulation mesh
/// \param node_mass is the node mass
/// \param corner_mass is the corner mass
///
/////////////////////////////////////////////////////////////////////////////
void init_corner_node_masses_zero(const mesh_t& mesh,
                                  const DCArrayKokkos<double>& node_mass,
                                  const DCArrayKokkos<double>& corner_mass){
                    
    // calculate the nodal mass
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        node_mass(node_gid) = 0.0;
    }); // end parallel over nodes

    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        corner_mass(corner_gid) = 0.0;
    });  // end parallel over corners

} // end setting masses equal to zero


/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_press_sspd_stress
///
/// \brief a function to initialize pressure, sound speed and stress
///
/// \param Materials holds the material models and global parameters
/// \param mesh is the simulation mesh
/// \param node_coords are the nodal coordinates of the mesh
/// \param node_mass is mass of the node
/// \param corner_mass is corner mass
/// \param MaterialPoints_mass is the mass at the material point for mat_id
/// \param num_mat_elems is the number of material elements for mat_id
///
/////////////////////////////////////////////////////////////////////////////
void calc_corner_node_masses(const Material_t& Materials,
                             const mesh_t& mesh,
                             const DCArrayKokkos<double>& node_coords,
                             const DCArrayKokkos<double>& node_mass,
                             const DCArrayKokkos<double>& corner_mass,
                             const DCArrayKokkos<double>& MaterialPoints_mass,
                             const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                             const size_t num_mat_elems){


    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid);  

        // calculate the fraction of matpt mass to scatter to each corner
        double corner_frac = 1.0/((double)mesh.num_nodes_in_elem);  // =1/8
        
        // partion the mass to the corners
        for(size_t corner_lid=0; corner_lid<mesh.num_nodes_in_elem; corner_lid++){
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
            corner_mass(corner_gid) += corner_frac*MaterialPoints_mass(mat_elem_lid);
        } // end for

    }); // end parallel for over mat elem local ids


    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {

            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            node_mass(node_gid) += corner_mass(corner_gid);
        } // end for elem_lid
    }); // end parallel loop over nodes in the mesh

} // end function calculate SGH mass