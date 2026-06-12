// End-to-end test: equispaced box -> GLL projection -> tg_vortex IC ->
// VTK_LAGRANGE_HEXAHEDRON .vtu writer.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include <mpi.h>

#include "ao_gll_converter.hpp"
#include "ao_visualization.hpp"

using namespace ao_sgh;
using namespace mtr;


static int  g_num_failures = 0;
static int  g_num_checks   = 0;
static bool g_verbose      = false;


static void check_close(const double actual,
                        const double expected,
                        const double tol,
                        const std::string& label)
{
    g_num_checks++;
    const double diff = std::fabs(actual - expected);
    if (diff > tol) {
        g_num_failures++;
        printf("  [FAIL] %s : actual=%.16e expected=%.16e diff=%.3e tol=%.3e\n",
               label.c_str(), actual, expected, diff, tol);
    }
    else if (g_verbose) {
        printf("  [ ok ] %s : actual=%.16e expected=%.16e diff=%.3e\n",
               label.c_str(), actual, expected, diff);
    }
}


static void check_true(const bool cond, const std::string& label)
{
    g_num_checks++;
    if (!cond) {
        g_num_failures++;
        printf("  [FAIL] %s\n", label.c_str());
    }
    else if (g_verbose) {
        printf("  [ ok ] %s\n", label.c_str());
    }
}


// Equispaced (nx, ny, nz) Pn box on [0,1]^3 with IJK-lex per-elem connectivity.
static void build_equispaced_box(DCArrayKokkos<size_t>& nodes_in_elem,
                                 DCArrayKokkos<double>& coords,
                                 const size_t nx,
                                 const size_t ny,
                                 const size_t nz,
                                 const size_t p_order)
{
    const size_t n_dofs_1d = p_order + 1;
    const size_t npe       = n_dofs_1d * n_dofs_1d * n_dofs_1d;
    const size_t num_elems = nx * ny * nz;

    const size_t npi = nx * p_order + 1;
    const size_t npj = ny * p_order + 1;
    const size_t npk = nz * p_order + 1;
    const size_t num_nodes = npi * npj * npk;

    const double dxe = 1.0 / (nx * p_order);
    const double dye = 1.0 / (ny * p_order);
    const double dze = 1.0 / (nz * p_order);

    coords        = DCArrayKokkos<double>(num_nodes, 3, "vis_test_coords");
    nodes_in_elem = DCArrayKokkos<size_t>(num_elems, npe, "vis_test_nodes_in_elem");

    for (size_t k = 0; k < npk; ++k) {
        for (size_t j = 0; j < npj; ++j) {
            for (size_t i = 0; i < npi; ++i) {
                const size_t gid = k * npi * npj + j * npi + i;
                coords.host(gid, 0) = i * dxe;
                coords.host(gid, 1) = j * dye;
                coords.host(gid, 2) = k * dze;
            }
        }
    }
    coords.update_device();

    for (size_t ek = 0; ek < nz; ++ek) {
        for (size_t ej = 0; ej < ny; ++ej) {
            for (size_t ei = 0; ei < nx; ++ei) {
                const size_t elem_gid = ek * ny * nx + ej * nx + ei;
                for (size_t k = 0; k <= p_order; ++k) {
                    for (size_t j = 0; j <= p_order; ++j) {
                        for (size_t i = 0; i <= p_order; ++i) {
                            const size_t lex = i + j * n_dofs_1d + k * n_dofs_1d * n_dofs_1d;
                            const size_t gi  = ei * p_order + i;
                            const size_t gj  = ej * p_order + j;
                            const size_t gk  = ek * p_order + k;
                            const size_t gid = gk * npi * npj + gj * npi + gi;
                            nodes_in_elem.host(elem_gid, lex) = gid;
                        }
                    }
                }
            }
        }
    }
    nodes_in_elem.update_device();
}


// Canonical tg_vortex velocity IC.
static void apply_tg_vortex_inline(const DCArrayKokkos<double>& coords,
                                   DCArrayKokkos<double>&       vel,
                                   const size_t                 num_nodes)
{
    const double PI = 3.14159265358979323846;
    for (size_t n = 0; n < num_nodes; ++n) {
        const double x = coords.host(n, 0);
        const double y = coords.host(n, 1);
        vel.host(n, 0) =  std::sin(PI * x) * std::cos(PI * y);
        vel.host(n, 1) = -std::cos(PI * x) * std::sin(PI * y);
        vel.host(n, 2) = 0.0;
    }
    vel.update_device();
}


static int count_in_file(const std::string& path, const std::string& needle)
{
    std::ifstream in(path);
    if (!in) return -1;
    std::stringstream ss; ss << in.rdbuf();
    const std::string s = ss.str();
    int n = 0;
    std::string::size_type pos = 0;
    while ((pos = s.find(needle, pos)) != std::string::npos) {
        ++n;
        pos += needle.size();
    }
    return n;
}


static void test_full_pipeline_2x2x2_p3()
{
    printf("\n=== test_full_pipeline_2x2x2_p3 ===\n");

    const size_t nx = 2, ny = 2, nz = 2;
    const size_t p  = 3;

    DCArrayKokkos<size_t> nodes_in_elem;
    DCArrayKokkos<double> coords;
    build_equispaced_box(nodes_in_elem, coords, nx, ny, nz, p);

    const size_t n1d = p + 1;
    const size_t npe = n1d * n1d * n1d;
    const size_t num_elems = nx * ny * nz;
    const size_t num_nodes = (nx * p + 1) * (ny * p + 1) * (nz * p + 1);

    equispaced_to_gll(nodes_in_elem, coords, num_elems, p);

    DCArrayKokkos<double> vel(num_nodes, 3, "vis_test_vel");
    apply_tg_vortex_inline(coords, vel, num_nodes);

    {
        const size_t origin_gid = nodes_in_elem.host(0, 0);
        check_close(vel.host(origin_gid, 0), 0.0, 1e-15, "tg_vortex_u(0,0,0)");
        check_close(vel.host(origin_gid, 1), 0.0, 1e-15, "tg_vortex_v(0,0,0)");
    }

    const std::string path = "ao_mesh_vis_test_box_2x2x2_p3.vtu";
    write_lagrange_hex_vtu(path,
                           nodes_in_elem,
                           coords,
                           vel,
                           num_elems,
                           num_nodes,
                           p);

    {
        std::ifstream in(path);
        check_true(in.good(), "vtu_file_opens");
    }

    check_true(count_in_file(path, "VTK_LAGRANGE_HEXAHEDRON") >= 0,
               "vtu_file_readable");

    {
        std::stringstream piece;
        piece << "NumberOfPoints=\"" << num_nodes
              << "\" NumberOfCells=\"" << num_elems << "\"";
        check_true(count_in_file(path, piece.str()) == 1,
                   "vtu_piece_header_matches");
    }

    check_true(count_in_file(path,
                             "<DataArray type=\"Float64\" NumberOfComponents=\"3\""
                             " Name=\"velocity\" format=\"ascii\"") == 1,
               "vtu_has_velocity_field");

    check_true(count_in_file(path,
                             "<DataArray type=\"Float64\" NumberOfComponents=\"3\""
                             " Name=\"Points\" format=\"ascii\"") == 1,
               "vtu_has_points_array");

    // Slots 0..7 of the first cell connectivity = the 8 element-0 corners
    // in VTK_HEXAHEDRON ordering.
    {
        std::ifstream in(path);
        std::string line;
        bool in_conn = false;
        std::string first_conn_line;
        while (std::getline(in, line)) {
            if (line.find("Name=\"connectivity\"") != std::string::npos) {
                in_conn = true;
                continue;
            }
            if (in_conn) {
                if (line.find("</DataArray>") != std::string::npos) {
                    break;
                }
                if (first_conn_line.empty() && line.find_first_not_of(" \t") != std::string::npos) {
                    first_conn_line = line;
                    break;
                }
            }
        }
        std::stringstream ss(first_conn_line);
        long long slots[8];
        for (int i = 0; i < 8; ++i) ss >> slots[i];

        const size_t n1d_local = n1d;
        auto lex = [&](size_t i, size_t j, size_t k) {
            return i + j * n1d_local + k * n1d_local * n1d_local;
        };
        check_true(static_cast<size_t>(slots[0]) == nodes_in_elem.host(0, lex(0, 0, 0)), "conn_slot0=(0,0,0)");
        check_true(static_cast<size_t>(slots[1]) == nodes_in_elem.host(0, lex(p, 0, 0)), "conn_slot1=(p,0,0)");
        check_true(static_cast<size_t>(slots[2]) == nodes_in_elem.host(0, lex(p, p, 0)), "conn_slot2=(p,p,0)");
        check_true(static_cast<size_t>(slots[3]) == nodes_in_elem.host(0, lex(0, p, 0)), "conn_slot3=(0,p,0)");
        check_true(static_cast<size_t>(slots[4]) == nodes_in_elem.host(0, lex(0, 0, p)), "conn_slot4=(0,0,p)");
        check_true(static_cast<size_t>(slots[5]) == nodes_in_elem.host(0, lex(p, 0, p)), "conn_slot5=(p,0,p)");
        check_true(static_cast<size_t>(slots[6]) == nodes_in_elem.host(0, lex(p, p, p)), "conn_slot6=(p,p,p)");
        check_true(static_cast<size_t>(slots[7]) == nodes_in_elem.host(0, lex(0, p, p)), "conn_slot7=(0,p,p)");
    }

    if (g_num_failures == 0) {
        printf("  wrote and validated %s (%zu cells, %zu DoFs)\n",
               path.c_str(), num_elems, num_nodes);
    }

    (void)npe;
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 1 && std::string(argv[1]) == "-v") {
        g_verbose = true;
    }

    MATAR_INITIALIZE(argc, argv);
    {
        if (rank == 0) {
            printf("ao_mesh_vis_test : equispaced -> GLL -> tg_vortex -> Lagrange-hex .vtu\n");
        }

        test_full_pipeline_2x2x2_p3();

        if (rank == 0) {
            printf("\n=== Summary ===\n");
            printf("checks : %d\n", g_num_checks);
            printf("failed : %d\n", g_num_failures);
            if (g_num_failures == 0) {
                printf("RESULT : PASSED\n");
            }
            else {
                printf("RESULT : FAILED\n");
            }
        }
    }
    MATAR_FINALIZE();

    MPI_Finalize();

    return (g_num_failures == 0) ? 0 : 1;
}
