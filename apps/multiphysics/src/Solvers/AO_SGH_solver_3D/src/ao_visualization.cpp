#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "ao_visualization.hpp"
#include "ao_gll_converter.hpp"

namespace ao_sgh
{

namespace {

// Tensor-product (i, j, k) -> VTK_LAGRANGE_HEXAHEDRON cell-local slot.
// Verbatim from VTK's reference implementation.
int PointIndexFromIJK(int i, int j, int k, const int* order)
{
    bool ibdy = (i == 0 || i == order[0]);
    bool jbdy = (j == 0 || j == order[1]);
    bool kbdy = (k == 0 || k == order[2]);
    int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

    if (nbdy == 3) {
        return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    }
    int offset = 8;
    if (nbdy == 2) {
        if (!ibdy) {
            return (i - 1)
                 + (j ? order[0] - 1 + order[1] - 1 : 0)
                 + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0)
                 + offset;
        }
        if (!jbdy) {
            return (j - 1)
                 + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1)
                 + (k ? 2 * (order[0] - 1 + order[1] - 1) : 0)
                 + offset;
        }
        offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
        return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
    }
    offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
    if (nbdy == 1) {
        if (ibdy) {
            return (j - 1)
                 + ((order[1] - 1) * (k - 1))
                 + (i ? (order[1] - 1) * (order[2] - 1) : 0)
                 + offset;
        }
        offset += 2 * (order[1] - 1) * (order[2] - 1);
        if (jbdy) {
            return (i - 1)
                 + ((order[0] - 1) * (k - 1))
                 + (j ? (order[2] - 1) * (order[0] - 1) : 0)
                 + offset;
        }
        offset += 2 * (order[2] - 1) * (order[0] - 1);
        return (i - 1)
             + ((order[0] - 1) * (j - 1))
             + (k ? (order[0] - 1) * (order[1] - 1) : 0)
             + offset;
    }
    offset += 2 * ( (order[1] - 1) * (order[2] - 1)
                  + (order[2] - 1) * (order[0] - 1)
                  + (order[0] - 1) * (order[1] - 1));
    return offset
         + (i - 1)
         + (order[0] - 1) * ((j - 1) + (order[1] - 1) * (k - 1));
}

} // anonymous namespace


void write_lagrange_hex_vtu(const std::string&              path,
                            const DCArrayKokkos<size_t>&    nodes_in_elem,
                            const DCArrayKokkos<double>&    coords,
                            const DCArrayKokkos<double>&    vel,
                            const size_t                    num_elems,
                            const size_t                    num_nodes,
                            const size_t                    p_order,
                            const std::vector<ScalarField>& extra_scalars)
{
    const int p     = static_cast<int>(p_order);
    const int n_1d  = p + 1;
    const int npe   = n_1d * n_1d * n_1d;

    // Type-72 assumes equispaced reference nodes; project coords for output.
    DCArrayKokkos<double> equi_coords(num_nodes, 3, "ao_sgh_vtu_equi_coords");
    gll_to_equispaced(nodes_in_elem, coords, equi_coords, num_elems, p_order);

    std::vector<size_t> slot_to_lex((size_t)npe, 0);
    int order[3] = {p, p, p};
    for (int k = 0; k <= p; ++k) {
        for (int j = 0; j <= p; ++j) {
            for (int i = 0; i <= p; ++i) {
                const int    slot = PointIndexFromIJK(i, j, k, order);
                const size_t lex  = (size_t)i
                                  + (size_t)j * n_1d
                                  + (size_t)k * n_1d * n_1d;
                slot_to_lex[(size_t)slot] = lex;
            }
        }
    }

    std::ofstream f(path);
    if (!f) {
        throw std::runtime_error("ao_sgh::write_lagrange_hex_vtu: cannot open " + path);
    }
    f << std::setprecision(17);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << num_nodes
      << "\" NumberOfCells=\"" << num_elems << "\">\n";

    {
        std::string pd_attrs = "Vectors=\"velocity\"";
        if (!extra_scalars.empty()) {
            pd_attrs += " Scalars=\"" + extra_scalars.front().name + "\"";
        }
        f << "      <PointData " << pd_attrs << ">\n";
    }
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
         " Name=\"velocity\" format=\"ascii\">\n";
    for (size_t n = 0; n < num_nodes; ++n) {
        f << "          "
          << vel.host(n, 0) << " "
          << vel.host(n, 1) << " "
          << vel.host(n, 2) << "\n";
    }
    f << "        </DataArray>\n";

    for (const ScalarField& sf : extra_scalars) {
        if (sf.data == nullptr) {
            throw std::runtime_error(
                "ao_sgh::write_lagrange_hex_vtu: ScalarField '" + sf.name +
                "' has null data pointer");
        }
        f << "        <DataArray type=\"Float64\" NumberOfComponents=\"1\""
             " Name=\"" << sf.name << "\" format=\"ascii\">\n";
        for (size_t n = 0; n < num_nodes; ++n) {
            f << "          " << sf.data->host(n) << "\n";
        }
        f << "        </DataArray>\n";
    }
    f << "      </PointData>\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
         " Name=\"Points\" format=\"ascii\">\n";
    for (size_t n = 0; n < num_nodes; ++n) {
        f << "          "
          << equi_coords.host(n, 0) << " "
          << equi_coords.host(n, 1) << " "
          << equi_coords.host(n, 2) << "\n";
    }
    f << "        </DataArray>\n"
      << "      </Points>\n";

    f << "      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        f << "         ";
        for (int s = 0; s < npe; ++s) {
            const size_t lex = slot_to_lex[(size_t)s];
            const size_t gid = nodes_in_elem.host(e, lex);
            f << " " << gid;
        }
        f << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e = 1; e <= num_elems; ++e) {
        f << "          " << (long long)e * npe << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        f << "          72\n";
    }
    f << "        </DataArray>\n"
      << "      </Cells>\n";

    f << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";
}


void sample_thermo_at_kine_gll(const DCArrayKokkos<double>& thermo_coef_per_elem,
                               const DCArrayKokkos<size_t>& nodes_in_elem,
                               const ref_elem_t&            kine_ref,
                               const ref_elem_t&            thermo_ref,
                               const size_t                 num_elems,
                               const size_t                 num_nodes,
                               DCArrayKokkos<double>&       out_kine_dof)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nt = thermo_ref.num_dofs_1d;

    // B_cross(i_kine, j_thermo) = psi_j(xi_i_kine_GLL).
    DCArrayKokkos<double> B_cross(nk, nt, "ao_sgh_vis_B_cross");

    std::vector<double> xi_kine(nk);
    std::vector<double> xi_thermo(nt);
    for (size_t i = 0; i < nk; ++i) xi_kine[i]   = kine_ref.dof_positions_1d(i);
    for (size_t j = 0; j < nt; ++j) xi_thermo[j] = thermo_ref.dof_positions_1d(j);

    for (size_t i = 0; i < nk; ++i) {
        for (size_t j = 0; j < nt; ++j) {
            double phi = 1.0;
            for (size_t m = 0; m < nt; ++m) {
                if (m != j) {
                    phi *= (xi_kine[i] - xi_thermo[m])
                         / (xi_thermo[j] - xi_thermo[m]);
                }
            }
            B_cross.host(i, j) = phi;
        }
    }
    B_cross.update_device();
    Kokkos::fence();

    CArrayKokkos<double> T1(num_elems, nt, nt, nk, "ao_sgh_vis_T1");
    CArrayKokkos<double> T2(num_elems, nt, nk, nk, "ao_sgh_vis_T2");

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t lz = 0; lz < nt; ++lz) {
            for (size_t ly = 0; ly < nt; ++ly) {
                for (size_t ix = 0; ix < nk; ++ix) {
                    double acc = 0.0;
                    for (size_t lx = 0; lx < nt; ++lx) {
                        const size_t lex = lx + ly * nt + lz * nt * nt;
                        acc += B_cross(ix, lx) * thermo_coef_per_elem(elem_gid, lex);
                    }
                    T1(elem_gid, lz, ly, ix) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t lz = 0; lz < nt; ++lz) {
            for (size_t iy = 0; iy < nk; ++iy) {
                for (size_t ix = 0; ix < nk; ++ix) {
                    double acc = 0.0;
                    for (size_t ly = 0; ly < nt; ++ly) {
                        acc += B_cross(iy, ly) * T1(elem_gid, lz, ly, ix);
                    }
                    T2(elem_gid, lz, iy, ix) = acc;
                }
            }
        }
    });
    Kokkos::fence();

    // Shared boundary DoFs: last-write-wins (DG jumps render as that side).
    FOR_ALL(elem_gid, 0, num_elems, {
        for (size_t iz = 0; iz < nk; ++iz) {
            for (size_t iy = 0; iy < nk; ++iy) {
                for (size_t ix = 0; ix < nk; ++ix) {
                    double acc = 0.0;
                    for (size_t lz = 0; lz < nt; ++lz) {
                        acc += B_cross(iz, lz) * T2(elem_gid, lz, iy, ix);
                    }
                    const size_t lex = ix + iy * nk + iz * nk * nk;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    out_kine_dof(gid) = acc;
                }
            }
        }
    });
    Kokkos::fence();
    out_kine_dof.update_host();
}


void sample_qpt_field_element_average(const DCArrayKokkos<double>& qpt_field,
                                      const DCArrayKokkos<double>& detj,
                                      const DCArrayKokkos<size_t>& nodes_in_elem,
                                      const quadrature_t&          quad,
                                      const ref_elem_t&            kine_ref,
                                      const size_t                 num_elems,
                                      const size_t                 num_nodes,
                                      DCArrayKokkos<double>&       out_kine_dof)
{
    const size_t nq = quad.num_qpts_1d;
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t n_qpts_3d = nq * nq * nq;

    const CArrayKokkos<double>& w1 = quad.qpt_weights_1d;

    FOR_ALL(elem_gid, 0, num_elems, {
        double num = 0.0;
        double den = 0.0;
        for (size_t qk = 0; qk < nq; ++qk) {
            for (size_t qj = 0; qj < nq; ++qj) {
                for (size_t qi = 0; qi < nq; ++qi) {
                    const size_t local_q = qi + qj * nq + qk * nq * nq;
                    const size_t g       = elem_gid * n_qpts_3d + local_q;
                    const double w_dJ    = w1(qi) * w1(qj) * w1(qk) * detj(g);
                    num += qpt_field(g) * w_dJ;
                    den += w_dJ;
                }
            }
        }
        const double avg = num / den;

        for (size_t iz = 0; iz < nk; ++iz) {
            for (size_t iy = 0; iy < nk; ++iy) {
                for (size_t ix = 0; ix < nk; ++ix) {
                    const size_t lex = ix + iy * nk + iz * nk * nk;
                    const size_t gid = nodes_in_elem(elem_gid, lex);
                    out_kine_dof(gid) = avg;
                }
            }
        }
    });
    Kokkos::fence();
    out_kine_dof.update_host();
}


namespace {

inline double lagrange_1d(const std::vector<double>& nodes,
                          const size_t              idx,
                          const double              x)
{
    double v = 1.0;
    const size_t n = nodes.size();
    for (size_t m = 0; m < n; ++m) {
        if (m != idx) {
            v *= (x - nodes[m]) / (nodes[idx] - nodes[m]);
        }
    }
    return v;
}

} // anonymous namespace


void write_lagrange_thermo_hex_vtu(const std::string&            path,
                                   const DCArrayKokkos<size_t>&  nodes_in_elem,
                                   const DCArrayKokkos<double>&  coords,
                                   const DCArrayKokkos<double>&  sie_coef,
                                   const DCArrayKokkos<double>&  pres_coef,
                                   const DCArrayKokkos<double>&  den_coef,
                                   const ref_elem_t&             kine_ref,
                                   const ref_elem_t&             thermo_ref,
                                   const size_t                  num_elems,
                                   const size_t                  p_order)
{
    const size_t nk = kine_ref.num_dofs_1d;
    const size_t nt = thermo_ref.num_dofs_1d;
    // Cells are written at the kine order so the geometry is curvilinear.
    // The field samples are evaluations of the degree-(p-1) thermo
    // polynomials, which the order-p cell interpolates exactly -- the
    // rendered fields remain Q_{p-1}.
    const size_t vis_order = p_order;
    const size_t nv_1d     = vis_order + 1;
    const size_t nv_3d     = nv_1d * nv_1d * nv_1d;

    std::vector<double> xi_kine(nk), xi_thermo(nt), xi_equi(nv_1d);
    for (size_t a = 0; a < nk; ++a) xi_kine[a]   = kine_ref.dof_positions_1d(a);
    for (size_t j = 0; j < nt; ++j) xi_thermo[j] = thermo_ref.dof_positions_1d(j);
    for (size_t i = 0; i < nv_1d; ++i) {
        xi_equi[i] = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(vis_order);
    }

    std::vector<double> B_kine_at_vis  (nv_1d * nk);
    std::vector<double> B_thermo_at_vis(nv_1d * nt);
    for (size_t i = 0; i < nv_1d; ++i) {
        for (size_t a = 0; a < nk; ++a) {
            B_kine_at_vis[i * nk + a] = lagrange_1d(xi_kine, a, xi_equi[i]);
        }
        for (size_t j = 0; j < nt; ++j) {
            B_thermo_at_vis[i * nt + j] = lagrange_1d(xi_thermo, j, xi_equi[i]);
        }
    }

    int order[3] = { static_cast<int>(vis_order),
                     static_cast<int>(vis_order),
                     static_cast<int>(vis_order) };
    std::vector<int> slot_for_lex(nv_3d);
    for (size_t kk = 0; kk < nv_1d; ++kk) {
        for (size_t jj = 0; jj < nv_1d; ++jj) {
            for (size_t ii = 0; ii < nv_1d; ++ii) {
                const int    slot = PointIndexFromIJK(
                    static_cast<int>(ii), static_cast<int>(jj),
                    static_cast<int>(kk), order);
                const size_t lex  = ii + jj * nv_1d + kk * nv_1d * nv_1d;
                slot_for_lex[lex] = slot;
            }
        }
    }

    std::ofstream f(path);
    if (!f) {
        throw std::runtime_error(
            "ao_sgh::write_lagrange_thermo_hex_vtu: cannot open " + path);
    }
    f << std::setprecision(17);

    const size_t num_total_points = num_elems * nv_3d;
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << num_total_points
      << "\" NumberOfCells=\"" << num_elems << "\">\n";

    f << "      <PointData Scalars=\"sie\">\n";

    auto eval_thermo = [&](const DCArrayKokkos<double>& coef,
                           size_t e, size_t ix, size_t iy, size_t iz) -> double
    {
        double v = 0.0;
        for (size_t jz = 0; jz < nt; ++jz) {
            const double bz = B_thermo_at_vis[iz * nt + jz];
            for (size_t jy = 0; jy < nt; ++jy) {
                const double by = B_thermo_at_vis[iy * nt + jy];
                for (size_t jx = 0; jx < nt; ++jx) {
                    const double bx = B_thermo_at_vis[ix * nt + jx];
                    const size_t lex = jx + jy * nt + jz * nt * nt;
                    v += bx * by * bz * coef.host(e, lex);
                }
            }
        }
        return v;
    };

    auto write_thermo_field = [&](const char* name,
                                  const DCArrayKokkos<double>& coef)
    {
        f << "        <DataArray type=\"Float64\" NumberOfComponents=\"1\""
             " Name=\"" << name << "\" format=\"ascii\">\n";
        for (size_t e = 0; e < num_elems; ++e) {
            for (size_t iz = 0; iz < nv_1d; ++iz) {
                for (size_t iy = 0; iy < nv_1d; ++iy) {
                    for (size_t ix = 0; ix < nv_1d; ++ix) {
                        f << "          " << eval_thermo(coef, e, ix, iy, iz) << "\n";
                    }
                }
            }
        }
        f << "        </DataArray>\n";
    };

    write_thermo_field("sie",      sie_coef);
    write_thermo_field("pressure", pres_coef);
    write_thermo_field("density",  den_coef);
    f << "      </PointData>\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
         " Name=\"Points\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        for (size_t iz = 0; iz < nv_1d; ++iz) {
            for (size_t iy = 0; iy < nv_1d; ++iy) {
                for (size_t ix = 0; ix < nv_1d; ++ix) {
                    double px = 0.0, py = 0.0, pz = 0.0;
                    for (size_t c = 0; c < nk; ++c) {
                        const double bz = B_kine_at_vis[iz * nk + c];
                        for (size_t b = 0; b < nk; ++b) {
                            const double by = B_kine_at_vis[iy * nk + b];
                            for (size_t a = 0; a < nk; ++a) {
                                const double bx  = B_kine_at_vis[ix * nk + a];
                                const double w   = bx * by * bz;
                                const size_t lex = a + b * nk + c * nk * nk;
                                const size_t gid = nodes_in_elem.host(e, lex);
                                px += w * coords.host(gid, 0);
                                py += w * coords.host(gid, 1);
                                pz += w * coords.host(gid, 2);
                            }
                        }
                    }
                    f << "          " << px << " " << py << " " << pz << "\n";
                }
            }
        }
    }
    f << "        </DataArray>\n"
      << "      </Points>\n";

    f << "      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        const long long base = static_cast<long long>(e * nv_3d);
        std::vector<long long> conn(nv_3d, 0);
        for (size_t lex = 0; lex < nv_3d; ++lex) {
            conn[static_cast<size_t>(slot_for_lex[lex])] = base + static_cast<long long>(lex);
        }
        f << "         ";
        for (size_t s = 0; s < nv_3d; ++s) f << " " << conn[s];
        f << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e = 1; e <= num_elems; ++e) {
        f << "          " << static_cast<long long>(e * nv_3d) << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        f << "          72\n";
    }
    f << "        </DataArray>\n"
      << "      </Cells>\n"
      << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";
}


void write_lagrange_kine_hex_vtu(const std::string&            path,
                                 const DCArrayKokkos<size_t>&  nodes_in_elem,
                                 const DCArrayKokkos<double>&  coords,
                                 const DCArrayKokkos<double>&  vel,
                                 const ref_elem_t&             kine_ref,
                                 const size_t                  num_elems,
                                 const size_t                  p_order)
{
    const size_t nk        = kine_ref.num_dofs_1d;
    const size_t vis_order = p_order;
    const size_t nv_1d     = vis_order + 1;
    const size_t nv_3d     = nv_1d * nv_1d * nv_1d;

    std::vector<double> xi_kine(nk), xi_equi(nv_1d);
    for (size_t a = 0; a < nk; ++a) xi_kine[a] = kine_ref.dof_positions_1d(a);
    for (size_t i = 0; i < nv_1d; ++i) {
        xi_equi[i] = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(vis_order);
    }

    std::vector<double> B_kine_at_vis(nv_1d * nk);
    for (size_t i = 0; i < nv_1d; ++i) {
        for (size_t a = 0; a < nk; ++a) {
            B_kine_at_vis[i * nk + a] = lagrange_1d(xi_kine, a, xi_equi[i]);
        }
    }

    int order[3] = { static_cast<int>(vis_order),
                     static_cast<int>(vis_order),
                     static_cast<int>(vis_order) };
    std::vector<int> slot_for_lex(nv_3d);
    for (size_t kk = 0; kk < nv_1d; ++kk) {
        for (size_t jj = 0; jj < nv_1d; ++jj) {
            for (size_t ii = 0; ii < nv_1d; ++ii) {
                const int    slot = PointIndexFromIJK(
                    static_cast<int>(ii), static_cast<int>(jj),
                    static_cast<int>(kk), order);
                const size_t lex  = ii + jj * nv_1d + kk * nv_1d * nv_1d;
                slot_for_lex[lex] = slot;
            }
        }
    }

    std::ofstream f(path);
    if (!f) {
        throw std::runtime_error(
            "ao_sgh::write_lagrange_kine_hex_vtu: cannot open " + path);
    }
    f << std::setprecision(17);

    const size_t num_total_points = num_elems * nv_3d;
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << num_total_points
      << "\" NumberOfCells=\"" << num_elems << "\">\n";

    f << "      <PointData Vectors=\"velocity\">\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
         " Name=\"velocity\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        for (size_t iz = 0; iz < nv_1d; ++iz) {
            for (size_t iy = 0; iy < nv_1d; ++iy) {
                for (size_t ix = 0; ix < nv_1d; ++ix) {
                    double vx = 0.0, vy = 0.0, vz = 0.0;
                    for (size_t c = 0; c < nk; ++c) {
                        const double bz = B_kine_at_vis[iz * nk + c];
                        for (size_t b = 0; b < nk; ++b) {
                            const double by = B_kine_at_vis[iy * nk + b];
                            for (size_t a = 0; a < nk; ++a) {
                                const double bx  = B_kine_at_vis[ix * nk + a];
                                const double w   = bx * by * bz;
                                const size_t lex = a + b * nk + c * nk * nk;
                                const size_t gid = nodes_in_elem.host(e, lex);
                                vx += w * vel.host(gid, 0);
                                vy += w * vel.host(gid, 1);
                                vz += w * vel.host(gid, 2);
                            }
                        }
                    }
                    f << "          " << vx << " " << vy << " " << vz << "\n";
                }
            }
        }
    }
    f << "        </DataArray>\n"
      << "      </PointData>\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
         " Name=\"Points\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        for (size_t iz = 0; iz < nv_1d; ++iz) {
            for (size_t iy = 0; iy < nv_1d; ++iy) {
                for (size_t ix = 0; ix < nv_1d; ++ix) {
                    double px = 0.0, py = 0.0, pz = 0.0;
                    for (size_t c = 0; c < nk; ++c) {
                        const double bz = B_kine_at_vis[iz * nk + c];
                        for (size_t b = 0; b < nk; ++b) {
                            const double by = B_kine_at_vis[iy * nk + b];
                            for (size_t a = 0; a < nk; ++a) {
                                const double bx  = B_kine_at_vis[ix * nk + a];
                                const double w   = bx * by * bz;
                                const size_t lex = a + b * nk + c * nk * nk;
                                const size_t gid = nodes_in_elem.host(e, lex);
                                px += w * coords.host(gid, 0);
                                py += w * coords.host(gid, 1);
                                pz += w * coords.host(gid, 2);
                            }
                        }
                    }
                    f << "          " << px << " " << py << " " << pz << "\n";
                }
            }
        }
    }
    f << "        </DataArray>\n"
      << "      </Points>\n";

    f << "      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        const long long base = static_cast<long long>(e * nv_3d);
        std::vector<long long> conn(nv_3d, 0);
        for (size_t lex = 0; lex < nv_3d; ++lex) {
            conn[static_cast<size_t>(slot_for_lex[lex])] = base + static_cast<long long>(lex);
        }
        f << "         ";
        for (size_t s = 0; s < nv_3d; ++s) f << " " << conn[s];
        f << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t e = 1; e <= num_elems; ++e) {
        f << "          " << static_cast<long long>(e * nv_3d) << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t e = 0; e < num_elems; ++e) {
        f << "          72\n";
    }
    f << "        </DataArray>\n"
      << "      </Cells>\n"
      << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";
}


void write_pvd_collection(const std::string&           path,
                          const std::vector<PvdEntry>& entries)
{
    std::ofstream f(path);
    if (!f) {
        throw std::runtime_error("ao_sgh::write_pvd_collection: cannot open " + path);
    }
    f << std::setprecision(17);
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <Collection>\n";
    for (const PvdEntry& e : entries) {
        f << "    <DataSet timestep=\"" << e.time
          << "\" part=\""                << e.part
          << "\" name=\""                << e.part_name
          << "\" file=\""                << e.file_path
          << "\"/>\n";
    }
    f << "  </Collection>\n"
      << "</VTKFile>\n";
}

} // end namespace ao_sgh
