// Standalone Pn box-mesh generator writing VTK_LAGRANGE_HEXAHEDRON (.vtu).
// Node positions are equispaced (VTK standard); the solver re-projects to
// GLL at init time. Freestanding: C++17, no Kokkos / MATAR / Fierro deps.
//
// Usage:
//   ao_sgh_box_vtu_gen --nx N --ny N --nz N --p P
//                      [--origin x,y,z] [--length lx,ly,lz]
//                      --out path.vtu

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


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


struct Args {
    int nx = 2, ny = 2, nz = 2;
    int p  = 1;
    double ox = 0.0, oy = 0.0, oz = 0.0;
    double lx = 1.0, ly = 1.0, lz = 1.0;
    std::string out = "box.vtu";
};


void parse_triple(const std::string& s, double& a, double& b, double& c)
{
    auto p1 = s.find(',');
    auto p2 = (p1 == std::string::npos) ? std::string::npos : s.find(',', p1 + 1);
    if (p1 == std::string::npos || p2 == std::string::npos) {
        std::cerr << "expected three comma-separated doubles, got: " << s << "\n";
        std::exit(1);
    }
    a = std::stod(s.substr(0, p1));
    b = std::stod(s.substr(p1 + 1, p2 - p1 - 1));
    c = std::stod(s.substr(p2 + 1));
}


Args parse(int argc, char** argv)
{
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << k << "\n";
                std::exit(1);
            }
            return std::string(argv[++i]);
        };

        if      (k == "--nx") a.nx = std::stoi(next());
        else if (k == "--ny") a.ny = std::stoi(next());
        else if (k == "--nz") a.nz = std::stoi(next());
        else if (k == "--p")  a.p  = std::stoi(next());
        else if (k == "--origin") parse_triple(next(), a.ox, a.oy, a.oz);
        else if (k == "--length") parse_triple(next(), a.lx, a.ly, a.lz);
        else if (k == "--out") a.out = next();
        else if (k == "-h" || k == "--help") {
            std::cout
                << "Usage:\n"
                << "  ao_sgh_box_vtu_gen --nx N --ny N --nz N --p P\n"
                << "                     [--origin x,y,z] [--length lx,ly,lz]\n"
                << "                     --out path.vtu\n";
            std::exit(0);
        }
        else {
            std::cerr << "unknown arg: " << k << "\n";
            std::exit(1);
        }
    }

    if (a.nx < 1 || a.ny < 1 || a.nz < 1) {
        std::cerr << "nx, ny, nz must each be >= 1\n";
        std::exit(1);
    }
    if (a.p < 1) {
        std::cerr << "p must be >= 1\n";
        std::exit(1);
    }
    return a;
}

}  // namespace


int main(int argc, char** argv)
{
    Args a = parse(argc, argv);

    const int p   = a.p;
    const int npi = a.nx * p + 1;
    const int npj = a.ny * p + 1;
    const int npk = a.nz * p + 1;

    const long long num_points = (long long)npi * npj * npk;
    const long long num_elems  = (long long)a.nx * a.ny * a.nz;
    const int npe = (p + 1) * (p + 1) * (p + 1);

    const double dxe = a.lx / (a.nx * p);
    const double dye = a.ly / (a.ny * p);
    const double dze = a.lz / (a.nz * p);

    std::cout
        << "ao_sgh_box_vtu_gen:\n"
        << "  elems   = " << a.nx << " x " << a.ny << " x " << a.nz
        << " = " << num_elems << "\n"
        << "  p_order = " << p << "  (nodes_per_elem = " << npe << ")\n"
        << "  points  = " << npi << " x " << npj << " x " << npk
        << " = " << num_points << "\n"
        << "  origin  = (" << a.ox << ", " << a.oy << ", " << a.oz << ")\n"
        << "  length  = (" << a.lx << ", " << a.ly << ", " << a.lz << ")\n"
        << "  out     = " << a.out << "\n";

    std::vector<long long> conn((size_t)num_elems * npe, -1);
    int order[3] = {p, p, p};

    for (int ek = 0; ek < a.nz; ++ek) {
        for (int ej = 0; ej < a.ny; ++ej) {
            for (int ei = 0; ei < a.nx; ++ei) {
                long long elem = ((long long)ek * a.ny + ej) * a.nx + ei;
                for (int k = 0; k <= p; ++k) {
                    for (int j = 0; j <= p; ++j) {
                        for (int i = 0; i <= p; ++i) {
                            int gi = ei * p + i;
                            int gj = ej * p + j;
                            int gk = ek * p + k;
                            long long pid = ((long long)gk * npj + gj) * npi + gi;
                            int slot = PointIndexFromIJK(i, j, k, order);
                            conn[(size_t)elem * npe + slot] = pid;
                        }
                    }
                }
            }
        }
    }

    std::ofstream f(a.out);
    if (!f) {
        std::cerr << "could not open " << a.out << " for writing\n";
        return 1;
    }
    f << std::setprecision(17);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << num_points
      << "\" NumberOfCells=\"" << num_elems << "\">\n";

    // Single-material ObjectId: reader requires this CellData field.
    f << "      <CellData>\n"
      << "        <DataArray type=\"Int32\" Name=\"ObjectId\" format=\"ascii\">\n";
    for (long long e = 0; e < num_elems; ++e) {
        f << "          0\n";
    }
    f << "        </DataArray>\n"
      << "      </CellData>\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\" format=\"ascii\">\n";
    for (int k = 0; k < npk; ++k) {
        for (int j = 0; j < npj; ++j) {
            for (int i = 0; i < npi; ++i) {
                double x = a.ox + i * dxe;
                double y = a.oy + j * dye;
                double z = a.oz + k * dze;
                f << "          " << x << " " << y << " " << z << "\n";
            }
        }
    }
    f << "        </DataArray>\n"
      << "      </Points>\n";

    f << "      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (long long e = 0; e < num_elems; ++e) {
        f << "         ";
        for (int n = 0; n < npe; ++n) {
            f << " " << conn[(size_t)e * npe + n];
        }
        f << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (long long e = 1; e <= num_elems; ++e) {
        f << "          " << e * npe << "\n";
    }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (long long e = 0; e < num_elems; ++e) {
        f << "          72\n";
    }
    f << "        </DataArray>\n"
      << "      </Cells>\n";

    f << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";

    std::cout << "wrote " << a.out << "\n";
    return 0;
}
