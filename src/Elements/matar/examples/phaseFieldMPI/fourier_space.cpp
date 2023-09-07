#include <iostream>

#include "fourier_space.h"
#include "mpi.h"

FourierSpace::FourierSpace(const std::array<int,3> & glob_nn_real,
                           const std::array<int,3> & loc_nn_cmplx, 
                           const std::array<int,3> & loc_start_index, 
                           const std::array<double,3> & delta) :
kx(loc_nn_cmplx[0]),
ky(loc_nn_cmplx[1]),
kz(loc_nn_cmplx[2])
{
    // set values of kx, ky, and kz
    set_kx_ky_kz(glob_nn_real, loc_nn_cmplx, loc_start_index, delta);
}


void FourierSpace::set_kx_ky_kz(const std::array<int,3> & glob_nn_real, 
                                const std::array<int,3> & loc_nn_cmplx, 
                                const std::array<int,3> & loc_start_index, 
                                const std::array<double,3> & delta)
{
    int nx = glob_nn_real[0];
    int ny = glob_nn_real[1];
    int nz = glob_nn_real[2];
    double dx = delta[0];
    double dy = delta[1];
    double dz = delta[2];
    int xstart = loc_start_index[0];
    int ystart = loc_start_index[1];
    int zstart = loc_start_index[2];

    // calculate kx
    FOR_ALL_CLASS(i, 0, kx.dims(0), {
            int ti;
            ti = i + xstart;
            if (ti > nx/2) ti = ti - nx;
            kx(i) = (double(ti) * twopi) / (double(nx) * dx);
    });

    // calculate ky
    FOR_ALL_CLASS(j, 0, ky.dims(0), {
            int tj;
            tj = j + ystart;
            if (tj > ny/2) tj = tj - ny;
            ky(j) = (double(tj) * twopi) / (double(ny) * dy);
    });

    // calculate kz
    FOR_ALL_CLASS(k, 0, kz.dims(0), {
            int tk;
            tk = k + zstart;
            if (tk > nz/2) tk = tk - nz;
            kz(k) = (double(tk) * twopi) / (double(nz) * dz);
    });
}
