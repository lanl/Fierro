#include "evpfft.h"
#include "euler.h"
#include <memory>
#include "Profiler.h"

void EVPFFT::write_texture()
{
  Profiler profiler(__FUNCTION__);

  if(my_rank == 0) {

    ag.update_host();

    MatrixTypeRealHost aux33(3,3);
    MatrixTypeRealHost fbar(3,3);
    int ig;
    real_t ph, th, om, one;

    std::unique_ptr <FILE, decltype(&fclose)> file_24 (fopen("tex.out", "w"), &fclose);

    one = 1.0;

    for (int j = 1; j <= 3; j++) {
      for (int i = 1; i <= 3; i++) {
        fbar(i,j) = 0.0;
      }
    }

    for (int i = 1; i <= 3; i++) {
      fbar(i,i) = delt(i);
    }

    fprintf(file_24.get(), "TEX.OUT from FFT\n");
    fprintf(file_24.get(), "Formatted to plot with POLE\n");
    for (int j = 1; j <= 3; j++) {
      for (int i = 1; i <= 3; i++) {
        fprintf(file_24.get(), "%7.3f", fbar(i,j));
      }
    }
    fprintf(file_24.get(), "\n");
    fprintf(file_24.get(), "B %10d\n", npts1*npts2*npts3);


    ig = 0;
    for (int k = 1; k <= npts3; k++) {
      for (int j = 1; j <= npts2; j++) {
        for (int i = 1; i <= npts1; i++) {

          ig += 1;

          for (int jj = 1; jj <= 3; jj++) {
            for (int ii = 1; ii <= 3; ii++) {
              aux33(ii,jj) = ag.host(jj,ii,i,j,k);
            }
          }

          euler(1,ph,th,om,aux33.pointer());

          fprintf(file_24.get(), 
                  "%9.3f%9.3f%9.3f%9.3f"
                  "%5d%5d%5d%5d%5d\n",
                  ph,th,om,one,
                  i,j,k,jgrain(i,j,k),jphase.host(i,j,k)
          );

        } // end for ii
      } // end for jj
    } // end for kk

  } // end if(iwtex .eq. 1)

}
