#include "evpfft.h"
#include "Profiler.h"

void EVPFFT::write_macro_state()
{
  Profiler profiler(__FUNCTION__);

  if (my_rank == 0) {
    fprintf(ofile_mgr.vm_file.get(), "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e\n",
            evm,evmp,dvm,dvmp,svm,svm1);

#if 0
    fprintf(str_str_file, 
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e\n",
            disgradmacro.host(1,1),disgradmacro.host(2,2),
            disgradmacro.host(3,3),disgradmacro.host(2,3),
            disgradmacro.host(1,3),disgradmacro.host(1,2),
            scauav(1,1),scauav(2,2),scauav(3,3),
            scauav(2,3),scauav(1,3),scauav(1,2),
            epav(1,1),epav(2,2),epav(3,3),
            epav(2,3),epav(1,3),epav(1,2),
            velgradmacro.host(1,1),velgradmacro.host(2,2),velgradmacro.host(3,3),
            velgradmacro.host(2,3),velgradmacro.host(1,3),velgradmacro.host(1,2),
            edotpav(1,1),edotpav(2,2),edotpav(3,3),
            edotpav(2,3),edotpav(1,3),edotpav(1,2),
            evm,evmp,dvm,dvmp,svm);
#endif
    // writing a line in str_str_file was split into multiple writes to avoid
    // some strange error I cant figure out.
    fprintf(ofile_mgr.str_str_file.get(),
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,",
            disgradmacro.host(1,1),disgradmacro.host(2,2),
            disgradmacro.host(3,3),disgradmacro.host(2,3),
            disgradmacro.host(1,3),disgradmacro.host(1,2),
            scauav(1,1),scauav(2,2),scauav(3,3),
            scauav(2,3),scauav(1,3),scauav(1,2) );

    fprintf(ofile_mgr.str_str_file.get(),
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,",
            epav(1,1),epav(2,2),epav(3,3),
            epav(2,3),epav(1,3),epav(1,2),
            velgradmacro.host(1,1),velgradmacro.host(2,2),velgradmacro.host(3,3),
            velgradmacro.host(2,3),velgradmacro.host(1,3),velgradmacro.host(1,2) );

    fprintf(ofile_mgr.str_str_file.get(),
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,%011.4e,"
            "%011.4e,%011.4e,%011.4e,%011.4e,%011.4e\n",
            edotpav(1,1),edotpav(2,2),edotpav(3,3),
            edotpav(2,3),edotpav(1,3),edotpav(1,2),
            evm,evmp,dvm,dvmp,svm );
  } // end if (my_rank == 0)
  return;
}
