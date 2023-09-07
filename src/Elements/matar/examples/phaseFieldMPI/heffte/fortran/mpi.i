%{
#include <mpi.h>
%}

/* -------------------------------------------------------------------------
 * \brief MPI Fortran compatibility datatypes
 *
 * Use MPI/Fortran compatibility code to wrap C MPI types, but use ISO-C
 * compatibility layer to pass the Fortran integers (and convert back to native
 * MPI type in the wrapper result code).
 */
%apply int { MPI_Comm };

%typemap(ftype) MPI_Comm
  "integer"
%typemap(fin) MPI_Comm
  "$1 = int($input, C_INT)"
%typemap(fout) MPI_Comm
  "$result = int($1)"

%typemap(in, noblock=1) MPI_Comm {
  $1 = MPI_Comm_f2c((MPI_Fint)*$input);
}
%typemap(out, noblock=1) MPI_Comm {
  $result = (int)MPI_Comm_c2f($1);
}
