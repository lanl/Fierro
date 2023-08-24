program HeffteFortranExample
    use heffte_fftw
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none
    type(heffte_fft3d_fftw) :: fft
    integer :: mpi_err, mpi_size, me
    integer :: i
    COMPLEX(C_DOUBLE_COMPLEX), dimension(:), allocatable :: input, output
    REAL(C_DOUBLE) :: err

! This example mimics heffte_example_fftw.cpp but using the Fortran wrappers
! The fftw backend is wrapped in module heffte_fftw which defines class heffte_fft3d_fftw
! The complex types defined by the iso_c_binding are compatible with native Fortran

! Initialize MPI and make sure that we are using two ranks
call MPI_Init(mpi_err)

call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)

call MPI_Comm_rank(MPI_COMM_WORLD, me, mpi_err)

if (mpi_size .ne. 2) then
    if (me == 0) then
        write(*,*) "This example is set for 2 MPI ranks; however, this is set with ", mpi_size
    endif
    call MPI_Finalize(mpi_err)
    call exit(1)
endif

! In Fortran, there is no way to create box3d instances inline without a memory leak
! Therefore, the constructor uses just integers corresponding to the low and high
! indexes of the inbox and outbox
if (me == 0) then
    fft = heffte_fft3d_fftw(0, 0, 0, 3, 3, 1, 0, 0, 0, 3, 3, 1, MPI_COMM_WORLD)
else
    fft = heffte_fft3d_fftw(0, 0, 2, 3, 3, 3, 0, 0, 2, 3, 3, 3, MPI_COMM_WORLD)
endif

! Allocate the input and output buffers
allocate(input(fft%size_inbox()))
allocate(output(fft%size_outbox()))

! Prepare the input
do i = 1, fft%size_inbox()
    input(i) = i - 1
enddo

! The Fortran wrappers use generic overloads to automatically handle the array types
! Thus, forward can be called with any valid combination of single and double precision
! of complex and real inputs. See the documentation of the fft3d class.
call fft%forward(input, output)

! Do a quick sanity check on the result.
if (me == 1) then
    if (abs(output(1) + 512) > 1.0D-14) then
        write(*,*) "Error detected after forward()"
    endif
endif

! Reset the input, will be overwritten with the result of the backend transform
do i = 1, fft%size_inbox()
    input(i) = 0
enddo

! Similar to the forward call, backward comes with overloads
! The scale_fftw_full contains the name of the backend to avoid conflicts in names
! Other modules have similarly scale_cufft_full and scale_rocfft_full
! and those can be used interchangeably
call fft%backward(output, input, scale_fftw_full)

! Compute the local error in the max norm
err = 0
do i = 1, fft%size_inbox()
    err = max( err, abs(input(i) - (i - 1)) )
enddo

! Write out the output one rank at a time
do i = 1, mpi_size
    if (me == i-1) then
        write(*,*) "MPI rank ", me, " observed error: ", err
    endif
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
enddo

! Memory has to be cleaned manually
! The fft%release() call is needed to invoke the C++ destructor
call fft%release()
deallocate(input, output)

call MPI_Finalize(mpi_err)

end program HeffteFortranExample
