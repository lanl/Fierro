program HeffteFortranTester
    use heffte_fftw
    use heffte_cufft
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none
    type(heffte_fft3d_fftw) :: fft_cpu
    type(heffte_fft3d_cufft) :: cufft
    integer :: mpi_err, me
    integer :: i
    REAL(C_DOUBLE), dimension(:), allocatable :: input
    COMPLEX(C_DOUBLE_COMPLEX), dimension(:), allocatable :: output

! The main purpose of the test is to make sure that we can load two module without a conflict
! Conflicts can happen in the class name, methods, enums, and helpers (error checking)
! The cufft variable checks for a conflict with declared type of that name
! Does not check for correctness of the input, that's not the point of the test

call MPI_Init(mpi_err)

call MPI_Comm_rank(MPI_COMM_WORLD, me, mpi_err)

fft_cpu = heffte_fft3d_fftw(0, 0, 2 * me, 3, 3, 2 * me + 1, 0, 0, 2 * me, 3, 3, 2 * me + 1, MPI_COMM_WORLD)
cufft = heffte_fft3d_cufft(0, 0, 2 * me, 3, 3, 2 * me + 1, 0, 0, 2 * me, 3, 3, 2 * me + 1, MPI_COMM_WORLD)

if (fft_cpu%size_inbox() .ne. 32) then
    write(*,*) "Wrong inbox size at rank: ", me
    error stop
endif
if (cufft%size_outbox() .ne. 32) then
    write(*,*) "Wrong outbox size at rank: ", me
    error stop
endif

allocate(input(fft_cpu%size_inbox()))
allocate(output(fft_cpu%size_outbox()))

do i = 1, fft_cpu%size_inbox()
    input(i) = i
enddo

call fft_cpu%forward(input, output, scale_fftw_symmetric)

call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

if (me == 0) then
    write(*,*) "FFTW-CUFFT mixed compile test: OK"
endif

call fft_cpu%release()
call cufft%release()
deallocate(input, output)

call MPI_Finalize(mpi_err)

end program HeffteFortranTester
