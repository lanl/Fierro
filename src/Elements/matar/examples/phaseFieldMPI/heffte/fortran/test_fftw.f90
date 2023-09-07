program HeffteFortranTester
    use heffte_fftw
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none
    type(heffte_fft3d_fftw) :: fft
    integer :: mpi_err, me
    integer :: i
    COMPLEX(C_DOUBLE_COMPLEX), dimension(:), allocatable :: input, output, reference

call MPI_Init(mpi_err)

call MPI_Comm_rank(MPI_COMM_WORLD, me, mpi_err)

if (me == 0) then
    fft = heffte_fft3d_fftw(0, 0, 0, 3, 3, 1, 0, 0, 0, 3, 3, 1, MPI_COMM_WORLD)
else
    fft = heffte_fft3d_fftw(0, 0, 2, 3, 3, 3, 0, 0, 2, 3, 3, 3, MPI_COMM_WORLD)
endif

allocate(input(fft%size_inbox()))
allocate(output(fft%size_outbox()))
allocate(reference(fft%size_outbox()))

do i = 1, fft%size_inbox()
    input(i) = i - 1
enddo
do i = 1, fft%size_outbox()
    reference(i) = 0.0
enddo
reference(1) = -512.0

call fft%forward(input, output)
if (me == 1) then
    do i = 1, fft%size_outbox()
        if (abs(output(i) - reference(i)) > 1.0D-14) then
            write(*,*) "Error in forward() exceeds tolerance of 1.E-14, error is: ", abs(input(i) - reference(i))
            error stop
        endif
    enddo
endif

do i = 1, fft%size_inbox()
    input(i) = 0
enddo

call fft%backward(output, input, scale_fftw_full)

do i = 1, fft%size_inbox()
    if (abs(input(i) - i + 1) > 1.0D-14) then
        write(*,*) "Error in backward() exceeds tolerance of 1.E-14, error is: ", abs(input(i) - i - 1)
        error stop
    endif
enddo

call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

if (me == 0) then
    write(*,*) "FFTW test: OK"
endif

call fft%release()
deallocate(input, output, reference)

call MPI_Finalize(mpi_err)

end program HeffteFortranTester
