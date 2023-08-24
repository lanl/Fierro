program main
    ! Note that the use of "using iso_c_binding" is not really needed if the 
    ! correct data types between fortran and c++ are defined
    ! for example fortran real = c++ float.
    ! However, to avoid unnecessary mistakes, specify the kind 
    ! for fortran data type that would be used as arguments in c++ subroutine calls. 
    ! For example real(kind=c_double) for c++ double

    use iso_c_binding

    implicit none

    integer(kind=c_int), parameter :: nx = 4
    integer(kind=c_int), parameter :: ny = 4

    real(kind=c_double) :: array_2D(nx,ny)
    real(kind=c_double) :: sum_of_elements
    integer :: i, j, k, n

    ! initialize kokkos
    write(*,*)'initializing kokkos'
    call kokkos_initialize()

    k = 0
    do j = 1, ny
        do i = 1, nx
          k = k+1
          array_2D(i,j) = k
        enddo
    enddo

    ! calling cpp function
    call square_array_elements(array_2D, nx, ny)

    print*, "printing squared array elements in fortran:"
    do j = 1, ny
        do i = 1, nx
          print*, array_2D(i,j)
        enddo
    enddo

    ! sum array elements and print result
    sum_of_elements = 0.0
    call sum_array_elements(array_2D, nx, ny, sum_of_elements)
    print*, "sum of elements in fortran = ", sum_of_elements
    n = nx*ny;
    print*, "sum of elements should be ", (n*(n+1)*(2*n+1)) / 6

    ! finalize kokkos
    write(*,*)'finalizing kokkos'
    call kokkos_finalize()

end program main
