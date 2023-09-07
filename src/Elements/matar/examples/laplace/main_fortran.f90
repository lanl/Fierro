module dimensions
    implicit none
    integer, parameter :: width = 1000
    integer, parameter :: height = 1000
end module dimensions

program main
    use dimensions
    implicit none
    
    double precision, parameter :: temp_tolerance = 0.01
    double precision            :: worst_dt = 100.0
    integer                     :: i, j
    integer                     :: iteration = 1
    real                        :: begin_time, end_time, elapsed

    double precision            :: temperature(0:height+1, 0:width+1)
    double precision            :: temperature_previous(0:height+1, 0:width+1) 

    ! Start measuring time
    call cpu_time(begin_time)

    ! initialize temperature profile
    call initialize(temperature_previous)

    do while (worst_dt > temp_tolerance)
        ! finite difference
        do j = 1, width
            do i = 1, height
                temperature(i,j) = 0.25 * (temperature_previous(i+1,j) &
                                        + temperature_previous(i-1,j)  &
                                        + temperature_previous(i,j+1)  &
                                        + temperature_previous(i,j-1));         
            enddo
        enddo
        
        ! calculate max difference between temperature and temperature_previous
        worst_dt = 0.0
        do j = 1, width
            do i = 1, height
                worst_dt = max(abs(temperature(i,j) - &
                                temperature_previous(i,j)), & 
                                worst_dt);
            enddo
        enddo

        ! update temperature_previous
        do j = 1, width
            do i = 1, height
                temperature_previous(i,j) = temperature(i,j)
            enddo
        enddo

        ! track progress
        if (mod(iteration, 100) == 0) then
            call track_progress(iteration, temperature)
        endif

        iteration = iteration + 1
    enddo

    ! Stop measuring time and calculate the elapsed time
    call cpu_time(end_time)
    elapsed = end_time - begin_time

    write(*, '("Total time was ", f9.6, " seconds.")') elapsed
    write(*, '("Max error at iteration ", i0, " was ", f9.6)') iteration-1, worst_dt

end program main

subroutine initialize(temperature_previous)
    use dimensions
    implicit none

    double precision :: temperature_previous(0:height+1, 0:width+1) 
    integer          :: i,j

    ! initialize temperature_previous to 0.0
    do j = 0, width+1
        do i = 0, height+1
            temperature_previous(i,j) = 0.0            
        enddo
    enddo

    ! setting the left and right boundary conditions
    do i = 0, height+1
        temperature_previous(i,0) = 0.0
        temperature_previous(i,width+1) = (100.0/height)*i
    enddo

    ! setting the top and bottom boundary condition
    do j = 0, width+1
        temperature_previous(0,j) = 0.0
        temperature_previous(height+1,j) = (100.0/width)*j 
    enddo
end subroutine initialize

subroutine track_progress(iteration, temperature)
    use dimensions
    implicit none

    integer          :: iteration
    double precision :: temperature(0:height+1, 0:width+1)
    integer          :: i

    write(*, '("---------- Iteration number: ", i0, " ----------")') iteration
    do i = height-5, height
        write(*, '("["i0,","i0,"]:",f6.2," ")', advance='no') &
                     i, i, temperature(i,i)
    enddo
    write(*,*)
end subroutine track_progress