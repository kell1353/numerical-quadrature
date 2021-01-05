!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The subroutines in this module are for reading the users input for basic
!! file parameters and checking if the inputs are valid reals and integers. (first subroutine) 
!!
!! Then writing the results of x_0, booles, large x_0, monte carlo, MC uncertainty 
!! into a file called results_basic.dat for analysis in Jupyter. (second subroutine)
!!
!! The subroutines in this module are for reading the users input for advanced
!! file parameters and checking if the inputs are valid reals and integers. (third subroutine) 
!!
!! Then writing the results of radius, box_booles, hollow_booles, hollow_mc, sigma_hollow 
!! into a file called results_advanced.dat for analysis in Jupyter. (fourth subroutine)
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use neutron_flux, only : box_flux_booles, large_x0_flux, box_flux_monte_carlo, total_flux_booles, hollow_box_flux_mc

implicit none

private
public :: read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine provides the user with information about the program and 
!! prompts the user to choose correct values for the program to run effectively.
!! 
!! It then reads each of the the users inputs and converts them into a string 
!! and checks if the the value is a greater than zero real or integer 
!! input depending on the input type. If not there is a message letting the user 
!! know and loops back to allow the user another entry.
!!
!! If there isn't an error the subroutine exits the loop and the program continues.
!!
!! If neither case is true the subroutine prints a message letting the user 
!! know that they entered a blank value and loops back to allow user another entry.
!!
!!----------------------------------------------------------------------
!! Output:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(out) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(out) ::  n_grid, m_samples

    print *, 'This program takes an individual real and integer inputs and calculates the neutron flux'
    print *, 'using the booles and monte carlo integration techniques over different volumes.'
    print *, 'This section you will provide inputs for the box approximations.'

    
    print *, '' ! Just to seperate the variables from description
    depth = read_real('depth, D')
    width = read_real('width, W')
    height = read_real('height, H')
    y_zero = read_real('intial y point, y_zero')
    x_min = read_real('minimum x value, x_min')
    x_max = read_real('maximum x value, x_max')
    x_step = read_real('value of step between x values, x_step')
    
    n_grid = read_integer('number of lattice points, N')
    m_samples = read_integer('number Monte Carlo samples, M')

end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! When reading real input from a user, checks have to be made to make 
!! sure that the user provided the correct type of input. 
!! 
!! We enclose the input reading inside an infinite loop that can only
!! be exited when a correct input is given.
!! 
!! The first check is to make sure that the user input is positive and non-zero.
!!
!! The second check is to make sure that the string is not empty 
!! (i.e. the user simply pressed the enter key)
!! 
!! The third check is made by using the 'read' statement to convert
!! the string into a number, if that is not possible iostat gives an
!! error code different from zero.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Variables:
!!
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    do
        read(*,'(a)',iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a valid number, please provide a positive real number'
        else
            print *, 'that was an empty input, please provide a positive real number'
        endif
    enddo
end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! When reading integer input from a user, checks have to be made to make 
!! sure that the user provided the correct type of input. 
!! 
!! We enclose the input reading inside an infinite loop that can only
!! be exited when a correct input is given.
!! 
!! The first check is to make sure that the user input is positive and non-zero.
!!
!! The second check is to make sure that the string is not empty 
!! (i.e. the user simply pressed the enter key)
!! 
!! The third check is made by using the 'read' statement to convert
!! the string into a number, if that is not possible iostat gives an
!! error code different from zero.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Variables:
!!
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a nonzero positive integer for the '//trim(name)//':'

    do
        read(*,'(a)',iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a valid number, please provide a positive real number'
        else
            print *, 'that was an empty input, please provide a positive real number'
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!Subroutine: write_neutron_flux
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in the users input for depth, width, height, 
!! y_zero, x_min, x_max, x_step, n_grid, m_samples as an input. 
!!
!! Then establishes the initial depth, width, height, y_zero, x_min, x_max, 
!! x_step, n_grid, m_samples. Defines x_zero, box_booles, box_mc, box_large_x0, 
!! sigma_box and the file to write to.
!!
!! Initiates a loop where it calculates the box_flux_booles, large_x0_flux
!! and calls the box_flux_monte_carlo using the users input values. Then 
!! starting at x_zero it begins a do loop incrementing by x_step and 
!! performs the calculations and writes the calculated values into the 
!! file until x_zero is greater than x_max which then exits the subroutine.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: x_zero, box_booles, box_mc, box_large_x0, sigma_box
    character(len=*), parameter :: file_name = 'results_basic.dat'
    integer :: unit

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'x_0', 'booles', 'large x_0', 'monte carlo', 'MC uncertainty'
    x_zero = x_min
    do 
        if(x_zero > x_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        box_large_x0 = large_x0_flux(depth, width, height, x_zero, y_zero)
        call box_flux_monte_carlo(depth, width, height, x_zero, y_zero, m_samples, box_mc, sigma_box)
        write(unit,'(5e28.16)') x_zero, box_booles, box_large_x0, box_mc, sigma_box
        x_zero = x_zero + x_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    Print *, ''
end subroutine write_neutron_flux


!-----------------------------------------------------------------------
!! Subroutine: read_advanced_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine provides the user with information about the program and 
!! prompts the user to choose correct values for the program to run effectively.
!! 
!! It then reads each of the the users inputs and converts them into a string 
!! and checks if the the value is a greater than zero real or integer 
!! input depending on the input type. If not there is a message letting the user 
!! know and loops back to allow the user another entry.
!!
!! If there isn't an error the subroutine exits the loop and the program continues.
!!
!! If neither case is true the subroutine prints a message letting the user 
!! know that they entered a blank value and loops back to allow user another entry.
!!
!---------------------------------------------------------------------
!! Output:
!!
!! x_zero       real        x coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!-----------------------------------------------------------------------
subroutine read_advanced_input(x_zero, r_min, r_max, r_step)
    implicit none
    real(dp), intent(out) :: x_zero, r_min, r_max, r_step
 
    print *, 'This section you will provide inputs for the hollow box approximations.'
    print *, '' ! To breakup inputs from description

    x_zero = read_real('x coordinate of the detector, x_zero')
    r_min = read_real('minimum r value, r_min')
    r_max = read_real('maximum r value, r_max')
    r_step = read_real('value of step between r values, r_step')
end subroutine read_advanced_input

!-----------------------------------------------------------------------
!Subroutine: write_advanced_flux
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This subroutine takes in the users input for ddepth, width, height, x_zero, 
!! y_zero, r_min, r_max, r_step, n_grid, m_samples as an input. 
!!
!! Then establishes the initial depth, width, height, x_zero, y_zero, r_min, 
!! r_max, r_step, n_grid, m_samples. Defines radius, box_booles, hollow_booles, 
!! hollow_mc, sigma_hollow and the file to write to.
!!
!! Calculates the box_flux_booles and it calculates the hollow_booles
!! and calls the hollow_box_flux_mc using the users input values. Initiates a loop 
!! starting at r_min it begins a do loop incrementing by r_step and 
!! performs the calculations and writes the calculated values into the 
!! file until radius is greater than r_max which then exits the subroutine.
!! 
!!----------------------------------------------------------------------
!! Input:
!! 
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero, r_min, r_max, r_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
    character(len=*), parameter :: file_name = 'results_advanced.dat'
    integer :: unit

    open(newunit=unit, file=file_name)
    write(unit,'(5a28)') 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'
    radius = r_min

    box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)

    do 
        if(radius > r_max) exit
        hollow_booles = total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid)
        call hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, m_samples, hollow_mc, sigma_hollow)
        write(unit,'(5e28.16)') radius, box_booles, hollow_booles, hollow_mc, sigma_hollow
        radius = radius + r_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
end subroutine write_advanced_flux

end module read_write
