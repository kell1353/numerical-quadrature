! Program: nuclear_reactor
! By: Austin Keller
!-----------------------------------------------------------------------------
!The program takes an individual real and integer inputs (chosen by the user) 
!and calculates the neutron flux using the booles and monte carlo integration 
!techniques. It calculates the flux in two different scerarios:
!
!The first one, the flux is calculated over a solid box reactor where the distance of the 
!detector is increasing away from the distance from the reactor. The flux is then calculated
!at each distance x. The program also calculates the large x_0 
!approximation of the neutron flux.
!
!The second one, the flux is calculated over a hollow box reactor where the radius of the 
!hollow sphere the hollow space inside the reactor. The flux is then calculated
!at each radius of the sphere. 
!
!It calculates the results of the Boooles and Monte Carlo approximations for both 
!scenarios and writes each scenarios results into seperate .dat to be evaluated in Jupyter.

! You're free to give different values when running the code, the ones 
! bellow are just a suggestion that worked fine for me.

! Initial Parameters
! depth = 40.0_dp
! width = 100.0_dp
! height = 60.0_dp
! y_zero = 30._dp
! x_min = 5_dp
! x_max = 200._dp
! x_step = 5._dp
! n_grid = 25
! m_samples = 10000

! Extra Parameters
! x_zero = 80._dp
! r_min = 1.0_dp
! r_max = 19._dp
! r_step = 1.0_dp
!-----------------------------------------------------------------------------
program nuclear_reactor
use types 

use read_write,  only : read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

implicit none
real(dp) :: depth, width, height, y_zero, x_min, x_max, x_step
integer :: n_grid, m_samples
real(dp) :: x_zero, r_min, r_max, r_step


! Basic part of the project
call read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
call write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)

! Advanced part of the project
call read_advanced_input(x_zero, r_min, r_max, r_step)
call write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)

end program nuclear_reactor