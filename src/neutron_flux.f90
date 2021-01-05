!-----------------------------------------------------------------------
!Module: neutron_flux
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The purpose of the subroutines and functions in this module are used to 
!! break up the calculations of the flux using the Booles and Monte Carlo
!! integration techniques for different scenarios. The flux of the solid 
!! box and flux of a hollow box.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! box_flux_monte_carlo
!! hollow_box_flux_mc
!!----------------------------------------------------------------------
!! Included functions:
!!
!! box_flux_booles
!! sphere_flux_booles
!! total_flux_booles
!! sphere_flux_kernel
!! flux_kernel
!! flux_kernel_vector
!! hollow_box_flux_kernel
!! large_x0_flux
!-----------------------------------------------------------------------
module neutron_flux
use types
use quadrature, only : booles_quadrature, monte_carlo_quad

implicit none

private
public :: box_flux_booles, large_x0_flux, box_flux_monte_carlo, hollow_box_flux_mc, total_flux_booles

contains

!-----------------------------------------------------------------------
!! Function: box_flux_booles
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function takes in user inputs and calculates the delta_x, delta_y,
!! delta_z using the depth, width, and height respectively. Then allocates
!! the arrays f_x, g_xy, h_xyz.
!!
!! Then initiates a nested do loop to run the triple integral using the 
!! booles_quadrature approach over a rectangular volume. It initially calculates the flux_kernel
!! from values 0 to height for each step of z. Then integrates that result
!! from values 0 to width for each step of y using the booles approach. Then integrates 
!! that result from values 0 to depth for each step of x using the booles approach. 
!! Then finally returns the value for flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_grid
    
    real(dp) :: delta_x, delta_y, delta_z
    real(dp), allocatable :: f_x(:), g_xy(:), h_xyz(:)
    integer :: n_bins, i_x, i_y, i_z
    real(dp) :: x, y, z

. 

    ! First we need to determine the distance between
    ! the lattice points at which the function to integrate
    ! will be evaluated $\Delta x$.


    ! bins:              1   2   3   4   5   6   7   8
    ! x interval:      |---|---|---|---|---|---|---|---|
    ! grid points:     1   2   3   4   5   6   7   8   9 
    
    ! interval length: |-------------------------------|
    !                  0                               depth
    ! delta x length:  |---|
    !                  0   delta_x

    n_bins = n_grid-1

    delta_x = depth/n_bins
    delta_y = width/n_bins
    delta_z = height/n_bins


    ! allocate memory for the arrays that will contain
    ! the evaluated function to integrate
    allocate(  f_x(1:n_grid))
    allocate( g_xy(1:n_grid))
    allocate(h_xyz(1:n_grid))

    do i_x = 1, n_grid
        x = (i_x-1)*delta_x

        do i_y = 1, n_grid
            y = (i_y-1)*delta_y

            do i_z = 1, n_grid
                z = (i_z-1)*delta_z
                h_xyz(i_z) = flux_kernel(x, y, z, x_zero, y_zero)
            enddo
     
            g_xy(i_y) = booles_quadrature(h_xyz, delta_z)
        enddo

        f_x(i_x) = booles_quadrature(g_xy, delta_y)
    enddo

    flux = booles_quadrature(f_x, delta_x)
end function box_flux_booles

!-----------------------------------------------------------------------
!! Function: flux_kernel
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function establishes the initial function to be integrated
!! over for a rectangular volume.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x            real        x coordinate of the small integration volume
!! y            real        y coordinate of the small integration volume
!! y            real        z coordinate of the small integration volume
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function flux_kernel(x, y, z, x0, y0) result(k)
    implicit none
    real(dp), intent(in) :: x, y, z, x0, y0
    
    k = 1/(4*pi*((x+x0)**2 + (y-y0)**2 + z**2))
end function flux_kernel

!-----------------------------------------------------------------------
!! Function: large_x0_flux
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function is an estimates the booles quadrature for large x0 of the 
!! detector, it approximates the reactor as a point source.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! d            real        Depth of the rectangular nuclear reactor
!! w            real        Width of the rectangular nuclear reactor
!! h            real        Height of the rectangular nuclear reactor
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!!----------------------------------------------------------------------
real(dp) function large_x0_flux(d, w, h, x0, y0) result(flux)
    implicit none
    real(dp), intent(in) :: d, w, h, x0, y0
    
    flux = (d*w*h)/((4*pi)*((x0+(d/2))**2+(y0-(w/2))**2+(h/2)**2))
end function large_x0_flux

!-----------------------------------------------------------------------
!! Subroutine: box_flux_monte_carlo
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Establishes the a, b, and data arrays. Sets the a array as the lower
!! integration limit of solid box, which in this case is 0. Sets the b array
!! as the upper integration limit, which in this case is depth, width, height
!! respectively. Sets the data array as the x_zero, and y_zero coordinate
!! of the detector to be passed on to the kernal.
!!
!! Then calls the monte_carlo_quad function using the flux_kernel_vector
!! as the procedure function, and the a, b, data arrays. Then returns the 
!! value for the flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine box_flux_monte_carlo(depth, width, height, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: flux, sigma_f
    
    real(dp) :: a(1:3), b(1:3), data(1:2)

    ! a is the lower integration limit in the x, y, z coordinates. 
    ! since the origin was placed at the corner of the nuclear reactor the
    ! lower limit is zero in all coordinates
    a = 0._dp

    ! b is the upper integration limit in the x, y, z coordinates.
    b(1) = depth
    b(2) = width
    b(3) = height

    ! This is the 'work array' contains parameters
    ! (other than the sample point) needed to evaluate the function to
    ! integrate
    data(1) = x_zero
    data(2) = y_zero
 
    call monte_carlo_quad(flux_kernel_vector, a, b, data, n_samples, flux, sigma_f)
end subroutine box_flux_monte_carlo

!-----------------------------------------------------------------------
!! Function: flux_kernel_vector
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Establishes the x_vector, data arrays. Sets the x, y, z as values that
!! are contained within x_vector array. Sets x0, y0 as values that are contained
!! within the data array.
!!
!! Then inputs these values into the flux_kernel and returns a 
!! result k, which is the function to be integrated over.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------

real(dp) function flux_kernel_vector(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0

    x = x_vector(1)
    y = x_vector(2)
    z = x_vector(3)
    x0 = data(1)
    y0 = data(2)

    k = flux_kernel(x, y, z, x0, y0)
end function flux_kernel_vector


!-----------------------------------------------------------------------
!! Function: sphere_flux_booles
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function takes in user inputs and calculates the delta_r, delta_theta
!! using the radius and pi respectively. Then allocate the arrays f_r, g_rtheta.
!!
!! Then initiates a nested do loop to run the double integral using the 
!! booles_quadrature approach over a spherical volume. It initially calculates the 
!! sphere_flux_kernel from values 0 to pi for each step of theta. Then integrates that result
!! from values 0 to radius for each step of r using the booles approach.
!! Then finally returns the value for flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! distance     real        Distance from the center of the reactor to the detector
!! radius       real        Radius of the spherical reactor
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function sphere_flux_booles(distance, radius, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: distance, radius
    integer, intent(in) :: n_grid

    real(dp) :: delta_r, delta_theta
    real(dp), allocatable :: f_r(:), g_rtheta(:)
    integer :: n_bins, i_r, i_theta
    real(dp) :: r, theta

    n_bins = n_grid-1

    delta_r = radius/n_bins
    delta_theta = pi/n_bins

    allocate( f_r(1:n_grid))
    allocate(g_rtheta(1:n_grid))


    do i_r = 1, n_grid
        r = (i_r-1)*delta_r

        do i_theta = 1, n_grid
            theta = (i_theta-1)*delta_theta
            g_rtheta(i_theta) = sphere_flux_kernel(r, theta, distance)
        enddo

        f_r(i_r) = booles_quadrature(g_rtheta, delta_r)
    enddo

end function sphere_flux_booles

!-----------------------------------------------------------------------
!! Function: sphere_flux_kernel
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function establishes the initial function to be integrated using 
!! the booles technique over a spherical volume.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! r_prime      real        r coordinate of the small integration volume
!! theta        real        theta coordinate of the small integration volume
!! big_r        integer     distance from the center of the sphere to the detector
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function sphere_flux_kernel(r_prime, theta, big_r) result(k)
    implicit none
    real(dp), intent(in) :: r_prime, theta, big_r
    
    k = (r_prime**2*(sin(theta)))/(2*(r_prime**2 + big_r**2 - 2*r_prime*big_r*cos(theta)))
end function sphere_flux_kernel

!-----------------------------------------------------------------------
!! Function: total_flux_booles
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! This function takes in the user inputs the volume of the box, radius
!! of the sphere, and position of the detector.
!! 
!! It then calculates the distance from the center of the box/sphere to the 
!! the location of the detector.
!!
!! This then uses the initial inputs and calculate the flux of the solid
!! box using the booles quadrature approach. Then uses the distance and 
!! radius to calculate the flux of the sphere using the also using the booles quadrature
!! approach. Then subtractes the flux of the sphere from the flux of the 
!! solid box to get the flux of the hollow box.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_grid

    real(dp) distance, box_flux, sphere_flux

    distance = sqrt(((depth/2) + x_zero)**2 + ((width/2) - y_zero)**2 + (height/2)**2)
    
    ! calculate box booles flux and spherical booles flux
    box_flux = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
    sphere_flux = sphere_flux_booles(distance, radius, n_grid)

    flux = box_flux - sphere_flux

end function total_flux_booles


!-----------------------------------------------------------------------
!! Function: hollow_box_flux_kernel
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Establishes the x_vector, data arrays. Sets the x, y, z as values that
!! are contained within x_vector array. Sets x0, y0, radius, depth, width, height
!! as values that are contained within the data array. 
!!
!! Then uses the depth, width, height to calculate the distance from
!! the sample point x, y, z to the center of the hollow box.
!!
!! Then uses that distance to check weather the samples point is inside 
!! the radius of the sphere (inside the box) if so it set k = 0.
!!
!! If not it then inputs these values into the flux_kernel and returns a 
!! result k, which is the function to be integrated over.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the sphere's radius and x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function hollow_box_flux_kernel(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0, radius, depth, width, height!, .... what other parameters do you need? 
    real(dp) :: distance_to_center

    x = x_vector(1)
    y = x_vector(2)
    z = x_vector(3)

    x0 = data(1)
    y0 = data(2)
    radius = data(3)
    depth = data(4)
    width = data(5)
    height = data(6)

    distance_to_center = sqrt(((depth/2)-x)**2 + ((width/2)-y)**2 + ((height/2)-z)**2)

    if (distance_to_center < radius) then
        k = 0._dp
    else
        k = flux_kernel(x, y, z, x0, y0)
    end if

end function hollow_box_flux_kernel


!-----------------------------------------------------------------------
!! Subroutine: hollow_box_flux_mc
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Establishes the a, b, and data arrays. Sets the a array as the lower
!! integration limit of solid box, which in this case is 0. Sets the b array
!! as the upper integration limit, which in this case is depth, width, height
!! respectively. Sets the data array as the x_zero, y_zero coordinates
!! of the detector, radius, depth, width, height which will be passed on to
!! the kernal.
!!
!! Then calls the monte_carlo_quad function using the hollow_box_flux_kernel_vector
!! as the procedure function, and the a, b, data arrays. Then returns the 
!! value for the flux.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) ::  flux, sigma_f

    real(dp) :: a(1:3), b(1:3), data(1:6)

    a = 0._dp
    ! b is the upper integration limit in the x, y, z coordinates.
    b(1) = depth
    b(2) = width
    b(3) = height

    data(1) = x_zero
    data(2) = y_zero
    data(3) = radius
    data(4) = depth
    data(5) = width
    data(6) = height
    
    call monte_carlo_quad(hollow_box_flux_kernel, a, b, data, n_samples, flux, sigma_f)
end subroutine hollow_box_flux_mc

end module neutron_flux