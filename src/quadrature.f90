!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! The purpose of the functions and subroutines in this module are for 
!! calculating the integrals of general functions using two different 
!! integaration approximation techniques Booles and Monte Carlo.
!!
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! monte_carlo_quad
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Takes in the array of the evaluated function checks to make sure the 
!! size of the array is divisible by four.
!!
!! Then calculates the integral, using the booles rule, on sections of 5 
!! and adds each section, which increases by delta_x for each iteration
!! reaches the end of the given array. 
!!
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i

    fx_size = size(fx)

    if (mod(fx_size-1, 4) /= 0) then 
       print *, 'fx array size in booles_quadrature has to be divisible by 4.'
       stop
    endif

    s = 0._dp

    do i = 5, fx_size, 4
       s = s + booles_rule(fx(i-4:i), delta_x)
    enddo
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! It takes in the array containing the evaluated function. The array is the
!! fx evaluted at each point and the given delta_x for the integral.
!! 
!! The function then gets the size of the array and makes sure that it is
!! exactly 5 points.
!!
!! Once that is cleared it takes each individual fx evaluation and assigns
!! them to a variable to be evaluated in using the Boole's Rule.
!!
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)
    ! additional test to make sure that the array
    ! received has 5 and only 5 points  
    if (fx_size /= 5) then
        print *, 'fx array size must be equal to 5.'
        stop
    endif
    
    fx0 = 7*fx(1)
    fx1 = 32*fx(2)
    fx2 = 12*fx(3)
    fx3 = 32*fx(4)
    fx4 = 7*fx(5)

    s = (delta_x*2*(fx4+fx3+fx2+fx1+fx0))/45
end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! It takes in the function to be integrated, the lower limit array (a), 
!! upper limit array (b) and the data array.
!!
!! Then checks if the arrays a and b are the same size. If not the program stops.
!!
!! Then initiates a loop that calls a random number between [0,1) rescales
!! it to be between [a,b) and inputs it into the fx(i) array and square it 
!! and inserts it into the gx(i) array. This loop continues n_samples are generated.
!!
!! Then computes the summation of all the fx(i)/num_samples. Then computes
!! the volume of the object given the x_vector array. It then mulitplies
!! the sum with the volmue to compute the result of the integral. 
!! 
!! Then computes the summation of the gx(i)/n_samples. Then subtracts the 
!! sum of the fx(i) in order to get the variance which is used for the Monte Carlo error calculation.
!! 
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:) ! a: n-dimensional, b: n-dimensional, data: n-dimensional
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s

    real(dp) :: f_x_squared, sum_fx, volume
    integer :: i, vector_size, vector_size_b
    real(dp), allocatable :: x_vector(:), fx(:), gx(:)! ...you might need to declare other arrays here


    vector_size = size(a)
    vector_size_b = size(b)

    if (vector_size /= vector_size_b) then
         print *, 'a and b arrays in monte_carlo_quad have to be the same size'
        stop       
    endif

    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples))
    allocate(gx(1:n_samples))

    do i=1,n_samples

        call random_number(x_vector) !generates an array with random numbers in the [0,1) interval
        x_vector = a + x_vector*(b-a) !rescaling to the integration volume [a,b)
        fx(i) = f(x_vector,data)
        gx(i) = f(x_vector,data)**2
    enddo
    ! Integral calculation

    !Sumation calculation
    sum_fx = 0._dp
    do i = 1, n_samples
        sum_fx = sum_fx + (fx(i)/n_samples)
    enddo

    !Volume calculation
    volume = 1._dp
    do i = 1, vector_size
        volume = volume*(b(i)-a(i))
    enddo
    s = volume*sum_fx

    !Error calculation
    do i = 1, n_samples
        f_x_squared = f_x_squared + (gx(i)/n_samples)
    enddo
    sigma_s = volume*sqrt((f_x_squared - sum_fx**2)/n_samples)

end subroutine monte_carlo_quad

end module quadrature