module module_constants

  implicit none
  
  public :: pi, nvmax, ntmax, nbmax, maxn, maxnz, rho_water, rho_ice, rho_air,  &
            C_w, C_a, C, e, delta_min, p_str, c_d, h_grid, alpha_EVP, beta_EVP, &
            N_total, L_2_accuracy, omega_e, num_time_steps, varepsilon, hour,   &
            T_period, alpha_Anderson, m_Anderson, num_anderson_iterations, N_Picard
  
  real*8, parameter    :: pi = 3.1415926535897931d0                                    ! Class-wide private constant
  integer, parameter   :: nvmax = 14000, ntmax = 30000, nbmax = 10000                  ! max number of verticies, triangles, boundary edges
  integer, parameter   :: maxn = 20000,  maxnz = 100000                                ! max size of systems, nonzero elements
  real*8, parameter    :: rho_water = 997d0                                            ! water density 
  real*8, parameter    :: rho_ice = 900d0                                              ! ice density
  real*8, parameter    :: rho_air = 1225d-3                                            ! air density
  real*8, parameter    :: C_w = 55d-4                                                  ! water-ice drag coefficient
  real*8, parameter    :: C_a = 225d-5                                                 ! air-ice drag coefficient
  real*8, parameter    :: C = 20d0                                                     ! constant in pressure definition
  real*8, parameter    :: e = 2d0                                                      ! ellipticity parameter
  real*8, parameter    :: delta_min = 2d-9                                             ! minimum value of coefficient delta
  real*8, parameter    :: p_str = 275d2                                                ! value of p* in pressure definition
  real*8, parameter    :: c_d = 50d-2                                                  ! value of c_d in numerical viscosity
  real*8, parameter    :: h_grid = 15000d0                                             ! parameter of quazi-uniform mesh
  real*8, parameter    :: alpha_EVP = 5d2, beta_EVP = 5d2                              ! iterative EVP-solver parameters
  integer, parameter   :: N_total = 500                                                ! max number of iterations in EVP-solver
  real*8, parameter    :: L_2_accuracy = 1d-7                                          ! L_2 accuracy in EVP-solver 
  real*8, parameter    :: omega_e = 7.292116d-5                                        ! angle velocity of Earh rotation 
  integer, parameter   :: num_time_steps = 400       
  real*8, parameter    :: varepsilon = 1d-30                                           ! machine epsilon
  real*8, parameter    :: hour = 36d2
  real*8, parameter    :: T_period = 3456d2
  real*8, parameter    :: alpha_Anderson = 5d-5     !1d-12                             ! alpha for anderson acceleration
  integer, parameter   :: m_Anderson = 5                                               ! degrees of freedom for Anderson acceleration
  integer, parameter   :: num_anderson_iterations = 500                                ! number of anderson iterations
  integer, parameter   :: N_Picard = 30                                                ! number of Picard iterations
  integer, parameter   :: strlen = 40                                                  ! length of strings  
  integer, parameter   :: max_number_of_diagonals  = 15                                ! max number of diagonals in sparse matricies    
  integer, parameter   :: N_evp = 300
  

end module module_constants
