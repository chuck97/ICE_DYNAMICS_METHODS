program Main
  
  use module_Classes
  use initialization_module
  use module_grid_values
  use module_numerical_integration
  use module_flux_correction
  use module_constants
  use module_assembling
  use module_dynamics
  use module_matricies
  use module_plotwrite
  
  use json_module
      
  
  implicit none
  
  !! global time step, current time
  real*8 :: time_step, time
  
  real*8 :: squ, max_u  
  
  integer :: i, j, k, num, s, r, ps
  
  real*8, allocatable                                  :: solution_anderson(:)
  real*8                                               :: previous_res


! Maximum size of matrix and the maximum number of non-zero entries
  integer :: ide, nonzero_raw_L, nonzero_L, nonzero_raw_C, nonzero_C, m_k
  real*8 , allocatable                                :: rhs_Anderson(:)

! Arrays for the matrix sparse row (CSR) format
  real*8 :: Right_hand_C(maxn), Right_hand_L(maxn), &
     solution_vec_L(maxn), solution_vec_C(maxn), solution_vec(maxn), sum_mass, sum_concentration
  type(Matrix) :: mass_matrix_high, mass_matrix_low, ND_mass_matrix, ND_mass_matrix_lum, &
    N_xx, N_xy, N_yx, N_yy 
    
  real*8, dimension(:,:), allocatable :: residuals_anderson, solutions_anderson, g_anderson
  real*8, allocatable                 :: init_resid(:)
  real*8, allocatable                 :: new_sol(:), new_g(:), new_residual(:)

! Work arrays for ILU factors and 8 BCG vectors
  
  integer, parameter :: MaxWr=maxnz+8*maxn, MaxWi=maxnz+10*maxn+1
  real*8 :: rW(MaxWr)
  integer :: iW(MaxWi)   
   
! ILU0 data
  integer :: ierr, ipaLU, ipjLU, ipjU, ipiw   
   
! Local variables
  integer :: ipBCG   
   
! BiCGStab data
  External :: matvec, prevec0
  Integer :: ITER, INFO, NUNIT
  Real*8 :: RESID
   
  integer, allocatable   :: work(:)
  integer                :: lwork

! External routines from the BLAS library
  real*8 :: ddot
  external :: ddot, dcopy  

! Local variables
  integer                        :: ipBCG_L, ipBCG_C
  integer                        :: imatvec_L(1), iprevec_L(1), imatvec_C(1), iprevec_C(1) 
  real*8                         :: resinit_L, resinit_C, prom, max_C_norm, prom_A(nvmax), x_coord, y_coord
  real*8                         :: u_new(2,nvmax), L2_err, L2_resud, resudal_mass(4000)
  integer                        :: num_iter, resudal_mass_index
  real*8, allocatable            :: prom_m(:)
  type(Matrix)                   :: A_matrix
  real*8, allocatable            :: Big_G_anderson(:, :), Big_F_anderson(:, :), f_anderson(:, :), &
                                     gamma_coefficients(:), residual_pr(:)
  real*8                         :: curr_res    
  integer                        :: nPic             
  Type(Matrix)                   :: Mass_VP, B_VP, lhs_VP    
   
  real*8                         :: lhs_a(maxnz)
  integer                        :: lhs_ia(maxn), lhs_ja(maxnz)  
  real*8                         :: rhs_VP_1(nvmax), rhs_VP_2(nvmax)   
   
   
  real*8                         :: solution_vec_1(maxn), solution_vec_2(maxn)
   
   
   
   ! ILU data
  real*8 :: tau1,tau2,partlur,partlurout
  integer :: verb, UsedWr, UsedWi

! Local variables
  integer :: imatvec_1(1), imatvec_2(1), iprevec_1(1), iprevec_2(1)
  real*8 :: resinit_1, resinit_2     
   
! BiCGStab data
  External :: prevec2
  Integer :: ITER_1, ITER_2, INFO_1, INFO_2, NUNIT_1, NUNIT_2
  Real*8 :: RESID_1, RESID_2
   
   
  !! json variables   
  logical                        :: is_found
  type(json_file)                :: json
  character(len = strlen)        :: path_to_json
    
    
  !! grid variables 
  integer                        :: Nbv, Nbl
  double precision, allocatable  :: bv(:, :)
  integer, allocatable           :: bl(:, :)
  double precision, allocatable  :: bltail(:,:)
   
  double precision, allocatable  :: bv_in(:)
  integer, allocatable           :: bl_in(:)
  double precision,allocatable   :: bltail_in(:)
   
  
  !! pathes to different folders
  character(len=:), allocatable  :: path_to_triangles
  character(len=:), allocatable  :: path_to_nodals
  character(len=:), allocatable  :: path_to_graphics
  
  
  !! boundary type 
  character(:), allocatable      :: boundary_type
  
  !! triangulation parameters
  real*8                  :: square_size
  integer                 :: number_of_vertex_per_square_edge
         
         
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
  
  !! read everything from .json file
  
  call json%initialize()
  
  call get_command_argument(1, path_to_json)
  if (len_trim(path_to_json) == 0) stop "no path to .json file!"
  
  call json%load(filename = trim(path_to_json))
  if (json%failed()) stop "failed to load .json!"
  
  call json%get('time step in hours', time_step, is_found)
  if (.not. is_found) stop "no value of time step!"
    
  call json%get('number of boundary nodes', Nbv, is_found)
  if (.not. is_found) stop "no value of number of boundary nodes!"
  
  call json%get('number of boundary edges', Nbl, is_found)
  if (.not. is_found) stop "no value of number of boundary edges!"
  
  allocate(bv_in(2*Nbv))
  allocate(bl_in(7*Nbl))
  allocate(bltail_in(24))
  
  call json%get('square size', square_size, is_found)
  if (.not. is_found) stop "no value of square size!"
  
  call json%get('number of grid verticies per square edge', number_of_vertex_per_square_edge, is_found)
  if (.not. is_found) stop "no value of number of grid verticies per square edge!"
  
  call json%get('boundary nodes', bv_in, is_found)
  if (.not. is_found) stop "no value of boundary nodes!"
  
  call json%get('boundary edges connectivity list', bl_in, is_found)
  if (.not. is_found) stop "no value of boundary edges connectivity list!"
  
  call json%get('curved boundary list', bltail_in, is_found)
  if (.not. is_found) stop "no value of curved boundary list!"
  
  call json%get('path to graphics folder', path_to_graphics, is_found)
  if (.not. is_found) stop "no value of path to graphics folder!"
  
  call json%get('path to triangles folder', path_to_triangles, is_found)
  if (.not. is_found) stop "no value of path to triangles folder!"
  
  call json%get('path to nodals folder', path_to_nodals, is_found)
  if (.not. is_found) stop "no value of path to nodals folder!"
  
  call json%get('velocity boundary condition', boundary_type, is_found)
  if (.not. is_found) stop "no value of velocity boundary type!"
  
  
  call json%destroy()
  
  print *, ".Json reading: done"
  
  !!!  grid initialization
  
  
  allocate(bv(2,Nbv))
  allocate(bl(7,Nbl))
  allocate(bltail(2,12))
  
  bv = reshape(bv_in, shape(bv))
  bl = reshape(bl_in, shape(bl))
  
  bltail = reshape(bltail_in, shape(bltail))
  
  do i = 1, 2
    do j = 1, Nbv
      bv(i,j) = bv(i,j)*square_size
    end do
  end do  
  
  call grid_initialization( &
  Nbv, Nbl, bv, bl, bltail, number_of_vertex_per_square_edge, square_size)
  
  deallocate(bv)
  deallocate(bl)
  deallocate(bltail)
  
  deallocate(bv_in)
  deallocate(bl_in)
  deallocate(bltail_in)
  
  print *, "Grid generating: done"
  
  !!! scalars, vectors, forcing initialization
  
  call vectors_initialization(boundary_type)
  call scalars_initialization()
  call wind_forcing_initialization(boundary_type)
  call water_forcing_initialization(boundary_type)
  
  print *, "External forcing setup: done"
  
  !!! setup time step in seconds
  
  time_step = time_step*hour
  
  print *, "Time step setup: done"
  
  !!! Assembling non-Direchlet matrix and non-Direchlet lumped matrix
  
  call non_Direchlet_mass_matrix_assembling()
  call non_Direchlet_mass_matrix_lumped_assembling()
   
  print *, "Assembling ND and lumped ND mass matrix: done"
  
  !!! Assembling  M_L and M_C in CSR format for transport equation
  
  call transport_mass_matrix_low_order_assembling(time_step)
  call transport_mass_matrix_high_order_assembling(time_step)
  
  print *, "Assembling mass matricies for transport: done"
  
  
  !!! Time stepping 
  
  do i = 1, number_of_triangles
  
    List_of_Triangles(i)%sigma1 = 0d0
    List_of_Triangles(i)%sigma2 = 0d0
    List_of_Triangles(i)%sigma12 = 0d0
  
  end do
  
  
  do while (num < num_time_steps)
  
    !!! wind and water setup
  
    do i = 1, number_of_non_Direchlet_elements
    
      x_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(1)
      y_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(2)
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1) = 5d0 + &
       (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*x_coord/square_size)*dsin(pi*y_coord/square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2) = 5d0 + &
       (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*y_coord/square_size)*dsin(pi*x_coord/square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1) = &
       (2d0*y_coord - 1d6)/(10d0*square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2) = &
       -(2d0*x_coord - 1d6)/(10d0*square_size)  
       
      List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = &
      List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = &
      List_of_non_Direchlet_Elements(i)%pointing_element%u(2) 
       
    
    end do
    
    !! compute initial residual each time step
    
    init_resid = init_residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, &
      List_of_Triangles, number_of_non_Direchlet_elements, number_of_triangles, time_step, 0)
  
  
    num_iter = 1
    
    call dot_epsilon_delta_recalculation(List_of_Triangles, number_of_triangles, 2)
    call P_0_recalculation(List_of_Triangles, number_of_triangles)
    
    !! mEVP-stepping
    
    do while((num_iter < (N_evp+1)))
    
      call dot_epsilon_delta_recalculation(List_of_Triangles, number_of_triangles, 1)
      call P_0_recalculation(List_of_Triangles, number_of_triangles)
      
      !!! sigma recalculation
    
      do i = 1, number_of_Triangles
      
        List_of_Triangles(i)%sigma1 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma1 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon1 - &
         List_of_Triangles(i)%delta)* &
         (List_of_Triangles(i)%P_0)/ &
         (List_of_Triangles(i)%delta + delta_min)
         
        List_of_Triangles(i)%sigma2 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma2 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon2)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2)
         
        List_of_Triangles(i)%sigma12 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma12 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon12)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2) 
         
      end do
      
      call velocity_recalculation(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
       ND_mass_matrix_lum, time_step)
      
      do i = 1, number_of_non_Direchlet_elements
  
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(1)
    
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(2)
       
        
        List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
    
        List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
        
      end do
      
      residual_array(num_iter) = L2_norm(residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
           number_of_non_Direchlet_elements, number_of_triangles, time_step, 0), 2*number_of_non_Direchlet_elements)/ &
           L2_norm(init_resid, 2*number_of_non_Direchlet_elements)
      
      print *, "iter", num_iter !, residual_array(num_iter)
      
      num_iter = num_iter + 1
        
    end do
      
      !call print_resudal(residual_array, N_evp, num)  
        
      
      print *, "max_u:", maximum_u(List_of_Elements, number_of_elements)
      
      do i = 1, number_of_non_Direchlet_elements
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
         
        List_of_non_Direchlet_Elements(i)%pointing_element%u(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) 
      
      end do
       
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         Mass transport         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      print *, "Mass transport:"
      
      do k = 1, number_of_elements
      
        List_of_Elements(k)%element_value = List_of_Elements(k)%m
      
      end do 
      
      do k = 1, number_of_elements
      
        solution_vec_L(k) = 0d0 
        solution_vec_C(k) = 0d0
        
      
      end do
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Right hand calculation for M_L (m) !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call transport_right_hand_low_order_assembling(Right_hand_L, mass_matrix_low, List_of_Elements, &
      number_of_elements, time_step)
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Right hand calculation for M_C (m) !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      call transport_right_hand_high_order_assembling(Right_hand_C, List_of_Elements, &
      number_of_elements, time_step)
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! set initial guess and compute initial residual for M_C !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
      !set initial guess to 0
      Call dcopy(number_of_elements,0d0,0,solution_vec_C,1)

     !compute initial residual norm
      resinit_C = ddot(number_of_elements,Right_hand_C,1,Right_hand_C,1)
      resinit_C = dsqrt(resinit_C)
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!           Solution for M_L   (m)      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      do k = 1, number_of_elements
      
        solution_vec_L(k) = Right_hand_L(k)/mass_matrix_low%a(k) 
      
      end do
          
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!     Iterative solution for M_C  (m)   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
          
      do k = 1, 3
      
        call sparce_matrix_vector_multiplication(mass_matrix_high, solution_vec_C, number_of_elements, prom_m) 
      
        do i = 1, number_of_elements
        
          solution_vec_C(i) = solution_vec_C(i) + (Right_hand_C(i) - prom_m(i))/mass_matrix_low%a(i)
        
        end do
      
         
      
      end do 
       
                   
      do k = 1, number_of_elements
      
        solution_vec_C(k) = solution_vec_C(k) + List_of_Elements(k)%element_value
        solution_vec_L(k) = solution_vec_L(k) + List_of_Elements(k)%element_value
      
      end do             
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Updating element values using flux correction method !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      call flux_correction_procedure(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements, &
      mass_matrix_low%a, solution_vec_C, solution_vec_L)
      
      
      do k = 1, number_of_elements
      
        
        List_of_Elements(k)%m = List_of_Elements(k)%element_value ! List_of_Elements(k)%element_value
        
      
      end do
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!        Concentration transport         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      print *, "Concentration transport:"
      
      do k = 1, number_of_elements
      
        solution_vec_L(k) = 0d0 
        solution_vec_C(k) = 0d0
        
        Right_hand_C(k) = 0d0
        Right_hand_L(k) = 0d0
      
      end do
      
      
      do k = 1, number_of_elements
      
        List_of_Elements(k)%element_value = List_of_Elements(k)%A
        
      end do 
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Right hand calculation for M_L (A) !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call transport_right_hand_low_order_assembling(Right_hand_L, mass_matrix_low, List_of_Elements, &
      number_of_elements, time_step)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Right hand calculation for M_C (A) !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      call transport_right_hand_high_order_assembling(Right_hand_C, List_of_Elements, &
      number_of_elements, time_step)
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! set initial guess and compute initial residual for M_C !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
      !set initial guess to 0
      Call dcopy(number_of_elements,0d0,0,solution_vec_C,1)

     !compute initial residual norm
      resinit_C = ddot(number_of_elements,Right_hand_C,1,Right_hand_C,1)
      resinit_C = dsqrt(resinit_C)
      
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!           Solution for M_L   (A)      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      do k = 1, number_of_elements
      
        solution_vec_L(k) = Right_hand_L(k)/mass_matrix_low%a(k) 
      
      end do
          
       
      do k = 1, 3
      
        call sparce_matrix_vector_multiplication(mass_matrix_high, solution_vec_C, number_of_elements, prom_A)
      
        do i = 1, number_of_elements
        
          solution_vec_C(i) = solution_vec_C(i) + (Right_hand_C(i) - prom_A(i))/mass_matrix_low%a(i)
        
        end do
      
         
      
      end do     
          
          
          
      do k = 1, number_of_elements
      
        solution_vec_C(k) = solution_vec_C(k) + List_of_Elements(k)%element_value
        solution_vec_L(k) = solution_vec_L(k) + List_of_Elements(k)%element_value
      
      end do     
                
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Updating element values using flux correction method !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      call flux_correction_procedure(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements, &
      mass_matrix_low%a, solution_vec_C, solution_vec_L)
      
      
      
      do k = 1, number_of_elements
             
        if (List_of_Elements(k)%element_value > 1d0) then
        
          List_of_Elements(k)%A = 1d0
           
        else
        
          List_of_Elements(k)%A = List_of_Elements(k)%element_value   
           
        end if
        
      end do
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!                     h recalculation                  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_elements
       
        List_of_Elements(i)%h = (List_of_Elements(i)%m)/(List_of_Elements(i)%A*rho_ice)
    
      end do
      
      
      
      !! next time step
      
      num = num + 1
      time = time + time_step
      
      print *, "time moment: ", time/hour , "hours"
      
      call plot_ParaView(num, path_to_graphics)
      call write_nodals(num, path_to_nodals)
      call write_triangles(num, path_to_triangles)
      
      
      print *, "iteration: ", num
      print *, "sum_mass:", L2_mass(List_of_Triangles, number_of_triangles, 1)
      print *, "sum_thickness:", L2_mass(List_of_Triangles, number_of_triangles, 2)
      print *, "sum_concentration:", L2_mass(List_of_Triangles, number_of_triangles, 3)
      
        
  end do
  
  print *, "Time stepping: done"
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!          Ploting solution                            !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  
  deallocate(prom_m)
  deallocate(solution_anderson)
  deallocate(residuals_anderson)
  deallocate(g_anderson)
  deallocate(f_anderson)
  deallocate(init_resid)
  deallocate(new_sol)
  deallocate(Big_G_anderson)
  deallocate(Big_F_anderson)
  deallocate(rhs_Anderson)
  deallocate(work)
  deallocate(new_residual)
  deallocate(new_g)
  deallocate(gamma_coefficients)
  deallocate(residual_pr)
      
  
end program Main
