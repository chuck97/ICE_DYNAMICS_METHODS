program Main
  
  use module_Classes
  use module_Grid
  use module_numerical_integration
  use module_flux_correction
  use module_constants
  use module_assembling
  use module_dynamics
  use module_matricies
      
  
  implicit none
  
  
  real*8 :: time_step, total_time, time, squ, max_u ! time 
    
  integer :: i, j, k, num, s, r, ps, number_of_triangles, number_of_elements, number_of_edges, &
  number_of_boundary_edges, number_of_non_Direchlet_elements
  
  

  type(Element), target, dimension(nvmax) :: List_of_Elements
  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements                     
  type(Edge), target, dimension((ntmax*3)/2 + nbmax) :: List_of_Edges                                     
  type(Triangle), target, dimension(ntmax) :: List_of_Triangles
  type(Matrix) :: Matr1, Matr2, Matr3, Matr4, Matr_res
    

  

!!!! BCG PARAMETERS


! Maximum size of matrix and the maximum number of non-zero entries
   integer :: ide, nonzero_raw_L, nonzero_L, nonzero_raw_C, nonzero_C 

! Arrays for the matrix sparse row (CSR) format
   real*8 :: Right_hand_C(maxn), Right_hand_L(maxn), &
     solution_vec_L(maxn), solution_vec_C(maxn), solution_vec(maxn), sum_mass, sum_concentration
   type(Matrix) :: mass_matrix_high, mass_matrix_low, ND_mass_matrix, ND_mass_matrix_lum, &
    N_xx, N_xy, N_yx, N_yy 

! Work arrays for ILU factors and 8 BCG vectors
  
   integer, parameter :: MaxWr=maxnz+8*maxn, MaxWi=maxnz+2*maxn+1
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
   real*8, allocatable            :: residual_array(:)
   real*8, allocatable            :: init_resid(:)
   
! Allocating all variables
   
  allocate(residual_array(N_evp))     
  
   
  time_step = hour/4d0
  time = 0d0  
  num = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Making list of Direchlet elements and initialization  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    
  call initialization(List_of_Triangles, List_of_Elements, List_of_Edges, &
     number_of_triangles, number_of_elements, number_of_edges, h_grid)
       
  print *, "Elements, triangles and edges initialization: done"
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !! Making list of non Direchlet elements !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  allocate(prom_m(number_of_elements))
  
  k = 1
   
  do i = 1, number_of_elements
   
    if(List_of_Elements(i)%on_boundary .eqv. .false.) then
     
      List_of_non_Direchlet_Elements(k)%pointing_element => List_of_Elements(i)
      List_of_non_Direchlet_Elements(k)%pointing_element%nd_identificator = k
      
      k = k + 1
     
    end if
   
  end do
   
  number_of_non_Direchlet_elements = k - 1
  
  allocate(init_resid(2*number_of_non_Direchlet_elements))
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!    Assembling non-Direchlet (phi, phi) matrix and lumped matrix     !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call non_Direchlet_mass_matrix_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
   ND_mass_matrix)
   
  call ND_mass_matrix_lum%init(ND_mass_matrix%nonzero, ND_mass_matrix%matr_size_x, ND_mass_matrix%matr_size_y)  
  call Matr1_eqv_Matr2(ND_mass_matrix_lum, ND_mass_matrix) 
  
  call lumping_square(ND_mass_matrix_lum)
  
  print *, "Assembling non-Direchlet and lumped matrix in CSR: done"
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!       Assembling  M_L in CSR format for transport equation           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  call mass_matrix_low_order_assembling(mass_matrix_low, List_of_Elements, &
   number_of_elements, time_step) 
   
  
  print *, "M_L assembling in CSR: done"
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!       Assembling  M_C in CSR format                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call mass_matrix_high_order_assembling(mass_matrix_high, List_of_Elements, &
   number_of_elements, time_step)
    
  
  print *, "M_C assembling in CSR: done"
  
  !!!!!!!!!!!!!!!!!!!
  !! Time stepping !!
  !!!!!!!!!!!!!!!!!!!
  
  do i = 1, number_of_triangles
  
    List_of_Triangles(i)%sigma1 = 0d0
    List_of_Triangles(i)%sigma2 = 0d0
    List_of_Triangles(i)%sigma12 = 0d0
  
  end do
  
  
  do while (num < num_time_steps)
  
      !! wind and water 
  
      do i = 1, number_of_non_Direchlet_elements
      
        x_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(1)
        y_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(2)
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1) = 5d0 + &
         (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*x_coord/L)*dsin(pi*y_coord/L)
        List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2) = 5d0 + &
         (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*y_coord/L)*dsin(pi*x_coord/L)
        List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1) = &
         (2d0*y_coord - 1d6)/(10d0*L)
        List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2) = &
         -(2d0*x_coord - 1d6)/(10d0*L)  
         
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
      
      call plot_ParaView(List_of_Elements, List_of_Triangles, number_of_elements, number_of_triangles, num)
      call write_nodals(List_of_Elements, number_of_elements, num)
      call write_triangles(List_of_Triangles, number_of_triangles, num)
      
      
      print *, "iteration: ", num
      print *, "sum_mass:", L2_mass(List_of_Triangles, number_of_triangles, 1)
      print *, "sum_thickness:", L2_mass(List_of_Triangles, number_of_triangles, 2)
      print *, "sum_concentration:", L2_mass(List_of_Triangles, number_of_triangles, 3)
      
        
  end do
  
  print *, "Time stepping: done"
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!          Ploting solution                            !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  call plot_GNU(List_of_Elements, number_of_elements)
  
  deallocate(prom_m)
  deallocate(residual_array)
  deallocate(init_resid)
      
  
end program Main
