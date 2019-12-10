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
    
  integer :: i, j, k, npic, nvp, num, s, r, ps, number_of_triangles, number_of_elements, number_of_edges, &
  number_of_boundary_edges, number_of_non_Direchlet_elements, ierr
  
  real*8 :: lhs_a(maxnz)
  integer :: lhs_ia(maxn), lhs_ja(maxnz)
  
  

  type(Element), target, dimension(nvmax) :: List_of_Elements
  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements                     
  type(Edge), target, dimension((ntmax*3)/2 + nbmax) :: List_of_Edges                                     
  type(Triangle), target, dimension(ntmax) :: List_of_Triangles
  type(Matrix) :: Matr1, Matr2, Matr3, Matr4, Matr_res
  real*8 :: rhs_VP_1(nvmax), rhs_VP_2(nvmax)
  real*8 :: res_array(500)
    

  

!!!! BCG PARAMETERS


! Maximum size of matrix and the maximum number of non-zero entries
   integer :: ide, nonzero_raw_L, nonzero_L, nonzero_raw_C, nonzero_C 

! Arrays for the matrix sparse row (CSR) format
   real*8 :: Right_hand_C(maxn), Right_hand_L(maxn), &
     solution_vec_L(maxn), solution_vec_C(maxn), solution_vec(maxn), sum_mass, sum_concentration, &
     solution_vec_1(maxn), solution_vec_2(maxn)
   type(Matrix) :: mass_matrix_high, mass_matrix_low, ND_mass_matrix, ND_mass_matrix_lum, &
    N_xx, N_xy, N_yx, N_yy, Mass_VP, B_VP, lhs_VP

! Work arrays for ILU factors and 8 BCG vectors
  
   integer, parameter :: MaxWr=maxnz+10*maxn, MaxWi=maxnz+7*maxn+1
   real*8 :: rW(MaxWr)
   integer :: iW(MaxWi)   
     
! ILU data
   real*8 :: tau1,tau2,partlur,partlurout
   integer :: verb, UsedWr, UsedWi

! Local variables
   integer :: ipBCG
   integer :: imatvec_1(1), imatvec_2(1), iprevec_1(1), iprevec_2(1)
   real*8 :: resinit_1, resinit_2, prom, max_C_norm, L2_error      
   
! BiCGStab data
   External :: matvec, prevec2
   Integer :: ITER_1, ITER_2, INFO_1, INFO_2, NUNIT_1, NUNIT_2
   Real*8 :: RESID_1, RESID_2

! External routines from the BLAS library
   real*8 :: ddot
   external :: ddot, dcopy  

! Local variables
   integer :: ipBCG_L, ipBCG_C
   integer :: imatvec_L(1), iprevec_L(1), imatvec_C(1), iprevec_C(1) 
   real*8 :: resinit_L, resinit_C, prom_m(nvmax), prom_A(nvmax), x_coord, y_coord
   real*8 :: u_new(2,nvmax), L2_err, L2_resud, resudal_mass(4000)
   integer :: num_iter, resudal_mass_index
   
  time_step = hour/4d0
  time = 0d0  
  num = 1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Making list of Direchlet elements and initialization  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    
  call initialization(List_of_Triangles, List_of_Elements, List_of_Edges, &
     number_of_triangles, number_of_elements, number_of_edges, h_grid)
       
  print *, "Elements, triangles and edges initialization: done"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !! Making list of non Direchlet elements !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  k = 1
   
  do i = 1, number_of_elements
   
    if (List_of_Elements(i)%on_boundary .eqv. .false.) then
     
      List_of_non_Direchlet_Elements(k)%pointing_element => List_of_Elements(i)
      List_of_non_Direchlet_Elements(k)%pointing_element%nd_identificator = k
      
      k = k + 1
     
    end if
   
  end do
  
  
   
  number_of_non_Direchlet_elements = k - 1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!    Assembling non-Direchlet (phi, phi) matrix and lumped matrix     !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call non_Direchlet_mass_matrix_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
   ND_mass_matrix)
   
  call ND_mass_matrix_lum%init(ND_mass_matrix%nonzero, ND_mass_matrix%matr_size)  
  call Matr1_eqv_Matr2(ND_mass_matrix_lum, ND_mass_matrix) 
  
  call lumping(ND_mass_matrix_lum)
  
  print *, "Assembling non-Direchlet and lumped matrix in CSR: done"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                    Assembling N_ii matricies                        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !call N_xx_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
  ! N_xx)
   
  !call N_xy_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
  ! N_xy) 
   
  !call N_yx_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
  ! N_yx) 
   
  !call N_yy_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
  ! N_yy) 
  
  
  
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
  
  do while (num < num_time_steps)
  
      
      do i = 1, number_of_non_Direchlet_elements
      
          List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = &
           List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
           
          List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = &
           List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
          
        
      end do
  
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
         
      
      end do
      
      
      do npic = 1, N_Picard
      
          call velocity_delta_str_P_recalculation(List_of_non_Direchlet_Elements, List_of_Triangles, &
            number_of_non_Direchlet_elements, number_of_triangles, npic)  
            
          call VP_Mass_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
            Mass_VP, ND_mass_matrix_lum, time_step)  
            
          call VP_B_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
            B_VP)  
           
          call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, &
          1, Mass_VP%a, Mass_VP%ja, Mass_VP%ia, B_VP%a , B_VP%ja, B_VP%ia, lhs_a, lhs_ja, lhs_ia, &
          maxnz, iw, ierr)
        
          
          call lhs_VP%init(lhs_ia(number_of_non_Direchlet_elements+1)-1, number_of_non_Direchlet_elements)
          lhs_VP%ia(1:(lhs_VP%matr_size+1)) = lhs_ia(1:(lhs_VP%matr_size+1))
          lhs_VP%ja(1:lhs_VP%nonzero) = lhs_ja(1:lhs_VP%nonzero)
          lhs_VP%a(1:lhs_VP%nonzero) = lhs_a(1:lhs_VP%nonzero) 
           
          call VP_rhs_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, rhs_VP_1, &
              rhs_VP_2, ND_mass_matrix_lum, time_step)
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!  Initialization of the preconditioner  !!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          verb = 0 ! verbose no
          tau1 = 1d-2
          tau2 = 1d-3
          partlur = 0.5
          ierr = 0

          call iluoo(number_of_non_Direchlet_elements, lhs_VP%ia, lhs_VP%ja, lhs_VP%a, tau1, tau2, verb, &
                rW, iW, MaxWr, MaxWi, partlur, partlurout, &
                UsedWr, UsedWi, ierr)
           
          ipBCG = UsedWr + 1 
           
          !set initial guess to 0
          call dcopy(number_of_non_Direchlet_elements,0d0,0,solution_vec_1,1)
          call dcopy(number_of_non_Direchlet_elements,0d0,0,solution_vec_2,1)

          !compute initial residual norm for 1 system
          resinit_1 = ddot(number_of_non_Direchlet_elements,rhs_VP_1,1,rhs_VP_1,1)
          resinit_1 = dsqrt(resinit_1)
           
           
          !compute initial residual norm for 2 system
          resinit_2 = ddot(number_of_non_Direchlet_elements,rhs_VP_2,1,rhs_VP_2,1)
          resinit_2 = dsqrt(resinit_2)
           
           
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!     Iterative solution    !!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
          
          ITER_1 = 1000             
          RESID_1 = 1d-20*resinit_1  
          INFO_1  = 0               
          NUNIT_1 = 6               
          iprevec_1(1) = number_of_non_Direchlet_elements           
          imatvec_1(1) = number_of_non_Direchlet_elements  
           
          ITER_2 = 1000             
          RESID_2 = 1d-20*resinit_2  
          INFO_2  = 0               
          NUNIT_2 = 6               
          iprevec_2(1) = number_of_non_Direchlet_elements           
          imatvec_2(1) = number_of_non_Direchlet_elements             
          
          call slpbcgs( &
              prevec2, iprevec_1, iW,rW, &
              matvec,  imatvec_1, lhs_VP%ia, lhs_VP%ja, lhs_VP%a, &
              rW(ipBCG), number_of_non_Direchlet_elements, 8, &
              number_of_non_Direchlet_elements, rhs_VP_1, solution_vec_1, &
              ITER_1, RESID_1, &
              INFO_1, NUNIT_1)
                    
          call slpbcgs( &
             prevec2, iprevec_2, iW,rW, &
             matvec,  imatvec_2, lhs_VP%ia, lhs_VP%ja, lhs_VP%a, &
             rW(ipBCG), number_of_non_Direchlet_elements, 8, &
             number_of_non_Direchlet_elements, rhs_VP_2, solution_vec_2, &
             ITER_2, RESID_2, &
             INFO_2, NUNIT_2)    
              
          do i = 1, number_of_non_Direchlet_elements
           
            List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = solution_vec_1(i)
            List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = solution_vec_2(i)
             
          end do   
          
          ! resuidal calculation !!!!!!!!!!!!!!!!!!!!!
          
          do i = 1, number_of_non_Direchlet_elements
          
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = &
             List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = &
             List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
          
          end do
          
          res_array(npic) =  resuidal(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
            number_of_non_Direchlet_elements, number_of_triangles, time_step)
          
          print *, "iteration", npic,"residual:",  res_array(npic)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!
                
          call B_VP%distr()
          call Mass_VP%distr()
          call lhs_VP%distr()
          
          ipBCG = 0
          UsedWr = 0
          
          do i = 1, number_of_non_Direchlet_elements
          
            rhs_VP_1(i) = 0d0
            rhs_VP_2(i) = 0d0
          
          end do
          
          do i = 1, MaxWr
          
            rW(i) = 0d0    
          
          end do
          
          do i = 1, MaxWi
          
            iW(i) = 0
          
          end do
         
      end do
      
      call velocity_delta_str_P_recalculation(List_of_non_Direchlet_Elements, List_of_Triangles, &
            number_of_non_Direchlet_elements, number_of_triangles, npic)
      
      call VP_coriolis_update(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, time_step)
      
      do i = 1, number_of_non_Direchlet_elements
          
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = &
             List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = &
             List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
          
      end do
          
      res_array(npic + 1) =  resuidal(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
            number_of_non_Direchlet_elements, number_of_triangles, time_step)
          
      print *, "final iteration residual:",  res_array(npic+1)
      
      call print_resudal(res_array, num, N_Picard)
      
      
      
      call triangles_values_calculations(List_of_Triangles, number_of_triangles)
      
      print *, "max_u:", maximum_u(List_of_Elements, number_of_elements)
      
      
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
      
        prom_m = sparce_matrix_vector_multiplication(mass_matrix_high, solution_vec_C)
      
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
      
        prom_A = sparce_matrix_vector_multiplication(mass_matrix_high, solution_vec_C)
      
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
      
  
end program Main
