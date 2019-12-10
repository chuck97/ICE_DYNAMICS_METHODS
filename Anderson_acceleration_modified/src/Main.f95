program Main
  
  use module_Classes
  use module_Grid
  use module_numerical_integration
  use module_flux_correction
  use module_constants
  use module_assembling
  use module_dynamics
  use module_matricies
  
  use json_module
      
  
  implicit none
  
  
  real*8 :: time_step, total_time, time, squ, max_u ! time 
  
  integer :: i, j, k, num, s, r, ps
  
  real*8, allocatable                                  :: solution_anderson(:)
  real*8                                               :: previous_res
    

  

!!!! BCG PARAMETERS


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
   
   logical                    :: is_found
   type(json_file)            :: json
   real*8                     :: js1, js2, js3, js4
   real*8       , allocatable :: js5(:)
   character(len=32)          :: path_to_json
    
   integer              :: Nbv, Nbl
   real*8, allocatable  :: bv(:, :)
   integer, allocatable :: bl(:, :)
   real*8               :: bltail(2,12)
   
   real*8, allocatable  :: bv_in(:)
   integer, allocatable :: bl_in(:)
   real*8, allocatable  :: bltail_in(:)   
    
   
  time_step = hour/4d0
  time = 0d0  
  num = 0
  
  
  !!! read grid information from json file
  
  call json%initialize()
  
  call get_command_argument(1, path_to_json)
  
  call json%load(filename = trim(path_to_json))
    
  call json%get('Nbv', Nbv, is_found); if (.not. is_found) stop 1
  call json%get('Nbl', Nbl, is_found); if (.not. is_found) stop 1
  allocate(bv_in(2*Nbv))
  allocate(bl_in(7*Nbl))
  allocate(bv(2,Nbv))
  allocate(bl(7,Nbl))
  call json%get('bv',  bv_in, is_found); if (.not. is_found) stop 1
  call json%get('bl',  bl_in, is_found); if (.not. is_found) stop 1
  call json%get('bltail',  bltail_in, is_found); if (.not. is_found) stop 1
  
  do i = 1, Nbv
  
    bv(:,i) = bv_in((i-1)*2 +1 : i*2)
  
  end do
  
  do i  = 1, Nbl
  
    
  
  end do
    
 ! bl = bl_in
 ! bltail = bltail_in
  
  call json%destroy()
  
  call initialization(Nbv, Nbl, square_size*bv, bl, bltail)
  
  deallocate(bv)
  deallocate(bl)
  deallocate(bv_in)
  deallocate(bl_in)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  allocate(prom_m(number_of_elements))
  
  
  
  allocate(solution_anderson(number_of_non_Direchlet_elements*2))
  allocate(residuals_anderson(m_Anderson + 1, number_of_non_Direchlet_elements*2))
  allocate(g_anderson(m_Anderson + 1, number_of_non_Direchlet_elements*2))
  allocate(f_anderson(m_Anderson + 1, number_of_non_Direchlet_elements*2))
  allocate(init_resid(number_of_non_Direchlet_elements*2))
  allocate(new_sol(number_of_non_Direchlet_elements*2))
  allocate(Big_G_anderson(m_Anderson, number_of_non_Direchlet_elements*2))
  allocate(Big_F_anderson(m_Anderson, number_of_non_Direchlet_elements*2))
  allocate(rhs_Anderson(2*number_of_non_Direchlet_elements))
  allocate(work(2*number_of_non_Direchlet_elements*m_Anderson))
  allocate(new_residual(2*number_of_non_Direchlet_elements))
  allocate(new_g(2*number_of_non_Direchlet_elements))
  allocate(gamma_coefficients(m_Anderson + 1))
  allocate(residual_pr(num_anderson_iterations*20))
  
  solution_anderson = 0d0
  residuals_anderson = 0d0
  g_anderson = 0d0
  init_resid = 0d0
  lwork = 2*number_of_non_Direchlet_elements
  
  
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
      
      num_iter = 0 
      
      init_resid = init_residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, &
        List_of_Triangles, number_of_non_Direchlet_elements, number_of_triangles, time_step, 0)
     
      !!!!!!!!!!!!!!!! update first m_Anderson+1 f and g using implicit VP Picard iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
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
          
         call lhs_VP%init(lhs_ia(number_of_non_Direchlet_elements+1)-1, number_of_non_Direchlet_elements, &
          number_of_non_Direchlet_elements)
          lhs_VP%ia(1:(lhs_VP%matr_size_x+1)) = lhs_ia(1:(lhs_VP%matr_size_x+1))
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
           
            solution_anderson(i) = solution_vec_1(i)
            solution_anderson(number_of_non_Direchlet_elements + i) = solution_vec_2(i)
            
            List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = solution_vec_1(i)
            List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = solution_vec_2(i)
            
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = solution_vec_1(i)
            List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = solution_vec_2(i)
            
             
          end do   
          
          !!!!!!!!!!! update f and g  !!!!!!!!!!!!!!!!
          
          
          if (npic .ge. m_Anderson) then
        
          do i = 1, m_Anderson
          
            residuals_anderson(i, :) = residuals_anderson(i+1, :)
            g_anderson(i, :) =  g_anderson(i+1, :) 
          
          end do
          
          
          residuals_anderson(m_Anderson + 1, :) = residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
              number_of_non_Direchlet_elements, number_of_triangles, time_step, 0)
              
          g_anderson(m_Anderson + 1, :) = solution_anderson - &
          alpha_Anderson*(1d0/L2_norm(init_resid, 2*number_of_non_Direchlet_elements))*residuals_anderson(m_Anderson + 1, :)   
          
        
        else
        
          residuals_anderson(npic, :) = residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
              number_of_non_Direchlet_elements, number_of_triangles, time_step, 0)
              
          g_anderson(npic, :) =  solution_anderson - &
          alpha_Anderson*(1d0/L2_norm(init_resid, 2*number_of_non_Direchlet_elements))*residuals_anderson(npic, :)
        
        end if
                  
          print *, npic, L2_norm(residuals_anderson(npic, :), 2*number_of_non_Direchlet_elements)/  &
          L2_norm(init_resid, 2*number_of_non_Direchlet_elements)   
            
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
            
            
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
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
   
      !!!!!     Anderson iterations !!!!!
      
      num_iter = 1
      
      curr_res = L2_norm(f_anderson(m_Anderson + 1, :), 2*number_of_non_Direchlet_elements)
      previous_res = L2_norm(f_anderson(m_Anderson, :), 2*number_of_non_Direchlet_elements)
   
      do while(num_iter < num_anderson_iterations)!((curr_res/previous_res) > 1d-3) .and. (num_iter < 500)) !num_iter < num_anderson_iterations)
        
        previous_res = curr_res
        
        m_k = m_Anderson!min(m_Anderson, num_iter)
        
        !! form matrix Big_G and Big_F
        
        do i = 1, m_k
        
          Big_G_anderson(i, :) =  g_anderson(i + 1, :) - g_anderson(i, :)
          Big_F_anderson(i, :) =  f_anderson(i + 1, :) - f_anderson(i, :)
        
        end do
        
        !! set rhs equals to f_k
        
        rhs_Anderson = residuals_anderson(m_k+1, :)
        
        !! Solve new least square problem
    
        !call dgels(	"N", &
        !   2*number_of_non_Direchlet_elements, &
        !   m_k, &
        !   1, &
        !   Big_F_anderson, &
        !   2*number_of_non_Direchlet_elements, &
        !   rhs_Anderson, &
        !   2*number_of_non_Direchlet_elements, &
        !   work, &
        !   lwork, &
        !   info)	
    
         !! Gamma coefficients calculation
         
         do i = 1, m_k
    
            gamma_coefficients(i) = rhs_Anderson(i)
    
         end do
        
        !! Update new solution
        
        solution_anderson = g_anderson(m_k+1, :)
        
        do i = 1, m_k
        
          solution_anderson = solution_anderson - gamma_coefficients(i)*Big_G_anderson(i, :)
        
        end do 
        
        !! claculate new f and g
        
        do i = 1, number_of_non_Direchlet_elements
      
              List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) = solution_anderson(i)
              List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) = &
               solution_anderson(i + number_of_non_Direchlet_elements)
      
        end do
        
        new_residual = residual(ND_mass_matrix_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
              number_of_non_Direchlet_elements, number_of_triangles, time_step, 0)
              
        new_g = solution_anderson - alpha_Anderson* &
        (1d0/(L2_norm(init_resid, 2*number_of_non_Direchlet_elements)))*new_residual
        
        !! Update f and g matrix
        
        !if (num_iter .ge. m_Anderson) then
        
          do i = 1, m_Anderson
          
            residuals_anderson(i, :) = residuals_anderson(i+1, :)
            g_anderson(i, :) =  g_anderson(i+1, :) 
          
          end do
          
          
          
        residuals_anderson(m_Anderson + 1, :) = new_residual
              
        g_anderson(m_Anderson + 1, :)  = new_g
        
        !! print residual 
        
        residual_pr(num_iter) = L2_norm(new_residual, 2*number_of_non_Direchlet_elements) &
         /L2_norm(init_resid, 2*number_of_non_Direchlet_elements)
         
        curr_res = L2_norm(new_residual, 2*number_of_non_Direchlet_elements)
        
        print *, "iter num", num_iter, "residual:", residual_pr(num_iter)
        !print *, "curr/prev", num_iter, curr_res/previous_res
         
         
          
        num_iter = num_iter + 1  
        
      end do
      
      
      call print_resudal(residual_pr, num_iter-1, num)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      
      do i = 1, number_of_non_Direchlet_elements
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u(1) = &
         solution_anderson(i)
         
        List_of_non_Direchlet_Elements(i)%pointing_element%u(2) = &
         solution_anderson(number_of_non_Direchlet_elements + i)
      
      end do
      
      print *, "max_u:", maximum_u(List_of_Elements, number_of_elements)
      
      call triangles_values_calculations(List_of_Triangles, number_of_triangles)
      
      
      
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
