module module_mEVP_dynamics

  use module_classes
  use module_constants
  use module_grid_values
  use module_assembling
  use module_numerical_integration, only : coefficients_calculation
  
  private :: alpha_mEVP,                            &
             beta_mEVP,                             &
             num_iter_mEVP,                         &
             P_0_recalculation,                     &
             dot_epsilon_delta_recalculation,       &
             velocity_recalculation,                &
             maximum_u,                             &
             sigma_grad_phi_scalar_multiplication
             
  public  :: mEVP_velocity_update
              
  real*8, parameter   :: alpha_mEVP = 5d2
  real*8, parameter   :: beta_mEVP = 5d2
  integer, parameter  :: num_iter_mEVP = 500
  
  contains
  
  subroutine dot_epsilon_delta_recalculation(ind)
  
    implicit none
    
    !! arguments:
    integer                                 , intent(in)    :: ind
    
    !! local variables:
    integer                                                 :: i
    real*8                                                  :: dot_epsilon11, dot_epsilon22, dot_epsilon12
    real*8                                                  :: u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2, du_dx, &
                                                               du_dy, dv_dx, dv_dy
    
    do i = 1, number_of_triangles
    
      if (ind == 1) then
      
       x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
       x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
       x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
       
       y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
       y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
       y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
       
       u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_old(1)
       u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_old(1)
       u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_old(1)
       
       v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_old(2)
       v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_old(2)
       v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_old(2)
       
       du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       
       dot_epsilon11 = du_dx
       dot_epsilon22 = dv_dy
       dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
       
       List_of_Triangles(i)%dot_epsilon1 = dot_epsilon11 + dot_epsilon22
       List_of_Triangles(i)%dot_epsilon2 = dot_epsilon11 - dot_epsilon22
       List_of_Triangles(i)%dot_epsilon12 = dot_epsilon12
       
       List_of_Triangles(i)%delta = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
        (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
          
       
       else
       
       x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
       x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
       x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
       
       y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
       y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
       y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
       
       u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_old(1)
       u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_old(1)
       u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_old(1)
       
       v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_old(2)
       v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_old(2)
       v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_old(2)
       
       du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
       
       dot_epsilon11 = du_dx
       dot_epsilon22 = dv_dy
       dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
       
       List_of_Triangles(i)%dot_epsilon1 = dot_epsilon11 + dot_epsilon22
       List_of_Triangles(i)%dot_epsilon2 = dot_epsilon11 - dot_epsilon22
       List_of_Triangles(i)%dot_epsilon12 = dot_epsilon12
       
       
       List_of_Triangles(i)%delta_old = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
          (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
          
       end if
     
    
    end do
    
  end subroutine dot_epsilon_delta_recalculation
  
  
  
  subroutine P_0_recalculation()
  
    implicit none
    
    !!local variables:
    real*8                                                  :: A_avg, h_avg, delt 
    integer                                                 :: i
    
    do i = 1, number_of_triangles
      A_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%A)
      
      h_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%h)
      
      delt = List_of_Triangles(i)%delta
      
      List_of_Triangles(i)%P_0 = h_avg*p_str*dexp(-C*(1d0 - A_avg))*delt/(delt + delta_min)
    end do
    
  end subroutine P_0_recalculation
  
  subroutine velocity_recalculation(time_step)
 
    implicit none
     
    !! arguments:
    real*8                                       , intent(in)    :: time_step
    
    !! local variables:
    integer                                                      :: i, j
    real*8, allocatable                                          :: m(:), a(:), u_a_1(:), u_a_2(:),     &
                                                                    u_a_abs(:), u_w_1(:), u_w_2(:),     &
                                                                    u_old_1(:), u_old_2(:), u_n_1(:),   &
                                                                    u_n_2(:), u_w_minus_u_old_abs(:),   &
                                                                    L(:), rhs_1(:), rhs_2(:)
    real*8                                                       :: Ee                                                                 
  
    allocate(m(number_of_non_Direchlet_elements))
    allocate(a(number_of_non_Direchlet_elements))
    allocate(u_a_1(number_of_non_Direchlet_elements))
    allocate(u_a_2(number_of_non_Direchlet_elements))
    allocate(u_a_abs(number_of_non_Direchlet_elements))
    allocate(u_w_1(number_of_non_Direchlet_elements))
    allocate(u_w_2(number_of_non_Direchlet_elements))
    allocate(u_old_1(number_of_non_Direchlet_elements))
    allocate(u_old_2(number_of_non_Direchlet_elements))
    allocate(u_n_1(number_of_non_Direchlet_elements))
    allocate(u_n_2(number_of_non_Direchlet_elements))
    allocate(u_w_minus_u_old_abs(number_of_non_Direchlet_elements))
    allocate(L(number_of_non_Direchlet_elements))
    allocate(rhs_1(number_of_non_Direchlet_elements))
    allocate(rhs_2(number_of_non_Direchlet_elements))
  
    do i = 1, number_of_non_Direchlet_elements
    
      m(i) = List_of_non_Direchlet_Elements(i)%pointing_element%m
      a(i) = List_of_non_Direchlet_Elements(i)%pointing_element%A
      u_a_1(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1) 
      u_a_2(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
      u_a_abs(i) = dsqrt(u_a_1(i)**2 + u_a_2(i)**2)
      u_w_1(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1)
      u_w_2(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2)
      u_old_1(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
      u_old_2(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
      u_n_1(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
      u_n_2(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
      u_w_minus_u_old_abs(i) = dsqrt((u_w_1(i) - u_old_1(i))**2 + (u_w_2(i) - u_old_2(i))**2)
      L(i) = (beta_mEVP + 1) + (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)
      Ee = 2d0*time_step*omega_e
    
    end do
    
    
    do i = 1, number_of_non_Direchlet_elements
    
      rhs_1(i) = beta_mEVP*u_old_1(i) + u_n_1(i) + (time_step/m(i))*C_a*a(i)*rho_air*u_a_abs(i)*u_a_1(i) + &
        (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)*u_w_1(i) - &
        (time_step/m(i))*(1d0/non_Direchlet_mass_matrix_lumped%a(i))*sigma_grad_phi_scalar_multiplication( &
        List_of_non_Direchlet_Elements(i)%pointing_element, 1) - Ee*u_w_2(i)
    
    end do
    
    do i = 1, number_of_non_Direchlet_elements
    
      rhs_2(i) = beta_mEVP*u_old_2(i) + u_n_2(i) + (time_step/m(i))*C_a*a(i)*rho_air*u_a_abs(i)*u_a_2(i) + &
        (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)*u_w_2(i) - &
        (time_step/m(i))*(1d0/non_Direchlet_mass_matrix_lumped%a(i))*sigma_grad_phi_scalar_multiplication( &
        List_of_non_Direchlet_Elements(i)%pointing_element, 2) + Ee*u_w_1(i)  
    
    end do
    
    do i = 1, number_of_non_Direchlet_elements
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_new(1) = &
      (L(i)*rhs_1(i) + Ee*rhs_2(i))/(L(i)**2 + Ee**2)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_new(2) = &
      (L(i)*rhs_2(i) - Ee*rhs_1(i))/(L(i)**2 + Ee**2)
    
    end do
    
    deallocate(m)
    deallocate(a)
    deallocate(u_a_1)
    deallocate(u_a_2)
    deallocate(u_a_abs)
    deallocate(u_w_1)
    deallocate(u_w_2)
    deallocate(u_old_1)
    deallocate(u_old_2)
    deallocate(u_n_1)
    deallocate(u_n_2)
    deallocate(u_w_minus_u_old_abs)
    deallocate(L)
    deallocate(rhs_1)
    deallocate(rhs_2)
  
  end subroutine velocity_recalculation
  
  
  function maximum_u() result(res_max)
  
    implicit none
    
    !!local variables:
    integer                                              :: i
    real*8                                               :: res_max, prom
    
    prom = -5d2
    
    do i = 1, number_of_elements
    
      if(dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2) > prom) then
      
        prom = dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2)
        
      end if  
    
    end do
    
    res_max = prom
  
  end function maximum_u
  
  
  
  subroutine mEVP_velocity_update(time_step)
  
    implicit none
    
    !!arguments:
    real*8, intent(in)   :: time_step
    
    !!local variables:
    integer              :: num_iter, i
    
    num_iter = 1
    
    do while((num_iter < (num_iter_mEVP+1)))
    
      call dot_epsilon_delta_recalculation(1)
      call P_0_recalculation()
     
      do i = 1, number_of_Triangles
        List_of_Triangles(i)%sigma1 = &
         (1d0 - (1d0/alpha_mEVP))*List_of_Triangles(i)%sigma1 + &
         (1d0/alpha_mEVP)*(List_of_Triangles(i)%dot_epsilon1 - &
         List_of_Triangles(i)%delta)* &
         (List_of_Triangles(i)%P_0)/ &
         (List_of_Triangles(i)%delta + delta_min)
         
        List_of_Triangles(i)%sigma2 = &
         (1d0 - (1d0/alpha_mEVP))*List_of_Triangles(i)%sigma2 + &
         (1d0/alpha_mEVP)*(List_of_Triangles(i)%dot_epsilon2)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2)
         
        List_of_Triangles(i)%sigma12 = &
         (1d0 - (1d0/alpha_mEVP))*List_of_Triangles(i)%sigma12 + &
         (1d0/alpha_mEVP)*(List_of_Triangles(i)%dot_epsilon12)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2) 
      end do
      
      call velocity_recalculation(time_step)
      
      do i = 1, number_of_non_Direchlet_elements
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(1)
    
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(2)
      end do
      
      num_iter = num_iter + 1
        
    end do
      
      do i = 1, number_of_non_Direchlet_elements
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
         
        List_of_non_Direchlet_Elements(i)%pointing_element%u(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) 
      
      end do
      
      print *, "max_u:", maximum_u()
    
  
  end subroutine mEVP_velocity_update
  
  
  
  function sigma_grad_phi_scalar_multiplication(elem, ind) result(mult_result) ! ind = 1 -- first component, ind = 2 -- second component

    implicit none
    
    !!arguments:
    type(Element), target, intent(in)   :: elem
    integer, intent(in)                 :: ind
    
    !!local variables:
    real*8                              :: sigma_11, sigma_12, sigma_22
    integer                             :: i
    real*8                              :: coefficients(3), mult_result
  
    mult_result = 0d0
  
    if (ind == 1) then
      do i = 1, elem%number_of_neighbour_triangles
        coefficients = coefficients_calculation(elem, &
         elem%neighbour_triangles_list(i)%pointing_triangle)
                
        sigma_11 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 + &
        elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
        sigma_22 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 - &
        elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
        sigma_12 = elem%neighbour_triangles_list(i)%pointing_triangle%sigma12
  
        mult_result = mult_result + &
        elem%neighbour_triangles_list(i)%pointing_triangle%size_of_triangle*( &
        sigma_11*coefficients(1) + sigma_12*coefficients(2))
      end do
      
      else
    
      do i = 1, elem%number_of_neighbour_triangles
        coefficients = coefficients_calculation(elem, &
         elem%neighbour_triangles_list(i)%pointing_triangle)
        sigma_11 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 + &
        elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
        sigma_22 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 - &
        elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
        sigma_12 = elem%neighbour_triangles_list(i)%pointing_triangle%sigma12
  
        mult_result = mult_result + &
        elem%neighbour_triangles_list(i)%pointing_triangle%size_of_triangle*( &
        sigma_12*coefficients(1) + sigma_22*coefficients(2))
      end do
    endif
    
  end function sigma_grad_phi_scalar_multiplication

  
end module module_mEVP_dynamics
