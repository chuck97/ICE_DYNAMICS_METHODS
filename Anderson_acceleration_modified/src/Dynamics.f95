module module_dynamics

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_matricies
  use module_grid_values
  use module_assembling
  
  implicit none
  
  private
  public :: P_0_recalculation, &
            dot_epsilon_delta_recalculation, &
            velocity_recalculation, &
            maximum_u, &
            L2_error, &
            residual, &
            N_assembling_resuid, &
            scalar_mult_P_dphi_resuid, &
            scalar_mult_resuid, &
            L2_norm, &
            init_residual, &
            triangles_values_calculations, &
            velocity_delta_str_P_recalculation, &
            VP_Mass_assembling, &
            VP_B_assembling, &
            scalar_mult_B, &
            VP_rhs_assembling, &
            scalar_mult_P_dphi
  
  contains
  
!! dot epsilon and delta recalculation  
  
subroutine dot_epsilon_delta_recalculation(ind)

  implicit none
  
  integer                                 , intent(in)    :: ind
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
  
  
!! P_0 recalculation

subroutine P_0_recalculation()

  implicit none
  
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


! velocity recalculation according dynamics equation

subroutine velocity_recalculation(time_step)
 
  implicit none

  real*8                                       , intent(in)    :: time_step
  integer                                                      :: i, j, k
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
    L(i) = (beta_EVP + 1) + (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)
    Ee = 2d0*time_step*omega_e
  
  end do
  
  
  do i = 1, number_of_non_Direchlet_elements
  
    rhs_1(i) = beta_EVP*u_old_1(i) + u_n_1(i) + (time_step/m(i))*C_a*a(i)*rho_air*u_a_abs(i)*u_a_1(i) + &
      (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)*u_w_1(i) - &
      (time_step/m(i))*(1d0/non_Direchlet_mass_matrix_lumped%a(i))*sigma_grad_phi_scalar_multiplication( &
      List_of_non_Direchlet_Elements(i)%pointing_element, 1) !- Ee*u_w_2(i)
  
  end do
  
  do i = 1, number_of_non_Direchlet_elements
  
    rhs_2(i) = beta_EVP*u_old_2(i) + u_n_2(i) + (time_step/m(i))*C_a*a(i)*rho_air*u_a_abs(i)*u_a_2(i) + &
      (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)*u_w_2(i) - &
      (time_step/m(i))*(1d0/non_Direchlet_mass_matrix_lumped%a(i))*sigma_grad_phi_scalar_multiplication( &
      List_of_non_Direchlet_Elements(i)%pointing_element, 2) !+ Ee*u_w_1(i)  
  
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


! maximum value of u 

function maximum_u() result(res_max)

  implicit none
 
  integer                                              ::  i, j, k
  real*8                                               :: res_max, prom
  
  prom = -5d2
  
  do i = 1, number_of_elements
  
    if(dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2) > prom) then
    
      prom = dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2)
      
    end if  
  
  end do
  
  res_max = prom

end function maximum_u


!! L2 error function

function L2_error() result(L2_res)

  implicit none 
  
  integer                                                   :: i, j, k
  real*8                                                    :: u_new(2,nvmax), L2_res, prom1, prom2, u_coefficients(3), &
                                                               v_coefficients(3), u_new_coefficients(3), &
                                                               v_new_coefficients(3), u_old_coefficients(3), &
                                                               v_old_coefficients(3), &
                                                               x1, x2, x3, y1, y2, y3
  
  prom1 = 0d0
  prom2 = 0d0
  
  do i = 1, number_of_Triangles
  
    u_new_coefficients = resud_u_coefficients_calculation(List_of_Triangles(i), 1, 1) 
    v_new_coefficients = resud_u_coefficients_calculation(List_of_Triangles(i), 2, 1)
  
    u_old_coefficients = resud_u_coefficients_calculation(List_of_Triangles(i), 1, 2) 
    v_old_coefficients = resud_u_coefficients_calculation(List_of_Triangles(i), 2, 2)
    
    do j = 1, 3
    
      u_coefficients(j) = u_new_coefficients(j) - u_old_coefficients(j)
      v_coefficients(j) = v_new_coefficients(j) - v_old_coefficients(j)
    
    end do
    
    x1 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
    x2 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
    x3 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
    
    y1 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
    y2 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
    y3 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
      
    prom1 = prom1 + integral_over_triangle(x1, x2, x3, y1, y2, y3, squared_linear_on_triangle, &
     u_coefficients(1), u_coefficients(2), u_coefficients(3), 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0)
    prom2 = prom2 + integral_over_triangle(x1, x2, x3, y1, y2, y3, squared_linear_on_triangle, &
     v_coefficients(1), v_coefficients(2), v_coefficients(3), 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0)
    
  
  end do
  
  L2_res = sqrt(prom1 + prom2)/1d12


end function L2_error


function residual(time_step, standart) result(resuidal_res)   ! standart = 1 -- unit residual vector
                                                                                                     ! standart = 0 -- not unit residual vector
  implicit none
  
  real*8                                                    :: resuidal_res(2*number_of_non_Direchlet_elements)
  type(Matrix)                                              :: M, M_cor, M_w, M_a, &
                                                               N_xx_11, N_xx_10, N_xy_01, N_xy_1m1, &
                                                               N_yx_1m1, N_yx_01, N_yy_11, N_yy_10
  real*8, allocatable                                       :: u_a(:), v_a(:), u_w(:), &
                                                               v_w(:), p_x(:), p_y(:), &
                                                               res1(:), res2(:), u_val(:), &
                                                               v_val(:), u_n(:), v_n(:)
  integer                                                   :: number_of_non_Direchlet_elements, number_of_triangles, i, j, k
  real*8, allocatable                                       :: xi_tr(:), eta_tr(:)
  real*8                                                    :: x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, &
                                                               du_dx, du_dy, dv_dx, dv_dy, dot_epsilon11, &
                                                               dot_epsilon22, dot_epsilon12, delta, A_avg, h_avg, &
                                                               P_0, u_minis_u, time_step, u_a_loc, v_a_loc, abs_vel_air
  real*8, allocatable                                       :: prom1(:), prom2(:), u_minus_u, prom_res
  real*8                                                    :: init_resuidal, norm
  integer                                                   :: standart
   
   
  allocate(u_a(number_of_non_Direchlet_elements)) 
  allocate(v_a(number_of_non_Direchlet_elements)) 
  allocate(u_w(number_of_non_Direchlet_elements)) 
  allocate(v_w(number_of_non_Direchlet_elements))
  allocate(p_x(number_of_non_Direchlet_elements))  
  allocate(p_y(number_of_non_Direchlet_elements))
  allocate(res1(number_of_non_Direchlet_elements))
  allocate(res2(number_of_non_Direchlet_elements))
  allocate(u_val(number_of_non_Direchlet_elements))
  allocate(v_val(number_of_non_Direchlet_elements))
  allocate(u_n(number_of_non_Direchlet_elements))
  allocate(v_n(number_of_non_Direchlet_elements))
  allocate(xi_tr(number_of_Triangles))
  allocate(eta_tr(number_of_Triangles))
  allocate(prom1(number_of_non_Direchlet_elements))
  allocate(prom2(number_of_non_Direchlet_elements))
 
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_val(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1)
    v_val(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2)
    u_n(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
    v_n(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
  
  end do
  
  do i = 1, number_of_triangles
  
     x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
     x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
     x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
     
     y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
     y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
     y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
     
     u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_resuid(1)
     u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_resuid(1)
     u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_resuid(1)
     
     v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_resuid(2)
     v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_resuid(2)
     v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_resuid(2)
     
     du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     
     dot_epsilon11 = du_dx
     dot_epsilon22 = dv_dy
     dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
     
     delta = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
      (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
      
     A_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%A)
    
     h_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%h)
    
     List_of_Triangles(i)%P_0_resuid = ((h_avg*p_str*dexp(-C*(1d0 - A_avg)))*delta)/ &
      (delta + delta_min) 
      
     List_of_Triangles(i)%xi_resuid = (h_avg*p_str*dexp(-C*(1d0 - A_avg)))/(2d0*(delta + delta_min)) 
       
     List_of_Triangles(i)%eta_resuid = List_of_Triangles(i)%xi_resuid/(e**2)
  
  end do
  
  call  N_assembling_resuid(N_xx_11, 1, 1, 1d0, 1d0) 
  call  N_assembling_resuid(N_xx_10, 1, 1, 1d0, 0d0)
  call  N_assembling_resuid(N_xy_01, 1, 2, 0d0, 1d0)
  call  N_assembling_resuid(N_xy_1m1, 1, 2, 1d0, -1d0)
  call  N_assembling_resuid(N_yx_1m1, 2, 1, 1d0, -1d0)
  call  N_assembling_resuid(N_yx_01, 2, 1, 0d0, 1d0)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = non_Direchlet_mass_matrix_lumped%a(i)/(time_step)
    M%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_w%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    M_w%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*C_w*rho_water*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*u_minus_u
    M_w%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_w%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! M_a
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_a%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    u_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
    abs_vel_air = dsqrt(u_a_loc**2 + v_a_loc**2) 
    M_a%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*C_a*rho_air*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*abs_vel_air
    M_a%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_a%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! M_cor
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_cor%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    M_cor%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*2d0*omega_e
    M_cor%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_cor%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! p_x, p_y
  
  do i = 1, number_of_non_Direchlet_elements
  
    p_x(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 1)/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)
    p_y(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 2)/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)
  
  end do
  
  ! u_w, v_w
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_w(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1)
    v_w(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2)
  
  end do
  
  ! u_a, v_a
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_a(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
  
  end do
  
  !!!!!!!!!!!!!!!!!!!!
  ! res1 calculation !
  !!!!!!!!!!!!!!!!!!!!
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = 0d0
    
  end do
  
  ! + M*u
    
  call amux(number_of_non_Direchlet_elements, u_val, prom1, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do 
  
  ! - M_cor*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! + M_w*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do
  
  ! + N_xx_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_yx_1m1*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_yx_1m1%a, N_yx_1m1%ja, N_yx_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m 
        
  end do
  
  ! - p_x
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - p_x(i)  
        
  end do
  
  ! + N_yy_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_xy_01*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_xy_01%a, N_xy_01%ja, N_xy_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - N_yy_10*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_10%a, N_yy_10%ja, N_yy_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - M*u^n
  
  call amux(number_of_non_Direchlet_elements, u_n, prom1, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! - M_a*u_a
  
  call amux(number_of_non_Direchlet_elements, u_a, prom1, M_a%a, M_a%ja, M_a%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! - M_w*u_w
  
  call amux(number_of_non_Direchlet_elements, u_w, prom1, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! + M_cor*v_w
  
  call amux(number_of_non_Direchlet_elements, v_w, prom1, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !!!!!!!!!!!!!!!!!!!!
  ! res2 calculation !
  !!!!!!!!!!!!!!!!!!!!
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = 0d0
    
  end do
  
  ! + M*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + M_cor*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + M_w*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + N_yy_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_xy_1m1*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_xy_1m1%a, N_xy_1m1%ja, N_xy_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - p_y
   
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - p_y(i)  
        
  end do
  
  
  ! + N_xx_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! + N_yx_01*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_yx_01%a, N_yx_01%ja, N_yx_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - N_xx_10*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_10%a, N_xx_10%ja, N_xx_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! - M*v_n
  
  call amux(number_of_non_Direchlet_elements, v_n, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! - M_a*v_a
  
  call amux(number_of_non_Direchlet_elements, v_a, prom2, M_a%a, M_a%ja, M_a%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
  ! - M_w*v_w
  
  call amux(number_of_non_Direchlet_elements, v_w, prom2, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
  ! - M_cor*u_w
  
  call amux(number_of_non_Direchlet_elements, u_w, prom2, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! res calculation
  
  if (standart == 1) then
  
    do i = 1, number_of_non_Direchlet_elements
  
      resuidal_res(i) = res1(i)
      resuidal_res(number_of_non_Direchlet_elements + i) = res2(i)
  
    end do
    
    norm = L2_norm(resuidal_res, 2*number_of_non_Direchlet_elements)
    
    do i = 1, 2*number_of_non_Direchlet_elements
  
    resuidal_res(i) = resuidal_res(i)/norm
  
    end do
  
  else
  
    do i = 1, number_of_non_Direchlet_elements
  
      resuidal_res(i) = res1(i)
      resuidal_res(number_of_non_Direchlet_elements + i) = res2(i)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%resid_u = res1(i)
      List_of_non_Direchlet_Elements(i)%pointing_element%resid_v = res2(i)
      
  
    end do
  
  end if
  
  
  
  deallocate(u_a) 
  deallocate(v_a) 
  deallocate(u_w) 
  deallocate(v_w)
  deallocate(p_x)  
  deallocate(p_y)
  deallocate(res1)
  deallocate(res2)
  deallocate(u_val)
  deallocate(v_val)
  deallocate(u_n)
  deallocate(v_n)
  deallocate(xi_tr)
  deallocate(eta_tr)
  deallocate(prom1)
  deallocate(prom2)
  
  call M_cor%distr()
  call M_w%distr()
  call M_a%distr()
  call N_xx_11%distr()
  call N_xx_10%distr()
  call N_xy_01%distr()
  call N_xy_1m1%distr()
  call N_yx_1m1%distr()
  call N_yx_01%distr()
  call N_yy_11%distr()
  call N_yy_10%distr()
  
  
end function residual 



! N matrix assembling for resuidal

subroutine N_assembling_resuid(Matr, x_y_fir, x_y_sec, xi_ind, eta_ind)

  implicit none

  type(Matrix), intent(inout)  :: Matr
  integer                      :: x_y_fir, x_y_sec
  real*8                       :: xi_ind, eta_ind
  real*8                       :: a(maxnz), prom
  integer                      :: ia(maxn), ja(maxnz)
  integer                      :: i, j, k, nonzero, nonzero_raw
  type(Element), pointer       :: elem_master, elem_side
  type(Triangle), pointer      :: trian
  
  ia(1) = 1
  k = 1
  
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
  
    nonzero_raw = 0
  
    elem_master => List_of_non_Direchlet_Elements(i)%pointing_element
  
    do j = 1, elem_master%number_of_neighbour_elements    
    
      elem_side => elem_master%neighbour_elements_list(j)%pointing_element
      
      if (elem_side%on_boundary .eqv. .false.) then
      
          prom = scalar_mult_resuid(elem_master, elem_side, x_y_fir, x_y_sec, xi_ind, eta_ind)
           
          if (abs(prom) > varepsilon) then
          
            a(k) = prom
            ja(k) = elem_side%nd_identificator
            k = k + 1
            nonzero_raw = nonzero_raw + 1
            
          end if   
      
      end if
      
    end do
    
    ia(i+1) = ia(i) + nonzero_raw
    nonzero = nonzero + nonzero_raw
   
  end do
  
  call Matr%init(nonzero, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  Matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  Matr%ja(1:nonzero) = ja(1:nonzero)
  Matr%a(1:nonzero) = a(1:nonzero)

end subroutine N_assembling_resuid



function scalar_mult_P_dphi_resuid(elem, d_phi_dxy) result(sc_res)

  implicit none

  type(Element), target   :: elem
  integer                 :: i, j, k, d_phi_dxy
  type(Triangle), pointer :: trian
  real*8                  :: sc_res
  real*8                  :: coeff(3)
  
  sc_res = 0d0
  
  do i = 1, elem%number_of_neighbour_triangles
        
     trian => elem%neighbour_triangles_list(i)%pointing_triangle
     coeff = coefficients_calculation(elem, trian)
     sc_res = sc_res + trian%size_of_triangle*coeff(d_phi_dxy)*trian%P_0_resuid
  
  end do

end function scalar_mult_P_dphi_resuid




function scalar_mult_resuid(elem1, elem2, d_phi1_dxy, d_phi2_dxy, xi_ind, eta_ind) result(B_res)

  implicit none

  type(Element), target, intent(in)   :: elem1, elem2
  integer                             :: i, j, k, d_phi1_dxy, d_phi2_dxy
  type(Triangle), pointer             :: trian
  real*8                              :: B_res
  real*8                              :: coeff1(3), coeff2(3)
  real*8                              :: xi_ind, eta_ind
  
  B_res = 0d0
  
  do i = 1, elem1%number_of_neighbour_triangles
  
    do j = 1, elem2%number_of_neighbour_triangles
    
      if (associated(elem1%neighbour_triangles_list(i)%pointing_triangle, &
        elem2%neighbour_triangles_list(j)%pointing_triangle)) then
        
          trian => elem2%neighbour_triangles_list(j)%pointing_triangle
        
          coeff1 = coefficients_calculation(elem1, trian)
          coeff2 = coefficients_calculation(elem2, trian)
        
          B_res = B_res + trian%size_of_triangle*coeff1(d_phi1_dxy)*coeff2(d_phi2_dxy)* &
           (xi_ind*trian%xi_resuid + eta_ind*trian%eta_resuid)  
        
      end if  
    
    end do
  
  end do

end function scalar_mult_resuid


!! L2 norm of vector

function L2_norm(vec, size_of_vec) result(L2_norm_res)

  implicit none

  real*8, intent(in)         :: vec(:)
  integer, intent(in)        :: size_of_vec
  real*8                     :: L2_norm_res
  integer                    :: i, j, k
  real*8                     :: prom
  
  prom = 0d0
  
  do i = 1, size_of_vec
  
    prom = prom + vec(i)**2
  
  end do
  
  L2_norm_res = dsqrt(prom)

end function L2_norm




function init_residual(time_step, standart) result(resuidal_res)

  implicit none

  real*8                                                    :: resuidal_res(2*number_of_non_Direchlet_elements)
  type(Matrix)                                              :: M, M_cor, M_w, M_a, &
                                                               N_xx_11, N_xx_10, N_xy_01, N_xy_1m1, &
                                                               N_yx_1m1, N_yx_01, N_yy_11, N_yy_10
  real*8, allocatable                                       :: u_a(:), v_a(:), u_w(:), &
                                                               v_w(:), p_x(:), p_y(:), &
                                                               res1(:), res2(:), u_val(:), &
                                                               v_val(:), u_n(:), v_n(:)
  integer                                                   :: i, j, k
  real*8, allocatable                                       :: xi_tr(:), eta_tr(:)
  real*8                                                    :: x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, &
                                                               du_dx, du_dy, dv_dx, dv_dy, dot_epsilon11, &
                                                               dot_epsilon22, dot_epsilon12, delta, A_avg, h_avg, &
                                                               P_0, u_minis_u, time_step, u_a_loc, v_a_loc, abs_vel_air
  real*8, allocatable                                       :: prom1(:), prom2(:), u_minus_u, prom_res
  real*8                                                    :: init_resuidal, norm
  integer                                                   :: standart
   
   
  allocate(u_a(number_of_non_Direchlet_elements)) 
  allocate(v_a(number_of_non_Direchlet_elements)) 
  allocate(u_w(number_of_non_Direchlet_elements)) 
  allocate(v_w(number_of_non_Direchlet_elements))
  allocate(p_x(number_of_non_Direchlet_elements))  
  allocate(p_y(number_of_non_Direchlet_elements))
  allocate(res1(number_of_non_Direchlet_elements))
  allocate(res2(number_of_non_Direchlet_elements))
  allocate(u_val(number_of_non_Direchlet_elements))
  allocate(v_val(number_of_non_Direchlet_elements))
  allocate(u_n(number_of_non_Direchlet_elements))
  allocate(v_n(number_of_non_Direchlet_elements))
  allocate(xi_tr(number_of_Triangles))
  allocate(eta_tr(number_of_Triangles))
  allocate(prom1(number_of_non_Direchlet_elements))
  allocate(prom2(number_of_non_Direchlet_elements))
 
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_val(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
    v_val(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
    u_n(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
    v_n(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
  
  end do
  
  do i = 1, number_of_triangles
  
     x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
     x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
     x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
     
     y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
     y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
     y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
     
     u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_resuid(1)
     u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_resuid(1)
     u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_resuid(1)
     
     v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_resuid(2)
     v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_resuid(2)
     v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_resuid(2)
     
     du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     
     dot_epsilon11 = du_dx
     dot_epsilon22 = dv_dy
     dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
     
     delta = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
      (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
      
     A_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%A)
    
     h_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%h)
    
     List_of_Triangles(i)%P_0_resuid = ((h_avg*p_str*dexp(-C*(1d0 - A_avg)))*delta)/ &
      (delta + delta_min) 
      
     List_of_Triangles(i)%xi_resuid = (h_avg*p_str*dexp(-C*(1d0 - A_avg)))/(2d0*(delta + delta_min)) 
       
     List_of_Triangles(i)%eta_resuid = List_of_Triangles(i)%xi_resuid/(e**2)
  
  end do
  
  call  N_assembling_resuid(N_xx_11, 1, 1, 1d0, 1d0) 
  call  N_assembling_resuid(N_xx_10, 1, 1, 1d0, 0d0)
  call  N_assembling_resuid(N_xy_01, 1, 2, 0d0, 1d0)
  call  N_assembling_resuid(N_xy_1m1, 1, 2, 1d0, -1d0)
  call  N_assembling_resuid(N_yx_1m1, 2, 1, 1d0, -1d0)
  call  N_assembling_resuid(N_yx_01, 2, 1, 0d0, 1d0)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = non_Direchlet_mass_matrix_lumped%a(i)/(time_step)
    M%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_w%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    M_w%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*C_w*rho_water*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*u_minus_u
    M_w%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_w%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! M_a
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_a%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    u_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
    abs_vel_air = dsqrt(u_a_loc**2 + v_a_loc**2) 
    M_a%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*C_a*rho_air*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*abs_vel_air
    M_a%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_a%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! M_cor
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_cor%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    M_cor%a(i) = non_Direchlet_mass_matrix_lumped%a(i)*2d0*omega_e
    M_cor%ia(i+1) = non_Direchlet_mass_matrix_lumped%ia(i+1)
    M_cor%ja(i) = non_Direchlet_mass_matrix_lumped%ja(i)
  
  end do
  
  ! p_x, p_y
  
  do i = 1, number_of_non_Direchlet_elements
  
    p_x(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 1)/ &
     (List_of_non_Direchlet_Elements(i)%pointing_element%m)
    p_y(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 2)/ &
     (List_of_non_Direchlet_Elements(i)%pointing_element%m)
  
  end do
  
  ! u_w, v_w
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_w(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1)
    v_w(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2)
  
  end do
  
  ! u_a, v_a
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_a(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
  
  end do
  
  !!!!!!!!!!!!!!!!!!!!
  ! res1 calculation !
  !!!!!!!!!!!!!!!!!!!!
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = 0d0
    
  end do
  
  ! + M*u
    
  call amux(number_of_non_Direchlet_elements, u_val, prom1, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do 
  
  ! - M_cor*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! + M_w*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do
  
  ! + N_xx_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_yx_1m1*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_yx_1m1%a, N_yx_1m1%ja, N_yx_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m 
        
  end do
  
  ! - p_x
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - p_x(i)  
        
  end do
  
  ! + N_yy_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_xy_01*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_xy_01%a, N_xy_01%ja, N_xy_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - N_yy_10*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_10%a, N_yy_10%ja, N_yy_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - M*u^n
  
  call amux(number_of_non_Direchlet_elements, u_n, prom1, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! - M_a*u_a
  
  call amux(number_of_non_Direchlet_elements, u_a, prom1, M_a%a, M_a%ja, M_a%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! - M_w*u_w
  
  call amux(number_of_non_Direchlet_elements, u_w, prom1, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)  
        
  end do
  
  ! + M_cor*v_w
  
  call amux(number_of_non_Direchlet_elements, v_w, prom1, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)  
        
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !!!!!!!!!!!!!!!!!!!!
  ! res2 calculation !
  !!!!!!!!!!!!!!!!!!!!
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = 0d0
    
  end do
  
  ! + M*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + M_cor*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + M_w*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  ! + N_yy_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! + N_xy_1m1*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_xy_1m1%a, N_xy_1m1%ja, N_xy_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - p_y
   
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - p_y(i)  
        
  end do
  
  
  ! + N_xx_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! + N_yx_01*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_yx_01%a, N_yx_01%ja, N_yx_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  ! - N_xx_10*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_10%a, N_xx_10%ja, N_xx_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! - M*v_n
  
  call amux(number_of_non_Direchlet_elements, v_n, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i) !/List_of_non_Direchlet_Elements(i)%pointing_element%m  
        
  end do
  
  
  ! - M_a*v_a
  
  call amux(number_of_non_Direchlet_elements, v_a, prom2, M_a%a, M_a%ja, M_a%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
  ! - M_w*v_w
  
  call amux(number_of_non_Direchlet_elements, v_w, prom2, M_w%a, M_w%ja, M_w%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
  ! - M_cor*u_w
  
  call amux(number_of_non_Direchlet_elements, u_w, prom2, M_cor%a, M_cor%ja, M_cor%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)  
        
  end do
  
  
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! res calculation
  
  if (standart == 1) then
  
    do i = 1, number_of_non_Direchlet_elements
  
      resuidal_res(i) = res1(i)
      resuidal_res(number_of_non_Direchlet_elements + i) = res2(i)
  
    end do
    
    norm = L2_norm(resuidal_res, 2*number_of_non_Direchlet_elements)
    
    do i = 1, number_of_non_Direchlet_elements
  
      resuidal_res(i) = res1(i)/norm
      resuidal_res(number_of_non_Direchlet_elements + i) = res2(i)/norm
  
    end do
  
  else
  
    do i = 1, number_of_non_Direchlet_elements
  
      resuidal_res(i) = res1(i)
      resuidal_res(number_of_non_Direchlet_elements + i) = res2(i)
  
    end do
  
  end if
  
  deallocate(u_a) 
  deallocate(v_a) 
  deallocate(u_w) 
  deallocate(v_w)
  deallocate(p_x)  
  deallocate(p_y)
  deallocate(res1)
  deallocate(res2)
  deallocate(u_val)
  deallocate(v_val)
  deallocate(u_n)
  deallocate(v_n)
  deallocate(xi_tr)
  deallocate(eta_tr)
  deallocate(prom1)
  deallocate(prom2)
  
  call M_cor%distr()
  call M_w%distr()
  call M_a%distr()
  call N_xx_11%distr()
  call N_xx_10%distr()
  call N_xy_01%distr()
  call N_xy_1m1%distr()
  call N_yx_1m1%distr()
  call N_yx_01%distr()
  call N_yy_11%distr()
  call N_yy_10%distr()
  
  
end function init_residual 

! calculate values on triangles (for comparisson)

subroutine triangles_values_calculations()

  implicit none

  integer                                   :: i, j, k
  real*8                                    :: sigma11, sigma22, sigma12, x0, x1, x2, &
                                               y0, y1, y2, u0, u1, u2, &
                                               v0, v1, v2, &
                                               du_dx, du_dy, dv_dx, dv_dy, &
                                               dot_epsilon11, dot_epsilon22, &
                                               dot_epsilon12, A_avg, h_avg
  
  do i = 1, number_of_triangles
  
     x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
     x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
     x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
     
     y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
     y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
     y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
     
     u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u(1)
     u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u(1)
     u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u(1)
     
     v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u(2)
     v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u(2)
     v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u(2)
      
     
     du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     
     dot_epsilon11 = du_dx
     dot_epsilon22 = dv_dy
     dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
     
     List_of_Triangles(i)%delta = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
      (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
      
     A_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%A)
    
     h_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%h)
    
     List_of_Triangles(i)%P_0 = ((h_avg*p_str*dexp(-C*(1d0 - A_avg)))*List_of_Triangles(i)%delta)/ &
      (List_of_Triangles(i)%delta + delta_min) 
      
     List_of_triangles(i)%xi = (h_avg*p_str*dexp(-C*(1d0 - A_avg)))/(2d0*( &
       List_of_Triangles(i)%delta + delta_min)) 
       
     List_of_triangles(i)%eta = List_of_triangles(i)%xi/(e**2)
     
     sigma11 = 2d0*List_of_triangles(i)%eta*(dot_epsilon11 - 5d-1*(du_dx + dv_dy)) + &
      List_of_triangles(i)%xi*(du_dx + dv_dy) - 5d-1*List_of_Triangles(i)%P_0
     sigma22 = 2d0*List_of_triangles(i)%eta*(dot_epsilon22 - 5d-1*(du_dx + dv_dy)) + &
      List_of_triangles(i)%xi*(du_dx + dv_dy) - 5d-1*List_of_Triangles(i)%P_0
     sigma12 = 2d0*List_of_triangles(i)%eta*dot_epsilon12
     
     List_of_Triangles(i)%sigma1 = sigma11 + sigma22
     List_of_Triangles(i)%sigma2 = sigma11 - sigma22
     List_of_Triangles(i)%sigma12 = sigma12
  
  end do

  

end subroutine triangles_values_calculations






subroutine velocity_delta_str_P_recalculation(npic)
  
  implicit none

  integer :: i, j, k
  real*8 :: x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, du_dx, du_dy, dv_dx, dv_dy, &
    dot_epsilon11, dot_epsilon12, dot_epsilon22, A_avg, h_avg
  integer :: npic  
  
  do i = 1, number_of_non_Direchlet_elements
  
    if (npic == 1) then
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
       List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
       List_of_non_Direchlet_Elements(i)%pointing_element%u(2) 
       
    end if
    
    if (npic == 2) then
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
       5d-1*(List_of_non_Direchlet_Elements(i)%pointing_element%u(1) + &
       List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1))
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
       5d-1*(List_of_non_Direchlet_Elements(i)%pointing_element%u(2) + &
       List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)) 
       
    else
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
       List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
       List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
       
    end if
    
  end do
  
  do i = 1, number_of_Triangles
  
     x0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(1)
     x1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(1)
     x2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(1)
     
     y0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%coordinates(2)
     y1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%coordinates(2)
     y2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%coordinates(2)
     
     u0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_str(1)
     u1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_str(1)
     u2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_str(1)
     
     v0 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%u_str(2)
     v1 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%u_str(2)
     v2 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%u_str(2)
      
     
     du_dx = ((y2 - y0)*(u1 - u0) - (y1 - y0)*(u2 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     du_dy = ((x1 - x0)*(u2 - u0) - (x2 - x0)*(u1 - u0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dx = ((y2 - y0)*(v1 - v0) - (y1 - y0)*(v2 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     dv_dy = ((x1 - x0)*(v2 - v0) - (x2 - x0)*(v1 - v0))/((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0))
     
     dot_epsilon11 = du_dx
     dot_epsilon22 = dv_dy
     dot_epsilon12 = (1d0/2d0)*(du_dy + dv_dx)
     
     List_of_Triangles(i)%delta_str = dsqrt((dot_epsilon11**2 + dot_epsilon22**2)*(1d0 + 1d0/(e**2)) + &
      (4d0/(e**2))*dot_epsilon12**2 + 2d0*dot_epsilon11*dot_epsilon22*(1d0 - 1d0/(e**2)))
      
     A_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%A + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%A)
    
     h_avg = (1d0/3d0)*(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%h + &
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%h)
    
     List_of_Triangles(i)%P_0 = ((h_avg*p_str*dexp(-C*(1d0 - A_avg)))*List_of_Triangles(i)%delta_str)/ &
      (List_of_Triangles(i)%delta_str + delta_min) 
      
     List_of_triangles(i)%xi = (h_avg*p_str*dexp(-C*(1d0 - A_avg)))/(2d0*( &
       List_of_Triangles(i)%delta_str + delta_min)) 
       
     List_of_triangles(i)%eta = List_of_triangles(i)%xi/(e**2)  
      
  end do
  
end subroutine velocity_delta_str_P_recalculation




subroutine VP_Mass_assembling(VP_mass_matrix, time_step)
  
  implicit none
  
  type(Matrix), intent(inout) :: VP_mass_matrix
  real*8, intent(in)          :: time_step
  real*8                      :: mas, conc, u_minus_u
  integer                     :: i
  real*8                      :: a(maxnz)
  integer                     :: ia(maxn), ja(maxnz)
  integer                     :: nonzero
  
    
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
  
    a(i) = non_Direchlet_mass_matrix_lumped%a(i)*(mas/time_step + C_w*conc*rho_water*u_minus_u)
    
  end do
  
  call VP_mass_matrix%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, &
    number_of_non_Direchlet_elements)
  VP_mass_matrix%ia = non_Direchlet_mass_matrix_lumped%ia
  VP_mass_matrix%ja = non_Direchlet_mass_matrix_lumped%ja
  VP_mass_matrix%a(1:number_of_non_Direchlet_elements) = a(1:number_of_non_Direchlet_elements)
  
  
end subroutine VP_Mass_assembling






! Assembling B matrix in CSR format

subroutine VP_B_assembling(B_VP)

  implicit none
  
  integer                                     :: i, j, k, r
  type(Matrix)                                :: B_VP
  type(Element), pointer                      :: elem_master, elem_side
  type(Triangle), pointer                     :: trian 
  real*8                                      :: a(maxnz), prom
  integer                                     :: ia(maxn), ja(maxnz)
  integer                                     :: nonzero, nonzero_raw
  
  ia(1) = 1
  k = 1
  
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
  
    nonzero_raw = 0
  
    elem_master => List_of_non_Direchlet_Elements(i)%pointing_element
  
    do j = 1, elem_master%number_of_neighbour_elements    
    
      elem_side => elem_master%neighbour_elements_list(j)%pointing_element
      
      if (elem_side%on_boundary .eqv. .false.) then
      
          prom = scalar_mult_B(elem_master, elem_side, 1, 1, 1d0, 1d0) + &
           scalar_mult_B(elem_master, elem_side, 2, 2, 1d0, 1d0)
           
          if (abs(prom) > varepsilon) then
          
            a(k) = prom
            ja(k) = elem_side%nd_identificator
            k = k + 1
            nonzero_raw = nonzero_raw + 1
            
          end if  
      
      end if
      
    end do
    
    ia(i+1) = ia(i) + nonzero_raw
    nonzero = nonzero + nonzero_raw
   
  end do
  
  call B_VP%init(nonzero, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  B_VP%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  B_VP%ja(1:nonzero) = ja(1:nonzero)
  B_VP%a(1:nonzero) = a(1:nonzero)

end subroutine VP_B_assembling



function scalar_mult_B(elem1, elem2, d_phi1_dxy, d_phi2_dxy, xi_ind, eta_ind) result(B_res)

  implicit none

  type(Element), target, intent(in)   :: elem1, elem2
  real*8, intent(in)                  :: xi_ind, eta_ind
  integer, intent(in)                 :: d_phi1_dxy, d_phi2_dxy
  type(Triangle), pointer             :: trian
  real*8                              :: B_res
  real*8                              :: coeff1(3), coeff2(3)
  integer                             :: i, j, k
  
  B_res = 0d0
  
  do i = 1, elem1%number_of_neighbour_triangles
  
    do j = 1, elem2%number_of_neighbour_triangles
    
      if (associated(elem1%neighbour_triangles_list(i)%pointing_triangle, &
        elem2%neighbour_triangles_list(j)%pointing_triangle)) then
        
        trian => elem2%neighbour_triangles_list(j)%pointing_triangle
        
        coeff1 = coefficients_calculation(elem1, trian)
        coeff2 = coefficients_calculation(elem2, trian)
        
        B_res = B_res + trian%size_of_triangle*coeff1(d_phi1_dxy)*coeff2(d_phi2_dxy)* &
         (xi_ind*trian%xi + eta_ind*trian%eta)  
        
      end if  
    
    end do
  
  end do

end function scalar_mult_B




subroutine VP_rhs_assembling(rhs_VP_1, rhs_VP_2, time_step)

  implicit none 
  
  integer                  :: i, j, k
  real*8                   :: rhs_VP_1(nvmax), rhs_VP_2(nvmax), mas, conc, u_minus_u, u_n, v_n, &
                              u_str, v_str, u_w, v_w, u_a, v_a, abs_vel_air, time_step
  type(Element), pointer   :: master_elem, side_elem
 
  
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    u_n = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
    v_n = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
    u_str = List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1)
    v_str = List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2)
    u_w = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1)
    v_w = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2)
    u_a = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
    abs_vel_air = dsqrt(u_a**2 + v_a**2) 
  
    rhs_VP_1(i) = non_Direchlet_mass_matrix_lumped%a(i)*(-2d0*mas*omega_e*v_str + (mas/time_step)*u_n + C_w*conc*rho_water*u_w*u_minus_u + &
     C_a*conc*rho_air*u_a*abs_vel_air) !- 2d0*mas*omega_e*v_w)
    rhs_VP_2(i) = non_Direchlet_mass_matrix_lumped%a(i)*(2d0*mas*omega_e*u_str + (mas/time_step)*v_n + C_w*conc*rho_water*v_w*u_minus_u + &
     C_a*conc*rho_air*v_a*abs_vel_air) !+ 2d0*mas*omega_e*u_w)
  
  end do
  
  do i = 1, number_of_non_Direchlet_elements
    
    master_elem => List_of_non_Direchlet_Elements(i)%pointing_element
    
    do j = 1, master_elem%number_of_neighbour_elements
    
      side_elem => master_elem%neighbour_elements_list(j)%pointing_element
      
      rhs_VP_1(i) = rhs_VP_1(i) - side_elem%u_str(2)*scalar_mult_B(side_elem, master_elem, 2, 1, 1d0, -1d0) - &
         side_elem%u_str(2)*scalar_mult_B(side_elem, master_elem, 1, 2, 0d0, 1d0) + &
         side_elem%u_str(1)*scalar_mult_B(side_elem, master_elem, 2, 2, 1d0, 0d0)
          
      rhs_VP_2(i) = rhs_VP_2(i) - side_elem%u_str(1)*scalar_mult_B(side_elem, master_elem, 1, 2, 1d0, -1d0) - &
         side_elem%u_str(1)*scalar_mult_B(side_elem, master_elem, 2, 1, 0d0, 1d0) + &
         side_elem%u_str(2)*scalar_mult_B(side_elem, master_elem, 1, 1, 1d0, 0d0)  
          
    end do 
    
      rhs_VP_1(i) = rhs_VP_1(i) + 5d-1*scalar_mult_P_dphi(master_elem, 1)
      rhs_VP_2(i) = rhs_VP_2(i) + 5d-1*scalar_mult_P_dphi(master_elem, 2)
  
  end do

end subroutine VP_rhs_assembling



function scalar_mult_P_dphi(elem, d_phi_dxy) result(sc_res)

  implicit none

  type(Element), target, intent(in) :: elem
  integer, intent(in)               :: d_phi_dxy
  integer                           :: i, j, k
  type(Triangle), pointer           :: trian
  real*8                            :: sc_res
  real*8                            :: coeff(3)
  
  sc_res = 0d0
  
  do i = 1, elem%number_of_neighbour_triangles
        
     trian => elem%neighbour_triangles_list(i)%pointing_triangle
     coeff = coefficients_calculation(elem, trian)
     sc_res = sc_res + trian%size_of_triangle*coeff(d_phi_dxy)*trian%P_0
  
  end do

end function scalar_mult_P_dphi

end module module_dynamics 
