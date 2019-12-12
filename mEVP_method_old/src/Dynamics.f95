module module_dynamics

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_matricies
  
  implicit none
  
  private
  public :: P_0_recalculation, dot_epsilon_delta_recalculation, &
   velocity_recalculation, maximum_u, L2_error, A_matrix_assembling, &
   R_assembling, N_assembling, N_assembling_resuid, scalar_mult_P_dphi_resuid, &
   scalar_mult_resuid, residual, init_residual, L2_norm
  
  contains
  
!! dot epsilon and delta recalculation  
  
subroutine dot_epsilon_delta_recalculation(List_of_Triangles, number_of_triangles, ind)

  implicit none
  
  type(Triangle), target, dimension(ntmax), intent(inout) :: List_of_Triangles 
  integer                                 , intent(in)    :: number_of_triangles, ind
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

subroutine P_0_recalculation(List_of_Triangles, number_of_triangles)

  implicit none
  
  type(Triangle), target, dimension(ntmax), intent(inout) :: List_of_Triangles
  integer                                 , intent(in)    :: number_of_triangles
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

subroutine velocity_recalculation(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 ND_mass_matrix_lum, time_step)

  type(container_element), dimension(nvmax)    , intent(inout) :: List_of_non_Direchlet_Elements
  integer                                      , intent(in)    :: number_of_non_Direchlet_elements
  type(Matrix)                                 , intent(in)    :: ND_mass_matrix_lum
  real*8                                       , intent(in)    :: time_step
  integer                                                      :: i, j, k
  real*8, allocatable                                          :: m(:), a(:), u_a_1(:), u_a_2(:),     &
                                                                  u_a_abs(:), u_w_1(:), u_w_2(:),     &
                                                                  u_old_1(:), u_old_2(:), u_n_1(:),   &
                                                                  u_n_2(:), u_w_minus_u_old_abs(:),      &
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
      (time_step/m(i))*(1d0/ND_mass_matrix_lum%a(i))*sigma_grad_phi_scalar_multiplication( &
      List_of_non_Direchlet_Elements(i)%pointing_element, 1) !- Ee*u_w_2(i)
  
  end do
  
  do i = 1, number_of_non_Direchlet_elements
  
    rhs_2(i) = beta_EVP*u_old_2(i) + u_n_2(i) + (time_step/m(i))*C_a*a(i)*rho_air*u_a_abs(i)*u_a_2(i) + &
      (time_step/m(i))*C_w*a(i)*rho_water*u_w_minus_u_old_abs(i)*u_w_2(i) - &
      (time_step/m(i))*(1d0/ND_mass_matrix_lum%a(i))*sigma_grad_phi_scalar_multiplication( &
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

function maximum_u(List_of_Elements, number_of_elements) result(res_max)

  type(Element), target,  dimension(nvmax), intent(in) :: List_of_Elements
  integer                                 , intent(in) :: number_of_elements  
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

function L2_error(List_of_Triangles, number_of_triangles) result(L2_res)

  type(Triangle), target, dimension(ntmax)      ,intent(in) :: List_of_Triangles
  integer                                       ,intent(in) :: number_of_triangles
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

! Assembling global A matrix (3*nt+2*nv , 3*nt+2*nv)

subroutine A_matrix_assembling(List_of_non_Direchlet_Elements, List_of_Triangles, number_of_non_Direchlet_elements, &
  number_of_triangles, nd_mass_lum, time_step, result_matrix)

  implicit none
  
  type(Triangle), target, dimension(ntmax)      ,intent(in) :: List_of_Triangles
  type(container_element), dimension(nvmax)     ,intent(in) :: List_of_non_Direchlet_Elements
  integer                                       ,intent(in) :: number_of_non_Direchlet_elements, number_of_triangles
  real*8                                        ,intent(in) :: time_step
  type(Matrix)                                  ,intent(in) :: nd_mass_lum
  type(Matrix)                               ,intent(inout) :: result_matrix
  integer                                                   :: i, j, k
  type(Matrix)                                              :: M1, M2, M3, M4, &
                                                               R1, R2, R3, R4, R5, R6, &
                                                               N1, N2, N3, N4, N5, N6, &
                                                               S, Z, R_x, R_y, N_x, N_y
  real*8, allocatable                                       :: mass(:), Cw(:)
  real*8, allocatable                                       :: a_w(:)
  integer, allocatable                                      :: ia_w(:), ja_w(:)
  real*8, allocatable                                       :: work(:) 
  logical                                                   :: values
  
  
  allocate(mass(number_of_non_Direchlet_elements))
  allocate(Cw(number_of_non_Direchlet_elements)) 
  allocate(work(number_of_non_Direchlet_elements))
  allocate(a_w(number_of_non_Direchlet_elements*15))
  allocate(ia_w(number_of_non_Direchlet_elements*2))
  allocate(ja_w(number_of_non_Direchlet_elements*15))
  
  do i = 1, number_of_non_Direchlet_elements
  
    mass(i) = List_of_non_Direchlet_Elements(i)%pointing_element%m
    Cw(i) = List_of_non_Direchlet_Elements(i)%pointing_element%A*C_w*rho_water* &
    dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u(1) - &
     List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
     (List_of_non_Direchlet_Elements(i)%pointing_element%u(2) - &
     List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    work(i) = mass(i) + time_step*Cw(i)
  
  end do     
  
  
  
  !! Assembling M1                                                       
    
  call diamua(number_of_non_Direchlet_elements, 1, nd_mass_lum%a, &
    nd_mass_lum%ja, nd_mass_lum%ia, work, a_w, ja_w, ia_w)
    
  call csort (number_of_non_Direchlet_elements, a_w, ja_w, ia_w, values)   
    
  call M1%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M1%ia(1:(number_of_non_Direchlet_elements+1)) = ia_w(1:(number_of_non_Direchlet_elements+1))                                       
  M1%ja(1:number_of_non_Direchlet_elements) = ja_w(1:number_of_non_Direchlet_elements)
  M1%a(1:number_of_non_Direchlet_elements) = a_w(1:number_of_non_Direchlet_elements)
                                                            
                                                               
  !! Assembling M4                                                             
   
  call M4%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M4%ia(1:(number_of_non_Direchlet_elements+1)) = ia_w(1:(number_of_non_Direchlet_elements+1))                                       
  M4%ja(1:number_of_non_Direchlet_elements) = ja_w(1:number_of_non_Direchlet_elements)
  M4%a(1:number_of_non_Direchlet_elements) = a_w(1:number_of_non_Direchlet_elements)  
  
  ia_w(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    ia_w(i+1) = 1
  
  end do         
  
  do i = 1, number_of_non_Direchlet_elements
  
    ja_w(i) = 0
    a_w(i) = 0d0
  
  end do                                                   
  
  
  !! Assembling M2
  
  do i = 1, number_of_non_Direchlet_elements
  
    work(i) = 2d0*time_step*omega_e*mass(i)
  
  end do
  
  call diamua(number_of_non_Direchlet_elements, 1, nd_mass_lum%a, &
    nd_mass_lum%ja, nd_mass_lum%ia, work, a_w, ja_w, ia_w)
    
  call csort (number_of_non_Direchlet_elements, a_w, ja_w, ia_w, values)  
  
  call M2%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M2%ia(1:(number_of_non_Direchlet_elements+1)) = ia_w(1:(number_of_non_Direchlet_elements+1))
  M2%ja(1:number_of_non_Direchlet_elements) = ja_w(1:number_of_non_Direchlet_elements)
  
  do i = 1, number_of_non_Direchlet_elements
  
    M2%a(i) = -1d0*a_w(i)  
    
  end do
  
  !! Assembling M3
  
  call M3%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M3%ia(1:(number_of_non_Direchlet_elements+1)) = ia_w(1:(number_of_non_Direchlet_elements+1))
  M3%ja(1:number_of_non_Direchlet_elements) = ja_w(1:number_of_non_Direchlet_elements)
  M3%a(1:number_of_non_Direchlet_elements) = a_w(1:number_of_non_Direchlet_elements)
  
  
  ia_w(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    ia_w(i+1) = 1
  
  end do         
  
  do i = 1, number_of_non_Direchlet_elements
  
    ja_w(i) = 0
    a_w(i) = 0d0
  
  end do  
  
  !! Assembling R_x, R_y
  
  call R_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, number_of_triangles, &
   1, R_x)
   
  call R_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, number_of_triangles, &
   2, R_y)
   
  
  !! Assembling R1
  
  call R1%init(R_x%nonzero, R_x%matr_size_x, R_x%matr_size_y)
  
  R1%ia(1:(R1%matr_size_y+1)) = R_x%ia(1:(R_x%matr_size_y+1)) 
  R1%ja(1:R1%nonzero) = R_x%ja(1:R_x%nonzero)
  
  do i = 1, R1%nonzero
  
    R1%a(i) = R_x%a(i)*(5d-1)*time_step
  
  end do
  
  !! Assembling R2
  
  call R2%init(R_x%nonzero, R_x%matr_size_x, R_x%matr_size_y)
  
  R2%ia(1:(R2%matr_size_y+1)) = R_x%ia(1:(R_x%matr_size_y+1)) 
  R2%ja(1:R2%nonzero) = R_x%ja(1:R_x%nonzero)
  
  do i = 1, R2%nonzero
  
    R2%a(i) = R_x%a(i)*(5d-1)*time_step
  
  end do
  
  !! Assembling R3
  
  call R3%init(R_y%nonzero, R_y%matr_size_x, R_y%matr_size_y)
  
  R3%ia(1:(R3%matr_size_y+1)) = R_y%ia(1:(R_x%matr_size_y+1)) 
  R3%ja(1:R3%nonzero) = R_y%ja(1:R_y%nonzero)
  
  do i = 1, R3%nonzero
  
    R3%a(i) = R_y%a(i)*time_step
  
  end do
  
  
  !! Assembling R4
  
  call R4%init(R_y%nonzero, R_y%matr_size_x, R_y%matr_size_y)
  
  R4%ia(1:(R4%matr_size_y+1)) = R_y%ia(1:(R_y%matr_size_y+1))
  R4%ja(1:R4%nonzero) = R_y%ja(1:R_y%nonzero)
  
  do i = 1, R4%nonzero
  
    R4%a(i) = R_y%a(i)*(5d-1)*time_step
  
  end do
  
  
  !! Assembling R5
  
  call R5%init(R_y%nonzero, R_y%matr_size_x, R_y%matr_size_y)
  
  R5%ia(1:(R5%matr_size_y+1)) = R_y%ia(1:(R_y%matr_size_y+1)) 
  R5%ja(1:R5%nonzero) = R_y%ja(1:R_y%nonzero)
  
  do i = 1, R5%nonzero
  
    R5%a(i) = R_y%a(i)*(-5d-1)*time_step
  
  end do
  
  !! Assembling R6
  
  call R6%init(R_x%nonzero, R_x%matr_size_x, R_x%matr_size_y)
  
  R6%ia(1:(R6%matr_size_y+1)) = R_x%ia(1:(R_x%matr_size_y+1)) 
  R6%ja(1:R6%nonzero) = R_x%ja(1:R_x%nonzero)
  
  do i = 1, R6%nonzero
  
    R6%a(i) = R_x%a(i)*time_step
  
  end do
  
  
  !! Assembling N_x, N_y
  
  call N_assembling(List_of_Triangles, number_of_triangles, number_of_non_Direchlet_elements, &
   1, N_x)
   
  call N_assembling(List_of_Triangles, number_of_triangles, number_of_non_Direchlet_elements, &
   2, N_y) 
  
  !! Assembling N1
  
  call N1%init(N_x%nonzero, N_x%matr_size_x, N_x%matr_size_y)
  
  N1%ia(1:(N1%matr_size_y + 1)) = N_x%ia(1:(N_x%matr_size_y + 1))
  N1%ja(1:N1%nonzero) = N_x%ja(1:N_x%nonzero)
  
  do i = 1, N1%nonzero
  
    N1%a(i) = (-1d0)*N_x%a(i)
  
  end do
  
  
  !! Assembling N2
  
  call N2%init(N_y%nonzero, N_y%matr_size_x, N_y%matr_size_y)
  
  N2%ia(1:(N2%matr_size_y + 1)) = N_y%ia(1:(N_y%matr_size_y + 1))
  N2%ja(1:N2%nonzero) = N_y%ja(1:N_y%nonzero)
  
  do i = 1, N2%nonzero
  
    N2%a(i) = (-1d0)*N_y%a(i)
  
  end do
  
  !! Assembling N3
  
  call N3%init(N_x%nonzero, N_x%matr_size_x, N_x%matr_size_y)
  
  N3%ia(1:(N3%matr_size_y + 1)) = N_x%ia(1:(N_x%matr_size_y + 1))
  N3%ja(1:N3%nonzero) = N_x%ja(1:N_x%nonzero)
  
  do i = 1, N3%nonzero
  
    N3%a(i) = (-1d0/(e**2))*N_x%a(i)
  
  end do
  
  !! Assembling N4
  
  call N4%init(N_y%nonzero, N_y%matr_size_x, N_y%matr_size_y)
  
  N4%ia(1:(N4%matr_size_y + 1)) = N_y%ia(1:(N_y%matr_size_y + 1))
  N4%ja(1:N4%nonzero) = N_y%ja(1:N_y%nonzero)
  
  do i = 1, N4%nonzero
  
    N4%a(i) = (1d0/(e**2))*N_y%a(i)
  
  end do
  
  !! Assembling N5
  
  call N5%init(N_y%nonzero, N_y%matr_size_x, N_y%matr_size_y)
  
  N5%ia(1:(N5%matr_size_y + 1)) = N_y%ia(1:(N_y%matr_size_y + 1))
  N5%ja(1:N5%nonzero) = N_y%ja(1:N_y%nonzero)
  
  do i = 1, N5%nonzero
  
    N5%a(i) = (-1d0/(2d0*(e**2)))*N_y%a(i)
  
  end do
  
  !! Assembling N6
  
  call N6%init(N_x%nonzero, N_x%matr_size_x, N_x%matr_size_y)
  
  N6%ia(1:(N6%matr_size_y + 1)) = N_x%ia(1:(N_x%matr_size_y + 1))
  N6%ja(1:N6%nonzero) = N_x%ja(1:N_x%nonzero)
  
  do i = 1, N6%nonzero
  
    N6%a(i) = (-1d0/(2d0*(e**2)))*N_x%a(i)
  
  end do
  
  
  !! Assembling S
  
  call S%init(number_of_triangles, number_of_triangles, number_of_triangles)
  
  S%ia(1) = 1
  
  do i = 1, number_of_triangles
  
    S%ia(i+1) = i+1
  
  end do
  
  do i = 1, number_of_triangles
  
    S%ja(i) = i
  
  end do
  
  do i = 1, number_of_triangles
  
    S%a(i) = 1d0 !List_of_Triangles(i)%size_of_triangle
  
  end do
  
  !! Assembzing Z
  
  call Z%init(0, number_of_triangles, number_of_triangles)
  
  Z%ia(1) = 1
  
  do i = 1, number_of_triangles
  
    Z%ia(i+1) = 1
  
  end do
  
  !! Assembling result matrix
   
  call big_block_matrix_assembling(M1, M2, R1, R2, R3, &
                                         M3, M4, R4, R5, R6, &
                                         N1, N2, S, Z, Z, &
                                         N3, N4, Z, S, Z, &
                                         N5, N6, Z, Z, S, result_matrix, number_of_triangles, &
                                         number_of_non_Direchlet_elements)
                                         
  call M1%distr()
  call M2%distr()
  call M3%distr()
  call M4%distr()
  call R1%distr()
  call R2%distr()
  call R3%distr()
  call R4%distr()
  call R5%distr()
  call R6%distr()
  call N1%distr()
  call N2%distr()
  call N3%distr()
  call N4%distr()
  call N5%distr()
  call N6%distr()
  call S%distr()
  call Z%distr()
  call R_x%distr()
  call R_y%distr()
  call N_x%distr()
  call N_y%distr()                        
  
  deallocate(mass)
  deallocate(Cw)
  deallocate(work)
  deallocate(a_w)
  deallocate(ia_w)
  deallocate(ja_w)

end subroutine A_matrix_assembling

!! R assemling subroutine

subroutine R_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, number_of_triangles, &
   ind, R_matr)

  implicit none 
  
  type(container_element), dimension(nvmax)     ,intent(in) :: List_of_non_Direchlet_Elements
  integer                                       ,intent(in) :: number_of_non_Direchlet_elements, ind, number_of_triangles
  type(Matrix)                               ,intent(inout) :: R_matr
  integer                                                   :: i, j, k, nonzero, nonzero_raw
  integer, allocatable                                      :: ia(:), ja(:)
  real*8, allocatable                                       :: a(:)
  type(Element), pointer                                    :: elem_master
  type(Triangle), pointer                                   :: trian
  real*8                                                    :: prom
  real*8                                                    :: coeff(3)
  logical                                                   :: values
  
  allocate(ia(number_of_non_Direchlet_elements*2))
  allocate(ja(number_of_non_Direchlet_elements*10))
  allocate(a(number_of_non_Direchlet_elements*10))
  
  ia(1) = 1
  k = 1
  
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
  
    nonzero_raw = 0
    elem_master => List_of_non_Direchlet_Elements(i)%pointing_element
    
    do j = 1, elem_master%number_of_neighbour_triangles
    
      trian => elem_master%neighbour_triangles_list(j)%pointing_triangle
      
      coeff = coefficients_calculation(elem_master, trian)
      
      if (ind == 1) then
      
        prom = coeff(1)*trian%size_of_triangle
      
      else
      
        prom = coeff(2)*trian%size_of_triangle
      
      end if
      
      if(abs(prom) > varepsilon) then
      
        a(k) = prom
        ja(k) = trian%identificator
        k = k + 1
        nonzero_raw = nonzero_raw + 1
        
      end if    
    
    end do
    
    ia(i+1) = ia(i) + nonzero_raw
    nonzero = nonzero + nonzero_raw
  
  end do
  
  call csort (number_of_non_Direchlet_elements, a, ja, ia, values)
  
  call R_matr%init(nonzero, number_of_triangles, number_of_non_Direchlet_elements)
  
  R_matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  R_matr%ja(1:nonzero) = ja(1:nonzero)
  R_matr%a(1:nonzero) = a(1:nonzero)
  
  deallocate(ia)
  deallocate(ja)
  deallocate(a)

end subroutine R_assembling



!! N assemling subroutine

subroutine N_assembling(List_of_Triangles, number_of_triangles, number_of_non_Direchlet_elements, &
   ind, N_matr)

  implicit none 
  
  type(Triangle), target, dimension(ntmax)      ,intent(in) :: List_of_Triangles
  integer                                       ,intent(in) :: number_of_non_Direchlet_elements, ind, number_of_triangles
  type(Matrix)                               ,intent(inout) :: N_matr
  integer                                                   :: i, j, k, nonzero, nonzero_raw
  integer, allocatable                                      :: ia(:), ja(:)
  real*8, allocatable                                       :: a(:)
  type(Element), pointer                                    :: elem
  type(Triangle), pointer                                   :: trian_master
  real*8                                                    :: prom
  real*8                                                    :: coeff(3)
  logical                                                   :: values
  
  allocate(ia(number_of_triangles*2))
  allocate(ja(number_of_triangles*10))
  allocate(a(number_of_triangles*10))
  
  ia(1) = 1
  k = 1
  
  nonzero = 0
  
  do i = 1, number_of_triangles
  
    nonzero_raw = 0
    trian_master => List_of_Triangles(i)
    
    do j = 1, 3
    
      elem => trian_master%neighbour_elements_list(j)%pointing_element
      
      if (elem%on_boundary .eqv. .false.) then
      
        coeff = coefficients_calculation(elem, trian_master)
      
        if (ind == 1) then
      
          prom = coeff(1)*(trian_master%P_0/(trian_master%delta + delta_min))
      
        else
      
          prom = coeff(2)*(trian_master%P_0/(trian_master%delta + delta_min))
      
        end if
      
        if(abs(prom) > varepsilon) then
      
          a(k) = prom
          ja(k) = elem%nd_identificator
          k = k + 1
          nonzero_raw = nonzero_raw + 1
        
        end if    
    
     end if
    
    end do
    
    ia(i+1) = ia(i) + nonzero_raw
    nonzero = nonzero + nonzero_raw
  
  end do
  
  call csort (number_of_triangles,a,ja,ia, values)
  
  call N_matr%init(nonzero, number_of_non_Direchlet_elements, number_of_triangles)
  
  N_matr%ia(1:(number_of_triangles+1)) = ia(1:(number_of_triangles+1))
  N_matr%ja(1:nonzero) = ja(1:nonzero)
  N_matr%a(1:nonzero) = a(1:nonzero)
  
  deallocate(ia)
  deallocate(ja)
  deallocate(a)

end subroutine N_assembling

function residual(nd_mass_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
  number_of_non_Direchlet_elements, number_of_triangles, time_step, standart) result(resuidal_res)   ! standart = 1 -- unit residual vector
                                                                                                     ! standart = 0 -- not unit residual vector
  implicit none

  type(Triangle), target, dimension(ntmax), intent(inout)   :: List_of_Triangles
  type(container_element), dimension(nvmax), intent(inout)  :: List_of_non_Direchlet_Elements
  real*8                                                    :: resuidal_res(2*number_of_non_Direchlet_elements)
  type(Matrix)                                              :: nd_mass_lum, M, M_cor, M_w, M_a, &
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
  
  call  N_assembling_resuid(N_xx_11, 1, 1, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements) 
  call  N_assembling_resuid(N_xx_10, 1, 1, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_01, 1, 2, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_1m1, 1, 2, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_1m1, 2, 1, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_01, 2, 1, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = nd_mass_lum%a(i)*(List_of_non_Direchlet_Elements(i)%pointing_element%m)/(time_step)
    M%ia(i+1) = nd_mass_lum%ia(i+1)
    M%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_w%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    M_w%a(i) = nd_mass_lum%a(i)*C_w*rho_water*List_of_non_Direchlet_Elements(i)%pointing_element%A*u_minus_u
    M_w%ia(i+1) = nd_mass_lum%ia(i+1)
    M_w%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! M_a
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_a%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    u_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
    abs_vel_air = dsqrt(u_a_loc**2 + v_a_loc**2) 
    M_a%a(i) = nd_mass_lum%a(i)*C_a*rho_air*List_of_non_Direchlet_Elements(i)%pointing_element%A*abs_vel_air
    M_a%ia(i+1) = nd_mass_lum%ia(i+1)
    M_a%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! M_cor
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_cor%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    M_cor%a(i) = nd_mass_lum%a(i)*2d0*omega_e*(List_of_non_Direchlet_Elements(i)%pointing_element%m)
    M_cor%ia(i+1) = nd_mass_lum%ia(i+1)
    M_cor%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! p_x, p_y
  
  do i = 1, number_of_non_Direchlet_elements
  
    p_x(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 1)
    p_y(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 2)
  
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
    
    res1(i) = res1(i) + prom1(i) 
        
  end do
  
  ! + N_yx_1m1*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_yx_1m1%a, N_yx_1m1%ja, N_yx_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! - p_x
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - p_x(i)  
        
  end do
  
  ! + N_yy_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! + N_xy_01*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_xy_01%a, N_xy_01%ja, N_xy_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! - N_yy_10*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_10%a, N_yy_10%ja, N_yy_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i)
        
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
    
    res2(i) = res2(i) + prom2(i)
        
  end do
  
  ! + N_xy_1m1*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_xy_1m1%a, N_xy_1m1%ja, N_xy_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)
        
  end do
  
  ! - p_y
   
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - p_y(i)  
        
  end do
  
  
  ! + N_xx_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) 
        
  end do
  
  
  ! + N_yx_01*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_yx_01%a, N_yx_01%ja, N_yx_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i) 
        
  end do
  
  ! - N_xx_10*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_10%a, N_xx_10%ja, N_xx_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)
        
  end do
  
  
  ! - M*v_n
  
  call amux(number_of_non_Direchlet_elements, v_n, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)
        
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

subroutine N_assembling_resuid(Matr, x_y_fir, x_y_sec, xi_ind, eta_ind, List_of_non_Direchlet_Elements, &
  number_of_non_Direchlet_elements)
  type(Matrix) :: Matr
  integer :: number_of_non_Direchlet_elements
  integer :: x_y_fir, x_y_sec
  real*8 :: xi_ind, eta_ind
  real*8 :: a(maxnz), prom
  integer :: ia(maxn), ja(maxnz)
  integer :: i, j, k, nonzero, nonzero_raw
  type(Element), pointer :: elem_master, elem_side
  type(Triangle), pointer :: trian
  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  
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

  type(Element), target :: elem
  integer :: i, j, k, d_phi_dxy
  type(Triangle), pointer :: trian
  real*8 :: sc_res
  real*8 :: coeff(3)
  
  sc_res = 0d0
  
  do i = 1, elem%number_of_neighbour_triangles
        
     trian => elem%neighbour_triangles_list(i)%pointing_triangle
     coeff = coefficients_calculation(elem, trian)
     sc_res = sc_res + trian%size_of_triangle*coeff(d_phi_dxy)*trian%P_0_resuid
  
  end do

end function scalar_mult_P_dphi_resuid




function scalar_mult_resuid(elem1, elem2, d_phi1_dxy, d_phi2_dxy, xi_ind, eta_ind) result(B_res)

  type(Element), target :: elem1, elem2
  integer :: i, j, k, d_phi1_dxy, d_phi2_dxy
  type(Triangle), pointer :: trian
  real*8 :: B_res
  real*8 :: coeff1(3), coeff2(3)
  real*8 :: xi_ind, eta_ind
  
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


function init_residual(nd_mass_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
  number_of_non_Direchlet_elements, number_of_triangles, time_step, standart) result(resuidal_res)

  implicit none

  type(Triangle), target, dimension(ntmax), intent(inout)   :: List_of_Triangles
  type(container_element), dimension(nvmax), intent(inout)  :: List_of_non_Direchlet_Elements
  real*8                                                    :: resuidal_res(2*number_of_non_Direchlet_elements)
  type(Matrix)                                              :: nd_mass_lum, M, M_cor, M_w, M_a, &
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
  
  call  N_assembling_resuid(N_xx_11, 1, 1, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements) 
  call  N_assembling_resuid(N_xx_10, 1, 1, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_01, 1, 2, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_1m1, 1, 2, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_1m1, 2, 1, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_01, 2, 1, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = nd_mass_lum%a(i)*(List_of_non_Direchlet_Elements(i)%pointing_element%m)/(time_step)
    M%ia(i+1) = nd_mass_lum%ia(i+1)
    M%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_w%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_resuid(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    M_w%a(i) = nd_mass_lum%a(i)*C_w*rho_water*List_of_non_Direchlet_Elements(i)%pointing_element%A*u_minus_u
    M_w%ia(i+1) = nd_mass_lum%ia(i+1)
    M_w%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! M_a
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_a%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    u_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
    v_a_loc = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
    abs_vel_air = dsqrt(u_a_loc**2 + v_a_loc**2) 
    M_a%a(i) = nd_mass_lum%a(i)*C_a*rho_air*List_of_non_Direchlet_Elements(i)%pointing_element%A*abs_vel_air
    M_a%ia(i+1) = nd_mass_lum%ia(i+1)
    M_a%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! M_cor
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_cor%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    M_cor%a(i) = nd_mass_lum%a(i)*2d0*omega_e*(List_of_non_Direchlet_Elements(i)%pointing_element%m)
    M_cor%ia(i+1) = nd_mass_lum%ia(i+1)
    M_cor%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! p_x, p_y
  
  do i = 1, number_of_non_Direchlet_elements
  
    p_x(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 1)
    p_y(i) = 5d-1*scalar_mult_P_dphi_resuid(List_of_non_Direchlet_Elements(i)%pointing_element, 2)
  
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
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! + N_yx_1m1*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_yx_1m1%a, N_yx_1m1%ja, N_yx_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! - p_x
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - p_x(i)  
        
  end do
  
  ! + N_yy_11*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_11%a, N_yy_11%ja, N_yy_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! + N_xy_01*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom1, N_xy_01%a, N_xy_01%ja, N_xy_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) + prom1(i)
        
  end do
  
  ! - N_yy_10*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom1, N_yy_10%a, N_yy_10%ja, N_yy_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res1(i) = res1(i) - prom1(i) 
        
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
    
    res2(i) = res2(i) + prom2(i)
        
  end do
  
  ! + N_xy_1m1*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_xy_1m1%a, N_xy_1m1%ja, N_xy_1m1%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)
        
  end do
  
  ! - p_y
   
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - p_y(i)  
        
  end do
  
  
  ! + N_xx_11*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_11%a, N_xx_11%ja, N_xx_11%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)  
        
  end do
  
  
  ! + N_yx_01*u
  
  call amux(number_of_non_Direchlet_elements, u_val, prom2, N_yx_01%a, N_yx_01%ja, N_yx_01%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) + prom2(i)
        
  end do
  
  ! - N_xx_10*v
  
  call amux(number_of_non_Direchlet_elements, v_val, prom2, N_xx_10%a, N_xx_10%ja, N_xx_10%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)
        
  end do
  
  
  ! - M*v_n
  
  call amux(number_of_non_Direchlet_elements, v_n, prom2, M%a, M%ja, M%ia)
  
  do i = 1, number_of_non_Direchlet_elements
    
    res2(i) = res2(i) - prom2(i)
        
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

!! L2 norm of vector

function L2_norm(vec, size_of_vec) result(L2_norm_res)

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


end module module_dynamics 
