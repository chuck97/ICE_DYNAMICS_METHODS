module module_dynamics

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_matricies
  
  implicit none
  
  private
  public :: VP_Mass_assembling, VP_B_assembling, &
    velocity_delta_str_P_recalculation, scalar_mult_B, VP_rhs_assembling, &
    scalar_mult_P_dphi, VP_coriolis_update, maximum_u, triangles_values_calculations, &
    resuidal, scalar_mult_resuid, N_assembling, scalar_mult_P_dphi_resuid, N_assembling_resuid
  
  contains

! VP Mass matrix assembling

subroutine VP_Mass_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
     VP_mass_matrix, lum_Mass_Matr, time_step)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements  
  type(Matrix) :: VP_mass_matrix, lum_Mass_Matr
  real*8 :: mas, conc, u_minus_u, time_step
  integer :: i, number_of_non_Direchlet_elements
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
    
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
  
    a(i) = lum_Mass_Matr%a(i)*(mas/time_step + C_w*conc*rho_water*u_minus_u)
    
  end do
  
  call VP_mass_matrix%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  VP_mass_matrix%ia = lum_Mass_Matr%ia
  VP_mass_matrix%ja = lum_Mass_Matr%ja
  VP_mass_matrix%a(1:number_of_non_Direchlet_elements) = a(1:number_of_non_Direchlet_elements)
  
  
end subroutine VP_Mass_assembling
  
! velocity and delta str recalculation

subroutine velocity_delta_str_P_recalculation(List_of_non_Direchlet_Elements, List_of_Triangles, &
 number_of_non_Direchlet_elements, number_of_triangles, npic)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements  
  type(Triangle), target, dimension(ntmax) :: List_of_Triangles
  integer :: i, j, k, number_of_non_Direchlet_elements, number_of_triangles
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


! scalar multiplication xi*(dphi_i/dx_j)

function scalar_mult_B(elem1, elem2, d_phi1_dxy, d_phi2_dxy, xi_ind, eta_ind) result(B_res)

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
         (xi_ind*trian%xi + eta_ind*trian%eta)  
        
      end if  
    
    end do
  
  end do

end function scalar_mult_B



! scalar multiplication P*(d(phi)/dx_i)

function scalar_mult_P_dphi(elem, d_phi_dxy) result(sc_res)

  type(Element), target :: elem
  integer :: i, j, k, d_phi_dxy
  type(Triangle), pointer :: trian
  real*8 :: sc_res
  real*8 :: coeff(3)
  
  sc_res = 0d0
  
  do i = 1, elem%number_of_neighbour_triangles
        
     trian => elem%neighbour_triangles_list(i)%pointing_triangle
     coeff = coefficients_calculation(elem, trian)
     sc_res = sc_res + trian%size_of_triangle*coeff(d_phi_dxy)*trian%P_0
  
  end do

end function scalar_mult_P_dphi


! Assembling B matrix in CSR format

subroutine VP_B_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, B_VP)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements, i, j, k, r
  type(Matrix) :: B_VP
  type(Element), pointer :: elem_master, elem_side
  type(Triangle), pointer :: trian 
  real*8 :: a(maxnz), prom
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero, nonzero_raw
  
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
  
  call B_VP%init(nonzero, number_of_non_Direchlet_elements)
  B_VP%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  B_VP%ja(1:nonzero) = ja(1:nonzero)
  B_VP%a(1:nonzero) = a(1:nonzero)

end subroutine VP_B_assembling



! rhs_B in VP method assembling

subroutine VP_rhs_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, rhs_VP_1, &
 rhs_VP_2, ND_mass_matrix_lum, time_step)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: i, j, k, number_of_non_Direchlet_elements
  real*8 :: rhs_VP_1(nvmax), rhs_VP_2(nvmax), mas, conc, u_minus_u, u_n, v_n, &
   u_str, v_str, u_w, v_w, u_a, v_a, abs_vel_air, time_step
  type(Matrix) :: ND_mass_matrix_lum
  type(Element), pointer :: master_elem, side_elem
 
  
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
  
    rhs_VP_1(i) = ND_mass_matrix_lum%a(i)*(-2d0*mas*omega_e*v_str + (mas/time_step)*u_n + C_w*conc*rho_water*u_w*u_minus_u + &
     C_a*conc*rho_air*u_a*abs_vel_air) !- 2d0*mas*omega_e*v_w)
    rhs_VP_2(i) = ND_mass_matrix_lum%a(i)*(2d0*mas*omega_e*u_str + (mas/time_step)*v_n + C_w*conc*rho_water*v_w*u_minus_u + &
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

! VP coriolis update

subroutine VP_coriolis_update(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, time_step)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: i, j, k, number_of_non_Direchlet_elements
  real*8 :: time_step
  real*8 :: rhs_1, rhs_2, mas, conc, u_minus_u, u_nvp, v_nvp, &
   u_str, v_str, u_w, v_w, u_a, v_a, abs_vel_air, fir, sec
   
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A 
    u_str = List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1)
    v_str = List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2)
    u_nvp = List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
    v_nvp = List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
    fir = mas/time_step + C_w*conc*rho_water*u_minus_u
    sec = 2d0*omega_e*mas
    
    rhs_1 = (mas/time_step)*u_nvp - 2d0*mas*omega_e*v_str + C_w*conc*rho_water*u_minus_u*u_nvp
    rhs_2 = (mas/time_step)*v_nvp + 2d0*mas*omega_e*u_str + C_w*conc*rho_water*u_minus_u*v_nvp
    
    List_of_non_Direchlet_Elements(i)%pointing_element%u(1) = (sec*rhs_2 + fir*rhs_1)/(fir**2 + sec**2)
    List_of_non_Direchlet_Elements(i)%pointing_element%u(2) = (fir*rhs_2 - sec*rhs_1)/(fir**2 + sec**2)
    
  
  end do   

end subroutine VP_coriolis_update

! maximum u function

function maximum_u(List_of_Elements, number_of_elements) result(res_max)

  type(Element), target,  dimension(nvmax) :: List_of_Elements
  integer :: number_of_elements, i, j, k
  real*8 :: res_max, prom
  
  prom = -5d2
  
  do i = 1, number_of_elements
  
    if(dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2) > prom) then
    
      prom = dsqrt(List_of_Elements(i)%u(1)**2 + List_of_Elements(i)%u(2)**2)
      
    end if  
  
  end do
  
  res_max = prom

end function maximum_u

! calculate values on triangles (for comparisson)

subroutine triangles_values_calculations(List_of_Triangles, number_of_triangles)

  type(Triangle), target, dimension(ntmax) :: List_of_Triangles
  integer :: number_of_triangles, i, j, k
  real*8 :: sigma11, sigma22, sigma12, x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, &
   du_dx, du_dy, dv_dx, dv_dy, dot_epsilon11, dot_epsilon22, dot_epsilon12, A_avg, h_avg
  
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

! resuidal calculation

! resuidal calculation

function resuidal(nd_mass_lum, List_of_non_Direchlet_Elements, List_of_Triangles, &
  number_of_non_Direchlet_elements, number_of_triangles, time_step) result(resuidal_res)

  type(Triangle), target, dimension(ntmax) :: List_of_Triangles
  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  real*8 :: resuidal_res
  type(Matrix) :: nd_mass_lum, M, M_cor, M_w, M_a, N_xx_11, N_xx_10, N_xy_01, N_xy_1m1, N_yx_1m1, N_yx_01, &
   N_yy_11, N_yy_10
  real*8 :: u_a(nvmax), v_a(nvmax), u_w(nvmax), v_w(nvmax), p_x(nvmax), p_y(nvmax), res1(nvmax), res2(nvmax), &
   u_val(nvmax), v_val(nvmax), u_n(nvmax), v_n(nvmax)
  integer :: number_of_non_Direchlet_elements, number_of_triangles, i, j, k
  real*8 :: xi_tr(ntmax), eta_tr(ntmax)
  real*8 :: x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, du_dx, du_dy, dv_dx, dv_dy, dot_epsilon11, &
   dot_epsilon22, dot_epsilon12, delta, A_avg, h_avg, P_0, u_minis_u, time_step, u_a_loc, v_a_loc, abs_vel_air
  real*8 :: prom1(nvmax), prom2(nvmax), u_minus_u, prom_res
  real*8 :: resuidal_init
   
  
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
  call  N_assembling_resuid(N_xy_01, 1, 1, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_1m1, 1, 1, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_1m1, 2, 2, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_01, 2, 2, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = nd_mass_lum%a(i)*(List_of_non_Direchlet_Elements(i)%pointing_element%m)/(time_step)
    M%ia(i+1) = nd_mass_lum%ia(i+1)
    M%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  prom_res = -1d10
  
  resuidal_res = 0d0
  
  do i = 1, number_of_non_Direchlet_elements
  
    resuidal_res = resuidal_res + res1(i)**2 
    resuidal_res = resuidal_res + res2(i)**2
  
  end do
  
  !resuidal_res = prom_res
  resuidal_res = dsqrt(resuidal_res)
  
  call M%distr()
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! initial resuidal calculation !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
  call  N_assembling_resuid(N_xy_01, 1, 1, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_xy_1m1, 1, 1, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_1m1, 2, 2, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yx_01, 2, 2, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_11, 2, 2, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling_resuid(N_yy_10, 2, 2, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
    
  ! M
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = nd_mass_lum%a(i)*(List_of_non_Direchlet_Elements(i)%pointing_element%m)/(time_step)
    M%ia(i+1) = nd_mass_lum%ia(i+1)
    M%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  
  ! M_w
  
  call M_w%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  call M_a%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
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
  
  prom_res = -1d10
  
  resuidal_init = 0d0
  
  do i = 1, number_of_non_Direchlet_elements
  
    resuidal_init = resuidal_init + res1(i)**2 
    resuidal_init = resuidal_init + res2(i)**2
  
  end do
  
  !resuidal_res = prom_res
  resuidal_init = dsqrt(resuidal_init)
  
  resuidal_res = resuidal_res/resuidal_init
  
  call M%distr()
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
  
  
end function resuidal



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
  
  call Matr%init(nonzero, number_of_non_Direchlet_elements)
  Matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  Matr%ja(1:nonzero) = ja(1:nonzero)
  Matr%a(1:nonzero) = a(1:nonzero)

end subroutine N_assembling_resuid


! N matrix assembling

subroutine N_assembling(Matr, x_y_fir, x_y_sec, xi_ind, eta_ind, List_of_non_Direchlet_Elements, &
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
  
  call Matr%init(nonzero, number_of_non_Direchlet_elements)
  Matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  Matr%ja(1:nonzero) = ja(1:nonzero)
  Matr%a(1:nonzero) = a(1:nonzero)

end subroutine N_assembling

! scalar multiplication xi*(dphi_i/dx_j) for resuidal

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

! scalar multiplication P*(d(phi)/dx_i) for resuidal

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




end module module_dynamics 
