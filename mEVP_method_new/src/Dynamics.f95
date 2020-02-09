module module_dynamics

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_matricies
  use module_grid_values
  use module_assembling
  
  implicit none
  
  private
  public :: L2_error,                              &
            N_assembling_resuid,                   &
            scalar_mult_P_dphi_resuid,             &
            scalar_mult_resuid,                    &
            L2_norm,                               &
            VP_Mass_assembling,                    &
            VP_B_assembling,                       &
            scalar_mult_B,                         &
            VP_rhs_assembling,                     &
            scalar_mult_P_dphi
  
  contains

!! L2 error function

function L2_error() result(L2_res)

  implicit none 
  
  integer                                                   :: i, j
  real*8                                                    :: L2_res, prom1, prom2, u_coefficients(3), &
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
  integer                 :: i, d_phi_dxy
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
  integer                             :: i, j, d_phi1_dxy, d_phi2_dxy
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
  integer                    :: i
  real*8                     :: prom
  
  prom = 0d0
  
  do i = 1, size_of_vec
  
    prom = prom + vec(i)**2
  
  end do
  
  L2_norm_res = dsqrt(prom)

end function L2_norm



subroutine VP_Mass_assembling(VP_mass_matrix, time_step)
  
  implicit none
  
  type(Matrix), intent(inout) :: VP_mass_matrix
  real*8, intent(in)          :: time_step
  real*8                      :: mas, conc, u_minus_u
  integer                     :: i
  real*8                      :: a(maxnz)
  
    
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
  
  integer                                     :: i, j, k
  type(Matrix)                                :: B_VP
  type(Element), pointer                      :: elem_master, elem_side
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
  integer                             :: i, j
  
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
  
  integer                  :: i, j
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
  integer                           :: i
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
