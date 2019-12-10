module module_dynamics

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_matricies
  
  implicit none
  
  private
  public ::   velocity_delta_str_P_recalculation, scalar_mult_B, VP_rhs_assembling, &
    scalar_mult_P_dphi, maximum_u, triangles_values_calculations, &
    residual, scalar_mult_resuid, N_assembling_resuid, scalar_mult_P_dphi_resuid, &
    N_assembling, VP_Mass_assembling, VP_B_assembling, plot_matr, frobenius_norm, &
    L2_norm
  
  contains
  
! velocity and delta str recalculation

subroutine velocity_delta_str_P_recalculation(List_of_non_Direchlet_Elements, List_of_Triangles, &
 number_of_non_Direchlet_elements, number_of_triangles, npic, xi_loc_ident)      !! xi_loc_ident = 1 -- gain viscosities and delta from previous local step
                                                                                 !! xi_loc_ident = 0 -- gain viscosities and delta from previous global step
  implicit none                
                                                                      
  type(container_element), dimension(nvmax), intent(inout)   :: List_of_non_Direchlet_Elements  
  type(Triangle), target, dimension(ntmax), intent(inout)    :: List_of_Triangles
  integer, intent(in)                                        :: number_of_non_Direchlet_elements, number_of_triangles, &
                                                                xi_loc_ident, npic
  integer                                                    :: i, j, k
  real*8                                                     :: x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2, du_dx, &
                                                                du_dy, dv_dx, dv_dy, dot_epsilon11, dot_epsilon12,     &
                                                                dot_epsilon22, A_avg, h_avg
  
  do i = 1, number_of_non_Direchlet_elements
  
    if (xi_loc_ident .eq. 1) then
    
      if (npic == 1) then
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
        
        List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u(2)  
         
      else
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
        
        List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
         
      end if
      
    else
    
       List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
        
        List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
      
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
      
  end do
  
end subroutine velocity_delta_str_P_recalculation


! scalar multiplication xi*(dphi_i/dx_j)

function scalar_mult_B(elem1, elem2, d_phi1_dxy, d_phi2_dxy, xi_ind, eta_ind) result(B_res)

  implicit none

  type(Element), target, intent(in)   :: elem1, elem2
  integer, intent(in)                 :: d_phi1_dxy, d_phi2_dxy
  real*8, intent(in)                  :: xi_ind, eta_ind
  integer                             :: i, j, k
  type(Triangle), pointer             :: trian
  real*8                              :: B_res
  real*8                              :: coeff1(3), coeff2(3)
  
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
     sc_res = sc_res + trian%size_of_triangle*coeff(d_phi_dxy)*(trian%P_0)
  
  end do

end function scalar_mult_P_dphi


! Assembling B matrix in CSR format

subroutine VP_B_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 ND_mass_matrix_lum, B_VP)
 
  implicit none

  type(container_element), dimension(nvmax), intent(in)  :: List_of_non_Direchlet_Elements
  integer                                  , intent(in)  :: number_of_non_Direchlet_elements
  integer                                                :: i, j, k, r
  type(Matrix)                             , intent(in)  :: ND_mass_matrix_lum 
  type(Matrix)                             , intent(out) :: B_VP
  type(Matrix)                                           :: N_xx_11, N_xx_10, N_xy_01, N_xy_1m1, N_yx_1m1, &
                                                            N_yx_01, N_yy_11, N_yy_10
  type(Matrix)                                           :: M1, M2, M3, M4
  integer                                                :: ierr
  real*8, allocatable                                    :: M1_1_a(:), M1_2_a(:), M2_a(:), M3_a(:), M4_1_a(:), &
                                                            M4_2_a(:)
  integer, allocatable                                   :: M1_1_ja(:), M1_1_ia(:), M1_2_ja(:), M1_2_ia(:), &
                                                            M2_ja(:), M2_ia(:), M3_ja(:), M3_ia(:), M4_1_ja(:), &
                                                            M4_1_ia(:), M4_2_ja(:), M4_2_ia(:)  
  integer                                                :: iW(number_of_non_Direchlet_elements)  
  integer, allocatable                                   :: M1_res_ia(:), M1_res_ja(:), M2_res_ia(:), M2_res_ja(:), &
                                                            M3_res_ia(:), M3_res_ja(:), &
                                                            M4_res_ia(:), M4_res_ja(:)
  real*8, allocatable                                    :: M1_res_a(:), M2_res_a(:), M3_res_a(:), M4_res_a(:)
  real*8, allocatable                                    :: diag(:)
  
  !! allocating arrays
  
  allocate(diag(number_of_non_Direchlet_elements + 100))
  
  allocate(M1_1_ia(number_of_non_Direchlet_elements + 1))
  allocate(M1_1_a(number_of_non_Direchlet_elements*20))
  allocate(M1_1_ja(number_of_non_Direchlet_elements*10))
  
  allocate(M1_2_ia(number_of_non_Direchlet_elements + 1))
  allocate(M1_2_a(number_of_non_Direchlet_elements*20))
  allocate(M1_2_ja(number_of_non_Direchlet_elements*10)) 
  
  allocate(M1_res_ia(number_of_non_Direchlet_elements + 1))
  allocate(M1_res_a(number_of_non_Direchlet_elements*20))
  allocate(M1_res_ja(number_of_non_Direchlet_elements*10))   
  
  allocate(M2_ia(number_of_non_Direchlet_elements + 1))
  allocate(M2_ja(number_of_non_Direchlet_elements*20))
  allocate(M2_a(number_of_non_Direchlet_elements*20))
  
  allocate(M2_res_ia(number_of_non_Direchlet_elements + 1))
  allocate(M2_res_ja(number_of_non_Direchlet_elements*20))
  allocate(M2_res_a(number_of_non_Direchlet_elements*20))
  
  allocate(M3_ia(number_of_non_Direchlet_elements*2 + 1))
  allocate(M3_ja(number_of_non_Direchlet_elements*20))
  allocate(M3_a(number_of_non_Direchlet_elements*20))
  
  allocate(M3_res_ia(number_of_non_Direchlet_elements + 1))
  allocate(M3_res_ja(number_of_non_Direchlet_elements*20))
  allocate(M3_res_a(number_of_non_Direchlet_elements*20)) 
  
  allocate(M4_1_ia(number_of_non_Direchlet_elements + 1))
  allocate(M4_1_a(number_of_non_Direchlet_elements*20))
  allocate(M4_1_ja(number_of_non_Direchlet_elements*20)) 
  
  allocate(M4_2_ia(number_of_non_Direchlet_elements + 1))
  allocate(M4_2_a(number_of_non_Direchlet_elements*20))
  allocate(M4_2_ja(number_of_non_Direchlet_elements*20)) 
  
  allocate(M4_res_ia(number_of_non_Direchlet_elements + 1))
  allocate(M4_res_a(number_of_non_Direchlet_elements*20))
  allocate(M4_res_ja(number_of_non_Direchlet_elements*20))
  
  !!! diagonal matrix m^-1 * M_L^-1 assembling
  
  do i = 1, number_of_non_Direchlet_elements
  
    diag(i) = 1d0/(ND_mass_matrix_lum%a(i)*List_of_non_Direchlet_Elements(i)%pointing_element%m)
    
  end do
  
  !!! N matricies assembling
  
  call  N_assembling(N_xx_11, 1, 1, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements) 
  call  N_assembling(N_xx_10, 1, 1, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_xy_01, 1, 2, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_xy_1m1, 1, 2, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_yx_1m1, 2, 1, 1d0, -1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_yx_01, 2, 1, 0d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_yy_11, 2, 2, 1d0, 1d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  call  N_assembling(N_yy_10, 2, 2, 1d0, 0d0, List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements)
  
  !!! reverse N_yy_10 and N_xx_10
  
  do i = 1, N_yy_10%nonzero
  
    N_yy_10%a(i) = -1d0*N_yy_10%a(i)
  
  end do
  
  do i = 1, N_xx_10%nonzero
  
    N_xx_10%a(i) = -1d0*N_xx_10%a(i)
  
  end do
  
  !!!!! M1 assembling
  
  M1_1_ia = 1
  
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, 1, &
          N_xx_11%a, N_xx_11%ja, N_xx_11%ia, N_yy_11%a, N_yy_11%ja, N_yy_11%ia, M1_1_a, M1_1_ja, M1_1_ia, &
          maxnz, iw, ierr)
          
  M1_2_ia = 1        
                          
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, 1, &
          M1_1_a, M1_1_ja, M1_1_ia, N_yy_10%a, N_yy_10%ja, N_yy_10%ia, M1_2_a, M1_2_ja, M1_2_ia, &
          maxnz, iw, ierr)
          
  M1_res_ia = 1        
          
  call diamua(number_of_non_Direchlet_elements, 1, M1_2_a, M1_2_ja, M1_2_ia, diag, M1_res_a, M1_res_ja, M1_res_ia)
  
  call M1%init(M1_res_ia(number_of_non_Direchlet_elements+1) - 1, number_of_non_Direchlet_elements)
  
  do i = 1, M1%nonzero
  
    M1%a(i) = M1_res_a(i)
    M1%ja(i) = M1_res_ja(i)
  
  end do
  
  do i = 1, (M1%matr_size + 1)
  
    M1%ia(i) = M1_res_ia(i)
  
  end do
  
  !!!!! M2 assembling
  
  M2_ia = 1
  
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, 1, &
          N_yx_1m1%a, N_yx_1m1%ja, N_yx_1m1%ia, N_xy_01%a, N_xy_01%ja, N_xy_01%ia, M2_a, M2_ja, M2_ia, &
          maxnz, iw, ierr)
  
  M2_res_ia = 1                
                  
  call diamua(number_of_non_Direchlet_elements, 1, M2_a, M2_ja, M2_ia, diag, M2_res_a, M2_res_ja, M2_res_ia)
  
  call M2%init(M2_res_ia(number_of_non_Direchlet_elements+1) - 1, number_of_non_Direchlet_elements)
  
  do i = 1, M2%nonzero
  
    M2%a(i) = M2_res_a(i)
    M2%ja(i) = M2_res_ja(i)
  
  end do
  
  do i = 1, (M2%matr_size + 1)
  
    M2%ia(i) = M2_res_ia(i)
  
  end do
  
  !!!!! M3 assembling
  
  M3_ia = 1
  
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, &
          1, N_xy_1m1%a, N_xy_1m1%ja, N_xy_1m1%ia, N_yx_01%a, N_yx_01%ja, N_yx_01%ia, M3_a, M3_ja, M3_ia, &
          maxnz, iw, ierr)
                  
  M3_res_ia = 1               
                  
  call diamua(number_of_non_Direchlet_elements, 1, M3_a, M3_ja, M3_ia, diag, M3_res_a, M3_res_ja, M3_res_ia)
  
  call M3%init(M3_res_ia(number_of_non_Direchlet_elements+1) - 1, number_of_non_Direchlet_elements)
  
  do i = 1, M3%nonzero
  
    M3%a(i) = M3_res_a(i)
    M3%ja(i) = M3_res_ja(i)
  
  end do
  
  do i = 1, (M3%matr_size + 1)
  
    M3%ia(i) = M3_res_ia(i)
  
  end do
  
  !!!!! M4 assembling
  
  M4_1_ia = 1
  
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, &
          1, N_yy_11%a, N_yy_11%ja, N_yy_11%ia, N_xx_11%a, N_xx_11%ja, N_xx_11%ia, M4_1_a, M4_1_ja, M4_1_ia, &
          maxnz, iw, ierr)
          
  M4_2_ia = 1       
                  
  call aplb(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements, &
          1, M4_1_a, M4_1_ja, M4_1_ia, N_xx_10%a, N_xx_10%ja, N_xx_10%ia, M4_2_a, M4_2_ja, M4_2_ia, &
          maxnz, iw, ierr)
          
  M4_res_ia = 1        
  
  call diamua(number_of_non_Direchlet_elements, 1, M4_2_a, M4_2_ja, M4_2_ia, diag, M4_res_a, M4_res_ja, M4_res_ia)
  
  call M4%init(M4_res_ia(number_of_non_Direchlet_elements+1) - 1, number_of_non_Direchlet_elements)
  
  do i = 1, M4%nonzero
  
    M4%a(i) = M4_res_a(i)
    M4%ja(i) = M4_res_ja(i)
  
  end do
  
  do i = 1, (M4%matr_size + 1)
  
    M4%ia(i) = M4_res_ia(i)
  
  end do
  
  
  !! B_VP matrix assembling
  
  !print *, "B_matrix"
   
  call four_block_matrix_assembling(M1, M2, M3, M4, B_VP)
  
  ! deallocating arrays
  
  deallocate(M1_1_ia)
  deallocate(M1_1_a)
  deallocate(M1_1_ja)
  
  deallocate(M1_2_ia)
  deallocate(M1_2_a)
  deallocate(M1_2_ja)
  
  deallocate(M1_res_ia)
  deallocate(M1_res_a)
  deallocate(M1_res_ja)
  
  deallocate(M2_ia)
  deallocate(M2_ja)
  deallocate(M2_a)
  
  deallocate(M2_res_ia)
  deallocate(M2_res_ja)
  deallocate(M2_res_a)
  
  deallocate(M3_ia)
  deallocate(M3_ja)
  deallocate(M3_a)
  
  deallocate(M3_res_ia)
  deallocate(M3_res_ja)
  deallocate(M3_res_a)
  
  deallocate(M4_1_ia)
  deallocate(M4_1_a)
  deallocate(M4_1_ja)
  
  deallocate(M4_2_ia)
  deallocate(M4_2_a)
  deallocate(M4_2_ja) 
  
  deallocate(M4_res_ia)
  deallocate(M4_res_a)
  deallocate(M4_res_ja)
  
  
  call N_xx_11%distr() 
  call N_xx_10%distr()
  call N_xy_01%distr()
  call N_xy_1m1%distr()
  call N_yx_1m1%distr()
  call N_yx_01%distr()
  call N_yy_11%distr()
  call N_yy_10%distr() 
  call M1%distr()
  call M2%distr()
  call M3%distr()
  call M4%distr()
  deallocate(diag)

end subroutine VP_B_assembling


! VP Mass matrix assembling

subroutine VP_Mass_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
     VP_mass_matrix, lum_Mass_Matr, time_step)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements  
  type(Matrix) :: VP_mass_matrix, lum_Mass_Matr
  real*8 :: mas, conc, u_minus_u, time_step
  integer :: i, j, k, number_of_non_Direchlet_elements
  type(Matrix) :: M11, M22, M33, M44
  
  !!!! M1 assembling
  
  call M11%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M11%ia(1) = 1
    
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
  
    M11%a(i) = 1d0/time_step + C_w*(conc/mas)*rho_water*u_minus_u
    M11%ia(i+1) = i+1
    M11%ja(i) = i
    
  end do
  
  !!!! M4 assembling
  
  call M44%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M44%ia(1) = 1
    
  do i = 1, number_of_non_Direchlet_elements
  
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    conc = List_of_non_Direchlet_Elements(i)%pointing_element%A
    u_minus_u = dsqrt((List_of_non_Direchlet_Elements(i)%pointing_element%u_str(1) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1))**2 + &
    (List_of_non_Direchlet_Elements(i)%pointing_element%u_str(2) - &
    List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2))**2)
  
    M44%a(i) = 1d0/time_step + C_w*(conc/mas)*rho_water*u_minus_u
    M44%ia(i+1) = i+1
    M44%ja(i) = i
    
  end do
  
  !!!! M2 assembling
  
  call M22%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M22%ia(1) = 1
    
  do i = 1, number_of_non_Direchlet_elements
  
    M22%a(i) = -2d0*omega_e
    M22%ia(i+1) = i+1
    M22%ja(i) = i
    
  end do
  
  !!!! M3 assembling
  
  call M33%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M33%ia(1) = 1
    
  do i = 1, number_of_non_Direchlet_elements
  
    M33%a(i) = 2d0*omega_e
    M33%ia(i+1) = i+1
    M33%ja(i) = i
    
  end do
  
  
  !!!! VP_mass_matrix assembling
  
  !print *, "mass matrix"
  
  call four_block_matrix_assembling(M11, M22, M33, M44, VP_mass_matrix)
  
  call M11%distr()
  call M22%distr()
  call M33%distr()
  call M44%distr()
  
end subroutine VP_Mass_assembling



! rhs_B in VP method assembling

subroutine VP_rhs_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, rhs_VP, &
  ND_mass_matrix_lum, time_step)

  type(container_element), dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: i, j, k, number_of_non_Direchlet_elements
  real*8 :: rhs_VP_1(nvmax), rhs_VP_2(nvmax), mas, conc, u_minus_u, u_n, v_n, &
   u_str, v_str, u_w, v_w, u_a, v_a, abs_vel_air, time_step, rhs_VP(2*nvmax)
  type(Matrix) :: ND_mass_matrix_lum
  type(Element), pointer :: master_elem
 
  
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
  
    rhs_VP_1(i) = (1d0/time_step)*u_n + C_w*(conc/mas)*rho_water*u_w*u_minus_u + &
     C_a*(conc/mas)*rho_air*u_a*abs_vel_air - 2d0*omega_e*v_w
    rhs_VP_2(i) = (1d0/time_step)*v_n + C_w*(conc/mas)*rho_water*v_w*u_minus_u + &
     C_a*(conc/mas)*rho_air*v_a*abs_vel_air + 2d0*omega_e*u_w
  
  end do
  
  do i = 1, number_of_non_Direchlet_elements
    
    master_elem => List_of_non_Direchlet_Elements(i)%pointing_element 
    mas = List_of_non_Direchlet_Elements(i)%pointing_element%m
    
    rhs_VP_1(i) = rhs_VP_1(i) + 5d-1*scalar_mult_P_dphi(master_elem, 1)/(mas*ND_mass_matrix_lum%a(i))
    rhs_VP_2(i) = rhs_VP_2(i) + 5d-1*scalar_mult_P_dphi(master_elem, 2)/(mas*ND_mass_matrix_lum%a(i))
  
  end do
  
  do i = 1, number_of_non_Direchlet_elements
  
    rhs_VP(i) = rhs_VP_1(i)
    rhs_VP(number_of_non_Direchlet_elements + i) = rhs_VP_2(i)
    
  end do

end subroutine VP_rhs_assembling

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
  call M%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
  
    M%a(i) = nd_mass_lum%a(i)/(time_step)
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
    M_w%a(i) = nd_mass_lum%a(i)*C_w*rho_water*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*u_minus_u
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
    M_a%a(i) = nd_mass_lum%a(i)*C_a*rho_air*List_of_non_Direchlet_Elements(i)%pointing_element%A/ &
    (List_of_non_Direchlet_Elements(i)%pointing_element%m)*abs_vel_air
    M_a%ia(i+1) = nd_mass_lum%ia(i+1)
    M_a%ja(i) = nd_mass_lum%ja(i)
  
  end do
  
  ! M_cor
  
  call M_cor%init(number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  
  M_cor%ia(1) = 1
  
  do i = 1, number_of_non_Direchlet_elements
   
    M_cor%a(i) = nd_mass_lum%a(i)*2d0*omega_e
    M_cor%ia(i+1) = nd_mass_lum%ia(i+1)
    M_cor%ja(i) = nd_mass_lum%ja(i)
  
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

subroutine N_assembling_resuid(N_Matr, x_y_fir, x_y_sec, xi_ind, eta_ind, List_of_non_Direchlet_Elements, &
  number_of_non_Direchlet_elements)
  
  implicit none
  
  type(Matrix), intent(inout)                           :: N_Matr
  integer, intent(in)                                   :: number_of_non_Direchlet_elements
  integer, intent(in)                                   :: x_y_fir, x_y_sec
  real*8,  intent(in)                                   :: xi_ind, eta_ind
  real*8                                                :: a(maxnz), prom
  integer                                               :: ia(maxn), ja(maxnz)
  integer                                               :: i, j, k, nonzero, nonzero_raw
  type(Element), pointer                                :: elem_master, elem_side
  type(Triangle), pointer                               :: trian
  type(container_element), dimension(nvmax), intent(in) :: List_of_non_Direchlet_Elements
  
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
  
  call N_Matr%init(nonzero, number_of_non_Direchlet_elements)
  N_Matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  N_Matr%ja(1:nonzero) = ja(1:nonzero)
  N_Matr%a(1:nonzero) = a(1:nonzero)

end subroutine N_assembling_resuid

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



! N matrix assembling

subroutine N_assembling(N_Matr, x_y_fir, x_y_sec, xi_ind, eta_ind, List_of_non_Direchlet_Elements, &
        number_of_non_Direchlet_elements)
  
  implicit none 
        
  type(Matrix)                             , intent(out)  :: N_Matr
  integer                                  , intent(in)   :: number_of_non_Direchlet_elements
  integer                                  , intent(in)   :: x_y_fir, x_y_sec
  real*8                                   , intent(in)   :: xi_ind, eta_ind
  type(container_element), dimension(nvmax), intent(in)   :: List_of_non_Direchlet_Elements
  real*8                                                  :: a(maxnz), prom
  integer                                                 :: ia(maxn), ja(maxnz)
  integer                                                 :: i, j, k, nonzero, nonzero_raw
  type(Element), pointer                                  :: elem_master, elem_side
  type(Triangle), pointer                                 :: trian
  
  ia(1) = 1
  k = 1
  
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
  
    nonzero_raw = 0
  
    elem_master => List_of_non_Direchlet_Elements(i)%pointing_element
  
    do j = 1, elem_master%number_of_neighbour_elements    
    
      elem_side => elem_master%neighbour_elements_list(j)%pointing_element
      
      if (elem_side%on_boundary .eqv. .false.) then
      
          prom = scalar_mult_B(elem_master, elem_side, x_y_fir, x_y_sec, xi_ind, eta_ind)
           
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
  
  call N_Matr%init(nonzero, number_of_non_Direchlet_elements)
  N_Matr%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  N_Matr%ja(1:nonzero) = ja(1:nonzero)
  N_Matr%a(1:nonzero) = a(1:nonzero)

end subroutine N_assembling



subroutine plot_matr(List_of_Elements, number_of_elements, mat)


   type(Element), target, dimension(nvmax) :: List_of_Elements
   integer :: i, j, number_of_elements
   type(Matrix) :: mat
   
   open(2,file='/home/chuck/Desktop/2D_ice_model/Article/VP_method_3/matrix.txt')
   
   write (2, *) mat%matr_size
   
   write (2, *)
   
   do i = 1, (mat%matr_size + 1)
           
     write (2, *) mat%ia(i)
        
   end do
   
   write (2, *)
   
   do i = 1, mat%nonzero
           
     write (2, *) mat%ja(i)
        
   end do
   
   write (2, *)
   
   do i = 1, mat%nonzero
           
     write (2, *) mat%a(i)
        
   end do
   
  
   close(2)
  

end subroutine plot_matr  



function frobenius_norm(mat) result(frob_result)

   integer :: i, j
   real*8 :: maxim, prom, frob_result
   type(Matrix) :: mat
   
   maxim = -1d10
   
   do i = 1, mat%matr_size
   
   prom = 0d0
   
   do j = mat%ia(i), (mat%ia(i+1)-1)
   
     prom = prom + abs(mat%a(j))
   
   end do
   
   if (prom > maxim) then
   
     maxim = prom
   
   end if
   
   end do
  
   frob_result = maxim

end function frobenius_norm  


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

end module module_dynamics 
