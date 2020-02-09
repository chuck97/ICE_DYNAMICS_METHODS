module module_residual

  use module_classes
  use module_grid_values
  use module_assembling
  use module_constants
  use module_numerical_integration
  use module_matricies
  use module_dynamics
  
  implicit none        
                     
  public ::          R_assembling,          &
                     N_assembling,          &
                     S_assembling,          &
                     R_x_matrix,            &
                     R_y_matrix,            &  
                     S_matrix,              &
                     N_x_old_matrix,        &
                     N_y_old_matrix,        &
                     N_new_distruction,     &
                     residual,              &
                     N_old_distruction,     &
                     N_x_new_matrix,        &
                     N_y_new_matrix 
                     
  
  type(Matrix)  ::   R_x_matrix, R_y_matrix
  type(Matrix)  ::   N_x_new_matrix, N_y_new_matrix
  type(Matrix)  ::   N_x_old_matrix, N_y_old_matrix
  type(Matrix)  ::   S_matrix
  
  contains
  
  
  subroutine S_assembling()
  
    !!local variables:
    integer  :: i
    
    call S_matrix%init(number_of_triangles, number_of_triangles, number_of_triangles)
    
    S_matrix%ia(1) = 1
    
    do i = 1, number_of_triangles
      S_matrix%ia(i+1) = S_matrix%ia(i) + 1
      S_matrix%ja(i) = i
      S_matrix%a(i) = List_of_Triangles(i)%size_of_triangle
    end do
  
  end subroutine S_assembling
  
  
  subroutine R_assembling(ind_xy)    !! ind_xy = 1 -- R_x, ind_xy = 2 -- R_y
    
    implicit none
    
    !!arguments:
    integer, intent(in)  :: ind_xy
    
    !!local variables:
    integer              :: i, j, r, nonzero_raw
    real*8, allocatable  :: a(:)
    integer, allocatable :: ia(:), ja(:)
    real*8               :: coefficients(3), siz
    
    
    allocate(a(max_number_of_diagonals*number_of_non_Direchlet_elements))
    allocate(ia(number_of_non_Direchlet_elements+1))
    allocate(ja(max_number_of_diagonals*number_of_non_Direchlet_elements))
    
    ia(1) = 1
    r = 1
    
    if (ind_xy == 1) then
    
      do i = 1, number_of_non_Direchlet_elements
        nonzero_raw = 0
        do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_triangles 
          coefficients = coefficients_calculation(List_of_non_Direchlet_Elements(i)%pointing_element, &
                List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle)
          siz = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%size_of_triangle
          a(r) = siz*coefficients(1)
          ja(r) = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%identificator
          r = r + 1
          nonzero_raw = nonzero_raw + 1
        end do
        ia(i+1) = ia(i) + nonzero_raw
      end do
      
      call R_x_matrix%init(r-1, number_of_triangles, number_of_non_Direchlet_elements)
      R_x_matrix%ia = ia
      R_x_matrix%ja = ja
      R_x_matrix%a = a
    
    else
      
      do i = 1, number_of_non_Direchlet_elements
        nonzero_raw = 0
        do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_triangles 
          coefficients = coefficients_calculation(List_of_non_Direchlet_Elements(i)%pointing_element, &
                List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle)
          siz = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%size_of_triangle
          a(r) = siz*coefficients(2)
          ja(r) = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%identificator
          r = r + 1
          nonzero_raw = nonzero_raw + 1
        end do
        ia(i+1) = ia(i) + nonzero_raw
      end do
      
      call R_y_matrix%init(r-1, number_of_triangles, number_of_non_Direchlet_elements)
      R_y_matrix%ia = ia
      R_y_matrix%ja = ja
      R_y_matrix%a = a
      
    end if
    
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
  
  end subroutine R_assembling
  
  
  subroutine N_assembling(ind_xy, ind_new_old)   !! ind_xy = 1 -- N_x, ind_xy = 2 -- N_y ; ind_new_old = 1 -- N_new, ind_new_old = 2 -- N_old 
  
    implicit none
    
    !!arguments:
    integer, intent(in)    :: ind_xy, ind_new_old
    
    !!local variables:
    integer              :: i, j, r, nonzero_raw
    real*8, allocatable  :: a(:)
    integer, allocatable :: ia(:), ja(:)
    real*8               :: coefficients(3), siz, P_0, delta, valu
    
    
    allocate(a(3*number_of_triangles))
    allocate(ia(number_of_triangles+1))
    allocate(ja(3*number_of_triangles))
    
    ia(1) = 1
    r = 1
    
    if(ind_xy == 1) then
      
      if(ind_new_old == 1) then
      
        do i = 1, number_of_Triangles
          nonzero_raw = 0
          do j = 1, 3 
            coefficients = coefficients_calculation(List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element, &
                  List_of_Triangles(i))
            siz = List_of_Triangles(i)%size_of_triangle
            P_0 = List_of_Triangles(i)%P_0
            delta = List_of_Triangles(i)%delta
            valu = siz*coefficients(1)/(delta_min + delta)
            if (abs(valu) < varepsilon) then
              a(r) = 0d0
            else
              a(r) = valu
            end if     
            ja(r) = List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%identificator
            r = r + 1
            nonzero_raw = nonzero_raw + 1
          end do
          ia(i+1) = ia(i) + nonzero_raw
        end do
      
        call N_x_new_matrix%init(r-1, number_of_non_Direchlet_elements, number_of_triangles)
        N_x_new_matrix%ia = ia
        N_x_new_matrix%ja = ja
        N_x_new_matrix%a = a
      
      else
      
        do i = 1, number_of_Triangles
          nonzero_raw = 0
          do j = 1, 3 
            coefficients = coefficients_calculation(List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element, &
                  List_of_Triangles(i))
            siz = List_of_Triangles(i)%size_of_triangle
            P_0 = List_of_Triangles(i)%P_0
            delta = List_of_Triangles(i)%delta_old
            valu = siz*coefficients(1)/(delta_min + delta)
            if (abs(valu) < varepsilon) then
              a(r) = 0d0
            else
              a(r) = valu
            end if   
            ja(r) = List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%identificator
            r = r + 1
            nonzero_raw = nonzero_raw + 1
          end do
          ia(i+1) = ia(i) + nonzero_raw
        end do
      
        call N_x_old_matrix%init(r-1, number_of_non_Direchlet_elements, number_of_triangles)
        N_x_old_matrix%ia = ia
        N_x_old_matrix%ja = ja
        N_x_old_matrix%a = a
      
      end if
      
    else
    
      if(ind_new_old == 1) then
      
        do i = 1, number_of_Triangles
          nonzero_raw = 0
          do j = 1, 3 
            coefficients = coefficients_calculation(List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element, &
                  List_of_Triangles(i))
            siz = List_of_Triangles(i)%size_of_triangle
            P_0 = List_of_Triangles(i)%P_0
            delta = List_of_Triangles(i)%delta
            valu = siz*coefficients(2)/(delta_min + delta)
            if (abs(valu) < varepsilon) then
              a(r) = 0d0
            else
              a(r) = valu
            end if   
            ja(r) = List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%identificator
            r = r + 1
            nonzero_raw = nonzero_raw + 1
          end do
          ia(i+1) = ia(i) + nonzero_raw
        end do
      
        call N_y_new_matrix%init(r-1, number_of_non_Direchlet_elements, number_of_triangles)
        N_y_new_matrix%ia = ia
        N_y_new_matrix%ja = ja
        N_y_new_matrix%a = a
      
      else
      
        do i = 1, number_of_Triangles
          nonzero_raw = 0
          do j = 1, 3 
            coefficients = coefficients_calculation(List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element, &
                  List_of_Triangles(i))
            siz = List_of_Triangles(i)%size_of_triangle
            P_0 = List_of_Triangles(i)%P_0
            delta = List_of_Triangles(i)%delta_old
            valu = siz*coefficients(2)/(delta_min + delta)
            if (abs(valu) < varepsilon) then
              a(r) = 0d0
            else
              a(r) = valu
            end if   
            ja(r) = List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%identificator
            r = r + 1
            nonzero_raw = nonzero_raw + 1
          end do
          ia(i+1) = ia(i) + nonzero_raw
        end do
      
        call N_y_old_matrix%init(r-1, number_of_non_Direchlet_elements, number_of_triangles)
        N_y_old_matrix%ia = ia
        N_y_old_matrix%ja = ja
        N_y_old_matrix%a = a
      
      end if
    
    end if
    
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
    
  
  end subroutine N_assembling
  
  
  function residual(old_new, is_rhs, time_step, square_size) result(residual_vector)  ! old_new = 1 - old, old_new = 2 - new, is_rhs = 1 - yes, is_rhs = 2 - no 
  
    implicit none
    
    !arguments:
    integer, intent(in)   :: old_new, is_rhs
    real*8, intent(in)    :: time_step, square_size
    real*8                :: residual_vector(2*number_of_non_Direchlet_elements + 3*number_of_triangles)
    
    !local variables:
    real*8, allocatable   :: u_vector(:), v_vector(:),                                                   &
                             sigma1_vector(:), sigma2_vector(:), sigma12_vector(:),                      &
                             u_res(:), v_res(:), sigma1_res(:), sigma2_res(:), sigma12_res(:),           &
                             prom_nd1(:), prom_nd2(:), prom_tr1(:), prom_tr2(:),                         &
                             m_vector(:), m_lumped_vector(:), a_vector(:), C_w_vector(:),                &
                             u_w_vector(:), v_w_vector(:), u_a_vector(:), v_a_vector(:),                 &
                             u_old_vector(:), v_old_vector(:),                                           &
                             u_rhs(:), v_rhs(:), sigma1_rhs(:), sigma2_rhs(:), sigma12_rhs(:),           &
                             dH_dx_list(:), dH_dy_list(:)
    real*8                :: dh_dx, dh_dy, size_t, prom_x, prom_y    
    integer               :: i, j                   
    
    allocate(u_vector(number_of_non_Direchlet_elements))
    allocate(v_vector(number_of_non_Direchlet_elements))
    allocate(u_res(number_of_non_Direchlet_elements))
    allocate(v_res(number_of_non_Direchlet_elements))
    allocate(sigma1_vector(number_of_triangles))
    allocate(sigma2_vector(number_of_triangles))
    allocate(sigma12_vector(number_of_triangles))
    allocate(sigma1_res(number_of_triangles))
    allocate(sigma2_res(number_of_triangles))
    allocate(sigma12_res(number_of_triangles))
    allocate(prom_nd1(number_of_non_Direchlet_elements))
    allocate(prom_nd2(number_of_non_Direchlet_elements))
    allocate(prom_tr1(number_of_triangles))
    allocate(prom_tr2(number_of_triangles))
    allocate(m_vector(number_of_non_Direchlet_elements))
    allocate(m_lumped_vector(number_of_non_Direchlet_elements))
    allocate(a_vector(number_of_non_Direchlet_elements))
    allocate(C_w_vector(number_of_non_Direchlet_elements))
    allocate(u_w_vector(number_of_non_Direchlet_elements))
    allocate(v_w_vector(number_of_non_Direchlet_elements))
    allocate(u_a_vector(number_of_non_Direchlet_elements))
    allocate(v_a_vector(number_of_non_Direchlet_elements))
    allocate(u_old_vector(number_of_non_Direchlet_elements))
    allocate(v_old_vector(number_of_non_Direchlet_elements))
    allocate(u_rhs(number_of_non_Direchlet_elements))
    allocate(v_rhs(number_of_non_Direchlet_elements))
    allocate(sigma1_rhs(number_of_triangles))
    allocate(sigma2_rhs(number_of_triangles))
    allocate(sigma12_rhs(number_of_triangles))
    allocate(dH_dx_list(number_of_non_Direchlet_elements))
    allocate(dH_dy_list(number_of_non_Direchlet_elements))
    
    
    do i = 1, number_of_non_Direchlet_elements
      u_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(1)
      v_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(2)
      u_res(i) = 0d0
      v_res(i) = 0d0
      prom_nd1(i) = 0d0
      prom_nd2(i) = 0d0
      m_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%m
      m_lumped_vector(i) = non_Direchlet_mass_matrix_lumped%a(i)
      a_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%a
      u_w_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1)
      v_w_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2)
      u_a_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1)
      v_a_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2)
      u_old_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
      v_old_vector(i) = List_of_non_Direchlet_Elements(i)%pointing_element%u(2)
      
      prom_x = 0d0
      prom_y = 0d0
      do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_triangles
        dh_dx = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%dH_dx
        dh_dy = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%dH_dy
        size_t = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_triangles_list(j)%pointing_triangle%size_of_triangle
        prom_x = prom_x + (1d0/3d0)*size_t*dh_dx
        prom_y = prom_y + (1d0/3d0)*size_t*dh_dy
      end do
      
      dH_dx_list(i) = prom_x
      dH_dy_list(i) = prom_y
     
    end do
    
    do i = 1, number_of_triangles
      sigma1_vector(i) = List_of_Triangles(i)%sigma_resid(1)/List_of_Triangles(i)%P_0
      sigma2_vector(i) = List_of_Triangles(i)%sigma_resid(2)/List_of_Triangles(i)%P_0
      sigma12_vector(i) = List_of_Triangles(i)%sigma_resid(3)/List_of_Triangles(i)%P_0 
      sigma1_res(i) = 0d0
      sigma2_res(i) = 0d0
      sigma12_res(i) = 0d0
      prom_tr1(i) = 0d0
      prom_tr2(i) = 0d0
    end do
    
    if (old_new == 1) then
    
      !assembling C_w vector 
      
      do i = 1, number_of_non_Direchlet_elements
        C_w_vector(i) = a_vector(i)*C_w*rho_water*dsqrt( &
        (List_of_non_Direchlet_Elements(i)%pointing_element%u(1) - u_w_vector(i))**2 + &
        (List_of_non_Direchlet_Elements(i)%pointing_element%u(2) - v_w_vector(i))**2)
      end do
      
      !!!!!!!!!!!!!!
      ! u_residual !
      !!!!!!!!!!!!!!
      
      !(m + dt*Cw)*M*u
      
      prom_nd1 = vector_mult_constant(C_w_vector, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(m_vector, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_vector, number_of_non_Direchlet_elements, u_res)
      
      !-2*dt*w_e*m*M*v
      prom_nd1 = vector_mult_constant(m_vector, number_of_non_Direchlet_elements, -2d0*time_step*omega_e)
      call vectors_elementwise_multiplication(prom_nd1,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2, v_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_sum(u_res, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
     
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd2(i)
      end do
    
      ! + (1/2)*dt*Rx*sigma1
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma1_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      ! + (1/2)*dt*Rx*sigma2
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma2_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      ! + dt*Ry*sigma12
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma12_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      !!!!!!!!!!!!!!
      ! v_residual !
      !!!!!!!!!!!!!!
      
      !2*dt*w_e*m*M*u
      prom_nd1 = vector_mult_constant(m_vector, number_of_non_Direchlet_elements, 2d0*time_step*omega_e)
      call vectors_elementwise_multiplication(prom_nd1,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2, u_vector, number_of_non_Direchlet_elements, prom_nd1)
     
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+(m + dt*Cw)*M*v
      
      prom_nd1 = vector_mult_constant(C_w_vector, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(m_vector, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+(1/2)*dt*Ry*sigma1
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma1_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !-(1/2)*dt*Ry*sigma2
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma2_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -5d-1*time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+dt*Rx*sigma12
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma12_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!
      ! sigma1_residual !
      !!!!!!!!!!!!!!!!!!!
      
      ! -N_x*u
      call sparce_matrix_vector_multiplication(N_x_old_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0)
      
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr2(i)
      end do
      
      ! - N_y*v
      call sparce_matrix_vector_multiplication(N_y_old_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0)
      call vectors_sum(prom_tr2, sigma1_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr1(i)
      end do
      
      ! + S*sigma1
      call sparce_matrix_vector_multiplication(S_matrix, sigma1_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma1_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr2(i)
      end do
      
      !!!!!!!!!!!!!!!!!!!
      ! sigma2_residual !
      !!!!!!!!!!!!!!!!!!!
      
      !-(1/e**2)*N_x*u
      call sparce_matrix_vector_multiplication(N_x_old_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/e**2)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr2(i)
      end do
      
      !+(1/e**2)*N_y*v
      call sparce_matrix_vector_multiplication(N_y_old_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, 1d0/e**2)
      call vectors_sum(prom_tr2, sigma2_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr1(i)
      end do
      
      !+S*sigma2
      call sparce_matrix_vector_multiplication(S_matrix, sigma2_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma2_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr2(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!!
      ! sigma12_residual !
      !!!!!!!!!!!!!!!!!!!!
      
      
      !-(1/2e**2)*N_y*u
      call sparce_matrix_vector_multiplication(N_y_old_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/(2d0*e**2))
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr2(i)
      end do
      
      !-(1/2e**2)*N_x*v
      call sparce_matrix_vector_multiplication(N_x_old_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/(2d0*e**2))
      call vectors_sum(prom_tr2, sigma12_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr1(i)
      end do
      
      !+S*sigma12
      call sparce_matrix_vector_multiplication(S_matrix, sigma12_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma12_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr2(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         u rhs            !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      !+m*M*u_old
      
      call vectors_elementwise_multiplication(m_lumped_vector, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_old_vector, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !+dt*a*C_a*M*u_a
      
      call vectors_elementwise_multiplication(m_lumped_vector, a_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_a_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step*C_a)
      call vectors_sum(prom_nd1, u_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !dt*C_w*M*u_w
      call  vectors_elementwise_multiplication(m_lumped_vector, C_w_vector, number_of_non_Direchlet_elements, prom_nd1)
      call  vectors_elementwise_multiplication(prom_nd1, u_w_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd1, u_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !-dt*m*g*(dh/dx)
      call vectors_elementwise_multiplication(dH_dx_list, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -1d0*time_step*g_gravity)
      call vectors_sum(prom_nd2, u_rhs, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd1(i)
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_non_Direchlet_elements
          u_rhs(i) = 0d0
        end do
      end if
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         v rhs            !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !+m*M*v_old
      
      call vectors_elementwise_multiplication(m_lumped_vector, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_old_vector, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !+dt*a*C_a*M*v_a
      
      call vectors_elementwise_multiplication(m_lumped_vector, a_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_a_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step*C_a)
      call vectors_sum(prom_nd1, v_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !dt*C_w*M*v_w
      call  vectors_elementwise_multiplication(m_lumped_vector, C_w_vector, number_of_non_Direchlet_elements, prom_nd1)
      call  vectors_elementwise_multiplication(prom_nd1, v_w_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd1, v_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !-dt*m*g*(dh/dy)
      call vectors_elementwise_multiplication(dH_dy_list, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -1d0*time_step*g_gravity)
      call vectors_sum(prom_nd2, v_rhs, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd1(i)
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_non_Direchlet_elements
          v_rhs(i) = 0d0
        end do
      end if
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         sigma1 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma1_rhs(i) = -(List_of_Triangles(i)%delta_old/ &
        (List_of_Triangles(i)%delta_old + delta_min))*List_of_Triangles(i)%size_of_triangle
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma1_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         sigma2 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma2_rhs(i) = 0d0
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma2_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !        sigma12 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma12_rhs(i) = 0d0
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma12_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Final old residual         !!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_non_Direchlet_elements
        residual_vector(i) = (u_res(i) - u_rhs(i))/(square_size**2)
      end do
      
      do i = 1, number_of_non_Direchlet_elements
        residual_vector(number_of_non_Direchlet_elements + i) = (v_res(i) - v_rhs(i))/(square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + i) = (sigma1_res(i) - sigma1_rhs(i))/ &
             (List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + number_of_triangles + i) = &
           (sigma2_res(i) - sigma2_rhs(i))/(List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + 2*number_of_triangles + i) = &
           (sigma12_res(i) - sigma12_rhs(i))/(List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
      
   
      
    else
    
      !assembling C_w vector 
      
      do i = 1, number_of_non_Direchlet_elements
        C_w_vector(i) = a_vector(i)*C_w*rho_water*dsqrt( &
        (u_vector(i) - u_w_vector(i))**2 + &
        (v_vector(i) - v_w_vector(i))**2)
      end do
      
      !!!!!!!!!!!!!!
      ! u_residual !
      !!!!!!!!!!!!!!
      
      !(m + dt*Cw)*M*u
      
      prom_nd1 = vector_mult_constant(C_w_vector, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(m_vector, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_vector, number_of_non_Direchlet_elements, u_res)
      
      !-2*dt*w_e*m*M*v
      prom_nd1 = vector_mult_constant(m_vector, number_of_non_Direchlet_elements, -2d0*time_step*omega_e)
      call vectors_elementwise_multiplication(prom_nd1,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2, v_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_sum(u_res, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
     
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd2(i)
      end do
    
      ! + (1/2)*dt*Rx*sigma1
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma1_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      ! + (1/2)*dt*Rx*sigma2
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma2_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      ! + dt*Ry*sigma12
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma12_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, u_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_res(i) = prom_nd1(i)
      end do
      
      !!!!!!!!!!!!!!
      ! v_residual !
      !!!!!!!!!!!!!!
      
      !2*dt*w_e*m*M*u
      prom_nd1 = vector_mult_constant(m_vector, number_of_non_Direchlet_elements, 2d0*time_step*omega_e)
      call vectors_elementwise_multiplication(prom_nd1,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2, u_vector, number_of_non_Direchlet_elements, prom_nd1)
     
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+(m + dt*Cw)*M*v
      
      prom_nd1 = vector_mult_constant(C_w_vector, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(m_vector, prom_nd1, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_elementwise_multiplication(prom_nd2,  m_lumped_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_vector, number_of_non_Direchlet_elements, prom_nd2)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+(1/2)*dt*Ry*sigma1
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma1_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, 5d-1*time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !-(1/2)*dt*Ry*sigma2
      call sparce_matrix_vector_multiplication(R_y_matrix, sigma2_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -5d-1*time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      !+dt*Rx*sigma12
      call sparce_matrix_vector_multiplication(R_x_matrix, sigma12_vector, number_of_triangles, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd2, v_res, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_res(i) = prom_nd1(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!
      ! sigma1_residual !
      !!!!!!!!!!!!!!!!!!!
      
      ! -N_x*u
      call sparce_matrix_vector_multiplication(N_x_new_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0)
 
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr2(i)
      end do
      
      ! - N_y*v
      call sparce_matrix_vector_multiplication(N_y_new_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0)
      call vectors_sum(prom_tr2, sigma1_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr1(i)
      end do
      
      ! + S*sigma1
      call sparce_matrix_vector_multiplication(S_matrix, sigma1_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma1_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma1_res(i) = prom_tr2(i)
      end do
      
      !!!!!!!!!!!!!!!!!!!
      ! sigma2_residual !
      !!!!!!!!!!!!!!!!!!!
      
      !-(1/e**2)*N_x*u
      call sparce_matrix_vector_multiplication(N_x_new_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/e**2)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr2(i)
      end do
      
      !+(1/e**2)*N_y*v
      call sparce_matrix_vector_multiplication(N_y_new_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, 1d0/e**2)
      call vectors_sum(prom_tr2, sigma2_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr1(i)
      end do
      
      !+S*sigma2
      call sparce_matrix_vector_multiplication(S_matrix, sigma2_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma2_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma2_res(i) = prom_tr2(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!!
      ! sigma12_residual !
      !!!!!!!!!!!!!!!!!!!!
      
      
      !-(1/2e**2)*N_y*u
      call sparce_matrix_vector_multiplication(N_y_new_matrix, u_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/(2d0*e**2))
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr2(i)
      end do
      
      !-(1/2e**2)*N_x*v
      call sparce_matrix_vector_multiplication(N_x_new_matrix, v_vector, number_of_non_Direchlet_elements, prom_tr1)
      prom_tr2 = vector_mult_constant(prom_tr1, number_of_triangles, -1d0/(2d0*e**2))
      call vectors_sum(prom_tr2, sigma12_res, number_of_triangles, prom_tr1)
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr1(i)
      end do
      
      !+S*sigma12
      call sparce_matrix_vector_multiplication(S_matrix, sigma12_vector, number_of_triangles, prom_tr1)
      call vectors_sum(prom_tr1, sigma12_res, number_of_triangles, prom_tr2)
      
      do i = 1, number_of_triangles
        sigma12_res(i) = prom_tr2(i)
      end do
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         u rhs            !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !+m*M*u_old
      
      call vectors_elementwise_multiplication(m_lumped_vector, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_old_vector, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !+dt*a*C_a*M*u_a
      
      call vectors_elementwise_multiplication(m_lumped_vector, a_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, u_a_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step*C_a)
      call vectors_sum(prom_nd1, u_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !dt*C_w*M*u_w
      call  vectors_elementwise_multiplication(m_lumped_vector, C_w_vector, number_of_non_Direchlet_elements, prom_nd1)
      call  vectors_elementwise_multiplication(prom_nd1, u_w_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd1, u_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd2(i)
      end do
      
      !-dt*m*g*(dh/dx)
      call vectors_elementwise_multiplication(dH_dx_list, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -1d0*time_step*g_gravity)
      call vectors_sum(prom_nd2, u_rhs, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        u_rhs(i) = prom_nd1(i)
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_non_Direchlet_elements
          u_rhs(i) = 0d0
        end do
      end if
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         v rhs            !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !+m*M*v_old
      
      call vectors_elementwise_multiplication(m_lumped_vector, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_old_vector, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !+dt*a*C_a*M*v_a
      
      call vectors_elementwise_multiplication(m_lumped_vector, a_vector, number_of_non_Direchlet_elements, prom_nd1)
      call vectors_elementwise_multiplication(prom_nd1, v_a_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step*C_a)
      call vectors_sum(prom_nd1, v_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !dt*C_w*M*v_w
      call  vectors_elementwise_multiplication(m_lumped_vector, C_w_vector, number_of_non_Direchlet_elements, prom_nd1)
      call  vectors_elementwise_multiplication(prom_nd1, v_w_vector, number_of_non_Direchlet_elements, prom_nd2)
      prom_nd1 = vector_mult_constant(prom_nd2, number_of_non_Direchlet_elements, time_step)
      call vectors_sum(prom_nd1, v_rhs, number_of_non_Direchlet_elements, prom_nd2)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd2(i)
      end do
      
      !-dt*m*g*(dh/dy)
      call vectors_elementwise_multiplication(dH_dy_list, m_vector, number_of_non_Direchlet_elements, prom_nd1)
      prom_nd2 = vector_mult_constant(prom_nd1, number_of_non_Direchlet_elements, -1d0*time_step*g_gravity)
      call vectors_sum(prom_nd2, u_rhs, number_of_non_Direchlet_elements, prom_nd1)
      
      do i = 1, number_of_non_Direchlet_elements
        v_rhs(i) = prom_nd1(i)
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_non_Direchlet_elements
          v_rhs(i) = 0d0
        end do
      end if
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         sigma1 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma1_rhs(i) = -(List_of_Triangles(i)%delta_old/ &
        (List_of_Triangles(i)%delta + delta_min))*List_of_Triangles(i)%size_of_triangle
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma1_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         sigma2 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma2_rhs(i) = 0d0
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma2_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !        sigma12 rhs       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_triangles
        sigma12_rhs(i) = 0d0
      end do
      
      if (is_rhs == 2) then
        do i = 1, number_of_triangles
          sigma12_rhs(i) = 0d0
        end do
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!          Final new residual         !!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_non_Direchlet_elements
        residual_vector(i) = (u_res(i) - u_rhs(i))/(square_size**2)
      end do
      
      do i = 1, number_of_non_Direchlet_elements 
        residual_vector(number_of_non_Direchlet_elements + i) = (v_res(i) - v_rhs(i))/(square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + i) = (sigma1_res(i) - sigma1_rhs(i))/ &
           (List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + number_of_triangles + i) = &
           (sigma2_res(i) - sigma2_rhs(i))/(List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
      
      do i = 1, number_of_triangles 
        residual_vector(2*number_of_non_Direchlet_elements + 2*number_of_triangles + i) = &
           (sigma12_res(i) - sigma12_rhs(i))/(List_of_triangles(i)%size_of_triangle*square_size**2)
      end do
    
    end if
    
    
    deallocate(u_vector)
    deallocate(v_vector)
    deallocate(u_res)
    deallocate(v_res)
    deallocate(sigma1_vector)
    deallocate(sigma2_vector)
    deallocate(sigma12_vector)
    deallocate(sigma1_res)
    deallocate(sigma2_res)
    deallocate(sigma12_res)
    deallocate(prom_nd1)
    deallocate(prom_nd2)
    deallocate(prom_tr1)
    deallocate(prom_tr2)
    deallocate(m_vector)
    deallocate(m_lumped_vector)
    deallocate(a_vector)
    deallocate(C_w_vector)
    deallocate(u_w_vector)
    deallocate(u_a_vector)
    deallocate(v_w_vector)
    deallocate(v_a_vector)
    deallocate(u_old_vector)
    deallocate(v_old_vector)
    deallocate(u_rhs)
    deallocate(v_rhs)
    deallocate(sigma1_rhs)
    deallocate(sigma2_rhs)
    deallocate(sigma12_rhs)
    deallocate(dH_dx_list)
    deallocate(dH_dy_list)
  
  end function residual
  
  
  
  subroutine N_new_distruction()
  
    implicit none
    
    call N_x_new_matrix%distr()
    call N_y_new_matrix%distr()
  
  end subroutine N_new_distruction
  
  
  
  subroutine N_old_distruction()
  
    implicit none
    
    call N_x_old_matrix%distr()
    call N_y_old_matrix%distr()
    
  end subroutine N_old_distruction
  
  
    

end module module_residual
