module module_flux_correction

  use module_Classes
  use module_grid_values
  use module_numerical_integration
  use module_constants

  implicit none
  
  private
  
  public :: flux_calculation, P_calculation, u_calculation, alpha_calculation, &
   flux_correction_procedure
  
  
  contains
  
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine flux_calculation(List_of_Triangles, number_of_triangles, solution_vec_C, mass_matrix_low_a)
  
    type(Triangle), target, dimension(ntmax) :: List_of_Triangles
    real*8 :: solution_vec_C(maxn)
    integer :: number_of_triangles, i, j, k, coef(2, 3)
    real*8 :: ux(3), numerator, denominator
    real*8 :: mass_matrix_low_a(maxnz)
    
    coef(1,1) = 2
    coef(2,1) = 3
    coef(1,2) = 1
    coef(2,2) = 3
    coef(1,3) = 1
    coef(2,3) = 2
    
    do k = 1, number_of_triangles
    
      do i = 1, 3
    
        ux(i) = ((c_d - 1d0)*List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element%element_value + &
          solution_vec_C(List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element%identificator))  
    
      end do
      
      
      
      do i = 1, 3
      
        numerator = 0d0
        denominator = 0d0
        
        denominator = mass_matrix_low_a(List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element%identificator)
        
        numerator = (scalar_multiplication_on_triangle(phi_times_phi, &
        List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element, &
        List_of_Triangles(k)%neighbour_elements_list(coef(1, i))%pointing_element, &
        List_of_Triangles(k)) + &
        scalar_multiplication_on_triangle(phi_times_phi, &
        List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element, &
        List_of_Triangles(k)%neighbour_elements_list(coef(2, i))%pointing_element, &
        List_of_Triangles(k)))*ux(i) - &
        scalar_multiplication_on_triangle(phi_times_phi, &
        List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element, &
        List_of_Triangles(k)%neighbour_elements_list(coef(1, i))%pointing_element, &
        List_of_Triangles(k))*ux(coef(1, i)) - &
        scalar_multiplication_on_triangle(phi_times_phi, &
        List_of_Triangles(k)%neighbour_elements_list(i)%pointing_element, &
        List_of_Triangles(k)%neighbour_elements_list(coef(2, i))%pointing_element, &
        List_of_Triangles(k))*ux(coef(2, i))
        
        
        
        List_of_Triangles(k)%flux_to_element(i) = numerator/denominator
        
      end do
       
      
    end do
    
  
  end subroutine flux_calculation
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine P_calculation(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements)
  
    type(Triangle), target, dimension(ntmax) :: List_of_Triangles
    type(Element), target, dimension(nvmax) :: List_of_Elements
    integer :: number_of_triangles, number_of_elements, i, j, k
    
    do i = 1, number_of_elements
    
      List_of_Elements(i)%P_pos = 0d0
      List_of_Elements(i)%P_neg = 0d0
    
    end do
    
    do i = 1, number_of_triangles
    
      do j = 1, 3
      
        if(List_of_Triangles(i)%flux_to_element(j) > varepsilon) then
        
          List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%P_pos = &
          List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%P_pos + &
          List_of_Triangles(i)%flux_to_element(j)
        
        else
        
          List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%P_neg = &
          List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%P_neg + &
          List_of_Triangles(i)%flux_to_element(j)
        
        end if
      
      end do
    
    end do
  
  end subroutine P_calculation
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 subroutine u_calculation(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements, solution_vec_L)
 
   type(Triangle), target, dimension(ntmax) :: List_of_Triangles
   type(Element), target, dimension(nvmax) :: List_of_Elements
   integer :: number_of_triangles, number_of_elements, i, j, k
   real*8 :: solution_vec_L(maxn)
   real*8 :: u_1_max, u_2_max, u_3_max, u_1_min, u_2_min, u_3_min
   real*8 :: u_max_el, u_min_el
   
   do i = 1, number_of_triangles
   
     !u_1_max = max(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator))
     !u_2_max = max(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator))
     !u_3_max = max(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator))
     
     !u_1_min = min(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator))
     !u_2_min = min(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator))
     !u_3_min = min(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%element_value, &
     !solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator))
     
     u_1_max = max(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator))
     u_2_max = max(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator))
     u_3_max = max(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator))
     
     u_1_min = min(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator))
     u_2_min = min(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator))
     u_3_min = min(solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator), &
     solution_vec_L(List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator))
     
   
     List_of_Triangles(i)%u_max_triangle = max(max(u_1_max, u_2_max), u_3_max)
     List_of_Triangles(i)%u_min_triangle = min(min(u_1_min, u_2_min), u_3_min)
   
   end do
   
   do i = 1, number_of_elements
   
     u_max_el = -1d10
     u_min_el = 1d10
   
     do j = 1, List_of_Elements(i)%number_of_neighbour_triangles
     
       if(List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle%u_max_triangle > &
        u_max_el)  then
       
         u_max_el = List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle%u_max_triangle
          
       end if 
       
       if(List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle%u_min_triangle < &
        u_min_el)  then
       
         u_min_el = List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle%u_min_triangle
          
       end if 
     
     end do
     
     List_of_Elements(i)%u_max_element = u_max_el  
     List_of_Elements(i)%u_min_element = u_min_el
     
     List_of_Elements(i)%Q_plus = List_of_Elements(i)%u_max_element - solution_vec_L(List_of_Elements(i)%identificator)
     List_of_Elements(i)%Q_minus = List_of_Elements(i)%u_min_element - solution_vec_L(List_of_Elements(i)%identificator)
   
   
     if(abs(List_of_Elements(i)%P_pos) > varepsilon) then
     
       List_of_Elements(i)%R_plus = min(1d0, List_of_Elements(i)%Q_plus/List_of_Elements(i)%P_pos)
       
     else
       
       List_of_Elements(i)%R_plus = 0d0
       
     end if  
     
     
     if(abs(List_of_Elements(i)%P_neg) > varepsilon) then
     
       List_of_Elements(i)%R_minus = min(1d0, List_of_Elements(i)%Q_minus/List_of_Elements(i)%P_neg)
       
     else
       
       List_of_Elements(i)%R_minus = 0d0
       
     end if  
     
     
   end do
   
   
   
 end subroutine u_calculation
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 subroutine alpha_calculation(List_of_Triangles, number_of_triangles)
 
   type(Triangle), target, dimension(ntmax) :: List_of_Triangles
   integer :: number_of_triangles, i, j, k
   real*8 :: R1, R2, R3
   
   do i = 1, number_of_triangles
   
     if(List_of_Triangles(i)%flux_to_element(1) > varepsilon) then
       
        R1 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%R_plus
           
     else
          
        R1 = List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%R_minus
    
     end if
     
     if(List_of_Triangles(i)%flux_to_element(2) > varepsilon) then
       
        R2 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%R_plus
           
     else
          
        R2 = List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%R_minus
    
     end if
          
     
     if(List_of_Triangles(i)%flux_to_element(3) > varepsilon) then
       
        R3 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%R_plus
           
     else
          
        R3 = List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%R_minus
    
     end if     
     
     List_of_Triangles(i)%alpha_correction = min(min(R1, R2), R3)
     
     if (List_of_Triangles(i)%alpha_correction < -varepsilon) then
     
       print *, "alpha: ", List_of_Triangles(i)%alpha_correction
     
     end if
             
   end do
   
 
 end subroutine alpha_calculation
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 subroutine flux_correction_procedure(List_of_Triangles, List_of_Elements, number_of_triangles, &
  number_of_elements, mass_matrix_low_a, solution_vec_C, solution_vec_L)
     
   type(Triangle), target, dimension(ntmax) :: List_of_Triangles
   type(Element), target, dimension(nvmax) :: List_of_Elements
   integer :: number_of_triangles, number_of_elements, i, j, k
   real*8 :: solution_vec_C(maxn), solution_vec_L(maxn), prom
   real*8 :: mass_matrix_low_a(maxnz)
   
   call flux_calculation(List_of_Triangles, number_of_triangles, solution_vec_C, mass_matrix_low_a)
 
   call P_calculation(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements)
      
   call u_calculation(List_of_Triangles, List_of_Elements, number_of_triangles, number_of_elements, solution_vec_L)
      
   call alpha_calculation(List_of_Triangles, number_of_triangles)
   
   do i = 1, number_of_elements
   
     List_of_Elements(i)%element_value = solution_vec_L(List_of_Elements(i)%identificator)
   
   end do
   
   do i = 1, number_of_triangles
   
     do j = 1, 3
     
       List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%element_value = &
       List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%element_value + &
       List_of_Triangles(i)%alpha_correction*List_of_Triangles(i)%flux_to_element(j)
       
     
     end do
   
   end do
 
 end subroutine flux_correction_procedure
 
 
  
  
  



end module module_flux_correction
