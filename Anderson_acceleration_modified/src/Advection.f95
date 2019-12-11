module module_advection

  use module_constants
  use module_classes
  use module_linear_system_solver
  use module_flux_correction
  use module_assembling
  use module_grid_values
  use module_dynamics
  use module_matricies
  
  implicit none
  
  private 
  public  :: scalar_advection
  
  contains 
  
  subroutine scalar_advection(time_step, ind) ! ind = 1 -- mass,  2 -- concentration, 3 -- height
  
    implicit none
    
    !! arguments:
    integer, intent(in)             :: ind
    real*8, intent(in)              :: time_step
    
    !! local variables:
    integer             :: i
    real*8, allocatable :: solution_vec_L(:), solution_vec_C(:)
    
    
    if (ind == 1) then
       do i = 1, number_of_elements  
         List_of_Elements(i)%element_value = List_of_Elements(i)%m
       end do 
    else
      if (ind == 2) then
        do i = 1, number_of_elements  
          List_of_Elements(i)%element_value = List_of_Elements(i)%A
        end do 
      else
        do i = 1, number_of_elements  
          List_of_Elements(i)%element_value = List_of_Elements(i)%h
        end do
      end if
    end if
    
    allocate(solution_vec_L(number_of_elements))
    allocate(solution_vec_C(number_of_elements))
    
    do i = 1, number_of_elements
      solution_vec_L(i) = 0d0 
      solution_vec_C(i) = 0d0
    end do
    
    call transport_right_hand_low_order_assembling(time_step)
    call transport_right_hand_high_order_assembling(time_step)
    
    
    
    !! low order solution
    do i = 1, number_of_elements
      solution_vec_L(i) = transport_right_hand_low_order(i)/ &
       transport_mass_matrix_low_order%a(i) 
    end do
    
    do i = 1, number_of_elements
      solution_vec_L(i) = solution_vec_L(i) + List_of_Elements(i)%element_value
    end do
    
    print *, "fr norm:", frobenius_norm(transport_mass_matrix_high_order)   
    
    
    
    !! high order solution
    call ILU_2_BICG_solver(transport_mass_matrix_high_order, &
                           transport_right_hand_high_order, &
                           number_of_elements, &
                           solution_vec_C)      
                           
    do i = 1, number_of_elements
      solution_vec_C(i) = solution_vec_C(i) + List_of_Elements(i)%element_value
    end do                       
                           
    !! flux correction update
    call flux_correction_procedure(solution_vec_C, solution_vec_L)
    
    !! value update
    if (ind == 1) then
       do i = 1, number_of_elements  
         List_of_Elements(i)%m = List_of_Elements(i)%element_value
       end do 
    else
      if (ind == 2) then
        do i = 1, number_of_elements  
          if (List_of_Elements(i)%element_value < 1d0) then
            List_of_Elements(i)%A = List_of_Elements(i)%element_value
          else
            List_of_Elements(i)%A = 1d0
          end if  
        end do 
      else
        do i = 1, number_of_elements  
          List_of_Elements(i)%h = List_of_Elements(i)%element_value
        end do
      end if
    end if                       
    
    
    !! deallocating arrays
    deallocate(solution_vec_L)
    deallocate(solution_vec_C)
    call right_hands_distruction()
  
  end subroutine scalar_advection


end module module_advection
