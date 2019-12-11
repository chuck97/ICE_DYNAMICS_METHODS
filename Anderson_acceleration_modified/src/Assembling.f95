module module_assembling

  use module_Classes
  use module_constants
  use module_numerical_integration
  use module_grid_values
  use module_matricies
  

  implicit none
  
  private
  
  public :: transport_mass_matrix_high_order_assembling, &    
            transport_mass_matrix_low_order_assembling , &    
            transport_right_hand_low_order_assembling, &      
            transport_right_hand_high_order_assembling, &     
            non_Direchlet_mass_matrix_assembling, &           
            non_Direchlet_mass_matrix_lumped_assembling, &    
            non_Direchlet_mass_matrix, &                      
            non_Direchlet_mass_matrix_lumped, &               
            transport_mass_matrix_high_order, &               
            transport_mass_matrix_low_order, &                
            transport_right_hand_high_order, &                
            transport_right_hand_low_order, &                 
            right_hands_distruction                           
            
            
  type(Matrix)        :: non_Direchlet_mass_matrix
  type(Matrix)        :: non_Direchlet_mass_matrix_lumped
  type(Matrix)        :: transport_mass_matrix_low_order
  type(Matrix)        :: transport_mass_matrix_high_order
  
  real*8, allocatable :: transport_right_hand_high_order(:)
  real*8, allocatable :: transport_right_hand_low_order(:)
            
  contains
  
  !Assembling non Direchlet mass matrix
  subroutine non_Direchlet_mass_matrix_assembling()
  
    implicit none
  
    !!local variables:
    integer                     :: i, j, k, r, nonzero_raw
    real*8, allocatable         :: a(:)
    integer, allocatable        :: ia(:), ja(:)
    integer                     :: nonzero
     
    allocate(a(max_number_of_diagonals*number_of_non_Direchlet_elements))
    allocate(ia(number_of_non_Direchlet_elements+1))
    allocate(ja(max_number_of_diagonals*number_of_non_Direchlet_elements))
     
    ia(1) = 1
    
    r = 1
    
    nonzero_raw = 0
    nonzero = 0
    
    do i = 1, number_of_non_Direchlet_elements
      do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
        if (List_of_non_Direchlet_Elements(i)%pointing_element% &
        neighbour_elements_list(j)%pointing_element%on_boundary .eqv. .false.) then
          a(r) = a(r) + & 
          scalar_multiplication(phi_times_phi, List_of_non_Direchlet_Elements(i)%pointing_element% &
          neighbour_elements_list(j)%pointing_element, List_of_non_Direchlet_Elements(i)%pointing_element) 
        end if
                    
        if (abs(a(r)) > varepsilon) then
          ja(r) = List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_elements_list(j)% &
          pointing_element%nd_identificator
          nonzero_raw = nonzero_raw + 1
          r = r + 1
        end if
      end do
      ia(i+1) = ia(i) + nonzero_raw
      nonzero = nonzero + nonzero_raw
      nonzero_raw = 0
    end do
     
    call non_Direchlet_mass_matrix%init(nonzero, &
           number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
    non_Direchlet_mass_matrix%ia(1:(non_Direchlet_mass_matrix%matr_size_y+1)) = &
            ia(1:(non_Direchlet_mass_matrix%matr_size_y+1))
    non_Direchlet_mass_matrix%ja(1:non_Direchlet_mass_matrix%nonzero) = &
            ja(1:non_Direchlet_mass_matrix%nonzero)
    non_Direchlet_mass_matrix%a(1:non_Direchlet_mass_matrix%nonzero) = &
            a(1:non_Direchlet_mass_matrix%nonzero)
            
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
     
  end subroutine non_Direchlet_mass_matrix_assembling

  
  !Assembling high order mass matrix for transport equation
  subroutine transport_mass_matrix_high_order_assembling(time_step)
   
    implicit none
    
    !!arguments:
    real*8,intent(in)       :: time_step
    
    !!local variables:
    integer                 :: i, j, k, r, nonzero_raw
    real*8, allocatable     :: a(:)
    integer, allocatable    :: ia(:), ja(:)
    integer, allocatable    :: nonzero
    
    allocate(a(max_number_of_diagonals*number_of_elements))
    allocate(ia(number_of_elements+1))
    allocate(ja(max_number_of_diagonals*number_of_elements))
    
    ia(1) = 1
    
    r = 1
    
    nonzero_raw = 0
    nonzero = 0
    
    do i = 1, number_of_elements
       
      do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        a(r) = a(r) + & 
        scalar_multiplication(phi_times_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
             List_of_Elements(i))
             
             
        if (abs(a(r)) > varepsilon) then
          ja(r) = List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%identificator
          nonzero_raw = nonzero_raw + 1
          a(r) = a(r)
          r = r + 1     
        end if  
      end do
      ia(i+1) = ia(i) + nonzero_raw
      nonzero = nonzero + nonzero_raw
      nonzero_raw = 0
    end do
    
    call transport_mass_matrix_high_order%init(nonzero, number_of_elements, number_of_elements)
    transport_mass_matrix_high_order%ia(1:(transport_mass_matrix_high_order%matr_size_y+1)) = &
     ia(1:(transport_mass_matrix_high_order%matr_size_y+1))
    transport_mass_matrix_high_order%ja(1:transport_mass_matrix_high_order%nonzero) = &
     ja(1:transport_mass_matrix_high_order%nonzero)
    transport_mass_matrix_high_order%a(1:transport_mass_matrix_high_order%nonzero) = &
     a(1:transport_mass_matrix_high_order%nonzero)
    
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
    
  end subroutine transport_mass_matrix_high_order_assembling


  !Assembling low order mass matrix for transport equation
  subroutine transport_mass_matrix_low_order_assembling(time_step)
  
    implicit none
  
    !!arguments:
    real*8,intent(in)      :: time_step
    
    !!local variables:
    integer                :: i, j, k, r, nonzero_raw
    real*8, allocatable    :: a(:)
    integer, allocatable   :: ia(:), ja(:)
    integer                :: nonzero
    
    allocate(a(max_number_of_diagonals*number_of_elements))
    allocate(ia(number_of_elements+1))
    allocate(ja(max_number_of_diagonals*number_of_elements))
    
    ia(1) = 1
    
    r = 1
    
    nonzero_raw = 0
    nonzero = 0
    
    do i = 1, number_of_elements 
      do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        a(r) = a(r) + & 
        scalar_multiplication(phi_times_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i)) 
      end do
      
      if (abs(a(r)) > varepsilon) then
           
        ja(r) = List_of_Elements(i)%identificator
        nonzero_raw = nonzero_raw + 1
        r = r + 1
             
      end if
         
      ia(i+1) = ia(i) + nonzero_raw
      nonzero = nonzero + nonzero_raw
      nonzero_raw = 0
      
    end do
    
    call transport_mass_matrix_low_order%init(nonzero, number_of_elements, number_of_elements)
    transport_mass_matrix_low_order%ia(1:(transport_mass_matrix_low_order%matr_size_y+1)) = &
     ia(1:(transport_mass_matrix_low_order%matr_size_y+1))
    transport_mass_matrix_low_order%ja(1:transport_mass_matrix_low_order%nonzero) = &
     ja(1:transport_mass_matrix_low_order%nonzero)
    transport_mass_matrix_low_order%a(1:transport_mass_matrix_low_order%nonzero) = &
      a(1:transport_mass_matrix_low_order%nonzero)
    
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
    
    
  end subroutine transport_mass_matrix_low_order_assembling


  !Assembling high order right hand for transport equation 
  subroutine transport_right_hand_high_order_assembling(time_step)
  
    implicit none
    
    !!arguments:
    real*8, intent(in)    :: time_step
   
    !!local variables:
    integer              :: i, j, k
    
    allocate(transport_right_hand_high_order(number_of_elements))
        
    do i = 1, number_of_elements
      transport_right_hand_high_order(i) = 0d0     
    end do
    
    do i = 1, number_of_elements
      do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        transport_right_hand_high_order(i) = transport_right_hand_high_order(i) + &
        time_step*List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value* &
        scalar_multiplication(u_phi_times_grad_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i)) - & 
        5d-1*(time_step**2)*List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value* &
        scalar_multiplication(u_nabla_u_phi_nabla_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i))
         
         transport_right_hand_high_order(i) = transport_right_hand_high_order(i)
         
      end do
    end do
    
    !print *,  time_step*scalar_multiplication(u_phi_times_grad_phi, List_of_Elements(6)%neighbour_elements_list(1)%pointing_element, &
    !     List_of_Elements(6))
    
  end subroutine transport_right_hand_high_order_assembling


  !Assembling low order right hand for transport equation 
  subroutine transport_right_hand_low_order_assembling(time_step)
   
    implicit none 
    
    !!arguments:
    real*8,intent(in)   :: time_step
    
    !!local variables:
    integer             :: i, j, k
    
    allocate(transport_right_hand_low_order(number_of_elements))
    
    do i = 1, number_of_elements
      transport_right_hand_low_order(i) = 0d0     
    end do
   
    do i = 1, number_of_elements  
      do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        transport_right_hand_low_order(i) = &
         transport_right_hand_low_order(i) + &
         List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value*( &
         c_d*scalar_multiplication(phi_times_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i))+ & 
         time_step*scalar_multiplication(u_phi_times_grad_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i)) - & 
         5d-1*(time_step**2)* &
         scalar_multiplication(u_nabla_u_phi_nabla_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
         List_of_Elements(i)))
      end do
      transport_right_hand_low_order(i) = transport_right_hand_low_order(i) - &
       List_of_Elements(i)%element_value*c_d* &
       transport_mass_matrix_low_order%a(List_of_Elements(i)%identificator)     
    end do    
  
  end subroutine transport_right_hand_low_order_assembling 


  !Right hands distruction
  subroutine right_hands_distruction()
    
    implicit none
    
    deallocate(transport_right_hand_low_order)
    deallocate(transport_right_hand_high_order)
  
  end subroutine right_hands_distruction


  !Assembling non Direchlet lumped mass matrix
  subroutine non_Direchlet_mass_matrix_lumped_assembling()
  
    implicit none
    
    call non_Direchlet_mass_matrix_lumped%init(non_Direchlet_mass_matrix%nonzero, &
     non_Direchlet_mass_matrix%matr_size_x, non_Direchlet_mass_matrix%matr_size_y)  
    call Matr1_eqv_Matr2(non_Direchlet_mass_matrix_lumped, non_Direchlet_mass_matrix)
    call lumping_square(non_Direchlet_mass_matrix_lumped)
     
  end subroutine non_Direchlet_mass_matrix_lumped_assembling
  
end module module_assembling
