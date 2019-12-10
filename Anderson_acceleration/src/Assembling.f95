module module_assembling

  use module_Classes
  use module_constants
  use module_numerical_integration
  

  implicit none
  
  private
  
  public :: mass_matrix_high_order_assembling , mass_matrix_low_order_assembling , &
  transport_right_hand_low_order_assembling, transport_right_hand_high_order_assembling, non_Direchlet_mass_matrix_assembling
  
  contains
  
  
  
  
!Assembling high order mass matrix for transport equation
  
subroutine mass_matrix_high_order_assembling(Matr, List_of_Elements, &
 number_of_elements, time_step)
 
  implicit none

  real*8                                  ,intent(in)    :: time_step
  integer                                                :: i, j, k, r, nonzero_raw
  type(Element), target, dimension(nvmax) ,intent(in)    :: List_of_Elements
  integer                                 ,intent(in)    :: number_of_elements
  type(Matrix)                            ,intent(inout) :: Matr
  real*8                                                 :: a(maxnz)
  integer                                                :: ia(maxn), ja(maxnz)
  integer                                                :: nonzero
  
  
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
        r = r + 1
           
      end if  
             
    end do
       
    ia(i+1) = ia(i) + nonzero_raw
    nonzero = nonzero + nonzero_raw
    nonzero_raw = 0
    
  end do
  
  call Matr%init(nonzero, number_of_elements, number_of_elements)
  Matr%ia = ia
  Matr%ja = ja
  Matr%a  =  a
  
end subroutine mass_matrix_high_order_assembling


!Assembling low order mass matrix for transport equation

subroutine mass_matrix_low_order_assembling(Matr, List_of_Elements, number_of_elements, time_step)

  implicit none

  real*8                                  ,intent(in)    :: time_step
  integer                                 ,intent(in)    :: number_of_elements
  integer                                                :: i, j, k, r, nonzero_raw
  type(Element), target, dimension(nvmax) ,intent(in)    :: List_of_Elements
  type(Matrix)                            ,intent(inout) :: Matr
  real*8                                                 :: a(maxnz)
  integer                                                :: ia(maxn), ja(maxnz)
  integer                                                :: nonzero
  
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
  
  call Matr%init(nonzero, number_of_elements, number_of_elements)
  Matr%ia = ia
  Matr%ja = ja
  Matr%a =  a
  
  
end subroutine mass_matrix_low_order_assembling




!Assembling high order right hand for transport equation 

subroutine transport_right_hand_high_order_assembling(Right_hand, List_of_Elements, number_of_elements, time_step)

  implicit none

  real*8                                   ,intent(inout) :: Right_hand(maxn)
  real*8                                   ,intent(in)    :: time_step
  integer                                  ,intent(in)    :: number_of_elements
  integer                                                 :: i, j, k
  type(Element), target, dimension(nvmax)  ,intent(in)    :: List_of_Elements
  

      
  do i = 1, maxn
    Right_hand(i) = 0d0     
  end do
  
      
      
  do i = 1, number_of_elements
      
    do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        
      Right_hand(i) = Right_hand(i) + &
      time_step*List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value* &
      scalar_multiplication(u_phi_times_grad_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
       List_of_Elements(i)) - & 
      5d-1*(time_step**2)*List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value* &
      scalar_multiplication(u_nabla_u_phi_nabla_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
       List_of_Elements(i))
        
    end do
       
  end do
  
end subroutine transport_right_hand_high_order_assembling



!Assembling low order right hand for transport equation 

subroutine transport_right_hand_low_order_assembling(Right_hand, Matr, List_of_Elements, &
 number_of_elements, time_step)
 
  implicit none 
      
  real*8                                  ,intent(inout) :: Right_hand(maxn)
  real*8                                  ,intent(in)    :: time_step
  integer                                 ,intent(in)    :: number_of_elements
  integer                                                :: i, j, k
  type(Element), target, dimension(nvmax) ,intent(in)    :: List_of_Elements
  type(Matrix)                            ,intent(in)    :: Matr
  
  
  do i = 1, maxn
    Right_hand(i) = 0d0     
  end do
  
      
      
  do i = 1, number_of_elements
      
    do j = 1, List_of_Elements(i)%number_of_neighbour_elements
        
      Right_hand(i) = Right_hand(i) + List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%element_value*( &
      c_d*scalar_multiplication(phi_times_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
      List_of_Elements(i))+ & 
      time_step*scalar_multiplication(u_phi_times_grad_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
      List_of_Elements(i)) - & 
      5d-1*(time_step**2)* &
      scalar_multiplication(u_nabla_u_phi_nabla_phi, List_of_Elements(i)%neighbour_elements_list(j)%pointing_element, &
      List_of_Elements(i)))
        
    end do
        
    Right_hand(i) = Right_hand(i) - List_of_Elements(i)%element_value*c_d* &
    Matr%a(List_of_Elements(i)%identificator)
        
        
  end do    

end subroutine transport_right_hand_low_order_assembling 



subroutine non_Direchlet_mass_matrix_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_mass)

  implicit none

  type(container_element), target, dimension(nvmax)  ,intent(in)    :: List_of_non_Direchlet_Elements
  integer                                            ,intent(in)    :: number_of_non_Direchlet_elements
  integer                                                           :: i, j, k, r, nonzero_raw
  type(Matrix)                                       ,intent(inout) :: Matr_mass
  real*8, allocatable                                               :: a(:)
  integer, allocatable                                              :: ia(:), ja(:)
  integer                                                           :: nonzero
  
  ! Assembling non Direchlet mass matrix
   
   allocate(a(15*number_of_non_Direchlet_elements))
   allocate(ia(number_of_non_Direchlet_elements+1))
   allocate(ja(15*number_of_non_Direchlet_elements))
   
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
   
  call Matr_mass%init(nonzero, number_of_non_Direchlet_elements, number_of_non_Direchlet_elements)
  Matr_mass%ia(1:(Matr_mass%matr_size_y+1)) = ia(1:(Matr_mass%matr_size_y+1))
  Matr_mass%ja(1:Matr_mass%nonzero) = ja(1:Matr_mass%nonzero)
  Matr_mass%a(1:Matr_mass%nonzero) =  a(1:Matr_mass%nonzero)
  
  deallocate(a)
  deallocate(ia)
  deallocate(ja)
  
   
end subroutine non_Direchlet_mass_matrix_assembling


end module module_assembling
