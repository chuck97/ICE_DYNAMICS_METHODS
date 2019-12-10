module module_assembling

  use module_Classes
  use module_constants
  use module_numerical_integration
  

  implicit none
  
  private
  
  public :: mass_matrix_high_order_assembling , mass_matrix_low_order_assembling , &
  transport_right_hand_low_order_assembling, transport_right_hand_high_order_assembling, N_xx_assembling, &
  N_xy_assembling, N_yx_assembling, N_yy_assembling, non_Direchlet_mass_matrix_assembling
  
  contains
  
  
  
  
!Assembling high order mass matrix for transport equation
  
subroutine mass_matrix_high_order_assembling(Matr, List_of_Elements, &
 number_of_elements, time_step)

  real*8 :: time_step
  integer :: i, j, k, r, nonzero_raw, number_of_elements
  type(Element), target, dimension(nvmax) :: List_of_Elements
  type(Matrix) :: Matr
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
  
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
  
  call Matr%init(nonzero, number_of_elements)
  Matr%ia = ia
  Matr%ja = ja
  Matr%a =  a
  
end subroutine mass_matrix_high_order_assembling


!Assembling low order mass matrix for transport equation

subroutine mass_matrix_low_order_assembling(Matr, List_of_Elements, number_of_elements, time_step)

  real*8 :: time_step
  integer :: i, j, k, r, nonzero_raw, number_of_elements
  type(Element), target, dimension(nvmax) :: List_of_Elements
  type(Matrix) :: Matr
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
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
  
  call Matr%init(nonzero, number_of_elements)
  Matr%ia = ia
  Matr%ja = ja
  Matr%a =  a
  
  
end subroutine mass_matrix_low_order_assembling




!Assembling high order right hand for transport equation 

subroutine transport_right_hand_high_order_assembling(Right_hand, List_of_Elements, number_of_elements, time_step)

  real*8 :: Right_hand(maxn), time_step
  integer :: i, j, k, number_of_elements
  type(Element), target, dimension(nvmax) :: List_of_Elements
  

      
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
      
  real*8 :: Right_hand(maxn), time_step
  integer :: i, j, k, number_of_elements
  type(Element), target, dimension(nvmax) :: List_of_Elements
  type(Matrix) :: Matr
  
  
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


! subroutine N_xx assembling
 
subroutine N_xx_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_res)
 
  type(container_element), target, dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements
  integer :: i, j, k, r
  type(Matrix) :: Matr_res
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero, nonzero_raw
  
  ia(1) = 1
  r = 1
  
  nonzero_raw = 0
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
    
    
    do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
       
      a(r) = a(r) + & 
      scalar_multiplication(phi_x_times_phi_x, List_of_non_Direchlet_Elements(i)%pointing_element, &
       List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_elements_list(j)%pointing_element)
         
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
  
  call Matr_res%init(nonzero, number_of_non_Direchlet_elements)
  Matr_res%ia = ia
  Matr_res%ja = ja
  Matr_res% a =  a
  
 
end subroutine N_xx_assembling      


! subroutine N_xy assembling
 
subroutine N_xy_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_res)
 
  type(container_element), target, dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements
  integer :: i, j, k, r, nonzero_raw
  type(Matrix) :: Matr_res
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
  ia(1) = 1
  
  r = 1
  
  nonzero_raw = 0
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
    
    
    do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
       
      a(r) = a(r) + & 
      scalar_multiplication(phi_x_times_phi_y, List_of_non_Direchlet_Elements(i)%pointing_element, &
       List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_elements_list(j)%pointing_element)
         
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
  
  call Matr_res%init(nonzero, number_of_non_Direchlet_elements)
  Matr_res%ia = ia
  Matr_res%ja = ja
  Matr_res% a =  a
  
 
end subroutine N_xy_assembling 


! subroutine N_yx assembling
 
subroutine N_yx_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_res)
 
  type(container_element), target, dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements
  integer ::  i, j, k, r, nonzero_raw
  type(Matrix) :: Matr_res
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
  ia(1) = 1
  
  r = 1
  
  nonzero_raw = 0
  Matr_res%nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
    
    
    do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
       
      a(r) = a(r) + & 
      scalar_multiplication(phi_y_times_phi_x, List_of_non_Direchlet_Elements(i)%pointing_element, &
       List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_elements_list(j)%pointing_element)
         
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
  
  call Matr_res%init(nonzero, number_of_non_Direchlet_elements)
  Matr_res%ia = ia
  Matr_res%ja = ja
  Matr_res% a =  a
 
end subroutine N_yx_assembling 

! subroutine N_yy assembling
 
subroutine N_yy_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_res)
 
  type(container_element), target, dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements
  integer :: i, j, k, r, nonzero_raw
  type(Matrix) :: Matr_res
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
  ia(1) = 1
  
  r = 1
  
  nonzero_raw = 0
  nonzero = 0
  
  do i = 1, number_of_non_Direchlet_elements
    
    
    do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
       
      a(r) = a(r) + & 
      scalar_multiplication(phi_y_times_phi_y, List_of_non_Direchlet_Elements(i)%pointing_element, &
       List_of_non_Direchlet_Elements(i)%pointing_element%neighbour_elements_list(j)%pointing_element)
         
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
  
  call Matr_res%init(nonzero, number_of_non_Direchlet_elements)
  Matr_res%ia = ia
  Matr_res%ja = ja
  Matr_res% a =  a
  
 
end subroutine N_yy_assembling 


subroutine non_Direchlet_mass_matrix_assembling(List_of_non_Direchlet_Elements, number_of_non_Direchlet_elements, &
 Matr_mass)

  type(container_element), target, dimension(nvmax) :: List_of_non_Direchlet_Elements
  integer :: number_of_non_Direchlet_elements, i, j, k, r, nonzero_raw
  type(Matrix) :: Matr_mass
  real*8 :: a(maxnz)
  integer :: ia(maxn), ja(maxnz)
  integer :: nonzero
  
  ! Assembling non Direchlet mass matrix
   
   ia(1) = 1
  
   r = 1
  
   nonzero_raw = 0
   nonzero = 0
  
   do i = 1, number_of_non_Direchlet_elements
    
     do j = 1, List_of_non_Direchlet_Elements(i)%pointing_element%number_of_neighbour_elements
       
       if (List_of_non_Direchlet_Elements(i)%pointing_element% &
       neighbour_elements_list(j)%pointing_element%on_boundary .eqv. .false.) then
       
         a(r) = scalar_multiplication(phi_times_phi, List_of_non_Direchlet_Elements(i)%pointing_element% &
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
   
  call Matr_mass%init(nonzero, number_of_non_Direchlet_elements)
  Matr_mass%ia(1:(number_of_non_Direchlet_elements+1)) = ia(1:(number_of_non_Direchlet_elements+1))
  Matr_mass%ja(1:nonzero) = ja(1:nonzero)
  Matr_mass%a(1:nonzero) =  a(1:nonzero)
  
end subroutine non_Direchlet_mass_matrix_assembling




end module module_assembling
