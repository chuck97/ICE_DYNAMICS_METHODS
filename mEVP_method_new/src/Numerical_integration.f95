module module_numerical_integration

  use module_classes
  use module_constants
  use module_grid_values
  
  implicit none
  
  private
  
  public :: from_baricentric_to_ordinary, integral_over_triangle, phi_times_phi, u_phi_times_grad_phi, &
  gradphi_times_gradphi, scalar_multiplication, coefficients_calculation, linear_on_triangle, &
  sigma_times_gradphi, sigma_grad_phi_scalar_multiplication, u_nabla_u_phi_nabla_phi, &
  L2_mass, linear_average_on_triangle, scalar_multiplication_on_triangle, u_coefficients_calculation, &
  squared_linear_on_triangle, resud_u_coefficients_calculation, phi_x_times_phi_x, phi_x_times_phi_y, &
  phi_y_times_phi_x, phi_y_times_phi_y, scalar_mult_two_vec

  
  ! Declare procedure interface
  interface
  
  
    function mult_phi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
      real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v 
      real*8  :: return_value
    end function mult_phi
  
  end interface
  
  contains 
  
  
  ! phi*phi
  
  function phi_times_phi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (a1*x + b1*y + c1)*(a2*x + b2*y + c2) 
    
  end function phi_times_phi
  
  !phi_x*phi_x
  
  function phi_x_times_phi_x(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (a1)*(a2) 
    
  end function phi_x_times_phi_x
  
  !phi_x*phi_y
  
  function phi_x_times_phi_y(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (a1)*(b2) 
    
  end function phi_x_times_phi_y
  
  !phi_y*phi_x
  
  function phi_y_times_phi_x(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (b1)*(a2) 
    
  end function phi_y_times_phi_x
  
  !phi_y*phi_y
  
  function phi_y_times_phi_y(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (b1)*(b2) 
    
  end function phi_y_times_phi_y
  
  ! linear function on triangle
  
  function linear_on_triangle(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = a1*x + b1*y + c1
    
  end function linear_on_triangle
  
  ! squared linear on triangle
  
  function squared_linear_on_triangle(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (a1*x + b1*y + c1)**2
    
  end function squared_linear_on_triangle
  
  
  ! u phi * grad(phi)    
  
  function u_phi_times_grad_phi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = (a_u*x + b_u*y + c_u)*(a1*x + b1*y + c1)*a2 + &
    (a_v*x + b_v*y + c_v)*(a1*x + b1*y + c1)*b2 
    
  end function u_phi_times_grad_phi
  
  
  ! grad(phi)*grad(phi)
  
  function gradphi_times_gradphi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result (return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8  :: return_value
    return_value = a1*a2 + b1*b2 
    
  end function gradphi_times_gradphi
  
  
  ! sigma*grad(phi)
  
  function sigma_times_gradphi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result(return_value)
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8 :: return_value
    
    return_value = a1*a2 + b1*b2
  
  end function sigma_times_gradphi
  
  ! ((u nabla*(u phi)), nabla phi) 
  
  function u_nabla_u_phi_nabla_phi(x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v) result(return_value)
  
    implicit none
  
    real*8, intent(in) :: x, y, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    real*8 :: return_value
    
    return_value = (a_u*x + b_u*y + c_u)*(a_u*x + b_u*y + c_u)*a1*a2 + &
     (a_u*x + b_u*y + c_u)*(a_v*x + b_v*y + c_v)*b1*a2 + &
     (a_u*x + b_u*y + c_u)*(a_v*x + b_v*y + c_v)*a1*b2 + &
     (a_v*x + b_v*y + c_v)*(a_v*x + b_v*y + c_v)*b1*b2 + &
     (a1*x + b1*y + c1)*(a_u + b_v)*(a2*(a_u*x + b_u*y + c_u) + b2*(a_v*x + b_v*y + c_v))
     
  
  end function u_nabla_u_phi_nabla_phi
  
  
  
  ! Calculation linear average of m, h or A on triangle

function linear_average_on_triangle(tri, ident) result(avg_tri_res) ! 1 -- m,  2 -- h, 3 -- A

  implicit none

  type(Triangle), target :: tri
  integer :: ident
  real*8 :: avg_tri_res
  real*8 :: x0, y0, h0
  real*8 :: x1, y1, h1
  real*8 :: x2, y2, h2
  real*8 :: a, b, c
  real*8 :: integr_tr
  
  x0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(1)
  y0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(2)
  
  
  x1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(1)
  y1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(2)
  
  
  x2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(1)
  y2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(2)
   
  
  if (ident == 1) then
  
    h0 = tri%neighbour_elements_list(1)%pointing_element%m
    h1 = tri%neighbour_elements_list(2)%pointing_element%m
    h2 = tri%neighbour_elements_list(3)%pointing_element%m
  
  end if
  
  
  if (ident == 2) then
  
    h0 = tri%neighbour_elements_list(1)%pointing_element%h
    h1 = tri%neighbour_elements_list(2)%pointing_element%h
    h2 = tri%neighbour_elements_list(3)%pointing_element%h
  
  end if
  
  if (ident == 3) then
  
    h0 = tri%neighbour_elements_list(1)%pointing_element%A
    h1 = tri%neighbour_elements_list(2)%pointing_element%A
    h2 = tri%neighbour_elements_list(3)%pointing_element%A
  
  end if
  

  a = ((h0 - h1)*(y0 - y2) - (y0 - y1)*(h0 - h2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  b = ((h0 - h2)*(x0 - x1) - (h0 - h1)*(x0 - x2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  c = (x0*(h2*y1 - h1*y2) + x1*(h0*y2 - h2*y0) + x2*(h1*y0 - h0*y1))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  
  avg_tri_res = integral_over_triangle(x0, x1, x2, y0, y1, y2, linear_on_triangle, a, b, c, 0d0, 0d0, 0d0, 0d0, 0d0, &
  0d0, 0d0, 0d0, 0d0)/tri%size_of_triangle

end function linear_average_on_triangle
  
  
  
  !! convertation from baricentric to ordinary coordinates a1, a2, a3 of triangle with verticies (x1, y1) , (x2, y2), (x3, y3)
  
  function from_baricentric_to_ordinary(a1, a2, a3, x1, x2, x3, y1, y2, y3) result(array_x_y)
  
    implicit none
    
    real*8, intent(in) :: a1, a2, a3, x1, x2, x3, y1, y2, y3
    real*8, dimension(2) :: array_x_y
    
    array_x_y(1) = x1*a1 + x2*a2 + x3*a3
    array_x_y(2) = y1*a1 + y2*a2 + y3*a3
    
  
  end function from_baricentric_to_ordinary 
  
  
  
  
  !!  cubature formula for function f over triangle with verticies (x1, y1), (x2, y2), (x3, y3)
  
  function integral_over_triangle(x1, x2, x3, y1, y2, y3, f, a1, b1, c1, a2, b2, c2, &
    a_u, b_u, c_u, a_v, b_v, c_v) &
    result(numeric_integral)
    
    implicit none
  
    real*8, intent(in) :: x1, x2, x3, y1, y2, y3, a1, b1, c1, a2, b2, c2, a_u, b_u, c_u, a_v, b_v, c_v
    procedure(mult_phi)  :: f
    real*8 :: numeric_integral
    real*8 :: intermediate_result
    real*8 :: three_ver_1_baricentric, three_ver_2_baricentric, three_weight
    real*8 :: three_ver_ordinary_1_x, three_ver_ordinary_1_y, three_ver_ordinary_2_x, three_ver_ordinary_2_y, &
         three_ver_ordinary_3_x, three_ver_ordinary_3_y
    real*8 :: six_ver_1_baricentric, six_ver_2_baricentric, six_ver_3_baricentric, &
     six_weight
    real*8 :: six_ver_ordinary_1_x, six_ver_ordinary_1_y, six_ver_ordinary_2_x, six_ver_ordinary_2_y,  &
     six_ver_ordinary_3_x, six_ver_ordinary_3_y
    real*8 :: six_ver_ordinary_4_x, six_ver_ordinary_4_y, six_ver_ordinary_5_x, six_ver_ordinary_5_y,  &
     six_ver_ordinary_6_x, six_ver_ordinary_6_y
    real*8 :: a, b, c, p, S
    real*8, dimension(2) :: three_ver_ordinary_1, three_ver_ordinary_2, three_ver_ordinary_3
    real*8, dimension(2) :: six_ver_ordinary_1, six_ver_ordinary_2, six_ver_ordinary_3, &
     six_ver_ordinary_4, six_ver_ordinary_5, six_ver_ordinary_6
    
    intermediate_result = 0d0
    numeric_integral = 0d0
    
    ! baricentric verticies with multiplicity 3 and their weight
    
    three_ver_2_baricentric = 0.124949503233232d0
    three_ver_1_baricentric = 0.437525248383384d0
    three_weight = 0.205950504760887d0
    
    ! convertation to ordinary coordinates
    
    
    
    three_ver_ordinary_1 = from_baricentric_to_ordinary(three_ver_2_baricentric, three_ver_1_baricentric, &
     three_ver_1_baricentric, x1, x2, x3, y1, y2, y3) 
    
    three_ver_ordinary_1_x = three_ver_ordinary_1(1)
    three_ver_ordinary_1_y = three_ver_ordinary_1(2)
    
    three_ver_ordinary_2 = from_baricentric_to_ordinary(three_ver_1_baricentric, three_ver_2_baricentric, &
     three_ver_1_baricentric, x1, x2, x3, y1, y2, y3)
     
    three_ver_ordinary_2_x = three_ver_ordinary_2(1)
    three_ver_ordinary_2_y = three_ver_ordinary_2(2)
    
    three_ver_ordinary_3 = from_baricentric_to_ordinary(three_ver_1_baricentric, three_ver_1_baricentric, &
     three_ver_2_baricentric, x1, x2, x3, y1, y2, y3) 
    
    three_ver_ordinary_3_x = three_ver_ordinary_3(1)
    three_ver_ordinary_3_y = three_ver_ordinary_3(2)
    
    ! calculation the contribution of 3-multiplicity points
    
    intermediate_result = intermediate_result + three_weight*f(three_ver_ordinary_1_x, &
     three_ver_ordinary_1_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + three_weight*f(three_ver_ordinary_2_x, &
     three_ver_ordinary_2_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + three_weight*f(three_ver_ordinary_3_x, &
     three_ver_ordinary_3_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    
    
    
    ! baricentric verticies with multiplicity 6 and their weight
    
    six_ver_1_baricentric = 0.797112651860071d0
    six_ver_2_baricentric = 0.165409927389841d0
    six_ver_3_baricentric = 0.037477420750088d0
    six_weight = 0.063691414286223d0
    
    ! convertation to ordinary coordinates
    
    
    
    six_ver_ordinary_1 = from_baricentric_to_ordinary(six_ver_1_baricentric, six_ver_2_baricentric, &
     six_ver_3_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_1_x = six_ver_ordinary_1(1)
    six_ver_ordinary_1_y = six_ver_ordinary_1(2)
    
    six_ver_ordinary_2 = from_baricentric_to_ordinary(six_ver_1_baricentric, six_ver_3_baricentric, &
     six_ver_2_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_2_x = six_ver_ordinary_2(1)
    six_ver_ordinary_2_y = six_ver_ordinary_2(2)
    
    
    six_ver_ordinary_3 = from_baricentric_to_ordinary(six_ver_2_baricentric, six_ver_1_baricentric, &
     six_ver_3_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_3_x = six_ver_ordinary_3(1)
    six_ver_ordinary_3_y = six_ver_ordinary_3(2)
    
    six_ver_ordinary_4 = from_baricentric_to_ordinary(six_ver_2_baricentric, six_ver_3_baricentric, &
     six_ver_1_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_4_x = six_ver_ordinary_4(1)
    six_ver_ordinary_4_y = six_ver_ordinary_4(2)
    
    six_ver_ordinary_5 = from_baricentric_to_ordinary(six_ver_3_baricentric, six_ver_1_baricentric, &
     six_ver_2_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_5_x = six_ver_ordinary_5(1)
    six_ver_ordinary_5_y = six_ver_ordinary_5(2)
    
    six_ver_ordinary_6 = from_baricentric_to_ordinary(six_ver_3_baricentric, six_ver_2_baricentric, &
     six_ver_1_baricentric, x1, x2, x3, y1, y2, y3)
    
    six_ver_ordinary_6_x = six_ver_ordinary_6(1)
    six_ver_ordinary_6_y = six_ver_ordinary_6(2)
      
    
     ! calculation the contribution of 6-multiplicity points
     
     
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_1_x, &
     six_ver_ordinary_1_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_2_x, &
     six_ver_ordinary_2_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_3_x, &
     six_ver_ordinary_3_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_4_x, &
     six_ver_ordinary_4_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_5_x, &
     six_ver_ordinary_5_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v)
    intermediate_result = intermediate_result + six_weight*f(six_ver_ordinary_6_x, &
     six_ver_ordinary_6_y, a1, b1, c1, a2, b2, c2,  a_u, b_u, c_u, a_v, b_v, c_v) 
    
    ! parameters of triangle
    
    
    
    a = dsqrt((x1 - x2)**2 + (y1 - y2)**2)
    b = dsqrt((x1 - x3)**2 + (y1 - y3)**2)
    c = dsqrt((x2 - x3)**2 + (y2 - y3)**2)
    p = (a+b+c)/2d0
    S = dsqrt(p*(p-a)*(p-b)*(p-c))
    
    ! finite result
    
    numeric_integral = S*intermediate_result
      
end function integral_over_triangle   




! scalar multiplication calculation

function scalar_multiplication(f, element1, element2) result(multiplication_result)

   implicit none
   
   type(Element), target, intent(in) :: element1, element2
   real*8                            :: multiplication_result
   logical                           :: not_assoc 
   integer                           :: i, j, k
   real*8                            :: coefficients_1(3), coefficients_2(3), &
                                        coefficients_u(3), coefficients_v(3)
   procedure(mult_phi)               :: f
   
   multiplication_result = 0d0
   not_assoc = .true.
   
   do i = 1, element1%number_of_neighbour_elements
    
      if (associated(element1%neighbour_elements_list(i)%pointing_element, element2)) then
      
        not_assoc = .false.
      
      end if
    
    end do
   
   if (not_assoc) then
       
      multiplication_result = 0d0
       
    else
      
      do i = 1, element1%number_of_neighbour_triangles
      
        do j = 1, element2%number_of_neighbour_triangles
        
          if (associated(element1%neighbour_triangles_list(i)%pointing_triangle, &
             element2%neighbour_triangles_list(j)%pointing_triangle)) then
             
               coefficients_1 = coefficients_calculation(element1, &
                element1%neighbour_triangles_list(i)%pointing_triangle)
            
               coefficients_2 = coefficients_calculation(element2, &
                element1%neighbour_triangles_list(i)%pointing_triangle)
                
               coefficients_u = u_coefficients_calculation(element1%neighbour_triangles_list(i)%pointing_triangle, 1)
               coefficients_v = u_coefficients_calculation(element1%neighbour_triangles_list(i)%pointing_triangle, 2) 
               
                
               multiplication_result = multiplication_result + integral_over_triangle( &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(1,1), &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(1,2), &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(1,3), &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(2,1), &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(2,2), &
               element1%neighbour_triangles_list(i)%pointing_triangle%coordinates(2,3), &
               f, coefficients_1(1), coefficients_1(2), coefficients_1(3), &
               coefficients_2(1), coefficients_2(2), coefficients_2(3), &
               coefficients_u(1), coefficients_u(2), coefficients_u(3), &
               coefficients_v(1), coefficients_v(2), coefficients_v(3)) 
               
               
          end if   
        
        end do
          
      
      end do
       
    end if   



end function scalar_multiplication


! scalar multiplication on given triangle


function scalar_multiplication_on_triangle(f, element1, element2, trian) result(multiplication_result)

   implicit none
   
   type(Element), target, intent(in)  :: element1, element2
   type(Triangle), target             :: trian
   real*8                             :: multiplication_result
   integer                            :: i, j, k
   real*8                             :: coefficients_1(3), coefficients_2(3), &
                                         coefficients_u(3), coefficients_v(3)
   procedure(mult_phi)                :: f
   
   multiplication_result = 0d0
   
   coefficients_1 = coefficients_calculation(element1, &
   trian)
   coefficients_2 = coefficients_calculation(element2, &
   trian)
   coefficients_u = u_coefficients_calculation(trian, 1)
   coefficients_v = u_coefficients_calculation(trian, 2) 
                
   multiplication_result = integral_over_triangle( &
   trian%coordinates(1,1), &
   trian%coordinates(1,2), &
   trian%coordinates(1,3), &
   trian%coordinates(2,1), &
   trian%coordinates(2,2), &
   trian%coordinates(2,3), &
   f, coefficients_1(1), coefficients_1(2), coefficients_1(3), &
   coefficients_2(1), coefficients_2(2), coefficients_2(3), &
   coefficients_u(1), coefficients_u(2), coefficients_u(3), &
   coefficients_v(1), coefficients_v(2), coefficients_v(3)) 
  

end function scalar_multiplication_on_triangle



! sigma grad(phi) scalar multiplication_function

function sigma_grad_phi_scalar_multiplication(elem, ind) result(mult_result)

  implicit none

  type(Element), target, intent(in)   :: elem
  integer, intent(in)                 :: ind
  real*8                              :: sigma_11, sigma_12, sigma_22
  integer                             :: i, j, k
  real*8                              :: coefficients(3), mult_result
  
  mult_result = 0d0
  
  if (ind == 1) then
  
    do i = 1, elem%number_of_neighbour_triangles
    
      coefficients = coefficients_calculation(elem, &
       elem%neighbour_triangles_list(i)%pointing_triangle)
                
      sigma_11 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 + &
      elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
      sigma_22 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 - &
      elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
      sigma_12 = elem%neighbour_triangles_list(i)%pointing_triangle%sigma12
  
      mult_result = mult_result + &
      elem%neighbour_triangles_list(i)%pointing_triangle%size_of_triangle*( &
      sigma_11*coefficients(1) + sigma_12*coefficients(2))
  
    end do
    
    else
    
    do i = 1, elem%number_of_neighbour_triangles
    
      coefficients = coefficients_calculation(elem, &
       elem%neighbour_triangles_list(i)%pointing_triangle)
  
      sigma_11 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 + &
      elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
      sigma_22 = 5d-1*(elem%neighbour_triangles_list(i)%pointing_triangle%sigma1 - &
      elem%neighbour_triangles_list(i)%pointing_triangle%sigma2)
      sigma_12 = elem%neighbour_triangles_list(i)%pointing_triangle%sigma12
  
      mult_result = mult_result + &
      elem%neighbour_triangles_list(i)%pointing_triangle%size_of_triangle*( &
      sigma_12*coefficients(1) + sigma_22*coefficients(2))
  
    end do
       
  endif

end function sigma_grad_phi_scalar_multiplication




!coefficients a, b, c calculation for multiplication of two functions (plane parametrization)

function coefficients_calculation(elem, trian) result(coeff)

  implicit none

  type(Element), target, intent(in)  :: elem
  type(Triangle), target, intent(in) :: trian
  real*8                             :: coeff(3)
  real*8                             :: x0, x1, x2, y0, y1, y2, a, b, c
  
  if (associated(trian%neighbour_elements_list(1)%pointing_element, elem)) then
  
    x0 = elem%coordinates(1)
    y0 = elem%coordinates(2)
    
    x1 = trian%neighbour_elements_list(2)%pointing_element%coordinates(1)
    y1 = trian%neighbour_elements_list(2)%pointing_element%coordinates(2)
    
    x2 = trian%neighbour_elements_list(3)%pointing_element%coordinates(1)
    y2 = trian%neighbour_elements_list(3)%pointing_element%coordinates(2)
  
  end if
  
  if (associated(trian%neighbour_elements_list(2)%pointing_element, elem)) then
  
    x0 = elem%coordinates(1)
    y0 = elem%coordinates(2)
    
    x1 = trian%neighbour_elements_list(1)%pointing_element%coordinates(1)
    y1 = trian%neighbour_elements_list(1)%pointing_element%coordinates(2)
    
    x2 = trian%neighbour_elements_list(3)%pointing_element%coordinates(1)
    y2 = trian%neighbour_elements_list(3)%pointing_element%coordinates(2)
  
  end if
  
  
  if (associated(trian%neighbour_elements_list(3)%pointing_element, elem)) then
  
    x0 = elem%coordinates(1)
    y0 = elem%coordinates(2)
    
    x1 = trian%neighbour_elements_list(1)%pointing_element%coordinates(1)
    y1 = trian%neighbour_elements_list(1)%pointing_element%coordinates(2)
    
    x2 = trian%neighbour_elements_list(2)%pointing_element%coordinates(1)
    y2 = trian%neighbour_elements_list(2)%pointing_element%coordinates(2)
  
  end if
  
  
  a = (y2 - y1)/((y0 - y1)*(x0 - x2) - (y0 - y2)*(x0 - x1))
  b = (x1 - x2)/((y0 - y1)*(x0 - x2) - (y0 - y2)*(x0 - x1))
  c = (x2*y1 - x1*y2)/((y0 - y1)*(x0 - x2) - (y0 - y2)*(x0 - x1))
  
  coeff(1) = a
  coeff(2) = b
  coeff(3) = c
  
end function coefficients_calculation


function u_coefficients_calculation(tri, ident) result(u_coeff_res)

  type(Triangle), target, intent(in) :: tri
  integer, intent(in)                :: ident
  real*8                             :: u_coeff_res(3)
  real*8                             :: x0, y0, h0
  real*8                             :: x1, y1, h1
  real*8                             :: x2, y2, h2
  real*8                             :: a, b, c
  
  
  x0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(1)
  y0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(2)
  
  x1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(1)
  y1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(2)
  
  
  x2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(1)
  y2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(2)
   
  
  if (ident == 1) then
  
    h0 = tri%neighbour_elements_list(1)%pointing_element%u(1)
    h1 = tri%neighbour_elements_list(2)%pointing_element%u(1)
    h2 = tri%neighbour_elements_list(3)%pointing_element%u(1)
  
  end if
  
  
  if (ident == 2) then
  
    h0 = tri%neighbour_elements_list(1)%pointing_element%u(2)
    h1 = tri%neighbour_elements_list(2)%pointing_element%u(2)
    h2 = tri%neighbour_elements_list(3)%pointing_element%u(2)
  
  end if
  
  
  a = ((h0 - h1)*(y0 - y2) - (y0 - y1)*(h0 - h2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  b = ((h0 - h2)*(x0 - x1) - (h0 - h1)*(x0 - x2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  c = (x0*(h2*y1 - h1*y2) + x1*(h0*y2 - h2*y0) + x2*(h1*y0 - h0*y1))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  
  u_coeff_res(1) = a
  u_coeff_res(2) = b
  u_coeff_res(3) = c
  
end function u_coefficients_calculation



function L2_mass(ident) result(L2_res) ! 1 -- m,  2 -- h, 3 -- A

  implicit none
  
  integer, intent(in)  :: ident
  integer              :: number_of_triangles, i, j
  real*8               :: prom, L2_res
  
  prom = 0d0
  
  do i = 1, number_of_triangles
  
    prom = prom + linear_average_on_triangle(List_of_Triangles(i), ident)*List_of_Triangles(i)%size_of_triangle
  
  end do
  
  L2_res = prom

end function L2_mass




function resud_u_coefficients_calculation(tri, ident, u_v) result(u_coeff_res) ! u_v = 1 -- u_new, 2 -- u_old

  implicit none

  type(Triangle), target, intent(in) :: tri
  integer, intent(in)                :: ident, u_v
  real*8                             :: u_coeff_res(3)
  real*8                             :: x0, y0, h0
  real*8                             :: x1, y1, h1
  real*8                             :: x2, y2, h2
  real*8                             :: a, b, c
  
  x0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(1)
  y0 = tri%neighbour_elements_list(1)%pointing_element%coordinates(2)
  
  
  x1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(1)
  y1 = tri%neighbour_elements_list(2)%pointing_element%coordinates(2)
  
  
  x2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(1)
  y2 = tri%neighbour_elements_list(3)%pointing_element%coordinates(2)
   
  
  if (ident == 1) then
  
    if(u_v == 1) then
  
      h0 = tri%neighbour_elements_list(1)%pointing_element%u_new(1)
      h1 = tri%neighbour_elements_list(2)%pointing_element%u_new(1)
      h2 = tri%neighbour_elements_list(3)%pointing_element%u_new(1)
       
    end if
    
    if(u_v == 2) then
  
      h0 = tri%neighbour_elements_list(1)%pointing_element%u_old(1)
      h1 = tri%neighbour_elements_list(2)%pointing_element%u_old(1)
      h2 = tri%neighbour_elements_list(3)%pointing_element%u_old(1)
       
    end if
  
  end if
  
  
  if (ident == 2) then
  
    if (u_v == 1) then
    
      h0 = tri%neighbour_elements_list(1)%pointing_element%u_new(2)
      h1 = tri%neighbour_elements_list(2)%pointing_element%u_new(2)
      h2 = tri%neighbour_elements_list(3)%pointing_element%u_new(2)
  
    end if
    
    if (u_v == 2) then
    
      h0 = tri%neighbour_elements_list(1)%pointing_element%u_old(2)
      h1 = tri%neighbour_elements_list(2)%pointing_element%u_old(2)
      h2 = tri%neighbour_elements_list(3)%pointing_element%u_old(2)
  
    end if
    
    
  end if
  
  
  a = ((h0 - h1)*(y0 - y2) - (y0 - y1)*(h0 - h2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  b = ((h0 - h2)*(x0 - x1) - (h0 - h1)*(x0 - x2))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  c = (x0*(h2*y1 - h1*y2) + x1*(h0*y2 - h2*y0) + x2*(h1*y0 - h0*y1))/((y0 - y2)*(x0 - x1) - (y0 - y1)*(x0 - x2))
  
  u_coeff_res(1) = a
  u_coeff_res(2) = b
  u_coeff_res(3) = c
  
end function resud_u_coefficients_calculation

!! scalar multiplication of two vectors

function scalar_mult_two_vec(vec1, vec2, size_of_vec) result(scalar_mult_res)

  implicit none

  real*8, intent(in)  :: vec1(:), vec2(:) 
  integer, intent(in) :: size_of_vec
  real*8              :: scalar_mult_res, prom
  integer             :: i, j, k
  
  prom = 0d0
  
  do i = 1, size_of_vec
  
    prom = prom + vec1(i)*vec2(i)
  
  end do
  
  scalar_mult_res = prom


end function scalar_mult_two_vec


end module module_numerical_integration
