module module_grid_values

  use module_Classes
  
  implicit none
  public ::  List_of_Triangles,                 &
             List_of_Elements,                  &
             List_of_Edges,                     &
             List_of_non_Direchlet_Elements,    &
             number_of_non_Direchlet_elements,  &
             number_of_elements,                &
             number_of_triangles,               &
             number_of_edges,                   &
             vectors_initialization,            &
             scalars_initialization 
  
  type(Triangle), target, dimension(ntmax)               :: List_of_Triangles                 !! List of Triangles
  type(Element),  target, dimension(nvmax)               :: List_of_Elements                  !! List of Elements
  type(Edge),     target, dimension((ntmax*3)/2 + nbmax) :: List_of_Edges                     !! List_of_Edges
  type(container_element), dimension(nvmax)              :: List_of_non_Direchlet_Elements    !! List of non-Direchlet elements
  
  integer                                                :: number_of_non_Direchlet_elements  !! Number of non-Direchlet elements
  integer                                                :: number_of_elements                !! Number of Elements
  integer                                                :: number_of_triangles               !! Number of triangles
  integer                                                :: number_of_edges                   !! Number of edges
  integer                                                :: number_of_boundary_edges          !! Number of boundary edges
  
  
  contains
  

  
  
  
!!  grid generation and triangles initialization



!! vectors initialization 
 
  subroutine vectors_initialization()
  
    implicit none
    
    !! local variables:
    
    integer :: i
  
    do i = 1, number_of_elements
    
      List_of_Elements(i)%u(1) = 0d0
      List_of_Elements(i)%u(2) = 0d0
      List_of_Elements(i)%u_old(1) = 0d0
      List_of_Elements(i)%u_old(2) = 0d0
    
    end do
  
  end subroutine vectors_initialization


!! scalars scalars initialization

  subroutine scalars_initialization()
  
    implicit none
    
    !! local variables:
    
    integer    :: i
    
    
    !! h initialization
    
    do i = 1, number_of_elements
    
      List_of_Elements(i)%h = 1d0
    
    end do 
    
    !! A initialization
    
    do i = 1, number_of_elements
    
      List_of_Elements(i)%A = 95d-2
    
    end do 
    
    
    !! m initialization
      
    do i = 1, number_of_elements
      
      List_of_Elements(i)%m = (List_of_Elements(i)%h)*(List_of_Elements(i)%A*rho_ice)
     
    end do
  
  end subroutine scalars_initialization


!! elements information

subroutine elements_information(List_of_Elements, number_of_elements)

  type(Element), target, dimension(nvmax) :: List_of_Elements        ! list of elements
  integer :: number_of_elements, i, j
  
  do i = 1, number_of_elements
  
    print *, "element", List_of_Elements(i)%identificator, ":"
    
    print *, "coordinates ", List_of_Elements(i)%coordinates(1), List_of_Elements(i)%coordinates(2)
    
    if (List_of_Elements(i)%on_boundary) then
      
      print *, "on boundary"
        
    end if
    
    print *, List_of_Elements(i)%number_of_neighbour_elements," neighbour elements:"
    
    do j = 1, List_of_Elements(i)%number_of_neighbour_elements
    
      write(*, fmt="(1x, I10)", advance="no") List_of_Elements(i)%neighbour_elements_list(j)%pointing_element%identificator
    
    end do
    
    print *
    
    print *, List_of_Elements(i)%number_of_neighbour_triangles," neighbour triangles:"
    
    do j = 1, List_of_Elements(i)%number_of_neighbour_triangles
    
      write(*, fmt="(1x, I10)", advance="no") List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle%identificator
    
    end do
    
    print *
    
    
    print *, List_of_Elements(i)%number_of_neighbour_edges," neighbour edges:"
    
    do j = 1, List_of_Elements(i)%number_of_neighbour_edges
    
      write(*, fmt="(1x, I10)" , advance="no") List_of_Elements(i)%neighbour_edges_list(j)%pointing_edge%identificator
    
    end do
    
    print *
    print *
    
    
  end do


end subroutine elements_information




!! triangles information

subroutine triangles_information(List_of_Triangles, number_of_triangles)

  type(Triangle), target, dimension(ntmax) :: List_of_Triangles        ! list of triangles
  integer :: number_of_triangles, i, j
  
  do i = 1, number_of_triangles
  
    print *, "triangle", List_of_Triangles(i)%identificator, ":"
    
    print *, "coordinates:"
    print *,  List_of_Triangles(i)%coordinates(1, 1), List_of_Triangles(i)%coordinates(2, 1)
    print *,  List_of_Triangles(i)%coordinates(1, 2), List_of_Triangles(i)%coordinates(2, 2)
    print *,  List_of_Triangles(i)%coordinates(1, 3), List_of_Triangles(i)%coordinates(2, 3)
    
    
    print *, List_of_Triangles(i)%number_of_neighbour_triangles," neighbour triangles:"
    
    do j = 1, List_of_Triangles(i)%number_of_neighbour_triangles
    
      write(*, fmt="(1x, I10)", advance="no") List_of_Triangles(i)%neighbour_triangles_list(j)%pointing_triangle%identificator
    
    end do
    
    print *
    
    print *, "3 neighbour elements:"
    
    do j = 1, 3
    
      write(*, fmt="(1x, I10)", advance="no") List_of_Triangles(i)%neighbour_elements_list(j)%pointing_element%identificator
    
    end do
    
    print *
    
    
    print *, "3 neighbour edges:"
    
    do j = 1, 3
    
      write(*, fmt="(1x, I10)" , advance="no") List_of_Triangles(i)%neighbour_edges_list(j)%pointing_edge%identificator
    
    end do
    
    print *
    print *
    
    
  end do


end subroutine triangles_information




!! edges information

subroutine edges_information(List_of_Edges, number_of_edges)

  type(Edge), target, dimension((ntmax*3)/2 + nbmax) :: List_of_Edges        ! list of edges
  integer :: number_of_edges, i, j
  
  do i = 1, number_of_edges
  
    print *, "edge", List_of_Edges(i)%identificator, ":"
    
    print *, "coordinates:"
    print *,  List_of_Edges(i)%coordinates(1, 1), List_of_Edges(i)%coordinates(2, 1)
    print *,  List_of_Edges(i)%coordinates(1, 2), List_of_Edges(i)%coordinates(2, 2)
    
    if (List_of_Edges(i)%on_boundary) then
    
      print *, "on boundary"
    
    end if
    
    
    print *, "2 neighbour elements:"
    
    do j = 1, 2
    
      write(*, fmt="(1x, I10)", advance="no") List_of_Edges(i)%neighbour_elements_list(j)%pointing_element%identificator
    
    end do
    
    
    print *
    print *
    
    
  end do


end subroutine edges_information




!! Plot subroutine  




! write triangles values to file by order of List_of_Triangles: delta, sigma1, sigma2, sigma12, P_0
    

end module module_grid_values
