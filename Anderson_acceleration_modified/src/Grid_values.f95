module module_grid_values

  use module_Classes
  use module_constants
  
  implicit none
  
  public :: List_of_Triangles,                   &
            List_of_Elements,                    &
            List_of_Edges,                       &
            List_of_non_Direchlet_Elements,      &
            number_of_non_Direchlet_elements,    &
            number_of_elements,                  &
            number_of_triangles,                 &
            number_of_edges,                     &
            number_of_boundary_edges             
  
  type(Triangle), target, dimension(ntmax)               :: List_of_Triangles                 
  type(Element),  target, dimension(nvmax)               :: List_of_Elements                  
  type(Edge),     target, dimension((ntmax*3)/2 + nbmax) :: List_of_Edges                     
  type(container_element), dimension(nvmax)              :: List_of_non_Direchlet_Elements    
  
  integer                                                :: number_of_non_Direchlet_elements  
  integer                                                :: number_of_elements                
  integer                                                :: number_of_triangles               
  integer                                                :: number_of_edges                   
  integer                                                :: number_of_boundary_edges          

end module module_grid_values
