module module_classes
  
  
  implicit none
  
  
  private
  public :: Triangle, Element, Edge, Boundary_Edge, Matrix, &
    container_element, container_triangle, container_edge


 !! Container for Elemets list of pointers

  type container_element
  
      type(Element), pointer :: pointing_element 
      
  end type container_element

 !! Container for Triangles lists of pointers   
  
  type container_triangle
  
      type(Triangle), pointer :: pointing_triangle
       
  end type container_triangle
  
  
  
 !! Container for Edges list of pointers
 
  type container_edge
  
      type(Edge), pointer     :: pointing_edge 
   
  end type container_edge
  
  
  !! Type for edges
 
  type Edge
  
     integer                               :: identificator               ! edge identificator
     real*8                                :: coordinates(2,2)            ! verticies coordinates
     logical                               :: on_boundary                 ! if edge on boundary?
     real*8                                :: length_of_edge              ! square of triangle
     type(container_element), dimension(2) :: neighbour_elements_list     ! list of elements for edge
     
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!! When edge is on boundary, the domain is on the left side while moving from first point to the second !!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
  end type Edge
  
  
  !! Type boundary edge
  
  type Boundary_Edge
  
  
     integer                               :: identificator                       
     real*8                                :: coordinates(2,2)                    
     integer                               :: domain_boundary                     
     type(container_element), dimension(2) :: neighbour_elements_list             
  
  
  end type Boundary_Edge

 !! Type for triangles
  
  type Triangle
  
     integer                                 :: identificator                                            
     real*8                                  :: dot_epsilon1, dot_epsilon2, dot_epsilon12                
     real*8                                  :: P_0                                                      
     real*8                                  :: coordinates(2,3)                                         
     real*8                                  :: sigma1, sigma2, sigma12
     real*8                                  :: sigma1_previous, sigma2_previous, sigma12_previous
     real*8                                  :: sigma_resid(3)                                   
     real*8                                  :: delta
     real*8                                  :: xi_resuid, eta_resuid, P_0_resuid
     real*8                                  :: delta_old         
     real*8                                  :: dH_dx, dH_dy                                          
     real*8                                  :: size_of_triangle                                         
     integer                                 :: number_of_neighbour_triangles                            
     type(container_element), dimension(3)   :: neighbour_elements_list                                  
     type(container_edge), dimension(3)      :: neighbour_edges_list                                     
     type(container_triangle), dimension(3)  :: neighbour_triangles_list                                  
     real*8                                  :: flux_to_element(3)                                       
     real*8                                  :: u_max_triangle, u_min_triangle                           
     real*8                                  :: alpha_correction                                         
     real*8                                  :: xi, eta                                                  
     real*8                                  :: delta_str
       
  end type Triangle

  
  
 !! Type for Elements (basis functions) (every element has different quantity of neighbour elements)  

  type Element
  
     integer                                      :: identificator                                         
     integer                                      :: nd_identificator                                      
     real*8                                       :: coordinates(2)                                        
     real*8                                       :: u(2)                                                  
     real*8                                       :: u_old(2)                                              
     real*8                                       :: u_water(2)
     real*8                                       :: u_resid(2)                                         
     real*8                                       :: u_air(2)
     real*8                                       :: u_new(2)
     real*8                                       :: water_level
     real*8                                       :: m                                                     
     real*8                                       :: old_m
     real*8                                       :: h                                                     
     real*8                                       :: A                                                     
     real*8                                       :: old_A
     real*8                                       :: element_value                                         
     logical                                      :: on_boundary                                           
     integer                                      :: number_of_neighbour_elements                          
     integer                                      :: number_of_neighbour_triangles                         
     integer                                      :: number_of_neighbour_edges                             
     type(container_element), dimension(20)       :: neighbour_elements_list                               
     type(container_triangle), dimension(20)      :: neighbour_triangles_list                              
     type(container_edge), dimension(20)          :: neighbour_edges_list                                  
     real*8                                       :: P_pos, P_neg                                          
     real*8                                       :: Q_plus, Q_minus                                       
     real*8                                       :: u_max_element, u_min_element                          
     real*8                                       :: R_plus, R_minus                                       
     real*8                                       :: resid_u, resid_v  
     real*8                                       :: u_str(2)
       
  end type Element
  
  
 !! Sparce matrix in CSR format type
  
  type Matrix
  
    real*8 , allocatable :: a (:)
    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    integer              :: nonzero
    integer              :: matr_size_x
    integer              :: matr_size_y
    
    contains 
    
     procedure, public :: init => init_matrix
     procedure, public :: distr => distruct_matrix
  
  end type Matrix
  
  contains
  
  ! Matrix constructor
  
  subroutine init_matrix(this, nonzero, matr_size_x, matr_size_y) 
  
    class(Matrix), intent(inout) :: this
    integer, intent(in)          :: nonzero
    integer, intent(in)          :: matr_size_x, matr_size_y
   
    this%nonzero = nonzero
    this%matr_size_x = matr_size_x
    this%matr_size_y = matr_size_y
    
    allocate(this%a(this%nonzero))
    allocate(this%ia(this%matr_size_y + 1))
    allocate(this%ja(this%nonzero))
  
  end subroutine init_matrix
  
  ! Matrix distruction
  
  subroutine distruct_matrix(this)
  
    class(Matrix), intent(inout) :: this
    
    this%nonzero = 0
    this%matr_size_x = 0
    this%matr_size_y = 0
    
    deallocate(this%ia)
    deallocate(this%ja)
    deallocate(this%a)
  
  end subroutine distruct_matrix
  
  
    
end module module_classes
