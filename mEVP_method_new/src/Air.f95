module module_air

  use module_classes
  use module_constants
  use module_grid_values
  
  implicit none
  
  private
  
  public :: air_forcing_update
  
  contains
  
  !! This routine initializes wind forcing
  subroutine air_forcing_update(boundary_type, time, square_size)
  
    implicit none
    
    !!arguments:
    character(len=*), intent(in) :: boundary_type
    real*8, intent(in)           :: time, square_size
    
    !!local variables
    integer     :: i
    real*8      :: x_coord, y_coord
    
    do i = 1, number_of_elements
      if (List_of_Elements(i)%on_boundary) then
        if (trim(boundary_type) == "sticking") then
          List_of_Elements(i)%u_air(1) = 0d0
          List_of_Elements(i)%u_air(2) = 0d0
        end if
      else
        !x_coord = List_of_elements(i)%coordinates(1)
        !y_coord = List_of_elements(i)%coordinates(2)
        List_of_Elements(i)%u_air(1) = 0d0 !5d0 + &
         !(dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*x_coord/square_size)*dsin(pi*y_coord/square_size)
        List_of_Elements(i)%u_air(2) = 15d0  !5d0 + &
         !(dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*y_coord/square_size)*dsin(pi*x_coord/square_size)
      end if
    end do
  
  end subroutine air_forcing_update

end module module_air
