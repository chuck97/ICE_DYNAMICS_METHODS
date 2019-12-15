module module_water

  use module_classes
  use module_constants
  use module_grid_values
  
  private
  
  public :: water_forcing_update
  
  contains
  
  !! This routine initializes water forcing
  subroutine water_forcing_update(boundary_type, time, square_size)
  
    implicit none
    
    !!arguments:
    character(len=*), intent(in) :: boundary_type
    real*8, intent(in)           :: time, square_size
    
    !!local variables
    integer     :: i
    real*8      :: x_coord, y_coord
    
    do i = 1, number_of_elements
      if (List_of_Elements(i)%on_boundary) then
        if (trim(boundary_type) == "cling") then
          List_of_Elements(i)%u_water(1) = 0d0
          List_of_Elements(i)%u_water(2) = 0d0
        end if
      else
        x_coord = List_of_elements(i)%coordinates(1)
        y_coord = List_of_elements(i)%coordinates(2)
        List_of_Elements(i)%u_water(1) = &
          (2d0*y_coord - 1d6)/(10d0*square_size)
        List_of_Elements(i)%u_water(2) = &
         -(2d0*x_coord - 1d6)/(10d0*square_size)  
      end if
    end do
  
  end subroutine water_forcing_update

end module module_water
