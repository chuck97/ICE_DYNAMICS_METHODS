module module_water

  use module_classes
  use module_constants
  use module_grid_values
  
  private
  
  public :: water_forcing_update, gradient_water_level_update
  
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
        List_of_Elements(i)%u_water(1) = 0d0 !&
          !(2d0*y_coord - 1d6)/(10d0*square_size)
        List_of_Elements(i)%u_water(2) = 0d0 !&
         !-(2d0*x_coord - 1d6)/(10d0*square_size)  
      end if
    end do
  
  end subroutine water_forcing_update
  
  
  subroutine gradient_water_level_update(time, square_size, is_geostrophic)
  
    implicit none
    
    !!arguments:
    real*8, intent(in)       :: time, square_size
    logical, intent(in)      :: is_geostrophic
    
    !!local arguments:
    integer                  :: i
    real*8                   :: u1, u2, u3, v1, v2, v3, u_ave, v_ave
    
    if (is_geostrophic) then
    
      do i = 1, number_of_triangles
        u1 = List_of_triangles(i)%neighbour_elements_list(1)%pointing_element%u_water(1)
        u2 = List_of_triangles(i)%neighbour_elements_list(2)%pointing_element%u_water(1)
        u3 = List_of_triangles(i)%neighbour_elements_list(3)%pointing_element%u_water(1)
        v1 = List_of_triangles(i)%neighbour_elements_list(1)%pointing_element%u_water(2)
        v2 = List_of_triangles(i)%neighbour_elements_list(2)%pointing_element%u_water(2)
        v3 = List_of_triangles(i)%neighbour_elements_list(3)%pointing_element%u_water(2)
        u_ave = (u1 + u2 + u3)/3d0
        v_ave = (v1 + v2 + v3)/3d0        
        List_of_triangles(i)%dH_dx = (2d0*omega_e/g_gravity)*v_ave
        List_of_triangles(i)%dH_dy = -(2d0*omega_e/g_gravity)*u_ave
      end do
    
    else
    
      stop "no input water level"
    
    end if
  
  end subroutine  gradient_water_level_update

end module module_water
