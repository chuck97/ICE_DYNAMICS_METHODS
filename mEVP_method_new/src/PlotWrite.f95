module module_plotwrite

  use module_classes
  use module_constants  
  use module_grid_values
  
  implicit none
  public :: plot_ParaView,   &
            itoa,            &
            write_nodals,    &
            write_triangles  
  
  contains
  
  !! itoa function   
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(I5.5)') i
    res = trim(tmp)
  end function itoa   
  
  
  !! This routine makes .plt file from current grid variables
  subroutine plot_ParaView(num, path_to_graphics)

    implicit none 
    
    !! arguments:
    integer, intent(in)               :: num
    character(len = *), intent(in)    :: path_to_graphics
    
    !! local variables:
    integer                           :: i, j, fid
    
    OPEN(1, FILE = trim(path_to_graphics)//trim('plt')//itoa(num)//trim('.vtk'))
    
    WRITE(1,"(a)") "# vtk DataFile Version 2.0"
    write (1, *) "Grid"
    write (1, *) "ASCII"
    write (1, *) "DATASET UNSTRUCTURED_GRID"
    write (1, *) "POINTS", number_of_elements, "FLOAT"
    
    do i = 1, number_of_elements
      write (1, *) List_of_Elements(i)%coordinates(1), List_of_Elements(i)%coordinates(2), &
           0d0
    end do
    
    write (1, *) "CELLS", number_of_triangles, 4*number_of_triangles
    
    do i = 1, number_of_triangles
      write (1, *) 3, List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%identificator - 1, &
       List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator - 1, &
       List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element%identificator - 1
    end do
    
    write (1, *) "CELL_TYPES", number_of_triangles
    
    do i = 1, number_of_triangles
    
    write (1, *) 5
    
    end do
    
    write (1, *)
    
    write (1, *) "POINT_DATA", number_of_elements
    write (1, *) "SCALARS mass float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_elements
      write (1, *) List_of_Elements(i)%m
    end do   
    write (1, *)
    write (1, *) "SCALARS concentration  float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_elements
      write (1, *) List_of_Elements(i)%A
    end do
    write (1, *)
    write (1, *) "SCALARS thickness  float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_elements
      write (1, *) List_of_Elements(i)%h
    end do
    write (1, *)
    write (1, *) "VECTORS velocity float"
    
    do i = 1, number_of_elements
      write (1, *) List_of_Elements(i)%u(1), List_of_Elements(i)%u(2), 0d0
    end do
    
    write (1, *)
    
    write (1, *) "CELL_DATA", number_of_triangles
    write (1, *) "SCALARS delta float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_triangles 
      write (1, *) List_of_Triangles(i)%delta
    end do
    
    write (1, *)
    write (1, *) "SCALARS sigma1 float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_triangles
      write (1, *) List_of_Triangles(i)%sigma1/(List_of_Triangles(i)%P_0 + varepsilon)
    end do
    
    write (1, *)
    write (1, *) "SCALARS sigma2 float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_triangles
      write (1, *) List_of_Triangles(i)%sigma2/(List_of_Triangles(i)%P_0 + varepsilon)
    end do
    
    write (1, *)
    write (1, *) "SCALARS sigma12 float 1"
    write (1, *) "LOOKUP_TABLE default"
    
    do i = 1, number_of_triangles
      write (1, *) List_of_Triangles(i)%sigma12/(List_of_Triangles(i)%P_0 + varepsilon)
    end do
    
    close(1)
    
  end subroutine plot_ParaView  
  
  
  !! This routine writes triangles values to file
  subroutine write_triangles(num, path_to_triangles)
  
    implicit none
    
    !! arguments:
    integer, intent(in)               :: num
    character(len=*), intent(in)      :: path_to_triangles
    
    !! local variables:
    integer                           :: i
  
    open(4, file = trim(path_to_triangles)//trim('triangles')//itoa(num)//trim('.txt'))
    write (4, *) "delta           sigma1           sigma2           sigma12           P_0"
    do i = 1, number_of_triangles
      write (4, *) List_of_Triangles(i)%delta, List_of_Triangles(i)%sigma1, List_of_Triangles(i)%sigma2, &
        List_of_Triangles(i)%sigma12, List_of_Triangles(i)%P_0
    end do

    close(4)

  end subroutine write_triangles
  
  !! This routine writes nodals values to file 

  subroutine write_nodals(num, path_to_nodals)
  
    implicit none
    
    !! arguments:
    integer, intent(in)               :: num
    character(len=*), intent(in)      :: path_to_nodals
    
    !!local variables
    integer                           :: i 
  
    open(3, file= &
    trim(path_to_nodals)//trim('nodals')//itoa(num)//trim('.txt'))
  
    write (3, *) "m                 A               h                 u                  v"
  
    do i = 1, number_of_elements
      write (3, *) List_of_Elements(i)%m, List_of_Elements(i)%A, List_of_Elements(i)%h, &
        List_of_Elements(i)%u(1), List_of_Elements(i)%u(2)
    end do
    
    close(3)
    
  end subroutine write_nodals

                   
end module module_plotwrite
