module module_VP_Picard_dynamics
  
  use module_classes
  use module_constants
  use module_grid_values
  use module_assembling
  
  implicit none
  
  private :: N_Picard_iterations, &
             
             
  public  :: VP_Picard_velocity_update
  
  integer, parameter      :: N_Picard_iterations
  
  contains
  
  subroutine VP_Picard_velocity_update(time_step, ind) ! ind = 1 -- implicit, ind = 2 -- half implicit
  
    implicit none
    
    !!arguments:
    real*8, intent(in)     :: time_step
    integer, intent(in)    :: ind
    
    !!local variables:
    integer                :: i
    
  
  end subroutine VP_Picard_velocity_update

end module module_VP_Picard_dynamics
