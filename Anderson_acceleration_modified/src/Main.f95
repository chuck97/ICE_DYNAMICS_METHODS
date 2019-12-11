program Main
  
  use module_classes
  use module_initialization
  use module_grid_values
  use module_constants
  use module_assembling
  use module_dynamics
  use module_plotwrite
  use module_advection
  use json_module
      
  
  implicit none
  
  !! global time step, current time:
  real*8                         :: time_step, time
  real*8                         :: squ, max_u
  integer                        :: num_iter, resudal_mass_index 
  real*8, allocatable            :: init_resid(:)
  
  
  !! mEVP variables
  
  real*8                         :: x_coord, y_coord
  
  !! local variables:
  integer                        :: i, j, k, num, s, r, ps

  !! Arrays for the matrix sparse row (CSR) format
  real*8                         :: Right_hand_C(maxn), Right_hand_L(maxn), &
                                   solution_vec_L(maxn), solution_vec_C(maxn), &
                                   solution_vec(maxn), sum_mass, sum_concentration
  
  !! json variables   
  logical                        :: is_found
  type(json_file)                :: json
  character(len = strlen)        :: path_to_json
     
  !! grid variables 
  integer                        :: Nbv, Nbl
  real*8,           allocatable  :: bv(:, :)
  integer, allocatable           :: bl(:, :)
  real*8, allocatable            :: bltail(:,:)
   
  real*8, allocatable            :: bv_in(:)
  integer, allocatable           :: bl_in(:)
  real*8, allocatable            :: bltail_in(:)
   
  !! pathes to different folders
  character(len=:), allocatable  :: path_to_triangles
  character(len=:), allocatable  :: path_to_nodals
  character(len=:), allocatable  :: path_to_graphics
  
  !! boundary type 
  character(:), allocatable      :: boundary_type
  
  !! triangulation parameters
  real*8                         :: square_size
  integer                        :: number_of_vertex_per_square_edge
 
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  START  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
         
  !! read everything from .json file
  
  call json%initialize()
  
  call get_command_argument(1, path_to_json)
  if (len_trim(path_to_json) == 0) stop "no path to .json file!"
  
  call json%load(filename = trim(path_to_json))
  if (json%failed()) stop "failed to load .json!"
  
  call json%get('time step in hours', time_step, is_found)
  if (.not. is_found) stop "no value of time step!"
    
  call json%get('number of boundary nodes', Nbv, is_found)
  if (.not. is_found) stop "no value of number of boundary nodes!"
  
  call json%get('number of boundary edges', Nbl, is_found)
  if (.not. is_found) stop "no value of number of boundary edges!"
  
  allocate(bv_in(2*Nbv))
  allocate(bl_in(7*Nbl))
  allocate(bltail_in(24))
  
  call json%get('square size in meters', square_size, is_found)
  if (.not. is_found) stop "no value of square size!"
  
  call json%get('number of grid verticies per square edge', number_of_vertex_per_square_edge, is_found)
  if (.not. is_found) stop "no value of number of grid verticies per square edge!"
  
  call json%get('boundary nodes coordinates unit size', bv_in, is_found)
  if (.not. is_found) stop "no value of boundary nodes!"
  
  call json%get('boundary edges connectivity list', bl_in, is_found)
  if (.not. is_found) stop "no value of boundary edges connectivity list!"
  
  call json%get('curved boundary list', bltail_in, is_found)
  if (.not. is_found) stop "no value of curved boundary list!"
  
  call json%get('path to graphics folder', path_to_graphics, is_found)
  if (.not. is_found) stop "no value of path to graphics folder!"
  
  call json%get('path to triangles folder', path_to_triangles, is_found)
  if (.not. is_found) stop "no value of path to triangles folder!"
  
  call json%get('path to nodals folder', path_to_nodals, is_found)
  if (.not. is_found) stop "no value of path to nodals folder!"
  
  call json%get('velocity boundary condition', boundary_type, is_found)
  if (.not. is_found) stop "no value of velocity boundary type!"
  
  
  call json%destroy()
  
  print *, ".Json reading: done"
  
  !!  grid initialization
  
  allocate(bv(2,Nbv))
  allocate(bl(7,Nbl))
  allocate(bltail(2,12))
  
  bv = reshape(bv_in, shape(bv))
  bl = reshape(bl_in, shape(bl))
  
  bltail = reshape(bltail_in, shape(bltail))
  
  do i = 1, 2
    do j = 1, Nbv
      bv(i,j) = bv(i,j)*square_size
    end do
  end do  
  
  call grid_initialization( &
  Nbv, Nbl, bv, bl, bltail, number_of_vertex_per_square_edge, square_size)
  
  deallocate(bv)
  deallocate(bl)
  deallocate(bltail)
  
  deallocate(bv_in)
  deallocate(bl_in)
  deallocate(bltail_in)
  
  print *, "Grid generating: done"
  
  !! scalars, vectors, forcing initialization
  
  call vectors_initialization(boundary_type)
  call scalars_initialization()
  call wind_forcing_initialization(boundary_type)
  call water_forcing_initialization(boundary_type)
  
  print *, "External forcing setup: done"
  
  !! setup time step in seconds
  
  time_step = time_step*hour
  
  print *, "Time step setup: done"
  
  !! Assembling non-Direchlet matrix and non-Direchlet lumped matrix
  
  call non_Direchlet_mass_matrix_assembling()
  call non_Direchlet_mass_matrix_lumped_assembling()
   
  print *, "Assembling ND and lumped ND mass matrix: done"
  
  !! Assembling  M_L and M_C in CSR format for transport equation
  
  call transport_mass_matrix_low_order_assembling(time_step)
  call transport_mass_matrix_high_order_assembling(time_step)
  
  print *, "Assembling mass matricies for transport: done"
  
  
  !! Time stepping 
  
  time = 0d0
  
  do i = 1, number_of_triangles
  
    List_of_Triangles(i)%sigma1 = 0d0
    List_of_Triangles(i)%sigma2 = 0d0
    List_of_Triangles(i)%sigma12 = 0d0
  
  end do
  
  allocate(init_resid(2*number_of_non_Direchlet_elements))
  
  do while (num < num_time_steps)
  
    !!! wind and water setup
  
    do i = 1, number_of_non_Direchlet_elements
    
      x_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(1)
      y_coord = List_of_non_Direchlet_Elements(i)%pointing_element%coordinates(2)
    
      List_of_non_Direchlet_Elements(i)%pointing_element%u_air(1) = 5d0 + &
       (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*x_coord/square_size)*dsin(pi*y_coord/square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_air(2) = 5d0 + &
       (dsin(2d0*pi*time/T_period) - 3d0)*dsin(2d0*pi*y_coord/square_size)*dsin(pi*x_coord/square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_water(1) = &
       (2d0*y_coord - 1d6)/(10d0*square_size)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_water(2) = &
       -(2d0*x_coord - 1d6)/(10d0*square_size)  
       
      List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(1) = &
      List_of_non_Direchlet_Elements(i)%pointing_element%u(1)
      
      List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(2) = &
      List_of_non_Direchlet_Elements(i)%pointing_element%u(2) 
       
    
    end do
    
    !! compute initial residual each time step
    
    init_resid = init_residual(time_step, 0)
  
  
    num_iter = 1
    
    call dot_epsilon_delta_recalculation(2)
    call P_0_recalculation()
    
    !! mEVP-stepping
    
    do while((num_iter < (N_evp+1)))
    
      call dot_epsilon_delta_recalculation(1)
      call P_0_recalculation()
      
      !!! sigma recalculation
    
      do i = 1, number_of_Triangles
      
        List_of_Triangles(i)%sigma1 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma1 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon1 - &
         List_of_Triangles(i)%delta)* &
         (List_of_Triangles(i)%P_0)/ &
         (List_of_Triangles(i)%delta + delta_min)
         
        List_of_Triangles(i)%sigma2 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma2 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon2)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2)
         
        List_of_Triangles(i)%sigma12 = &
         (1d0 - (1d0/alpha_EVP))*List_of_Triangles(i)%sigma12 + &
         (1d0/alpha_EVP)*(List_of_Triangles(i)%dot_epsilon12)* &
         (List_of_Triangles(i)%P_0)/ &
         ((List_of_Triangles(i)%delta + delta_min)*e**2) 
         
      end do
      
      call velocity_recalculation(time_step)
      
      do i = 1, number_of_non_Direchlet_elements
  
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(1)
    
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_new(2)
       
        
        List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(1) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
    
        List_of_non_Direchlet_Elements(i)%pointing_element%u_resid(2) = &
        List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2)
        
      end do
      
      num_iter = num_iter + 1
        
    end do
      
      do i = 1, number_of_non_Direchlet_elements
      
        List_of_non_Direchlet_Elements(i)%pointing_element%u(1) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1)
         
        List_of_non_Direchlet_Elements(i)%pointing_element%u(2) = &
         List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) 
      
      end do
      
      print *, "max_u:", maximum_u()
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         Mass transport         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call scalar_advection(time_step, 1)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!     Concentration transport    !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call scalar_advection(time_step, 2)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!        h recalculation         !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, number_of_elements
        List_of_Elements(i)%h = (List_of_Elements(i)%m)/(List_of_Elements(i)%A*rho_ice)
      end do
      
      
      !! next time step
      
      num = num + 1
      time = time + time_step
      
      print *, "time moment: ", time/hour , "hours"
      
      call plot_ParaView(num, path_to_graphics)
      call write_nodals(num, path_to_nodals)
      call write_triangles(num, path_to_triangles)
      
        
  end do
  
  print *, "Time stepping: done"
  
  deallocate(init_resid)    
      
end program Main
