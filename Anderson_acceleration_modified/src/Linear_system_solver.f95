module linear_system_solver_module

  use module_classes
  use module_constants
  
  implicit none
  private ::  
  public  :: ILU_0_BICG
  
  contains
  
  subroutine ILU_0_BICG(Matrix, right_hand, system_size, result_vector)
  
    implicit none
    
    !! arguments:
    type(Matrix), intent(in)  :: Matrix
    integer, intent(in)       :: system_size
    real*8, intent(in)        :: right_hand(system_size)
    real*8, intent(in)        :: result_vector(system_size)
    
    !! local variables:
    integer                   :: i, j, k, ierr
    integer                   :: ipBCG
    
    ! ILU data
    real*8                    :: tau1, tau2, partlur, partlurout
    integer                   :: verb, UsedWr, UsedWi
    
    ! Work arrays for ILU factors and 8 BCG vectors
    integer, parameter        :: MaxWr=maxnz+10*maxn, MaxWi=maxnz+7*maxn+1
    real*8                    :: rW(MaxWr)
    integer                   :: iW(MaxWi)
    
    ! residual vector
    integer                   :: resinit(system_size)
      
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Initialization of the preconditioner  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    verb = 0 ! verbose no
    tau1 = 1d-2
    tau2 = 1d-3
    partlur = 0.5
    ierr = 0
    
    call iluoo(system_size, Matr%ia, Matr%ja, Matr%a, tau1, tau2, verb, &
          rW, iW, MaxWr, MaxWi, partlur, partlurout, &
          UsedWr, UsedWi, ierr)
     
    ipBCG = UsedWr + 1 
     
    !set initial guess to 0
    call dcopy(system_size, 0d0, 0, result_vector, 1)
    
    !compute initial residual norm for 1 system
    resinit = ddot(number_of_non_Direchlet_elements,rhs_VP_1,1,rhs_VP_1,1)
      
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!     Iterative solution    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    ITER_1 = 1000             
    RESID_1 = 1d-20*resinit_1  
    INFO_1  = 0               
    NUNIT_1 = 6               
    iprevec_1(1) = number_of_non_Direchlet_elements           
    imatvec_1(1) = number_of_non_Direchlet_elements  
     
    ITER_2 = 1000             
    RESID_2 = 1d-20*resinit_2  
    INFO_2  = 0               
    NUNIT_2 = 6               
    iprevec_2(1) = number_of_non_Direchlet_elements           
    imatvec_2(1) = number_of_non_Direchlet_elements             
    
    call slpbcgs( &
        prevec2, iprevec_1, iW,rW, &
        matvec,  imatvec_1, lhs_VP%ia, lhs_VP%ja, lhs_VP%a, &
        rW(ipBCG), number_of_non_Direchlet_elements, 8, &
        number_of_non_Direchlet_elements, rhs_VP_1, solution_vec_1, &
        ITER_1, RESID_1, &
        INFO_1, NUNIT_1)
              
    call slpbcgs( &
       prevec2, iprevec_2, iW,rW, &
       matvec,  imatvec_2, lhs_VP%ia, lhs_VP%ja, lhs_VP%a, &
       rW(ipBCG), number_of_non_Direchlet_elements, 8, &
       number_of_non_Direchlet_elements, rhs_VP_2, solution_vec_2, &
       ITER_2, RESID_2, &
       INFO_2, NUNIT_2)    
        
    do i = 1, number_of_non_Direchlet_elements
     
      List_of_non_Direchlet_Elements(i)%pointing_element%u_old(1) = solution_vec_1(i)
      List_of_non_Direchlet_Elements(i)%pointing_element%u_old(2) = solution_vec_2(i)
       
    end do
  
  end subroutine ILU_0_BICG

end module linear_system_solver_module
