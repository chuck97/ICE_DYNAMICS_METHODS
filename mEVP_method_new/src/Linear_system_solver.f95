module module_linear_system_solver

  use module_classes
  use module_constants
  
  implicit none
  private :: residual_tolerance
  public  :: ILU_2_BICG_solver
  
  real*8, parameter  :: residual_tolerance = 1d-20
  
  contains
  
  subroutine ILU_2_BICG_solver(Matr, right_hand, system_size, result_vector)
  
    implicit none
    
    !! arguments:
    type(Matrix), intent(in)  :: Matr
    integer, intent(in)       :: system_size
    real*8, intent(in)        :: right_hand(system_size)
    real*8, intent(inout)     :: result_vector(system_size)
    
    !! local variables:
    integer                   :: i, j, k, ierr
    integer                   :: ipBCG
    
    !! ILU data:
    real*8                    :: tau1, tau2, partlur, partlurout
    integer                   :: verb, UsedWr, UsedWi
    
    !! Work arrays for ILU factors and 8 BCG vectors:
    integer, parameter        :: MaxWr=5*maxnz+10*maxn, MaxWi=5*maxnz+7*maxn+1
    real*8                    :: rW(MaxWr)
    integer                   :: iW(MaxWi)
    integer                   :: imatvec(1), iprevec(1) 
    
    !! residual vector:
    real*8                    :: resinit
    
    !! BiCGStab data:
    external                  :: matvec, prevec2
    integer                   :: ITER, INFO, NUNIT
    real*8                    :: RESID
    
    !! External routines from the BLAS library:
    real*8                    :: ddot
    external                  :: ddot, dcopy  
      
    
    
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
          
    if(ierr.ne.0) then
       write(*,'(A,I7)') 'Initialization of iluoo failed, ierr=', ierr
       stop 911
    end if      
     
    ipBCG = UsedWr + 1 
     
    !set initial guess to 0
    call dcopy(system_size, 0d0, 0, result_vector, 1)
    
    !compute initial residual norm
    resinit = ddot(system_size, right_hand, 1, right_hand, 1)
    resinit = dsqrt(resinit)
    
    if(resinit .le. varepsilon) then 
      write(*,'(A)') 'rhs=0, nothing to solve!'
      stop 911
    end if
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!     Iterative solution    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    ITER = 1000             
    RESID = residual_tolerance*resinit
    INFO  = 0               
    NUNIT = 6               
    iprevec(1) = system_size           
    imatvec(1) = system_size  
     
    call slpbcgs( &
        prevec2, iprevec, iW,rW, &
        matvec,  imatvec, Matr%ia(1:(system_size+1)), Matr%ja(1:(Matr%ia(system_size+1)-1)), &
        Matr%a(1:(Matr%ia(system_size+1)-1)), &
        rW(ipBCG), system_size, 8, &
        system_size, right_hand, result_vector, &
        ITER, RESID, &
        INFO, NUNIT)
        
    if(INFO.ne.0) then
      write(*,'(A)') 'BiCGStab failed'
      stop 911
    end if   
                  
  end subroutine ILU_2_BICG_solver

end module module_linear_system_solver
