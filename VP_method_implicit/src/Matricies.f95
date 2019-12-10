module module_matricies

  use module_constants
  use module_numerical_integration
  use module_Classes

  implicit none
  
  private
  
  public :: lumping, sparce_matrix_vector_multiplication, four_block_matrix_assembling, &
  print_matrix, Matr1_eqv_Matr2
  
  contains
  
  ! Mass lumping in CSR
  
  subroutine lumping(Matr)
  
    type(Matrix) :: Matr
    integer :: i, j, k, m_size
    real*8 :: prom_value
    real*8 :: a(maxnz)
    integer :: ia(maxn), ja(maxnz)
    integer :: nonzero
    
    ia(1) = 1
    m_size = Matr%matr_size
    
    do i = 1, m_size
    
      ia(i+1) = ia(i) + 1
      ja(i) = i
      prom_value = 0d0
      
      do  k = Matr%ia(i), (Matr%ia(i+1)-1)
      
        prom_value = prom_value + Matr%a(k)
      
      end do
      
      a(i) = prom_value
    
    end do
    
    call Matr%distr()
    call Matr%init(m_size, m_size)
    
    Matr%a = a
    Matr%ia = ia
    Matr%ja = ja
    
   
  end subroutine lumping
  
  
  ! CSR Matrix - vector multiplication function 
  
  function sparce_matrix_vector_multiplication(Matr, vector) result(mult_res)
  
    real*8 ::vector(nvmax), mult_res(nvmax)
    type(Matrix) :: Matr
    integer :: i, j, k
    
    k = 1
    
    do i = 1, Matr%matr_size
    
      mult_res(i) = 0d0
      
      do j = 1, (Matr%ia(i+1) - Matr%ia(i))
      
         mult_res(i) = mult_res(i) + vector(Matr%ja(k))*Matr%a(k)
         k = k + 1
      
      end do
    
    end do
  
  end function sparce_matrix_vector_multiplication 
  
  !! 4-block matrix assembling 
  
  subroutine four_block_matrix_assembling(M1, M2, M3, M4, M)
  
  implicit none
  
  type(Matrix), intent(in)     :: M1, M2, M3, M4
  type(Matrix), intent(inout)  :: M
  integer                      :: i, j, k, n, nz, ierr, nztot, nrowc, ncolc
							   
  integer                      :: nz1, nz2, nz3, nz4
							   
  integer, allocatable         :: ia(:), ja(:)
  real*8, allocatable          :: a(:)
							   
  integer, allocatable         :: ib(:), jb(:)
  real*8, allocatable          :: b(:)

  n = M1%matr_size
  
  nz1 = M1%nonzero
  nz2 = M2%nonzero
  nz3 = M3%nonzero
  nz4 = M4%nonzero
  
  nztot = (nz1 + nz2 + nz3 + nz4)*2
    
  allocate(ia(2*n + 3))
  allocate(ja(nztot))
  allocate(a(nztot))
  
  allocate(ib(2*n + 3))
  allocate(jb(nztot))
  allocate(b(nztot))
  
 
  !! set a equals to M1
  ia(1:(n+1)) = M1%ia(1:(M1%matr_size+1))
  ja(1: M1%nonzero) = M1%ja(1: M1%nonzero)
  a(1: M1%nonzero) = M1%a(1: M1%nonzero)
  
  !! add M2 n*n
        
  call addblk(n, n, a, ja, ia, 1, n + 1, 1, &
      n, n, M2%a(1:M2%nonzero), M2%ja(1:M2%nonzero), &
      M2%ia(1:(M2%matr_size+1)), nrowc, ncolc, b, jb, ib, nztot, ierr)
          
  !! add M3 n*n
  
  call addblk(n, 2*n, b, jb, ib, n + 1, 1, 1, &
      n, n, M3%a(1:M3%nonzero), M3%ja(1:M3%nonzero), &
      M3%ia(1:(M3%matr_size+1)), nrowc, ncolc, a, ja, ia, nztot, ierr)
         
  !! add M4 n*n
  
  call addblk(2*n, 2*n, a, ja, ia, n + 1, n + 1, 1, &
      n, n, M4%a(1:M4%nonzero), M4%ja(1:M4%nonzero), &
      M4%ia(1:(M4%matr_size+1)), nrowc, ncolc, b, jb, ib, nztot, ierr)    
      
  if ((M1%nonzero + M2%nonzero + M3%nonzero + M4%nonzero) .ne. (ib(nrowc+1)-1)) then
  
    print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print *, "Wrong 4-block matrix assembling"
    print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    stop
  
  end if
 
  call M%init((ib(nrowc+1)-1), 2*n)
  
  M%ia(1:(2*n+1)) = ib(1:(2*n+1))
  M%ja(1:M%nonzero) = jb(1:M%nonzero)
  M%a(1:M%nonzero) = b(1:M%nonzero)
          
  
  deallocate(ia)
  deallocate(ja)
  deallocate(a)
  
  deallocate(ib)
  deallocate(jb)
  deallocate(b) 
  
  end subroutine four_block_matrix_assembling
  
  
  !! print CSR-matrix subroutine
  
  subroutine print_matrix(Matr)
  
  type(Matrix) :: Matr
  integer :: i,j,k
  
  print *, "Matrix size:", Matr%matr_size
  print *, "Nonzero:", Matr%nonzero
  print *, "a:", Matr%a(1:Matr%nonzero)
  print *, "ia:", Matr%ia(1:(Matr%matr_size+1))
  print *, "ja:", Matr%ja(1:Matr%nonzero)
  
   
  end subroutine print_matrix
  
!! Matr1 eqv Matr2 subroutine
  
  subroutine Matr1_eqv_Matr2(Matr1, Matr2)
  
    type(Matrix) :: Matr1, Matr2
    
    Matr1%nonzero = Matr2%nonzero
    Matr1%matr_size = Matr2%matr_size
    Matr1%ia = Matr2%ia
    Matr1%a = Matr2%a
    Matr1%ja = Matr2%ja
    
  
  end subroutine Matr1_eqv_Matr2
                                                                                         
end module module_matricies


