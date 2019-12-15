module module_matricies

  use module_constants
  use module_Classes

  implicit none
  
  private
  
  public :: lumping_square,                        &
            Matr1_eqv_Matr2,                       &
            sparce_matrix_vector_multiplication,   &
            big_block_matrix_assembling,           &
            frobenius_norm
  
  contains
  
  ! Mass lumping in CSR
  
  subroutine lumping_square(Matr)
  
    implicit none
  
    type(Matrix),        intent(inout) :: Matr
    integer                            :: i, k, m_size
    real*8                             :: prom_value
    integer,   allocatable             :: ia(:), ja(:)
    real*8,    allocatable             :: a(:)
    
    allocate(a(Matr%nonzero))
    allocate(ia(Matr%matr_size_y + 1))
    allocate(ja(Matr%nonzero))
    
    ia(1) = 1
    m_size = Matr%matr_size_x
    
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
    
    call Matr%init(m_size, m_size, m_size)
    
    Matr%a(1:m_size) = a(1:m_size)
    Matr%ia(1:(m_size+1)) = ia(1:(m_size+1))
    Matr%ja(1:m_size) = ja(1:m_size)
    
    deallocate(a)
    deallocate(ia)
    deallocate(ja)
   
  end subroutine lumping_square

!! Matr1 eqv Matr2 subroutine
  
  subroutine Matr1_eqv_Matr2(Matr1, Matr2)
  
    implicit none
  
    type(Matrix)    ,intent(inout) :: Matr1
    type(Matrix)    ,intent(in)    :: Matr2
    
    
    Matr1%nonzero = Matr2%nonzero
    Matr1%matr_size_x = Matr2%matr_size_x
    Matr1%matr_size_y = Matr2%matr_size_y
    Matr1%ia = Matr2%ia
    Matr1%a = Matr2%a
    Matr1%ja = Matr2%ja
    
  
  end subroutine Matr1_eqv_Matr2
  
  
  subroutine sparce_matrix_vector_multiplication(matr, vector, size_of_vec, result_vec) 
  
    implicit none 
  
    real*8,                intent(in)    :: vector(:)
    real*8,                intent(inout) :: result_vec(:)
    integer,               intent(in)    :: size_of_vec 
    type(Matrix),          intent(in)    :: matr
  
    if (size_of_vec .ne. matr%matr_size_x) then
  
      print *, "WRONG MATRIX OR VECTOR SIZE"
  
    end if
    
    call amux(matr%matr_size_y, vector(1:size_of_vec), result_vec, matr%a, matr%ja, matr%ia)
  
  end subroutine sparce_matrix_vector_multiplication
  
  ! 5*5 block matrix
  
  subroutine big_block_matrix_assembling(M1, M2, M3, M4, M5, &
                                         M6, M7, M8, M9, M10, &
                                         M11, M12, M13, M14, M15, &
                                         M16, M17, M18, M19, M20, &
                                         M21, M22, M23, M24, M25, Matr_res, notr, nov)
    implicit none                                     
                                         
    type(Matrix), intent(in)    :: M1, M2, M3, M4, M5, &
                                         M6, M7, M8, M9, M10, &
                                         M11, M12, M13, M14, M15, &
                                         M16, M17, M18, M19, M20, &
                                         M21, M22, M23, M24, M25
    type(Matrix), intent(inout) :: Matr_res                                     
    integer, intent(in)         :: notr, nov
    integer                     :: nzmx, ierr, nrowc, ncolc, matr_size_tot, nonzero_tot
    integer, allocatable        :: ia(:), ib(:), ja(:), jb(:)
    real*8,  allocatable        :: a(:), b(:)                                           
    
    nonzero_tot = M1%nonzero + M2%nonzero + M3%nonzero + M4%nonzero + M5%nonzero + &
                  M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
                  M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
                  M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
                  M21%nonzero + M22%nonzero + M23%nonzero + M24%nonzero + M25%nonzero
                  
    matr_size_tot = 3*notr + 2*nov
    
    !nzmx = nonzero_tot
    
    allocate(a(nonzero_tot))
    allocate(b(nonzero_tot))
    allocate(ia(matr_size_tot+1))
    allocate(ib(matr_size_tot+1))
    allocate(ja(nonzero_tot))
    allocate(jb(nonzero_tot))
    
    ! set a and b to zero
    
    ia(1:(nov+1)) = M1%ia(1:(M1%matr_size_y+1))
    ja(1: M1%nonzero) = M1%ja(1: M1%nonzero)
    a(1: M1%nonzero) = M1%a(1: M1%nonzero)
    
    nzmx = 2*nonzero_tot
          
    ! add M2 nov*nov
        
    call addblk(nov, nov, a, ja, ia, 1, nov + 1, 1, &
      M2%matr_size_y, M2%matr_size_x, M2%a(1:M2%nonzero), M2%ja(1:M2%nonzero), &
      M2%ia(1:(M2%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc - nov, "ncolc", ncolc - 2*nov
    print *, "ierr", ierr
    print *, "check 1", (M1%nonzero + M2%nonzero) - (ib(nrowc+1)-1)
    
    ! add M3 nov*not
    
    call addblk(nov, 2*nov, b, jb, ib, 1, 2*nov+1, 1, &
      M3%matr_size_y, M3%matr_size_x, M3%a(1:M3%nonzero), M3%ja(1:M3%nonzero), &
      M3%ia(1:(M3%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 2", (M1%nonzero + M2%nonzero+ M3%nonzero) - (ia(nrowc+1)-1)
    
    ! add M4 nov*not
    
    call addblk(nov, 2*nov + notr, a, ja, ia, 1, 2*nov+ notr + 1, 1, &
      M4%matr_size_y, M4%matr_size_x, M4%a(1:M4%nonzero), M4%ja(1:M4%nonzero), &
      M4%ia(1:(M4%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 3", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero) - (ib(nrowc+1)-1)
    
    ! add M5 nov*not
    
    call addblk(nov, 2*nov + 2*notr, b, jb, ib, 1, 2*nov+ 2*notr + 1, 1, &
      M5%matr_size_y, M5%matr_size_x, M5%a(1:M5%nonzero), M5%ja(1:M5%nonzero), &
      M5%ia(1:(M5%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 4", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero) - (ia(nrowc+1)-1)
    
    ! add M6 nov*nov
    
    call addblk(nov, 2*nov + 3*notr, a, ja, ia, nov + 1, 1, 1, &
      M6%matr_size_y, M6%matr_size_x, M6%a(1:M6%nonzero), M6%ja(1:M6%nonzero), &
      M6%ia(1:(M6%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 5", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero) - (ib(nrowc+1)-1)
     
    ! add M7 nov*nov
    
    call addblk(2*nov, 2*nov + 3*notr, b, jb, ib, nov + 1, nov + 1, 1, &
      M7%matr_size_y, M7%matr_size_x, M7%a(1:M7%nonzero), M7%ja(1:M7%nonzero), &
      M7%ia(1:(M7%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 6", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero) - (ia(nrowc+1)-1) 
     
   ! add M8 nov*not
    
    call addblk(2*nov, 2*nov + 3*notr, a, ja, ia, nov + 1, 2*nov + 1, 1, &
      M8%matr_size_y, M8%matr_size_x, M8%a(1:M8%nonzero), M8%ja(1:M8%nonzero), &
      M8%ia(1:(M8%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 7", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero) - (ib(nrowc+1)-1)
     
   ! add M9 nov*not
    
    call addblk(2*nov, 2*nov + 3*notr, b, jb, ib, nov + 1, 2*nov + notr + 1, 1, &
      M9%matr_size_y, M9%matr_size_x, M9%a(1:M9%nonzero), M9%ja(1:M9%nonzero), &
      M9%ia(1:(M9%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 8", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero) - (ia(nrowc+1)-1)     
     
     
   ! add M10 nov*not
    
    call addblk(2*nov, 2*nov + 3*notr, a, ja, ia, nov + 1, 2*nov + 2*notr + 1, 1, &
      M10%matr_size_y, M10%matr_size_x, M10%a(1:M10%nonzero), M10%ja(1:M10%nonzero), &
      M10%ia(1:(M10%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 9", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero) - (ib(nrowc+1)-1)   
     
   ! add M11 not*nov
    
    call addblk(2*nov, 2*nov + 3*notr, b, jb, ib, 2*nov + 1, 1, 1, &
      M11%matr_size_y, M11%matr_size_x, M11%a(1:M11%nonzero), M11%ja(1:M11%nonzero), &
      M11%ia(1:(M11%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 10", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero) - (ia(nrowc+1)-1)  
     
   ! add M12 not*nov
    
    call addblk(2*nov + notr, 2*nov + 3*notr, a, ja, ia, 2*nov + 1, nov + 1, 1, &
      M12%matr_size_y, M12%matr_size_x, M12%a(1:M12%nonzero), M12%ja(1:M12%nonzero), &
      M12%ia(1:(M12%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 11", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero) - (ib(nrowc+1)-1)    
     
   ! add M13 not*not
    
    call addblk(2*nov + notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 1, 2*nov + 1, 1, &
      M13%matr_size_y, M13%matr_size_x, M13%a(1:M13%nonzero), M13%ja(1:M13%nonzero), &
      M13%ia(1:(M13%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 12", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero) - (ia(nrowc+1)-1)   
     
   ! add M14 not*not
    
    call addblk(2*nov + notr, 2*nov + 3*notr, a, ja, ia, 2*nov + 1, 2*nov + notr + 1, 1, &
      M14%matr_size_y, M14%matr_size_x, M14%a(1:M14%nonzero), M14%ja(1:M14%nonzero), &
      M14%ia(1:(M14%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 13", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero) - (ib(nrowc+1)-1)    
     
   ! add M15 not*not
    
    call addblk(2*nov + notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 1, 2*nov + 2*notr + 1, 1, &
      M15%matr_size_y, M15%matr_size_x, M15%a(1:M15%nonzero), M15%ja(1:M15%nonzero), &
      M15%ia(1:(M15%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 14", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero) - (ia(nrowc+1)-1)  
     
   ! add M16 not*nov
    
    call addblk(2*nov + notr, 2*nov + 3*notr, a, ja, ia, 2*nov + notr+ 1, 1, 1, &
      M16%matr_size_y, M16%matr_size_x, M16%a(1:M16%nonzero), M16%ja(1:M16%nonzero), &
      M16%ia(1:(M16%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 15", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero) - (ib(nrowc+1)-1)  
     
    ! add M17 not*nov
    
    call addblk(2*nov + 2*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + notr + 1, nov +  1, 1, &
      M17%matr_size_y, M17%matr_size_x, M17%a(1:M17%nonzero), M17%ja(1:M17%nonzero), &
      M17%ia(1:(M17%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 16", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero) - (ia(nrowc+1)-1)   
     
    ! add M18 not*not
    
    call addblk(2*nov + 2*notr, 2*nov + 3*notr, a, ja, ia, 2*nov + notr + 1, 2*nov +  1, 1, &
      M18%matr_size_y, M18%matr_size_x, M18%a(1:M18%nonzero), M18%ja(1:M18%nonzero), &
      M18%ia(1:(M18%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 17", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero) - (ib(nrowc+1)-1)  
     
   ! add M19 not*not
    
    call addblk(2*nov + 2*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + notr + 1, 2*nov + notr + 1, 1, &
      M19%matr_size_y, M19%matr_size_x, M19%a(1:M19%nonzero), M19%ja(1:M19%nonzero), &
      M19%ia(1:(M19%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 18", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero) - (ia(nrowc+1)-1) 
     
   ! add M20 not*not
    
    call addblk(2*nov + 2*notr, 2*nov + 3*notr, a, ja, ia, 2*nov + notr + 1, 2*nov + 2*notr + 1, 1, &
      M20%matr_size_y, M20%matr_size_x, M20%a(1:M20%nonzero), M20%ja(1:M20%nonzero), &
      M20%ia(1:(M20%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 19", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero) - (ib(nrowc+1)-1)    
    
    ! add M21 not*nov
    
    call addblk(2*nov + 2*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 2*notr + 1, 1, 1, &
      M21%matr_size_y, M21%matr_size_x, M21%a(1:M21%nonzero), M21%ja(1:M21%nonzero), &
      M21%ia(1:(M21%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 20", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
     M21%nonzero) - (ia(nrowc+1)-1)  
     
    ! add M22 not*nov
    
    call addblk(2*nov + 3*notr, 2*nov + 3*notr, a, ja, ia, 2*nov + 2*notr + 1, nov +  1, 1, &
      M22%matr_size_y, M22%matr_size_x, M22%a(1:M22%nonzero), M22%ja(1:M22%nonzero), &
      M22%ia(1:(M22%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 21", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
     M21%nonzero + M22%nonzero) - (ib(nrowc+1)-1)   
     
     
    ! add M23 not*not
    
    call addblk(2*nov + 3*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 2*notr + 1, 2*nov +  1, 1, &
      M23%matr_size_y, M23%matr_size_x, M23%a(1:M23%nonzero), M23%ja(1:M23%nonzero), &
      M23%ia(1:(M23%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 22", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
     M21%nonzero + M22%nonzero + M23%nonzero) - (ia(nrowc+1)-1)  
     
     
    ! add M24 not*not
    
    call addblk(2*nov + 3*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 2*notr + 1, 2*nov + notr + 1, 1, &
      M24%matr_size_y, M24%matr_size_x, M24%a(1:M24%nonzero), M24%ja(1:M24%nonzero), &
      M24%ia(1:(M24%matr_size_y+1)), nrowc, ncolc, b, jb, ib, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 23", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
     M21%nonzero + M22%nonzero + M23%nonzero + M24%nonzero) - (ib(nrowc+1)-1)  
     
     
    ! add M25 not*not
    
    call addblk(2*nov + 3*notr, 2*nov + 3*notr, b, jb, ib, 2*nov + 2*notr + 1, 2*nov + 2*notr + 1, 1, &
      M25%matr_size_y, M25%matr_size_x, M25%a(1:M25%nonzero), M25%ja(1:M25%nonzero), &
      M25%ia(1:(M25%matr_size_y+1)), nrowc, ncolc, a, ja, ia, nzmx, ierr)
    
    print *, "nrowc", nrowc, "ncolc", ncolc
    print *, "ierr", ierr
    print *, "check 24", (M1%nonzero + M2%nonzero+ M3%nonzero + M4%nonzero + M5%nonzero + &
     M6%nonzero + M7%nonzero + M8%nonzero + M9%nonzero + M10%nonzero + &
     M11%nonzero + M12%nonzero + M13%nonzero + M14%nonzero + M15%nonzero + &
     M16%nonzero + M17%nonzero + M18%nonzero + M19%nonzero + M20%nonzero + &
     M21%nonzero + M22%nonzero + M23%nonzero + M24%nonzero + M25%nonzero) - (ia(nrowc+1)-1)  
     
    call Matr_res%init(ia(nrowc+1)-1, nrowc, ncolc)
    
    Matr_res%ia(1:(Matr_res%matr_size_y+1)) = ia(1:(nrowc+1))
    Matr_res%ja(1:Matr_res%nonzero) = ja(1:(ia(nrowc+1)-1))
    Matr_res%a(1:Matr_res%nonzero) = a(1:(ia(nrowc+1)-1))
    
    
    ! frobenius norm of blocks
    
    print *
    print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print *
    print *, "M1 frob:", frobenius_norm(M1)
    print *, "M2 frob:", frobenius_norm(M2)
    print *, "M3 frob:", frobenius_norm(M6)
    print *, "M4 frob:", frobenius_norm(M7)
    print *
    print *, "R1 frob:", frobenius_norm(M3)
    print *, "R2 frob:", frobenius_norm(M4)
    print *, "R3 frob:", frobenius_norm(M5)
    print *, "R4 frob:", frobenius_norm(M8)
    print *, "R5 frob:", frobenius_norm(M9)
    print *, "R6 frob:", frobenius_norm(M10)
    print *
    print *, "N1 frob:", frobenius_norm(M11)
    print *, "N2 frob:", frobenius_norm(M12)
    print *, "N3 frob:", frobenius_norm(M16)
    print *, "N4 frob:", frobenius_norm(M17)
    print *, "N5 frob:", frobenius_norm(M21)
    print *, "N6 frob:", frobenius_norm(M22)
    print * 
    print *, "S frob:", frobenius_norm(M13)
    print *
    print *, "Z frob:", frobenius_norm(M14)
    print *
    print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print *
    print *, "result block Matrix frob:", frobenius_norm(Matr_res)
    print *
    print *, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print *
    
    
    deallocate(a)
    deallocate(b)
    deallocate(ia)
    deallocate(ib)
    deallocate(ja)
    deallocate(jb)
                  
  
  end subroutine big_block_matrix_assembling      
  
  ! frobenius norm of matrix
  
  function frobenius_norm(Matr) result(frb)   
  
    implicit none
  
    type(Matrix), intent(in)  :: Matr
    real*8                    :: frb
    integer                   :: i, j
    real*8                    :: prom, max_row
    
    max_row = -1d10
    
    do i = 1, Matr%matr_size_y
    
      prom = 0d0
    
      do j = Matr%ia(i), (Matr%ia(i+1)-1)
      
        prom = prom + abs(Matr%a(j)) 
      
      end do
      
      if (prom > max_row) then
      
        max_row = prom
      
      end if
    
    end do
    
    frb = max_row
  
  end function frobenius_norm             
  
  
  
                                                                                         
end module module_matricies


