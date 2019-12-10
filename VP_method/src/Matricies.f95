module module_matricies

  use module_constants
  use module_numerical_integration
  use module_Classes

  implicit none
  
  private
  
  public :: lumping, sparce_matrix_vector_multiplication, sparce_matrix_mult_number, &
  sparse_matrix_plus_sparse_matrix, comparisson, sum_rows, four_block_matrix_assembling, &
  print_matrix, sparce_matrix_mult_diagonal_matrix, &
  Matr1_eqv_Matr2
  
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
  
  
  ! CSR Matrix - number multiplication subroutine
  
  subroutine sparce_matrix_mult_number(Matr, numb) 
  
    type(Matrix) :: Matr
    real*8 :: numb
    integer :: i, j, k 
    
    do i = 1, Matr%nonzero
    
      Matr%a(i) = Matr%a(i)*numb
    
    end do
  
  end subroutine sparce_matrix_mult_number
  
  
  
  ! sparse_matrix_plus_sparse_matrix subroutine
  
  subroutine sparse_matrix_plus_sparse_matrix(Matr_first, Matr_second, Matr_res)
  
    integer :: pos_1(maxn), pos_2(maxn), pos_sum(maxn), len_1, len_2, len_sum, i, j, k, sum_ind, &
     ind1, ind2, indsum, indicator, end1, end2, system_size
    real*8 :: val_1(maxn), val_2(maxn), val_sum(maxn)
    real*8 :: a(maxnz)
    integer :: ia(maxn), ja(maxnz)
    integer :: nonzero
    type(Matrix) :: Matr_first, Matr_second, Matr_res
    
    ia(1) = 1
    nonzero = 0
    
    do i = 1, Matr_first%matr_size
      
      call sum_rows(Matr_first%ja(Matr_first%ia(i):(Matr_first%ia(i+1)-1)), &
       Matr_second%ja(Matr_second%ia(i):(Matr_second%ia(i+1)-1)), &
       pos_sum, Matr_first%a(Matr_first%ia(i):(Matr_first%ia(i+1)-1)), &
       Matr_second%a(Matr_second%ia(i):(Matr_second%ia(i+1)-1)), &
       val_sum, (Matr_first%ia(i+1) - Matr_first%ia(i)), (Matr_second%ia(i+1) - Matr_second%ia(i)), len_sum)
        
      nonzero = nonzero + len_sum
      ia(i+1) = ia(i) + len_sum
      ja(ia(i):(ia(i+1)-1)) = pos_sum(1:len_sum)
      a(ia(i):(ia(i+1)-1)) = val_sum(1:len_sum)
      
    end do
    
    call Matr_res%init(nonzero, Matr_first%matr_size)
    Matr_res%a(1:nonzero) = a(1:nonzero)
    Matr_res%ia(1:(Matr_first%matr_size+1)) = ia(1:(Matr_first%matr_size+1))
    Matr_res%ja(1:nonzero) = ja(1:nonzero)
    
  end subroutine sparse_matrix_plus_sparse_matrix
  
  
  
  ! comparisson two arrays elements subroutine
  
  subroutine comparisson(position_1, position_2, value_1, value_2, value_res, indicator) ! indicator = 1, if position_1 < position_2
                                                                                         ! indicator = 2, if position_2 < position_1
    integer :: position_1, position_2, indicator                                         ! indicator = 0, if position_1 = position_2
    real*8 :: value_1, value_2, value_res  
    
    if (position_1 < position_2) then
    
      indicator = 1
      value_res = value_1
      
    end if
    
    if (position_2 < position_1) then
    
      indicator = 2
      value_res = value_2
     
    end if
    
    if (position_1 == position_2) then
    
      indicator = 0
      value_res = value_1 + value_2
    
    end if
  
  end subroutine comparisson
  
  
  ! sum of two arrays (rows)
  
  subroutine sum_rows(pos_1, pos_2, pos_sum, val_1, val_2, val_sum, len_1, len_2, len_sum)
  
    integer :: pos_1(maxn), pos_2(maxn), pos_sum(maxn), len_1, len_2, len_sum, i, j, k, sum_ind, &
     ind1, ind2, indsum, indicator, end1, end2, bound1, bound2
    real*8 :: val_1(maxn), val_2(maxn), val_sum(maxn), value_res
    
    ind1 = 1
    ind2 = 1
    indsum = 1
    end1 = 0
    end2 = 0
    bound1 = pos_1(ind1)
    bound2 = pos_2(ind2)
    
    do while ((end1 == 0) .or. (end2 == 0))
    
      call comparisson(bound1, bound2, val_1(ind1), val_2(ind2), value_res, indicator)
      
      if (indicator == 1) then
      
        val_sum(indsum) = value_res
        pos_sum(indsum) = pos_1(ind1)
        ind1 = ind1 + 1
        indsum = indsum + 1
        
      end if
      
      if (indicator == 2) then
      
        val_sum(indsum) = value_res
        pos_sum(indsum) = pos_2(ind2)
        ind2 = ind2 + 1
        indsum = indsum + 1
        
      end if
      
      if (indicator == 0) then
      
        val_sum(indsum) = value_res
        pos_sum(indsum) = pos_1(ind1)
        ind1 = ind1 + 1
        ind2 = ind2 + 1
        indsum = indsum + 1
        
      end if
      
      if (ind1 >= (len_1 + 1)) then
      
        end1 = 1
        bound1 = maxn + 10
        
      else
      
        bound1 = pos_1(ind1)  
      
      end if
    
      if (ind2 >= (len_2 + 1)) then
      
        end2 = 1
        bound2 = maxn + 10
        
      else
      
        bound2 = pos_2(ind2)    
      
      end if
    
    end do
    
    len_sum = indsum - 1
  
  end subroutine sum_rows
  
  
  
  !! 4-block matrix assembling 
  
  subroutine four_block_matrix_assembling(Matr1, Matr2, Matr3, Matr4, Matr_res)
  
  type(Matrix) :: Matr1, Matr2, Matr3, Matr4, Matr_res
  integer :: i, j, k, len1, len2
  real*8 :: a_sum(maxnz)
  integer :: ja_sum(maxnz)
  integer :: point, n
  
  point = 1
  Matr_res%matr_size = 2*Matr1%matr_size
  Matr_res%nonzero = Matr1%nonzero + Matr2%nonzero + Matr3%nonzero + Matr4%nonzero
  Matr_res%ia(1) = 1
  n = Matr1%matr_size
  
  do i = 1, n
  
    len1 = Matr1%ia(i+1) - Matr1%ia(i)
    len2 = Matr2%ia(i+1) - Matr2%ia(i)
    Matr_res%a(point:(point+len1-1)) = Matr1%a(Matr1%ia(i):(Matr1%ia(i+1)-1))
    Matr_res%ja(point:(point+len1-1)) = Matr1%ja(Matr1%ia(i):(Matr1%ia(i+1)-1))
    point = point + len1
    Matr_res%a(point:(point+len2-1)) = Matr2%a(Matr2%ia(i):(Matr2%ia(i+1)-1))
    do k = 0, (len2-1)
    
      Matr_res%ja(point+k) = n + Matr2%ja(Matr2%ia(i)+k)
    
    end do
    point = point + len2
    Matr_res%ia(i+1) = Matr_res%ia(i) + len1 + len2
    
  end do
  
  do i = 1, n
  
    len1 = Matr3%ia(i+1) - Matr3%ia(i)
    len2 = Matr4%ia(i+1) - Matr4%ia(i)
    Matr_res%a(point:(point+len1-1)) = Matr3%a(Matr1%ia(i):(Matr3%ia(i+1)-1))
    Matr_res%ja(point:(point+len1-1)) = Matr3%ja(Matr3%ia(i):(Matr3%ia(i+1)-1))
    point = point + len1
    Matr_res%a(point:(point+len2-1)) = Matr4%a(Matr4%ia(i):(Matr4%ia(i+1)-1))
    do k = 0, (len2-1)
    
      Matr_res%ja(point+k) = n + Matr4%ja(Matr4%ia(i)+k)
    
    end do
    point = point + len2
    Matr_res%ia(n+i+1) = Matr_res%ia(n+i) + len1 + len2
    
  end do
  
  
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
  
  !! sparce_matrix_mult_diagonal_matrix subroutine
  
  subroutine sparce_matrix_mult_diagonal_matrix(Matr, vector)
  
    type(Matrix) :: Matr
    real*8 :: vector(nvmax)
    integer :: i, j, k
    
    do i = 1, Matr%matr_size
    
      do k = Matr%ia(i), (Matr%ia(i+1) - 1)
      
        Matr%a(k) = Matr%a(k)*vector(i)
      
      end do 
    
    end do
  
  end subroutine sparce_matrix_mult_diagonal_matrix
  
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


