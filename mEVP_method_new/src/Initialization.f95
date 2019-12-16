module module_initialization

  use module_Classes
  use module_constants
  use module_grid_values
  
  implicit none
  private :: append 
  public  :: grid_initialization,           &
             vectors_initialization,        &
             scalars_initialization
  
  contains
  
  !! append element to array
  subroutine append(edg, first, second, current_size)

    implicit none

    !! arguments:
    integer, intent(inout) :: edg(2, ((2*nvmax)*3)/2 + nbmax)
    integer, intent(inout) :: current_size
    integer, intent(in) :: first, second
    
    !! local variables:
    integer             :: i
    
    do i = 1, current_size
      if((first == edg(1,i)) .and. (second == edg(2, i))) then
        return
      end if
    end do
    edg(1, current_size+1) = first
    edg(2, current_size+1) = second
    current_size = current_size + 1
    
  end subroutine append 
  
  !! This routine generates mesh using Ani-2D Library and
  !! and assembles lists of elements/non Direchlet elements/edges/triangles 
  subroutine grid_initialization(Nbv, Nbl, bv, bl, bltail, &
                   number_of_vertex_per_square_edge, square_size)

  implicit none
  
  !! arguments:
  integer, intent(in)           :: Nbv, Nbl, number_of_vertex_per_square_edge
  real*8, intent(in)            :: bv(2,Nbv),bltail(2,12)
  integer, intent(in)           :: bl(7,Nbl)
  real*8, intent(in)            :: square_size
  
  
  !! local variables:
  integer            :: i, j, k, l
  integer            :: ierr
  integer            :: edg(2, (ntmax*3)/2 + nbmax)
  real*8             :: h, a, b, c, p, S, squ            
  real*8             :: v1, v2, R, x1, x2, dist
  integer            :: iW(3*ntmax + 2*nbmax + nvmax + 10), nPP(nvmax*2), IPP(nvmax*10), IEE(3, ntmax) 
  integer            :: nEP(nvmax), IEP(3*ntmax), IRE(3, ntmax)
  integer            :: IPE(3, ntmax), map_tri(3, ntmax), on_bond(nvmax), neighb_tr, neighb_edg
  integer            :: nv,nt,nb,nc
  real*8             :: vrt(2,nvmax)
  integer            ::  tri(3,ntmax), labelT(ntmax)
  integer            ::  bnd(2,nbmax), labelB(nbmax)
  real*8             :: crv(2,nbmax)
  integer            :: iFNC(nbmax)
  integer            :: aft2dboundary
  external           :: aft2dboundary
  external              listE2E        
  external              orientBoundary  
  external              listP2P         
  external              backReferences   
  external              listConv       
  
  h = square_size/number_of_vertex_per_square_edge
                                                         
  !! ====================  Building quazi-uniform triangular mesh  ====================== !!
  
  ierr = aft2dboundary(Nbv, bv, Nbl, bl, bltail, h, &  ! geometry
                        nv, vrt, &                     ! coordinates of nodes
                        nt, tri, labelT, &             ! triangles and their material
                        nb, bnd, labelB, &             ! boundary edges and their labels
                        nc, crv, iFNC)                 ! curved edges and parametrization (dummy)
  if (ierr.ne.0) stop ' error in function aft2dboundary'

  write(*,*) 'mesh: number of triangles/vertices ', nt, nv
  
  number_of_triangles = nt
  number_of_elements = nv
   
  if(nv.gt.nvmax) stop ' too many nodes'
  if(nt.gt.ntmax) stop ' too many triangles'
  if(nb.gt.nbmax) stop ' too many boundary edges'
  
  !! ==================================================================================== !!
  
  !! =====================  Drawing the grid  =========================================== !!
  
  call graph(nv,vrt, nt,tri, 'mesh_final.ps') ! drawing the grid
  
  !! ==================================================================================== !!
  
  !! .....................  Making bondary edges orientation  ........................... !!
  
  call orientBoundary(nv, nb, nt, vrt, bnd, tri, iW, 3*ntmax + 2*nbmax + nvmax + 10)
  
  !! .................................................................................... !!
  
  !! ............................ Making edges list .................................. !!
 
  k = 0
 
  do i = 1, (ntmax*3)/2 + nbmax
    edg(1, i) = 0
    edg(2, i) = 0
  end do
  
  do i = 1, nt
     if (tri(1,i) < tri(2,i)) then
       call append(edg, tri(1,i), tri(2,i),k)
     else
       call append(edg, tri(2,i), tri(1,i),k)  
     end if
     if (tri(1,i) < tri(3,i)) then
       call append(edg, tri(1,i), tri(3,i),k)
     else
       call append(edg, tri(3,i), tri(1,i),k)  
     end if
     if (tri(2,i) < tri(3,i)) then
       call append(edg, tri(2,i), tri(3,i),k)
     else
       call append(edg, tri(3,i), tri(2,i),k)  
     end if
  end do   
  
  number_of_edges = k
  k = 0
  
  !! ................................................................................. !! 
  
  !!  .................. Making identification list of boundary points .................  !!
  
  do i = 1, nv
      on_bond(i) = 0
  end do
  
  do i = 1, nb
    on_bond(bnd(1,i)) = 1
    on_bond(bnd(2,i)) = 1  
  end do
  
  !!  ..................................................................................  !!
  
  ! ................................. Making points connectivity list .............. !!
    
  call listP2P(nv, nt, nvmax*10, tri, nPP, IPP, iW)
    
  ! ................................................................................ !!
  ! ................ Making triangles connectivity list for all points .............. !!
    
  call backReferences(nv, nt, 3, 3, tri, nEP, IEP)  
  ! ................................................................................. !!
  !! ===================  Assembling elements list ====================================  !!
  
  do i = 1, number_of_elements
    List_of_Elements(i)%identificator = i
    List_of_Elements(i)%coordinates(1) = vrt(1,i)    
    List_of_Elements(i)%coordinates(2) = vrt(2,i)  
    if (on_bond(i) == 1) then                   
      List_of_Elements(i)%on_boundary = .true.
    else
      List_of_Elements(i)%on_boundary = .false. 
    end if 
     
     
    if(i == 1) then
      List_of_Elements(i)%number_of_neighbour_elements = nPP(1) 
      do j = 1, List_of_Elements(1)%number_of_neighbour_elements
        List_of_Elements(1)%neighbour_elements_list(j)%pointing_element => List_of_Elements(IPP(j)) 
      end do
    end if
    
    ! other points
    if (i /= 1) then
      List_of_Elements(i)%number_of_neighbour_elements = nPP(i) - nPP(i-1)
      do j = 1, List_of_Elements(i)%number_of_neighbour_elements
          List_of_Elements(i)%neighbour_elements_list(j)%pointing_element => List_of_Elements(IPP(nPP(i-1) + j)) 
      end do  
    end if  
   
    ! first point 
    if (i == 1) then
      List_of_Elements(i)%number_of_neighbour_triangles = nEP(1) 
      do j = 1, List_of_Elements(1)%number_of_neighbour_triangles
        List_of_Elements(1)%neighbour_triangles_list(j)%pointing_triangle => List_of_Triangles(IEP(j)) 
      end do
    end if
    
    ! other points
    if (i /= 1) then
      List_of_Elements(i)%number_of_neighbour_triangles = nEP(i) - nEP(i-1)
      do j = 1, List_of_Elements(i)%number_of_neighbour_triangles
          List_of_Elements(i)%neighbour_triangles_list(j)%pointing_triangle => List_of_Triangles(IEP(nEP(i-1)+j)) 
      end do  
    end if  
    
    
    ! neighbour edges
    neighb_edg = 0
    do j = 1, number_of_edges
      if ((edg(1, j) == i) .or. (edg(2, j) == i)) then
        neighb_edg = neighb_edg + 1
        List_of_Elements(i)%neighbour_edges_list(neighb_edg)%pointing_edge => List_of_Edges(j)
      end if
    end do
    List_of_Elements(i)%number_of_neighbour_edges = neighb_edg
  end do
  
  !! ====================================================================================== !!
  
  !! ===================  Assembling edges list ====================================  !! 
  
 
  do i = 1, number_of_edges
    List_of_Edges(i)%identificator = i
    List_of_Edges(i)%coordinates(1,1) = List_of_Elements(edg(1,i))%coordinates(1) 
    List_of_Edges(i)%coordinates(2,1) = List_of_Elements(edg(1,i))%coordinates(2)
    List_of_Edges(i)%coordinates(1,2) = List_of_Elements(edg(2,i))%coordinates(1)
    List_of_Edges(i)%coordinates(2,2) = List_of_Elements(edg(2,i))%coordinates(2)
   
   if ((on_bond(edg(1,i)) == 1) .and. (on_bond(edg(2,i)) == 1)) then                   
      List_of_Edges(i)%on_boundary = .true.   
    else                                         
      List_of_Edges(i)%on_boundary = .false.    
    end if
    
   List_of_Edges(i)%length_of_edge = dsqrt((List_of_Edges(i)%coordinates(1,1) - List_of_Edges(i)%coordinates(1,2))**2 + &
      (List_of_Edges(i)%coordinates(2,1) - List_of_Edges(i)%coordinates(2,2))**2) 
   List_of_Edges(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(edg(1,i)) 
   List_of_Edges(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(edg(2,i))   
 end do
  
  
 !! ===============================================================================  !!
  
  !! .......................... Making triangles connectivity list for all triangles ............. !!
  call listE2E(nv, nT, tri, IEE, nEP, IEP)
  !! ............................................................................................. !!
  
  !! ===========================  Assembling triangles list ===================================== !!
  
  do i = 1, number_of_triangles
    List_of_Triangles(i)%identificator = i
    List_of_Triangles(i)%coordinates(1,1) = vrt(1,tri(1,i)) !x1
    List_of_Triangles(i)%coordinates(2,1) = vrt(2,tri(1,i)) !y1
    List_of_Triangles(i)%coordinates(1,2) = vrt(1,tri(2,i)) !x2
    List_of_Triangles(i)%coordinates(2,2) = vrt(2,tri(2,i)) !y2
    List_of_Triangles(i)%coordinates(1,3) = vrt(1,tri(3,i)) !x3
    List_of_Triangles(i)%coordinates(2,3) = vrt(2,tri(3,i)) !y3
    ! Square calculation
    a = dsqrt((List_of_Triangles(i)%coordinates(1,1) - List_of_Triangles(i)%coordinates(1,2))**2 + &
      (List_of_Triangles(i)%coordinates(2,1) - List_of_Triangles(i)%coordinates(2,2))**2)
    b = dsqrt((List_of_Triangles(i)%coordinates(1,1) - List_of_Triangles(i)%coordinates(1,3))**2 + &
      (List_of_Triangles(i)%coordinates(2,1) - List_of_Triangles(i)%coordinates(2,3))**2)
    c = dsqrt((List_of_Triangles(i)%coordinates(1,2) - List_of_Triangles(i)%coordinates(1,3))**2 + &
      (List_of_Triangles(i)%coordinates(2,2) - List_of_Triangles(i)%coordinates(2,3))**2)
    p = (a+b+c)/2d0
    S = dsqrt(p*(p-a)*(p-b)*(p-c))
    List_of_Triangles(i)%size_of_triangle = S
    ! Neighbour elements for triangle
    if ((tri(1,i) < tri(2,i)) .and. (tri(2,i) < tri(3,i)) ) then
      List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(1,i))
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(2,i))
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(3,i))
    end if

    if ((tri(1,i) < tri(3,i)) .and. (tri(3,i) < tri(2,i)) ) then    
     List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(1,i))
     List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(3,i))
     List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(2,i))
    end if
    
    
    if ((tri(2,i) < tri(1,i)) .and. (tri(1,i) < tri(3,i)) ) then
      List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(2,i))
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(1,i))
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(3,i))
    end if
    
    if ((tri(2,i) < tri(3,i)) .and. (tri(3,i) < tri(1,i)) ) then
      List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(2,i))
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(3,i))
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(1,i))
    end if
    
    if ((tri(3,i) < tri(1,i)) .and. (tri(1,i) < tri(2,i)) ) then
      List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(3,i))
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(1,i))
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(2,i))
    end if
   
    if ((tri(3,i) < tri(2,i)) .and. (tri(2,i) < tri(1,i)) ) then 
      List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element => List_of_Elements(tri(3,i))
      List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element => List_of_Elements(tri(2,i))
      List_of_Triangles(i)%neighbour_elements_list(3)%pointing_element => List_of_Elements(tri(1,i))
    end if
    
    ! Neighbour triangles for triangle
    
    neighb_tr = 0
    
    if (IEE(1,i) /= 0) then
      neighb_tr = neighb_tr + 1
      List_of_Triangles(i)%neighbour_triangles_list(neighb_tr)%pointing_triangle => List_of_Triangles(IEE(1,i)) 
    end if
    
    if (IEE(2,i) /= 0) then
      neighb_tr = neighb_tr + 1
      List_of_Triangles(i)%neighbour_triangles_list(neighb_tr)%pointing_triangle => List_of_Triangles(IEE(2,i)) 
    end if
    
    if (IEE(3,i) /= 0) then
      neighb_tr = neighb_tr + 1
      List_of_Triangles(i)%neighbour_triangles_list(neighb_tr)%pointing_triangle => List_of_Triangles(IEE(3,i)) 
    end if
    
    List_of_Triangles(i)%number_of_neighbour_triangles = neighb_tr
    
    ! Neighbour edges for triangle
    do j = 1, List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%number_of_neighbour_edges
      if (List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%neighbour_edges_list(j)% &
       pointing_edge%neighbour_elements_list(2)%pointing_element%identificator ==  &
       List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%identificator) then
        List_of_Triangles(i)%neighbour_edges_list(1)%pointing_edge => List_of_Edges(List_of_Triangles(i)% &
        neighbour_elements_list(1)%pointing_element%neighbour_edges_list(j)%pointing_edge%identificator)
      end if
   
      if (List_of_Triangles(i)%neighbour_elements_list(1)%pointing_element%neighbour_edges_list(j)% &
      pointing_edge%neighbour_elements_list(2)%pointing_element%identificator == List_of_Triangles(i)% &
      neighbour_elements_list(3)%pointing_element%identificator) then
        List_of_Triangles(i)%neighbour_edges_list(2)%pointing_edge => List_of_Edges(List_of_Triangles(i)% &
        neighbour_elements_list(1)%pointing_element%neighbour_edges_list(j)%pointing_edge%identificator)     
      end if  
    end do
    
    do j = 1, List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%number_of_neighbour_edges
    
      if (List_of_Triangles(i)%neighbour_elements_list(2)%pointing_element%neighbour_edges_list(j)% &
      pointing_edge%neighbour_elements_list(2)%pointing_element%identificator == List_of_Triangles(i)% &
      neighbour_elements_list(3)%pointing_element%identificator) then
        List_of_Triangles(i)%neighbour_edges_list(3)%pointing_edge => List_of_Edges(List_of_Triangles(i)% &
        neighbour_elements_list(2)%pointing_element%neighbour_edges_list(j)%pointing_edge%identificator)    
      end if
    end do
  end do
  
  
  !! ============================================================================================ !!
  
  !! Making list of non Dirichlet elements
  
  k = 1
   
  do i = 1, number_of_elements
    if(List_of_Elements(i)%on_boundary .eqv. .false.) then
      List_of_non_Direchlet_Elements(k)%pointing_element => List_of_Elements(i)
      List_of_non_Direchlet_Elements(k)%pointing_element%nd_identificator = k
      k = k + 1
    end if
  end do
   
  number_of_non_Direchlet_elements = k - 1
  
end subroutine grid_initialization
  
  !! This routine initializes vector variables u(1), u(2)
  subroutine vectors_initialization(boundary_type)
  
    implicit none
    
    !!arguments:
    character(len=*), intent(in) :: boundary_type
    
    !! local variables
    integer                        :: i
    
    do i = 1, number_of_elements
      if (List_of_Elements(i)%on_boundary) then
        if (trim(boundary_type) == "cling") then
          List_of_Elements(i)%u(1) = 0d0
          List_of_Elements(i)%u(2) = 0d0
        end if
      else
        List_of_Elements(i)%u(1) = 0d0
        List_of_Elements(i)%u(2) = 0d0  
      end if
    end do

  end subroutine vectors_initialization
  
  !! This routine initializes scalar variables m, A, h
  subroutine scalars_initialization(square_size)
  
    implicit none
    
    !! arguments:
    real*8, intent(in)             :: square_size
    
    !! local variables
    integer                        :: i
       
    do i = 1, number_of_elements
      if ((List_of_Elements(i)%coordinates(2) < 25d-2*square_size) .and. &
          (List_of_Elements(i)%coordinates(2) > 5d-2*square_size) .and.  &
          (List_of_Elements(i)%coordinates(1) > 5d-2*square_size) .and. &
          (List_of_Elements(i)%coordinates(1) < 95d-2*square_size)) then
        List_of_Elements(i)%h = 2d0
        List_of_Elements(i)%A = 80d-2
      else
        List_of_Elements(i)%h = 10d-2
        List_of_Elements(i)%A = 10d-2
      end if
      List_of_Elements(i)%m = List_of_Elements(i)%A*List_of_Elements(i)%h*rho_ice
    end do

  end subroutine scalars_initialization
  
end module module_initialization
