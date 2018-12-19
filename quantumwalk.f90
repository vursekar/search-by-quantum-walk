!---------------------------------------------
! Module
!---------------------------------------------
! nsites :: number of sites in lattice
! tsteps :: number of discrete time steps
! p      :: probability that link is broken
!---------------------------------------------
module constants
  integer, parameter :: nsites=10, tsteps=500
  double precision   :: p=0.d0
end module constants

!---------------------------------------------
! Program
!---------------------------------------------
program quantumwalk
  use constants
  implicit none
  double complex     :: lattice(0:1,-nsites:nsites)
  integer            :: t, i

  open(unit=21,file='walk.txt')

  lattice = 0.d0
  lattice(0,0) =  cmplx(1.d0,0.d0)
  lattice(1,0) =  cmplx(0.d0,0.d0)
  call normalize(lattice)

  do t=1, tsteps
    call shift(lattice)
    call normalize(lattice)
  enddo

  do i=-nsites,nsites
    write(21,*) i, abs(lattice(0,i))**2+abs(lattice(1,i))**2
  enddo


end program quantumwalk

subroutine operate_coin(m)
  use constants
  implicit none
  double complex :: m(0:1), coin(2,2)

! Define Hadamard coin
!  coin = 1.d0
  coin = 1.d0/sqrt(2.d0)
  coin(2,2) = -coin(2,2)

  m = matmul(coin,m)
end subroutine operate_coin

subroutine shift(lattice)
  use constants
  implicit none
  double complex   :: lattice(0:1,-nsites:nsites), work(0:1,-nsites:nsites)
  integer          :: j
  double precision :: r
  logical          :: periodic, grover, measure

  periodic = .false.
  grover   = .true.
  measure = .false.

  work = lattice

  do j=-nsites,nsites
    call operate_coin(work(:,j))
  enddo

  do j=-nsites+1 ,nsites-1
    call random_number(r)
    if(r.lt.p) cycle
    if(grover) then
      lattice(1,j) = work(0,j+1)
      lattice(0,j) = work(1,j-1)
    else
      lattice(0,j) = work(0,j+1)
      lattice(1,j) = work(1,j-1)
    endif

    if(measure) then
      call random_number(r)
      if(r.lt.abs(lattice(1,j))**2) then
        lattice(1,j) = 1.d0
        lattice(0,j) = 0.d0
      else
        lattice(1,j) = 0.d0
        lattice(0,j) = 1.d0
      endif
    endif
  enddo

  if(periodic) then
    lattice(0,-nsites) = work(0,-nsites+1)
    lattice(1,-nsites) = work(1,nsites)
    lattice(0,nsites) = work(0,-nsites)
    lattice(1,nsites) = work(1,nsites-1)
  endif

end subroutine shift

subroutine normalize(lattice)
  use constants
  implicit none
  double complex :: lattice(0:1,-nsites:nsites)
  double precision :: A
  integer          :: i

  A = 0.d0
  do i=-nsites,nsites
    A = A + abs(lattice(0,i))**2 + abs(lattice(1,i))**2
  enddo

  lattice = lattice/sqrt(A)

end subroutine normalize
