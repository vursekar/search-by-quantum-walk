!------------------------------------------------------
! MODULE with all constants
!------------------------------------------------------
! nn     :: the width/length of the 2d lattice
! nsites :: nn/2 or half the width/length of 2d lattice
! tsteps :: number of discrete time steps in each walk
!------------------------------------------------------
module constants
  integer, parameter :: nn=10, nsites=nn/2, tsteps=1000
  double precision, parameter :: pi=4.d0*datan(1.d0)
end module constants

!------------------------------------------------------
! MAIN PROGRAM
!------------------------------------------------------
program quantumwalk2d
  use constants
  implicit none
!-----------------------------------------------------------------------
! Variables
! lattice  :: 4x(2n+1)^2 array with a ampltitudes of all lattice states
! coin     :: 4x4 coin operator
! grover   :: if true, walk uses Grover coin; else uses Hadamard coin
! inverted :: if true, sends coin --> -coin
! single   :: if true, runs a single walk; final state in walk_2d.txt
! marked   :: if true, runs search algorithm; else runs regular walk
! phase_error :: if true, introduces unitary noise to each coin operation
! delta      :: sets the scale of unitary noise
! j0,k0    :: location of marked vertex
! nwalks   :: for multi-walk case, determines how many walks to run
! nrep     :: determines how many times a walk should be repeated
!-----------------------------------------------------------------------
  double complex                :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double complex                :: coin(1:4,1:4)
  integer                       :: t, i, j, iwalks, irep
  logical, parameter            :: grover=.true., inverted=.false.
  logical, parameter            :: single=.true., marked=.true., phase_error=.true.
  double precision, parameter   :: p=0.d0, eps=0.d0
  integer, parameter            :: j0=5, k0=5, nwalks=1, nrep=1
  double precision              :: max_amp_tot(2), max_amp(2), delta

! Set initial noise to 0
  delta = 0.d0

! Set coin operator
  call set_coin(coin,grover,inverted)

! Code for just a single run

  if(single) then
    call run_walk(single,lattice,coin,p,eps,marked,delta,phase_error,j0,k0,max_amp)

! Write final probability amplitudes into walk_2d.txt
    open(unit=21,file='walk_2d.txt')
    do i=-nsites,nsites
      do j=-nsites,nsites
        write(21,*) i, j, abs(lattice(1,i,j))**2+abs(lattice(2,i,j))**2+abs(lattice(3,i,j))**2+abs(lattice(4,i,j))**2
      enddo
      write(21,*)
    enddo

! Code for multiple runs
  else

! Records the maximum amplitude reached for different delta
    open(unit=25,file='maxamp_delta.txt')

! nwalks walks are performed
    do iwalks=1, nwalks
      print*, iwalks, delta
      max_amp_tot = 0.d0

! For each delta, take the average of nrep runs
      do irep=1,nrep
        call run_walk(single,lattice,coin,p,eps,marked,delta,phase_error,j0,k0,max_amp)
        max_amp_tot = max_amp_tot + max_amp
      enddo

      max_amp_tot = max_amp_tot/dble(nrep)
      write(25,*) delta, max_amp_tot

! Increment delta
      delta = delta + 1.d-3
    enddo

  endif

end program quantumwalk2d


!---------------------------------------------
! SUBROUTINES
!---------------------------------------------

!---------------------------------------------------------
! run_walk :: Runs a quantum walk on an nsites+1xnsites+1
! lattice for t=tsteps with coin=coin and unitary noise
! with scale delta
!---------------------------------------------------------
subroutine run_walk(single,lattice,coin,p,eps,marked,delta,phase_error,j0,k0,max_amp)
  use constants
  implicit none
  double complex                :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double complex                :: coin(1:4,1:4)
  integer                       :: t, i, j
  logical                       :: single,marked, phase_error
  double precision              :: p, eps, delta, max_amp(2), amp
  integer                       :: j0, k0

! If single run, stores probability of marked site and origin at each time step
  if(single) then
    open(unit=22,file='prob_marked.txt')
    open(unit=23,file='prob_origin.txt')
    open(unit=24,file='maxamp.txt')
  endif

! Initialize lattice in uniform superposition state

!    lattice = 0.d0
!    lattice(1,0,0) =  cmplx( 1.d0,0.d0)
!    lattice(2,0,0) =  cmplx( 1.d0,0.d0)
!    lattice(3,0,0) =  cmplx( 1.d0,0.d0)
!    lattice(4,0,0) =  cmplx(-1.d0,0.d0)

    lattice = cmplx(1.d0,0.d0)
    call normalize(lattice)

    max_amp(1) = abs(lattice(1,j0,k0))**2+abs(lattice(2,j0,k0))**2+abs(lattice(3,j0,k0))**2+abs(lattice(4,j0,k0))**2
    max_amp(2) = 0

! Run walk for tsteps
    do t=1, tsteps

      if(single.and.mod(t,100).eq.0) print*, t

! Choose which shift operator to use: default is flip-flop
!      call shift(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
!      call shift_shenvim(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
      call shift_shenviff(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
      call normalize(lattice)

      amp = abs(lattice(1,j0,k0))**2+abs(lattice(2,j0,k0))**2+abs(lattice(3,j0,k0))**2+abs(lattice(4,j0,k0))**2

      if(amp.gt.max_amp(1)) then
        max_amp(1) = amp
        max_amp(2) = t
        write(24,*) t, amp
      endif

      if(single) then
        write(22,*) t, amp
        write(23,*) t, abs(lattice(1,0,0))**2+abs(lattice(2,0,0))**2+abs(lattice(3,0,0))**2+abs(lattice(4,0,0))**2
      endif

    enddo

end subroutine run_walk

!---------------------------------------------------------
! set_coin :: Defines coin operator
!---------------------------------------------------------
subroutine set_coin(coin,grover,inverted)
  use constants
  implicit none
  double complex :: coin(4,4)
  logical        :: grover, inverted

  coin = cmplx(1.d0/2.d0,0.d0)
  if(inverted) coin = -coin

  if(grover) then
    coin(1,1) = -coin(1,1)
    coin(2,2) = -coin(2,2)
    coin(3,3) = -coin(3,3)
    coin(4,4) = -coin(4,4)
  else
    coin(2,2) = -coin(2,2)
    coin(2,4) = -coin(2,4)
    coin(4,2) = -coin(4,2)
    coin(3,3) = -coin(3,3)
    coin(3,4) = -coin(3,4)
    coin(4,3) = -coin(4,3)
  endif

end subroutine set_coin


!---------------------------------------------------------
! normalize :: Normalizes wavefunction
!---------------------------------------------------------
subroutine normalize(lattice)
  use constants
  implicit none
  double complex   :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double precision :: a
  integer          :: i, j

  a = 0.d0
  do i=-nsites,nsites
    do j=-nsites,nsites
      a = a + abs(lattice(1,i,j))**2 + abs(lattice(2,i,j))**2 + abs(lattice(3,i,j))**2 + abs(lattice(4,i,j))**2
    enddo
  enddo
  lattice = lattice/sqrt(a)

end subroutine normalize


!---------------------------------------------------------
! operate_coin :: Operates coin on single lattice site
!---------------------------------------------------------
subroutine operate_coin(m,coin,eps,delta,phase_error)
  use constants
  use random
  implicit none
  integer, parameter :: n=4
  double complex     :: m(1:4), coin(4,4), id(n,n)
  double complex     :: noise(n,n), phase(n,n), work(n,n)
  double precision   :: eps, r, rp(3*n*(n-1)/2+1)
  integer            :: i, j, icount
  double precision   :: delta, phi, psi, xi
  logical            :: phase_error

! Generates unitary perturbation to coin operator
! according to algorithm by Kus

  if(phase_error) then
    id = 0.d0
    do i=1,n
      id(i,i) = 1.d0
    enddo
    phase = id

    call random_number(rp)
    rp = rp*delta
    rp(6) = 0.d0
    rp(12) = 0.d0
    rp(15) = 0.d0

    icount=0
    do i=2, n
      do j=i-1,1,-1
        phi = abs(asin(rp(icount+1)))
        psi = pi*2.d0*rp(icount+2)
        xi  = pi*2.d0*rp(icount+3)

        work = id
        work(j,j) = +cmplx(cos(psi),+sin(psi))*cos(phi)
        work(j,i) = +cmplx(cos(xi), +sin(xi))*sin(phi)
        work(i,j) = -cmplx(cos(xi), -sin(xi))*sin(phi)
        work(i,i) = +cmplx(cos(psi),-sin(psi))*cos(phi)

        phase = matmul(phase,work)
        icount = icount+3
      enddo
    enddo
    phase = cmplx(cos(rp(3*n*(n-1)/2+1)),+sin(rp(3*n*(n-1)/2+1)))*phase
  endif

! Generates non-unitary noise matrix added to coin
  noise = 0.d0
  if(eps.gt.0.d0) then
    do i=1,4
      do j=1,4
        call random_number(r)
  !      r = random_normal()
        noise(i,j) = r*eps
      enddo
    enddo
  endif

! Operate with coin+noise
  m = matmul(coin,m)
  if(phase_error) m = matmul(phase,m)

end subroutine operate_coin

!---------------------------------------------------------
! shift :: shifts lattice according to 00 --> down, left;
! 01 --> down, right ; 10 --> up, left; 11 --> up, right.
!---------------------------------------------------------

subroutine shift(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
  use constants
  implicit none
  double complex   :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: work(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: coin(4,4)
  integer          :: j, k, j0, k0, jp, jm, kp, km
  double precision :: r, p, eps, delta
  logical          :: marked, phase_error

  work = lattice

! Generate copy of lattice post-coin operation
  do j=-nsites,nsites
    do k=-nsites,nsites
      call operate_coin(work(:,j,k),coin,eps,delta,phase_error)
    enddo
  enddo

  if(marked) work(:,j0,k0) = -lattice(:,j0,k0)

! Update lattice according to post-coin copy
  do j=-nsites ,nsites
    do k=-nsites ,nsites

      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if(j.eq.-nsites) jm=+nsites
      if(k.eq.-nsites) km=+nsites
      if(j.eq.+nsites) jp=-nsites
      if(k.eq.+nsites) kp=-nsites

      lattice(1,j,k) = work(1,jm,km)
      lattice(2,j,k) = work(2,jp,km)
      lattice(3,j,k) = work(3,jm,kp)
      lattice(4,j,k) = work(4,jp,kp)

    enddo
  enddo

end subroutine shift


!---------------------------------------------------------
! shift :: shifts lattice according to 00 --> down;
! 01 --> up ; 10 --> left; 11 --> right. Spin is flipped
! at each step.
!---------------------------------------------------------
subroutine shift_shenviff(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
  use constants
  implicit none
  double complex   :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: work(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: coin(4,4), f
  integer          :: j, k, j0, k0
  integer          :: jp, jm, kp, km
  double precision :: r, p, eps, delta
  logical          :: marked, phase_error

  work = lattice

  f = cmplx(-1.d0,0.d0)
  f = f/abs(f)

! Generate copy of lattice post-coin operation
  do j=-nsites,nsites
    do k=-nsites,nsites
      if(marked.and.j.eq.j0.and.k.eq.k0) then
        call operate_coin(work(:,j,k),f*coin,eps,delta,phase_error)
      else
        call operate_coin(work(1:4,j,k),coin,eps,delta,phase_error)
      endif
    enddo
  enddo

! Update lattice according to post-coin copy
  do j=-nsites ,nsites
    do k=-nsites ,nsites

      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if(j.eq.-nsites) jm=+nsites
      if(k.eq.-nsites) km=+nsites
      if(j.eq.+nsites) jp=-nsites
      if(k.eq.+nsites) kp=-nsites

      lattice(1,j,k) = work(2,jp,k)
      lattice(2,j,k) = work(1,jm,k)
      lattice(3,j,k) = work(4,j,kp)
      lattice(4,j,k) = work(3,j,km)

    enddo
  enddo


end subroutine shift_shenviff

!---------------------------------------------------------
! shift :: shifts lattice according to 00 --> down;
! 01 --> up ; 10 --> left; 11 --> right. Spin remains the
! same after each step.
!---------------------------------------------------------

subroutine shift_shenvim(lattice,coin,p,eps,marked,delta,phase_error,j0,k0)
  use constants
  implicit none
  double complex   :: lattice(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: work(1:4,-nsites:nsites,-nsites:nsites)
  double complex   :: coin(4,4), f
  integer          :: j, k, j0, k0, jp, jm, kp, km
  double precision :: r, p, eps, delta
  logical          :: marked, phase_error

  work = lattice
  f = cmplx(-1.d0,0.d0)
  f = f/abs(f)

! Generate copy of lattice post-coin operation
  do j=-nsites,nsites
    do k=-nsites,nsites
      if(marked.and.j.eq.j0.and.k.eq.k0) then
        call operate_coin(work(:,j,k),f*coin,eps,delta,phase_error)
      else
        call operate_coin(work(:,j,k),coin,eps,delta,phase_error)
      endif
    enddo
  enddo


! Update lattice according to post-coin copy
  do j=-nsites ,nsites
    do k=-nsites ,nsites

      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if(j.eq.-nsites) jm=+nsites
      if(k.eq.-nsites) km=+nsites
      if(j.eq.+nsites) jp=-nsites
      if(k.eq.+nsites) kp=-nsites

      lattice(1,j,k) = work(1,jp,k)
      lattice(2,j,k) = work(2,jm,k)
      lattice(3,j,k) = work(3,j,kp)
      lattice(4,j,k) = work(4,j,km)

    enddo
  enddo

end subroutine shift_shenvim
