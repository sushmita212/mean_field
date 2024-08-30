subroutine hybrd_solve
use zeeman_mod
implicit none

integer :: Nh,INFOh,lwa, ij,il
double precision,dimension(:),allocatable :: FVEC,lm1,wa
double precision :: TOL,fnorm,enorm
external :: FCN


Nh=6
lwa=Nh*(3*Nh+13)/2+5
allocate(lm1(Nh),FVEC(Nh))
allocate(wa(lwa))

TOL=1d-13
!print *, "error tolerance=",tol
call hybrd1(fcn,nh,lm1,fvec,tol,infoh,wa,lwa)
fnorm=enorm(nh,fvec)

!print *,w
!print *,"the lagrange multipliers are"
!print *, lm1
!print *,"the constraints are"
!print *, fvec

lm=lm1
!print *,"the lagrange multipliers are"
!print *, lm
!do ij=1,8
 !       do il=1,8
  !              print *, ij,il,hcons(ij,il)
   !     end do
!end do       
end subroutine hybrd_solve


