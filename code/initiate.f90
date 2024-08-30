subroutine initiate()

use zeeman_mod
implicit none
integer :: p,q

!input K,hx,hy,hz,np,nq(grid size)
open(unit=12,file='para.in')
read(12,*)
read(12,*)K,h_max,h_div,np,nq,zeeman_info,kappa_info,n_iter,iter_max,stability,cutoff,conv_cut

hx=0.d0;        hy=0.d0;        hz=0.d0
h_pert=hx*hy*hz/(0.262**2.d0)

! switch zeeman, kappa term on/off based on zeeman_info, kappa_info being
! 1.d0/0.d0
hx=zeeman_info*hx;      hy=zeeman_info*hy;      hz=zeeman_info*hz
h_pert=kappa_info*h_pert

!print *, h_pert,hx,hy,hz
!close(12)

print '(a30,f9.6)', "Kitaev interaction strength=", K
!print '(a12,3f9.6)', "hx,hy,hz=",hx,hy,hz
!print *, "kappa interaction strengt=",h_pert



!K=1.d0
!hx=0.d0
!hy=0.d0
!hz=0.d0
!np=10
!nq=20
!print *, K,hx,hy,hz,np,nq

pi=4.d0*atan(1.d0)

allocate(U(np,nq,8,8))

JOBZ='V'
UPLO='U'
N=8
LDA=8
LWORK=20


!LATTICE SETUP
!lattice vectors (n1,n2)
n1=sqrt(3.d0)*(/sin(pi/6.d0),cos(pi/6.d0)/)
n2=sqrt(3.d0)*(/sin(-pi/6.d0),cos(-pi/6.d0)/)

!nearest neighbor distances
del_x=n2
del_y=-n1
del_z=n1-n2

!print *, n1,n2
!print *, "delx=",del_x
!print *, "dely=",del_y
!print *, "delz=",del_z

!reciprocal lattice vectors
f1=(/2.d0*pi/sqrt(3.d0),2.d0*pi/3.d0/)
f2=(/-2.d0*pi/sqrt(3.d0),2.d0*pi/3.d0/)


kvec=transpose(reshape((/f1,f2/),(/2,2/)))
        k1_p=(/2.d0*pi/(3.d0*sqrt(3.d0)),2.d0*pi/3.d0/)
        k2_p=(/-2.d0*pi/(3.d0*sqrt(3.d0)),2.d0*pi/3.d0/)
        m_p=(/pi/(sqrt(3.d0)),pi/3.d0/)
        g_p=(/0.d0,0.d0/)
        print *,'K-point=', k1_p
        print *,'M-point=', m_p
        print *,'Gamma-point=',g_p 
do p=1,4
        do q=1,4
                zero4(p,q)=0.d0
        end do
end do

om_v1(1,1:4)=(/(0.5d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,-0.5d0)/)
om_v1(2,1:4)=(/(0.d0,0.d0),(0.d0,-0.5d0),(0.5d0,0.d0),(0.d0,0.d0)/)
om_v1(3,1:4)=(/(0.5d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.5d0)/)
om_v1(4,1:4)=(/(0.d0,0.d0),(0.d0,0.5d0),(0.5d0,0.d0),(0.d0,0.d0)/)

v1(1:4,1:4)=om_v1;      v1(1:4,5:8)=zero4
v1(5:8,1:4)=zero4;      v1(5:8,5:8)=om_v1

open(21,file='kappa0.in')
open(22,file='kappax.in')
open(23,file='kappay.in')
open(24,file='kappaz.in')

read(21,*)
read(21,*) x(1),x(2),x(3)
read(22,*)
read(22,*) x(4),x(8),x(12)
read(23,*)
read(23,*) x(13),x(17),x(21)
read(24,*)
read(24,*) x(22),x(26),x(30)

x(31:48)=(0.d0,0.d0)

open(25,file='ma.in')
open(26,file='mb.in')

read(25,*)
read(25,*) m(1),m(2),m(3)
read(26,*)
read(26,*) m(4),m(5),m(6)

open(31,file='lm.in')
read(31,*)
read(31,*) lm(1), lm(2), lm(3),lm(4),lm(5),lm(6)

open(101,file='chiA.in')
open(102,file='chiB.in')
read(101,*)
read(101,*) x(49),x(50),x(51)
read(102,*)
read(102,*) x(52),x(53),x(54)

end subroutine initiate
