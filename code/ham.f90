subroutine ham(kg)
use zeeman_mod
implicit none
integer ::l,r
double precision :: kg(2),dt1,dt2,delx,dely,delz

        dt1=dot_product(kg,n1)
        dt2=dot_product(kg,n2)
        
        delx=dot_product(kg,del_x)
        dely=dot_product(kg,del_y)
        delz=dot_product(kg,del_z)

!       print *, dt1,dt2

!bond terms

        ABcc=x(4)*exp(i*dt1)+x(17)*exp(i*dt2)+x(30)
        ABxx=x(1)*exp(i*dt1);   AByy=x(2)*exp(i*dt2);   ABzz=x(3)
        ABcx=-x(31)*exp(i*dt1); ABcy=-x(35)*exp(i*dt2); ABcz=-x(39)
        ABxc=-x(40)*exp(i*dt1); AByc=-x(44)*exp(i*dt2); ABzc=-x(48)

        MAB=(0.d0,0.d0)

        MAB(1,1)=ABcc
        MAB(2,2)=ABxx;  MAB(3,3)=AByy;  MAB(4,4)=ABzz
        MAB(1,2)=ABcx;  MAB(1,3)=ABcy;  MAB(1,4)=ABcz
        MAB(2,1)=ABxc;  MAB(3,1)=AByc;  MAB(4,1)=ABzc
        MAB=-i*K*MAB/2.d0
       
        !do l=1,4
         !       do r=1,4
          !              print *, l,r,MAB(l,r)
           !     end do
        !end do

        hb(1:4,1:4)=zero4;      hb(1:4,5:8)=MAB;        hb(5:8,1:4)=conjg(transpose(MAB));      h(5:8,5:8)=zero4
!on-site terms
        AAcx=-hx-4.d0*K*m(4);  AAcy=-hy-4.d0*K*m(5);  AAcz=-hz-4.d0*K*m(6)
        AAxc=hx+4.d0*K*m(4);  AAyc=hy+4.d0*K*m(5);  AAzc=hz+4.d0*K*m(6)
        AAxy=-hz;  AAxz=hy
        AAyx=hz;  AAyz=-hx
        AAzx=-hy;  AAzy=hx
        
        MAA=(0.d0,0.d0)
        MAA(1,2)=AAcx;  MAA(1,3)=AAcy;  MAA(1,4)=AAcz
        MAA(2,1)=AAxc;  MAA(3,1)=AAyc;  MAA(4,1)=AAzc
        MAA(2,3)=AAxy;  MAA(2,4)=AAxz
        MAA(3,2)=AAyx;  MAA(3,4)=AAyz
        MAA(4,2)=AAzx;  MAA(4,3)=AAzy
        MAA=i*MAA/4.d0

        BBcx=-hx-4.d0*K*m(1);  BBcy=-hy-4.d0*K*m(2);  BBcz=-hz-4.d0*K*m(3)
        BBxc=hx+4.d0*K*m(1);  BByc=hy+4.d0*K*m(2);  BBzc=hz+4.d0*K*m(3)
        BBxy=-hz;  BBxz=hy
        BByx=hz;  BByz=-hx
        BBzx=-hy;  BBzy=hx
               
        MBB=(0.d0,0.d0)
        MBB(1,2)=BBcx;  MBB(1,3)=BBcy;  MBB(1,4)=BBcz
        MBB(2,1)=BBxc;  MBB(3,1)=BByc;  MBB(4,1)=BBzc
        MBB(2,3)=BBxy;  MBB(2,4)=BBxz
        MBB(3,2)=BByx;  MBB(3,4)=BByz
        MBB(4,2)=BBzx;  MBB(4,3)=BBzy
        MBB=i*MBB/4.d0

        hos(1:4,1:4)=MAA;       hos(1:4,5:8)=zero4;     hos(5:8,1:4)=zero4;     hos(5:8,5:8)=MBB

!inter-sublattice terms from kappa
        kappaAB=zero4      
        kappaAB(1,1)=i*h_pert*(x(4)*x(17)*x(30))*(exp(i*dt1)*(x(49)+x(52))+exp(i*dt2)*(x(50)+x(53))+(x(51)+x(54)))/2.d0

        kappaBA=conjg(transpose(kappaAB))

!intra-sublattice terms from kappa
        kappaAA=zero4
        kappa2AA=i*h_pert*((exp(-i*delx)-exp(i*delx))*x(17)*x(30)+(exp(-i*dely)-exp(i*dely))*x(4)*x(30)+(exp(-i*delz)-exp(i*delz))*x(4)*x(17))/2.d0
        kappaiAA=-i*h_pert*(x(4)*x(17)*x(30))*((exp(-i*delx)-exp(i*delx))*x(1)+(exp(-i*dely)-exp(i*dely))*x(2)+(exp(-i*delz)-exp(i*delz))*x(3))/2.d0
        kappaAA(1,1)=kappa2AA+kappaiAA               
        kappaBB=-kappaAA
        
        hkappa(1:4,1:4)=kappaAA;        hkappa(1:4,5:8)=kappaAB;        hkappa(5:8,1:4)=kappaBA;        hkappa(5:8,5:8)=kappaBB

end subroutine

