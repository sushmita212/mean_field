subroutine sc_eq()
use zeeman_mod
implicit none
external :: zheev
integer :: ip, iq, il,im, band
double complex :: fvec(6)
double precision::k_ind(2,1),kg1(2,1),kg(2),dt1,dt2,err1_real(6),err1_imag(6),delx,dely,delz

        !print *, "the lagrange multipliers are"
        !print *, lm
        MCA=(0.d0,0.d0)
        MCA(1,2)=-lm(1);     MCA(1,3)=-lm(2);     MCA(1,4)=-lm(3)
        MCA(2,1)=lm(1);      MCA(3,1)=lm(2);      MCA(4,1)=lm(3)
        MCA(2,3)=lm(3);      MCA(2,4)=-lm(2)   
        MCA(3,2)=-lm(3);     MCA(3,4)=lm(1)
        MCA(4,2)=lm(2);      MCA(4,3)=-lm(1)
        MCA=MCA*i/4.d0

        MCB=(0.d0,0.d0)
        MCB(1,2)=-lm(4);     MCB(1,3)=-lm(5);     MCB(1,4)=-lm(6)
        MCB(2,1)=lm(4);      MCB(3,1)=lm(5);      MCB(4,1)=lm(6)
        MCB(2,3)=lm(6);      MCB(2,4)=-lm(5)   
        MCB(3,2)=-lm(6);     MCB(3,4)=lm(4)
        MCB(4,2)=lm(5);      MCB(4,3)=-lm(4)
        MCB=MCB*i/4.d0

        hcons(1:4,1:4)=MCA;     hcons(1:4,5:8)=zero4;   hcons(5:8,1:4)=zero4;   hcons(5:8,5:8)=MCB
        
        hmaj=hb+hos+hcons+hkappa
        h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        
        do ip=1,np
         do iq=1,nq
                
                k_ind=reshape((/ip*1.d0/np,iq*1.d0/nq/),(/2,1/))
                kg1=matmul(transpose(kvec),k_ind)
                kg=kg1(:,1)
                
                call ham(kg)
                hmaj=hb+hos+hcons+hkappa
                
                h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
                call zheev(JOBZ, UPLO, n, h, LDA, W, WORK, LWORK, RWORK,INFO)
                U(ip,iq,:,:)=h
        
         end do
        end do
        !print*, conjg(U(1,2,3,4))
        fvec=0.d0
        x=(0.d0,0.d0)
        m=(0.d0,0.d0)
        do ip=1,np
                do iq=1,nq
                        k_ind=reshape((/ip*1.d0/np,iq*1.d0/nq/),(/2,1/))
                        kg1=matmul(transpose(kvec),k_ind)
                        kg=kg1(:,1)
                        dt1=dot_product(kg,n1)
                        dt2=dot_product(kg,n2)
                        delx=dot_product(kg,del_x)
                        dely=dot_product(kg,del_y)
                        delz=dot_product(kg,del_z)

                        do band=5,8
                        !kappa0=ic_ic_j
                        
                        x(1)=x(1)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(1)=x(1)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                                  
                        x(2)=x(2)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(2)=x(2)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
 
                        x(3)=x(3)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))
                        x(3)=x(3)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))
                        
                        !kappax
                        x(4)=x(4)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(4)=x(4)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        
                        x(5)=x(5)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(5)=x(5)+(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        
                        x(6)=x(6)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(6)=x(6)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)

                        x(7)=x(7)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(7)=x(7)-(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
 
                        x(8)=x(8)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(8)=x(8)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        
                        x(9)=x(9)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(9)=x(9)-(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        
                        x(10)=x(10)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(10)=x(10)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)

                        x(11)=x(11)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(11)=x(11)+(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)

                        x(12)=x(12)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(12)=x(12)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        
                        !kappay
                        x(13)=x(13)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(13)=x(13)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        
                        x(14)=x(14)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(14)=x(14)+(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        
                        x(15)=x(15)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(15)=x(15)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)

                        x(16)=x(16)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(16)=x(16)-(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
 
                        x(17)=x(17)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(17)=x(17)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        
                        x(18)=x(18)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(18)=x(18)-(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        
                        x(19)=x(19)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(19)=x(19)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)

                        x(20)=x(20)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(20)=x(20)+(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)

                        x(21)=x(21)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(21)=x(21)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
 
                        !kappaz
                        x(22)=x(22)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))
                        x(22)=x(22)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))
                        
                        x(23)=x(23)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))
                        x(23)=x(23)+(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))
                        
                        x(24)=x(24)-i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))
                        x(24)=x(24)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))

                        x(25)=x(25)-(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))-U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))
                        x(25)=x(25)-(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))-U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))
 
                        x(26)=x(26)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,8,band))+U(ip,iq,2,band)*conjg(U(ip,iq,6,band)))
                        x(26)=x(26)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,8,band))+U(ip,iq,4,band)*conjg(U(ip,iq,6,band)))
                        
                        x(27)=x(27)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))-U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))
                        x(27)=x(27)-(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))-U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))
                        
                        x(28)=x(28)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))
                        x(28)=x(28)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))

                        x(29)=x(29)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))
                        x(29)=x(29)+(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))

                        x(30)=x(30)-i*(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))
                        x(30)=x(30)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))
                        
                        !eta
                        x(31)=x(31)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(31)=x(31)+(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        
                        x(32)=x(32)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(32)=x(32)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        
                        x(33)=x(33)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(33)=x(33)+(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)

                        x(34)=x(34)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(34)=x(34)+(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
 
                        x(35)=x(35)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(35)=x(35)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        
                        x(36)=x(36)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(36)=x(36)+(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        
                        x(37)=x(37)-(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))
                        x(37)=x(37)+(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))

                        x(38)=x(38)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,7,band))+U(ip,iq,2,band)*conjg(U(ip,iq,5,band)))
                        x(38)=x(38)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,7,band))+U(ip,iq,4,band)*conjg(U(ip,iq,5,band)))

                        x(39)=x(39)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))+U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))
                        x(39)=x(39)+(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))+U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))

                        !eta(tilde)
                        x(40)=x(40)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(40)=x(40)-(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        
                        x(41)=x(41)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        x(41)=x(41)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt1)
                        
                        x(42)=x(42)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)
                        x(42)=x(42)-(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt1)

                        x(43)=x(43)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(43)=x(43)-(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        
                        x(44)=x(44)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        x(44)=x(44)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))*exp(-i*dt2)
                        
                        x(45)=x(45)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)
                        x(45)=x(45)-(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))*exp(-i*dt2)

                        x(46)=x(46)-(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))-U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))
                        x(46)=x(46)-(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))-U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))
                        
                        x(47)=x(47)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,8,band))+U(ip,iq,1,band)*conjg(U(ip,iq,6,band)))
                        x(47)=x(47)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,8,band))+U(ip,iq,3,band)*conjg(U(ip,iq,6,band)))
                        
                        x(48)=x(48)-(U(ip,iq,1,band)*conjg(U(ip,iq,7,band))-U(ip,iq,1,band)*conjg(U(ip,iq,5,band)))
                        x(48)=x(48)-(U(ip,iq,3,band)*conjg(U(ip,iq,7,band))-U(ip,iq,3,band)*conjg(U(ip,iq,5,band)))

                        !magnetization
                        m(1)=m(1)-(U(ip,iq,2,band)*conjg(U(ip,iq,3,band))+U(ip,iq,2,band)*conjg(U(ip,iq,1,band)))
                        m(1)=m(1)+(U(ip,iq,4,band)*conjg(U(ip,iq,3,band))+U(ip,iq,4,band)*conjg(U(ip,iq,1,band)))

                        m(2)=m(2)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,3,band))+U(ip,iq,2,band)*conjg(U(ip,iq,1,band)))
                        m(2)=m(2)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,3,band))+U(ip,iq,4,band)*conjg(U(ip,iq,1,band)))

                        m(3)=m(3)-(U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band)))
                        m(3)=m(3)+(U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band)))

                        m(4)=m(4)-(U(ip,iq,6,band)*conjg(U(ip,iq,7,band))+U(ip,iq,6,band)*conjg(U(ip,iq,5,band)))
                        m(4)=m(4)+(U(ip,iq,8,band)*conjg(U(ip,iq,7,band))+U(ip,iq,8,band)*conjg(U(ip,iq,5,band)))

                        m(5)=m(5)+i*(U(ip,iq,6,band)*conjg(U(ip,iq,7,band))+U(ip,iq,6,band)*conjg(U(ip,iq,5,band)))
                        m(5)=m(5)+i*(U(ip,iq,8,band)*conjg(U(ip,iq,7,band))+U(ip,iq,8,band)*conjg(U(ip,iq,5,band)))

                        m(6)=m(6)-(U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band)))
                        m(6)=m(6)+(U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band)))

                     
                        !constraints
                        !on A sublattice
                        !(ib^xc+ib^yb^z)==0 
                        fvec(1)=fvec(1)-2.d0*(U(ip,iq,2,band)*conjg(U(ip,iq,3,band))-U(ip,iq,4,band)*conjg(U(ip,iq,1,band)))
                        
                        !(ib^yc+ib^zb^x)==0 
                        fvec(2)=fvec(2)+i*(U(ip,iq,2,band)*conjg(U(ip,iq,3,band))+U(ip,iq,2,band)*conjg(U(ip,iq,1,band)))
                        fvec(2)=fvec(2)+i*(U(ip,iq,4,band)*conjg(U(ip,iq,3,band))+U(ip,iq,4,band)*conjg(U(ip,iq,1,band)))
                        fvec(2)=fvec(2)+i*(-U(ip,iq,1,band)*conjg(U(ip,iq,4,band))+U(ip,iq,1,band)*conjg(U(ip,iq,2,band)))
                        fvec(2)=fvec(2)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,4,band))-U(ip,iq,3,band)*conjg(U(ip,iq,2,band)))
                        
                        !(ib^zc+ib^xb^y)==0
                        fvec(3)=fvec(3)-(U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band)))
                        fvec(3)=fvec(3)+(U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band)))
                        fvec(3)=fvec(3)-(U(ip,iq,2,band)*conjg(U(ip,iq,4,band))+U(ip,iq,2,band)*conjg(U(ip,iq,2,band)))
                        fvec(3)=fvec(3)+(U(ip,iq,4,band)*conjg(U(ip,iq,4,band))+U(ip,iq,4,band)*conjg(U(ip,iq,2,band)))
                         
                        !on B sublattice
                        !(ib^xc+ib^yb^z)==0
                        fvec(4)=fvec(4)-2.d0*(U(ip,iq,6,band)*conjg(U(ip,iq,7,band))-U(ip,iq,8,band)*conjg(U(ip,iq,5,band)))
                        
                        !(ib^yc+ib^zb^x)==0
                        fvec(5)=fvec(5)+i*(U(ip,iq,6,band)*conjg(U(ip,iq,7,band))+U(ip,iq,6,band)*conjg(U(ip,iq,5,band)))
                        fvec(5)=fvec(5)+i*(U(ip,iq,8,band)*conjg(U(ip,iq,7,band))+U(ip,iq,8,band)*conjg(U(ip,iq,5,band)))
                        fvec(5)=fvec(5)+i*(-U(ip,iq,5,band)*conjg(U(ip,iq,8,band))+U(ip,iq,5,band)*conjg(U(ip,iq,6,band)))
                        fvec(5)=fvec(5)+i*(U(ip,iq,7,band)*conjg(U(ip,iq,8,band))-U(ip,iq,7,band)*conjg(U(ip,iq,6,band)))
                       
                        !(ib^zc+ib^xb^y)==0
                        fvec(6)=fvec(6)-(U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band)))
                        fvec(6)=fvec(6)+(U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band)))
                        fvec(6)=fvec(6)-(U(ip,iq,6,band)*conjg(U(ip,iq,8,band))+U(ip,iq,6,band)*conjg(U(ip,iq,6,band)))
                        fvec(6)=fvec(6)+(U(ip,iq,8,band)*conjg(U(ip,iq,8,band))+U(ip,iq,8,band)*conjg(U(ip,iq,6,band)))
 
                        !NNN parameters
                        x(49)=x(49)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band)))*exp(-i*delx)
                        x(49)=x(49)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band)))*exp(-i*delx)
 
                        x(50)=x(50)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band)))*exp(-i*dely)
                        x(50)=x(50)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band)))*exp(-i*dely)
 
                        x(51)=x(51)+i*(U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band)))*exp(-i*delz)
                        x(51)=x(51)+i*(U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band)))*exp(-i*delz)

                        x(52)=x(52)+i*(U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band)))*exp(i*delx)
                        x(52)=x(52)+i*(U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band)))*exp(i*delx)
 
                        x(53)=x(53)+i*(U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band)))*exp(i*dely)
                        x(53)=x(53)+i*(U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band)))*exp(i*dely)
 
                        x(54)=x(54)+i*(U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band)))*exp(i*delz)
                        x(54)=x(54)+i*(U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band)))*exp(i*delz)


                        
                        end do
                end do
        end do
        x=x/(np*nq)
        m=m/(2.d0*np*nq)
        fvec(1)=fvec(1)/(np*nq)
        fvec(2)=fvec(2)/(np*nq)
        fvec(3)=fvec(3)/(np*nq)
        fvec(4)=fvec(4)/(np*nq)
        fvec(5)=fvec(5)/(np*nq)
        fvec(6)=fvec(6)/(np*nq)
        
        !print *, fvec(1),fvec(2),fvec(3),fvec(4),fvec(5),fvec(6)

        !print *, kg

        !checking the constraints
        do ip=1,6
                err_real(ip)=real(fvec(ip))
                err_imag(ip)=real(fvec(ip))
                err1_real(ip)=abs(real(fvec(ip)))
                err1_imag(ip)=abs(aimag(fvec(ip)))
        end do 
        
        !print *, "the constraint values are"
        !print *, "real part"               
        !print *, err_real
        !print *, "imag part"
        !print *, err_imag

        constr_real=maxval(err1_real)
        constr_imag=maxval(err1_imag)
        !print *, "the smallest constraints arre: real,imag"
        !print *, constr_real, constr_imag



end subroutine sc_eq
