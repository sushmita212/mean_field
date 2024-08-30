!this subroutine contains the mf constraint equations to be solved simultaeously
subroutine fcn(nh,lm1,fvec,iflag)
use zeeman_mod
external :: zheev
double precision :: fvec(nh),lm1(6)
integer :: iflag,nh,ip,iq,band
double precision::k_ind(2,1),kg1(2,1),kg(2)
       nh=6
       ! print *, shape(fvec)
       !contraints
        
        MCA=(0.d0,0.d0)
        MCA(1,2)=-lm1(1);     MCA(1,3)=-lm1(2);     MCA(1,4)=-lm1(3)
        MCA(2,1)=lm1(1);      MCA(3,1)=lm1(2);      MCA(4,1)=lm1(3)
        MCA(2,3)=lm1(3);      MCA(2,4)=-lm1(2)   
        MCA(3,2)=-lm1(3);     MCA(3,4)=lm1(1)
        MCA(4,2)=lm1(2);      MCA(4,3)=-lm1(1)
        MCA=MCA*i/4.d0

        MCB=(0.d0,0.d0)
        MCB(1,2)=-lm1(4);     MCB(1,3)=-lm1(5);     MCB(1,4)=-lm1(6)
        MCB(2,1)=lm1(4);      MCB(3,1)=lm1(5);      MCB(4,1)=lm1(6)
        MCB(2,3)=lm1(6);      MCB(2,4)=-lm1(5)   
        MCB(3,2)=-lm1(6);     MCB(3,4)=lm1(4)
        MCB(4,2)=lm1(5);      MCB(4,3)=-lm1(4)
        MCB=MCB*i/4.d0

        hcons(1:4,1:4)=MCA;     hcons(1:4,5:8)=zero4;   hcons(5:8,1:4)=zero4;   hcons(5:8,5:8)=MCB
        
       ! hmaj=hb+hos+hcons+hkappa
       ! h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        
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
        do ip=1,np
                do iq=1,nq
                        do band=5,8
                        fvec(1)=fvec(1)-2.d0*real((U(ip,iq,2,band)*conjg(U(ip,iq,3,band))-U(ip,iq,4,band)*conjg(U(ip,iq,1,band))))
                        
                        fvec(2)=fvec(2)+real(i*(U(ip,iq,2,band)*conjg(U(ip,iq,3,band))+U(ip,iq,2,band)*conjg(U(ip,iq,1,band))))
                        fvec(2)=fvec(2)+real(i*(U(ip,iq,4,band)*conjg(U(ip,iq,3,band))+U(ip,iq,4,band)*conjg(U(ip,iq,1,band))))
                        fvec(2)=fvec(2)+real(i*(-U(ip,iq,1,band)*conjg(U(ip,iq,4,band))+U(ip,iq,1,band)*conjg(U(ip,iq,2,band))))
                        fvec(2)=fvec(2)+real(i*(U(ip,iq,3,band)*conjg(U(ip,iq,4,band))-U(ip,iq,3,band)*conjg(U(ip,iq,2,band))))
                        
                        fvec(3)=fvec(3)-real((U(ip,iq,1,band)*conjg(U(ip,iq,3,band))+U(ip,iq,1,band)*conjg(U(ip,iq,1,band))))
                        fvec(3)=fvec(3)+real((U(ip,iq,3,band)*conjg(U(ip,iq,3,band))+U(ip,iq,3,band)*conjg(U(ip,iq,1,band))))
                        fvec(3)=fvec(3)-real((U(ip,iq,2,band)*conjg(U(ip,iq,4,band))+U(ip,iq,2,band)*conjg(U(ip,iq,2,band))))
                        fvec(3)=fvec(3)+real((U(ip,iq,4,band)*conjg(U(ip,iq,4,band))+U(ip,iq,4,band)*conjg(U(ip,iq,2,band))))
                        
                        fvec(4)=fvec(4)-2.d0*real((U(ip,iq,6,band)*conjg(U(ip,iq,7,band))-U(ip,iq,8,band)*conjg(U(ip,iq,5,band))))
                         
                        fvec(5)=fvec(5)+real(i*(U(ip,iq,6,band)*conjg(U(ip,iq,7,band))+U(ip,iq,6,band)*conjg(U(ip,iq,5,band))))
                        fvec(5)=fvec(5)+real(i*(U(ip,iq,8,band)*conjg(U(ip,iq,7,band))+U(ip,iq,8,band)*conjg(U(ip,iq,5,band))))
                        fvec(5)=fvec(5)+real(i*(-U(ip,iq,5,band)*conjg(U(ip,iq,8,band))+U(ip,iq,5,band)*conjg(U(ip,iq,6,band))))
                        fvec(5)=fvec(5)+real(i*(U(ip,iq,7,band)*conjg(U(ip,iq,8,band))-U(ip,iq,7,band)*conjg(U(ip,iq,6,band))))
                        
                        fvec(6)=fvec(6)-real((U(ip,iq,5,band)*conjg(U(ip,iq,7,band))+U(ip,iq,5,band)*conjg(U(ip,iq,5,band))))
                        fvec(6)=fvec(6)+real((U(ip,iq,7,band)*conjg(U(ip,iq,7,band))+U(ip,iq,7,band)*conjg(U(ip,iq,5,band))))
                        fvec(6)=fvec(6)-real((U(ip,iq,6,band)*conjg(U(ip,iq,8,band))+U(ip,iq,6,band)*conjg(U(ip,iq,6,band))))
                        fvec(6)=fvec(6)+real((U(ip,iq,8,band)*conjg(U(ip,iq,8,band))+U(ip,iq,8,band)*conjg(U(ip,iq,6,band))))
                        
                     
                        end do
                end do
        end do

        fvec(1)=fvec(1)/(np*nq)
        fvec(2)=fvec(2)/(np*nq)
        fvec(3)=fvec(3)/(np*nq)
        fvec(4)=fvec(4)/(np*nq)
        fvec(5)=fvec(5)/(np*nq)
        fvec(6)=fvec(6)/(np*nq)
      
return
end subroutine fcn

