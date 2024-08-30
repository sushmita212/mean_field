subroutine gap()
use zeeman_mod
implicit none
       

        
        
        ek1=0.d0;       ek2=0.d0;       em=0.d0;        eg=0.d0
  
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
        
       
        call ham(k1_p)
        hmaj=hb+hos+hcons+hkappa
                
        h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        call zheev(JOBZ, UPLO, n, h, LDA, W, WORK, LWORK, RWORK,INFO)
        ek1=w(5)
        
        call ham(k2_p)
        hmaj=hb+hos+hcons+hkappa
                
        h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        call zheev(JOBZ, UPLO, n, h, LDA, W, WORK, LWORK, RWORK,INFO)
        ek2=w(5)

        call ham(m_p)
        hmaj=hb+hos+hcons+hkappa
                
        h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        call zheev(JOBZ, UPLO, n, h, LDA, W, WORK, LWORK, RWORK,INFO)
        em=w(5)

        call ham(g_p)
        hmaj=hb+hos+hcons+hkappa
                
        h=4.d0*matmul(v1,matmul(hmaj,conjg(transpose(v1))))
        call zheev(JOBZ, UPLO, n, h, LDA, W, WORK, LWORK, RWORK,INFO)
        eg=w(5)

        


       
end subroutine
