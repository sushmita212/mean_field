program test
implicit none
integer :: p,q,r
double precision :: pi
character(len=10) :: w_fmt
        
         
        open(23,file='pipi.txt',action='write',position='append')

        w_fmt='(3es20.12)'        
        p=61
        q=5
        r=2.d0*p/3.d0
        print *, r

        pi=4.d0*atan(1.d0)      

        print *, pi,pi,pi,pi
        print(w_fmt), pi,pi,pi,pi

        write(23,*) pi,pi,pi
        write(23,w_fmt)pi,pi,pi

        close(23)
end program
