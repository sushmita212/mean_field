program comp_to_r
implicit none
double complex :: x(3),x1,x2(2)
double precision :: xr(3), xrp(3)

        x=(/(0.5d0,0.1d0),(9.3d0,9.3d0),(0.4d0,4.5d0)/)
        x1=(-0.5d0,0.1d0)
        x2=(/-(9.3d0,9.3d0),(0.4d0,4.5d0)/)
        x(1)=x1
        x(2:3)=x2
        
        x=reshape((/x1,x2/),(/3/))
        xr=real(x)
        
        print *, x
        print *,xr
        
        xrp=abs(xr)
        
        print *, xr

        print *, xrp-xr
        print *, maxval(xrp-xr)
end program
