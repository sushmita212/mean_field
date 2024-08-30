program main
use zeeman_mod
implicit none
double precision :: start, finish
integer :: iter,h_iter
       
        call initiate()

        call open_files()
        call open_final()
        call cpu_time(start)
        do h_iter=1,h_div

                h_mag=0.d0+h_max*(h_iter-1)/(h_div-1)
                hx=zeeman_info*h_mag/sqrt(3.d0)
                hy=zeeman_info*h_mag/sqrt(3.d0)
                hz=zeeman_info*h_mag/sqrt(3.d0)
                print *, 'h_mag=',h_mag 
                print *, 'hx,hy,hz=', hx,hy,hz     
        do iter=1,n_iter
                print *, "iteration:",iter
                temp_m=m
                temp_x=x
                call sc_eq()
               if(constr_real>cutoff) then
                        print *, "calling hybrd"        
                        !if(iter>150)then
                         !       n_iter=n_iter+30
                        !end if        
                        m=temp_m
                        x=temp_x
                        call hybrd_solve()
                        call sc_eq()
                end if
                call mfe()
                call plaquette()
                !call gap() 
                call write_data()              
        end do
                print *, 'checking convergence'
                call conv()
                call gap()
                call write_final()
 
              
        
        end do

        call cpu_time(finish)
        print *, "time taken=",(finish-start)/60.d0,"mins" 
                
       
        
        
end program main
