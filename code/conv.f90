subroutine conv()
use zeeman_mod
implicit none

integer :: iter,stability_count
double complex :: par_old(60), par_new(60)
double precision :: par_old_r(60), par_new_r(60), diff

                open(354,file='iter.txt',action='write',position='append')

                par_old=reshape((/temp_x,temp_m/),(/60/))
                par_new=reshape((/x,m/),(/60/))

                par_old_r=real(par_old)
                par_new_r=real(par_new)

                diff=maxval(abs(par_new_r-par_old_r))

                print *, 'diff=',diff
              
                iter=n_iter
                stability_count=0
                do while ((diff>conv_cut.or.stability_count<stability).and.iter<iter_max)
                        !print *, 'iteration=',iter+n_iter
                        temp_m=m
                        temp_x=x
                        call sc_eq()
                        if(constr_real>cutoff) then
                                print *, "calling hybrd"        
                                      
                                m=temp_m
                                x=temp_x
                                call hybrd_solve()
                                call sc_eq()
                        end if
                        call mfe()
                        call plaquette()
                        !call gap() 
                        call write_data()              
                
                        par_old=reshape((/temp_x,temp_m/),(/60/))
                        par_new=reshape((/x,m/),(/60/))

                        par_old_r=real(par_old)
                        par_new_r=real(par_new)

                        diff=maxval(abs(par_new_r-par_old_r))
                        iter=iter+1
                        print *, 'iteration=',iter        
                        print *, 'diff=',diff
                        
                        if(diff<conv_cut) then
                                stability_count=stability_count+1
                                print *, 'stable for',  stability_count, 'iterations'
                        end if
                        
                end do

                write(354,*) iter

end subroutine
