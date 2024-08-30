subroutine mfe()
use zeeman_mod
implicit none

        e_k=-K*((real(x(1))*real(x(4))+real(x(2))*real(x(17))+real(x(3))*real(x(30)))-(real(x(31))*real(x(40))+real(x(35))*real(x(44))+real(x(39))*real(x(48)))&
-4.d0*(real(m(1))*real(m(4))+real(m(2))*real(m(5))+real(m(3))*real(m(6))))

        e_z=2.d0*(hx*(real(m(1))+real(m(4)))+hy*(real(m(2))+real(m(5)))+hz*(real(m(3))+real(m(6)))) 

        e_tot=e_k+e_z



end subroutine
