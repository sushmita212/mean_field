subroutine write_data()
use zeeman_mod
implicit none
character(len=10) :: w_fmt
        
                !w_fmt='(3es20.12)'        

                !kappa0
                write(32,*) real(x(1)),real(x(2)),real(x(3))
                !kappax
                write(33,*) real(x(4)),real(x(8)),real(x(12))
                write(34,*) real(x(5)),real(x(7))
                write(35,*) real(x(9)),real(x(11))
                write(36,*) real(x(6)),real(x(10))
                !kappay
                write(37,*) real(x(13)),real(x(17)),real(x(21))
                write(38,*) real(x(14)),real(x(16))
                write(39,*) real(x(18)),real(x(20))
                write(40,*) real(x(15)),real(x(19))
                !kappaz
                write(41,*) real(x(22)),real(x(26)),real(x(30))
                write(42,*) real(x(23)),real(x(25))
                write(43,*) real(x(27)),real(x(29))
                write(44,*) real(x(24)),real(x(28))
                !etax
                write(45,*) real(x(31)),real(x(32)),real(x(33))
                !etay
                write(46,*) real(x(34)),real(x(35)),real(x(36))
                !etaz
                write(47,*) real(x(37)),real(x(38)),real(x(39))
                !etatilx
                write(48,*) real(x(40)),real(x(41)),real(x(42))
                !etatily
                write(49,*) real(x(43)),real(x(44)),real(x(45))
                !etatilz
                write(50,*) real(x(46)),real(x(47)),real(x(48))
                !ma
                write(51,*) real(m(1)),real(m(2)),real(m(3))
                !mb
                write(52,*) real(m(4)),real(m(5)),real(m(6))
                !constraint
                write(53,*) err_real(1),err_real(2),err_real(3)
                write(54,*) err_real(4),err_real(5),err_real(6)
                !lm
                write(55,*) lm(1),lm(2),lm(3)
                write(56,*) lm(4),lm(5),lm(6)
                !NNN parameters
                write(57,*) real(x(49)),real(x(50)),real(x(51))
                write(58,*) real(x(52)),real(x(53)),real(x(54))
                
                write(101,*) e_k,e_z,e_tot
                write(102,*)  w_avg
 end subroutine
