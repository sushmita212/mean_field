subroutine open_files()
use zeeman_mod
implicit none

       
        open(32,file='mid_data/kappa0_mid.txt',action='write',position='append')
        open(33,file='mid_data/kappaxd_mid.txt',action='write',position='append')
        open(34,file='mid_data/kappaxodxy_mid.txt',action='write',position='append')
        open(35,file='mid_data/kappaxodyz_mid.txt',action='write',position='append')
        open(36,file='mid_data/kappaxodxz_mid.txt',action='write',position='append')
        open(37,file='mid_data/kappayd_mid.txt',action='write',position='append')
        open(38,file='mid_data/kappayodxy_mid.txt',action='write',position='append')
        open(39,file='mid_data/kappayodyz_mid.txt',action='write',position='append')
        open(40,file='mid_data/kappayodxz_mid.txt',action='write',position='append')
        open(41,file='mid_data/kappazd_mid.txt',action='write',position='append')
        open(42,file='mid_data/kappazodxy_mid.txt',action='write',position='append')
        open(43,file='mid_data/kappazodyz_mid.txt',action='write',position='append')
        open(44,file='mid_data/kappazodxz_mid.txt',action='write',position='append')
        open(45,file='mid_data/etax_mid.txt',action='write',position='append')
        open(46,file='mid_data/etay_mid.txt',action='write',position='append')
        open(47,file='mid_data/etaz_mid.txt',action='write',position='append')
        open(48,file='mid_data/etatilx_mid.txt',action='write',position='append')
        open(49,file='mid_data/etatily_mid.txt',action='write',position='append')
        open(50,file='mid_data/etatilz_mid.txt',action='write',position='append')
        open(51,file='mid_data/ma_mid.txt',action='write',position='append')
        open(52,file='mid_data/mb_mid.txt',action='write',position='append')
        open(53,file='mid_data/constra_mid.txt',action='write',position='append')
        open(54,file='mid_data/constrb_mid.txt',action='write',position='append')
        open(55,file='mid_data/lama_mid.txt',action='write',position='append')
        open(56,file='mid_data/lamb_mid.txt',action='write',position='append') 
        open(57,file='mid_data/chiA_mid.txt',action='write',position='append')
        open(58,file='mid_data/chiB_mid.txt',action='write',position='append')
        open(101,file='mid_data/mfe_mid.txt',action='write',position='append')
        open(102,file='mid_data/plaq_mid.txt',action='write',position='append')



end subroutine open_files

