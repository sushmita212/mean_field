subroutine open_final()
use zeeman_mod
implicit none

        open(59,file='h_mag.txt',action='write',position='append')
        open(60,file='kappa0.txt',action='write',position='append')
        open(61,file='kappaxd.txt',action='write',position='append')
        open(62,file='kappaxodxy.txt',action='write',position='append')
        open(63,file='kappaxodyz.txt',action='write',position='append')
        open(64,file='kappaxodxz.txt',action='write',position='append')
        open(65,file='kappayd.txt',action='write',position='append')
        open(66,file='kappayodxy.txt',action='write',position='append')
        open(67,file='kappayodyz.txt',action='write',position='append')
        open(68,file='kappayodxz.txt',action='write',position='append')
        open(69,file='kappazd.txt',action='write',position='append')
        open(70,file='kappazodxy.txt',action='write',position='append')
        open(71,file='kappazodyz.txt',action='write',position='append')
        open(72,file='kappazodxz.txt',action='write',position='append')
        open(73,file='etax.txt',action='write',position='append')
        open(74,file='etay.txt',action='write',position='append')
        open(75,file='etaz.txt',action='write',position='append')
        open(76,file='etatilx.txt',action='write',position='append')
        open(77,file='etatily.txt',action='write',position='append')
        open(78,file='etatilz.txt',action='write',position='append')
        open(79,file='ma.txt',action='write',position='append')
        open(80,file='mb.txt',action='write',position='append')
        open(81,file='constra.txt',action='write',position='append')
        open(82,file='constrb.txt',action='write',position='append')
        open(83,file='lama.txt',action='write',position='append')
        open(84,file='lamb.txt',action='write',position='append') 
        open(85,file='chiA.txt',action='write',position='append')
        open(86,file='chiB.txt',action='write',position='append')
        open(103,file='mfe.txt',action='write',position='append')
        open(104,file='plaq.txt',action='write',position='append')
        open(105,file='gap.txt',action='write',position='append')



end subroutine

