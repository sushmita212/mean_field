module zeeman_mod

implicit none
double precision ::K,hx,hy,hz,n1(2),n2(2),f1(2),f2(2),kvec(2,2),pi,lm(6),constr_real,constr_imag,err_real(6),err_imag(6),del_x(2),del_y(2),del_z(2),h_pert,h_mag,h_max
double precision ::zeeman_info,kappa_info
double precision ::e_k,e_z,e_tot,w_avg,ek1,ek2,em,eg
integer :: np,nq,n_iter,h_div
double complex, dimension(:,:,:,:), allocatable ::U
double complex, parameter :: i=(0.d0,1.d0)
!double complex :: MAA(4,4), MAB(4,4), MBB(4,4)
double complex :: x(54), m(6)
double complex :: temp_x(54),temp_m(6)
double complex :: zero4(4,4),om_v1(4,4),v1(8,8)
double complex :: MAB(4,4),ABcc,ABcx,ABcy,ABcz,ABxc,AByc,ABzc,ABxx,AByy,ABzz
double complex :: MAA(4,4),AAcx,AAcy,AAcz,AAxc,AAyc,AAzc,AAxy,AAxz,AAyx,AAyz,AAzx,AAzy
double complex :: MBB(4,4),BBcx,BBcy,BBcz,BBxc,BByc,BBzc,BBxy,BBxz,BByx,BByz,BBzx,BBzy
double complex :: MCA(4,4),MCB(4,4)
double complex :: kappaAA(4,4),kappaBB(4,4),kappaAB(4,4),kappaBA(4,4),kappa2AA,kappaiAA
double complex :: hb(8,8),hos(8,8),hcons(8,8),hmaj(8,8),hkappa(8,8),h(8,8)
integer :: iter_max,stability
double precision :: cutoff,conv_cut
double precision :: k1_p(2),k2_p(2),m_p(2),g_p(2)
 !LAPACK PARAMETERS
!!NEED TO CHECK DIMENSIONS
character :: JOBZ,UPLO
integer :: N,LDA,LWORK,INFO
double precision :: W(8),RWORK(40)
double complex ::WORK(20) 


end module zeeman_mod
