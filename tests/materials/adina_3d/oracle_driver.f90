! Test harness only. Original routines are extracted at run time from the
! user's ADINA tree; no copy of those routines is distributed with this test.
program oracle
  implicit double precision (a-h,o-z)
  double precision input
  common /EL/ IND,ICOUNT,NPAR(20),NUMEG,NEGL,NEGNL,IMASS,IDAMP,ISTAT,NDOF,KLIN,IEIG,IMASSN,IDAMPN
  common /VAR/ NG,MODEX,IUPDT,KSTEP,ITEMAX,IEQREF,ITE,KPRI,IREF,IEQUIT,IPRI,KPLOTN,KPLOTE
  common /CRACK/ GAMMA,STIFAC,SHEFAC,SIGMAT,SIGMAC,EPSC,SIGMAU,EPSU,RAM5,RBM5,RCM5,ICRACK,MOD45
  common /CONCRT/ BETA,GAMA,RKAPA,ALFA,SIGP(6),TEP(6),EP(6),YP(3),E,VNU,RK,G,E12,E13,E23,EPSCP,SIGCP,FALSTR,ILFSET
  common /GAUSS/ XG(6,6),WGT(6,6),EVAL2(9,2),EVAL3(27,3),F1,F2,F3
  common /MTMD3D/ D(6,6),STRESS(6),STRAIN(6),IPT,NEL,IPS,ISVE
  dimension ekk(39),sp1(6),sp31(6),sp32(6),sp33(6),sig(6),eps(6),dca(4),crkstr(3)
  dimension dummy(100),nods(27),nod9(27),input(6)
  npar=0; npar(15)=5; npar(10)=1; npar(11)=1
  icount=1; kpri=1; ieqref=1; ipri=1
  stifac=1d-4; shefac=.5d0; sigmat=3; sigmac=-30; epsc=-.002d0; sigmau=-6; epsu=-.006d0
  beta=.5d0; gama=1; rkapa=.7d0; alfa=0
  ekk=0; ekk(1)=30000; ekk(2)=.2d0
  xg=0; wgt=0; wgt(1,1)=2
  xg(1,3)=-sqrt(.6d0); xg(3,3)=sqrt(.6d0)
  wgt(1,3)=5d0/9; wgt(2,3)=8d0/9; wgt(3,3)=5d0/9
  sp1=[0d0,.25d0,.5d0,.75d0,1d0,1.2d0]
  sp31=[1d0,1.4d0,1.7d0,2.2d0,2.5d0,2.8d0]
  sp32=[1.3d0,1.5d0,2d0,2.3d0,2.7d0,3.2d0]
  sp33=[1.25d0,1.45d0,1.95d0,2.25d0,2.65d0,3.15d0]
  sig=0; eps=0; evmax=0; evgrav=0; pgrav=0; dca=0; dca(4)=100; crkstr=0
  nel=1; ipt=1; ips=0; isve=0; dummy=0; nods=0; nod9=0
  do
    read(*,*,iostat=ios) input
    if(ios/=0) exit
    ! OpenSees xy,yz,xz -> ADINA xy,xz,yz. Reorder only at the boundary.
    strain=input; strain(5)=input(6); strain(6)=input(5)
    call CMOD3D(nel,ekk,sp1,sp31,sp32,sp33,sig,eps,evmax,evgrav,pgrav,dca,crkstr, &
      stress,strain,d,ipt,ips,isve,nods,dummy,dummy,dummy,nod9,dummy)
    write(*,'(6(es25.16,1x))') stress(1:4),stress(6),stress(5)
  end do
end program
subroutine THERM3()
  error stop 'thermal branch outside oracle scope'
end subroutine
subroutine CAUCH3()
  error stop 'finite deformation branch outside oracle scope'
end subroutine
