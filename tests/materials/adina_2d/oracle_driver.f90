! Oracle harness: run the ORIGINAL native A06 routines without editing them.
! The runner extracts them into its temporary build directory at test time.
program oracle2d
  implicit double precision (a-h,o-z)
  double precision input
  character(len=512) parameter_file
  common /EL/ IND,ICOUNT,NPAR(20),NUMEG,NEGL,NEGNL,IMASS,IDAMP,ISTAT,NDOF,KLIN,IEIG,IMASSN,IDAMPN
  common /VAR/ NG,MODEX,IUPDT,KSTEP,ITEMAX,IEQREF,ITE,KPRI,IREF,IEQUIT,IPRI,KPLOTN,KPLOTE
  common /CRACK/ GAMMA,STIFAC,SHEFAC,SIGMAT,SIGMAC,EPSC,SIGMAU,EPSU,RAM5,RBM5,RCM5,ICRACK,MOD45
  common /CONCR2/ BETA,GAMA,RKAPA,ALFA,SIGP(4),TEP(4),EP(4),YP(3),E,VNU,RK,G,E12,E14,E24,EPSCP,SIGCP,FALSTR,ILFSET
  common /GAUSS/ XG(6,6),WGT(6,6),EVAL2(9,2),EVAL3(27,3),F1,F2,F3
  common /MATMOD/ STRESS(4),STRAIN(4),D(4,4),IPT,NEL,IPS,ISVE
  dimension ekk(39),sp1(6),sp31(6),sp32(6),sp33(6),sig(4),eps(4),crkstr(3)
  dimension dummy(100),nods(9),nod5(9),input(3)
  dimension specimen_parameters(37)
  npar=0; npar(5)=2; npar(15)=5; npar(10)=1
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
  ! Optional actual specimen parameters in the same explicit order as the
  ! native factory. No original ADINA routine or COMMON layout is changed.
  call get_command_argument(1,parameter_file)
  if(len_trim(parameter_file)>0) then
    open(unit=91,file=trim(parameter_file),status='old',action='read')
    read(91,*) specimen_parameters
    close(91)
    ekk(1:2)=specimen_parameters(1:2)
    sigmat=specimen_parameters(3);sigmac=specimen_parameters(4)
    epsc=specimen_parameters(5);sigmau=specimen_parameters(6);epsu=specimen_parameters(7)
    stifac=specimen_parameters(8);shefac=specimen_parameters(9)
    beta=specimen_parameters(10);gama=specimen_parameters(11)
    rkapa=specimen_parameters(12);alfa=specimen_parameters(13)
    sp1=specimen_parameters(14:19);sp31=specimen_parameters(20:25)
    sp32=specimen_parameters(26:31);sp33=specimen_parameters(32:37)
  end if
  sig=0; eps=0; evmax=0; evgrav=0; pgrav=0; angle=1000; crkstr=0
  nel=1; ipt=1; ips=0; isve=0; dummy=0; nods=0; nod5=0
  do
    read(*,*,iostat=ios) input
    if(ios/=0) exit
    strain(1:3)=input; strain(4)=0
    call CDMOD(nel,ekk,sp1,sp31,sp32,sp33,sig,eps,evmax,evgrav,pgrav,angle,crkstr, &
      stress,strain,d,ipt,ips,isve,nods,dummy,dummy,dummy,nod5,dummy)
    write(*,'(3(es25.16,1x))') stress(1:3)
  end do
end program
subroutine THERM2()
  error stop 'thermal branch outside oracle scope'
end subroutine
subroutine CAUCHY()
  error stop 'finite deformation branch outside oracle scope'
end subroutine
