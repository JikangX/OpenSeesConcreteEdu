! Material-point driver for the original JiangJJ2007.for routine Get_DS
! (compiled unchanged). Emulates MSC.Marc HYPELA2 calls for one integration
! point with one converged call per increment: e = strain at increment start,
! de = increment, s = stress at increment start; the returned stress becomes
! the next start (Marc's converged state).
! stdin: fc E0 nu epsCrush, then one total strain row (6, engineering shear)
! per increment.
! stdout per increment: stress(6) NCrack Beta(max) D11 D12 D44 0
program oracle
  use TypeDef
  use My_Mod
  implicit none
  type(typ_Concrete) :: C
  real*8 :: e(6), de(6), s(6), d(6,6), enew(6), cptim, timinc
  integer :: ios
  read(*,*) C%fc, C%E0, C%ENU, C%EPS_Crush
  C%NCrack=0; C%Beta=0.; C%SIG=0.; C%EPS=0.; C%Stress=0.; C%Strain=0.
  e=0.; s=0.; cptim=0.; timinc=1.
  do
    read(*,*,iostat=ios) enew
    if (ios/=0) exit
    de=enew-e
    call Get_DS(C, d, e, de, s, 6, cptim, timinc)
    write(*,'(12es25.16e3)') s, dble(C%NCrack), C%Beta, d(1,1), d(1,2), d(4,4), 0.d0
    e=enew
    cptim=cptim+timinc
  end do
end program
