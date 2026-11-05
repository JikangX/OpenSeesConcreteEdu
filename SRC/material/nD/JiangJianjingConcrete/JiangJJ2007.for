C	=====================================================
C	清华大学研究生课程《钢筋混凝土有限元》教学程序
C    非线性弹性全量模型(江见鲸模型)
C	利用MSC. Marc软件提供的 HYPELA2 子程序
C	主要参考文献：
C	1. MARC Volumn D user subroutine and Special Routines"
C	3. 过镇海 "钢筋混凝土原理"
C	=====================================================
	module TypeDef !定义混凝土材性模块
		type :: typ_Concrete
			real*8 fc, E0, ENU, EPS_Crush;  !抗压强度+，初始弹性模量，初始泊松比，压碎应变-
			real*8 Es, ENUs, Et             !割线模量，割线泊松比,  切线模量
			integer NCrack	                !判断是否破坏
			real*8 SIG(6), EPS(6),dEPS(6);  !开始时应力，应变， 应变增量
			real*8 Stress(6), Strain(6)     !结束时应力，应变
			real*8 Beta                     !非线性指标
			real*8 Dela(6,6),  Ds(6,6)      !刚度矩阵，弹性刚度矩阵，割线刚度矩阵
		end type typ_Concrete
	end module TypeDef

	module My_MOD !开辟公共变量空间
		use TypeDef
		integer FisrtStep                                              !判断是否是第一步计算
		real*8 time,dtime                                              !当前时间，时间增量
		type(typ_Concrete) :: My_ConA(1000,8)                          !定义混凝土数组
		type(typ_Concrete) :: My_ConB(size(My_ConA,1),size(My_ConA,2)) !定义数据备份
	end module


	subroutine Get_DS(C, D, E,DE,S, NGENS,cptim, timinc)
!     C      混凝土结构体(in & out)
!	D(6x6) 迭代本构矩阵(out)
!	E(6)   开始时刻的应变(in)
!	DE(6)  应变增量(in)
!	S(6)   开始时刻的应变(in & out)
!	NGENS  应变维数(in)
!	cptim  当前计算时间(in)
!	timinc 当前时间增量(in)

		use typedef 
		use My_Mod
		implicit none

		type(typ_concrete) :: C  ! 混凝土结构体
		integer :: n(2), ngens
		real*8 :: e(ngens),de(ngens),d(ngens,ngens),s(ngens) 		
		real*8 :: cptim, timinc

		real*8 Beta1,strain_m    ! 非线性指标临时变量
		real*8 s_m,J2,J3,r,sita,TempA, TempB, TempC;

		open(77,file='debug.txt',position='append') ! 打开程序调试文件
		C%SIG=s; C%EPS=e; C%dEPS=de;                ! 应变赋值
c	------------------------------------
		call Get_D(C%Dela, C%E0, C%ENU,ngens)       !得到弹性矩阵
c	--------------------------------------
C		判断是否压坏
		strain_m=(C%EPS(1)+C%EPS(2)+C%EPS(3))/3.
		if(Strain_m>0.) Strain_m=0.
		if(minval(C%EPS(1:3)-Strain_m)<C%EPS_Crush) then
			C%NCrack=100 !彻底破坏
		end if

c	如果已经压碎，应力清零，刚度为很小值，结束计算
		if(C%NCrack==100) then
			D=1E-6*C%Dela
			C%SIG=0.;
			C%Stress=0.;
			s=0.
			return
		end if

C	计算主应力和割线刚度
		C%Stress=C%SIG
		C%Strain=C%EPS+C%dEPS;
		s_m=(C%Stress(1)+C%Stress(2)+C%Stress(3))/3. !计算平均应力
		s(1:3)=C%Stress(1:3)-s_m                     !计算应力偏量
		s(4:6)=C%Stress(4:6)                         !计算应力偏量
		J2=-s(1)*s(2)-s(2)*s(3)-s(3)*s(1)+s(4)**2+s(5)**2+s(6)**2 !计算J2
		J3=s(1)*s(2)*s(3)+2.*s(4)*s(5)*s(6) ! 计算J3
	1		-s(1)*s(5)**2-s(2)*s(6)**2-s(3)*s(4)**2
		r=sqrt(4.*J2/3.)
		if(r.ne.0.) then
			sita=acos(4.*J3/r**3)/3.                 !计算sita
			if(4.*J3/r**3> 1.d0) sita=acos( 1.d0)/3. !计算sita
			if(4.*J3/r**3<-1.d0) sita=acos(-1.d0)/3. !计算sita
		else
			sita=0.                                  !计算sita
		end if

		TempA=1.2856/C%fc**2;                        !求解过程参看课本例4.1
		TempB=(1.4268+10.2551*cos(sita))/C%fc;
		TempC=3.2128*s_m*3./C%fc-1.;
		if(TempB**2-4.*TempA*TempC>=0.) then
			Beta1=-TempB+sqrt(TempB**2-4.*TempA*TempC)
			Beta1=Beta1/2./TempA
			Beta1=sqrt(J2)/Beta1 !根据江见鲸模型求解Beta
		else
			write(77,*) "Wrong"                     ! 超出
			!stop
			Beta1=1.;
		end if
		if(Beta1>C%Beta) then !Beta应该始终增大(对于全量模型)
			C%Beta=Beta1
		end if
		! 计算割线刚度和泊松比
		if(Beta1>1.) Beta1=1.
		
		C%Es=C%E0/2.*(1.+sqrt(1.-Beta1))             !计算割线模型
		C%Et=sqrt(1.-Beta1)*C%E0                     !计算切线模型
c		if(Beta1<0.8) C%ENUs=C%ENU
c		if(Beta1.ge.0.8) C%ENUs=0.42-(0.42-C%ENU)*
c	1			sqrt(1.-((Beta1-.8)/.2)**2)
		C%ENUs=C%ENU
		
		call Get_D(C%Ds, C%Es, C%ENUs, ngens)        !计算割线刚度矩阵
		
		if(Beta1<0.99999d0) then                     !如果没有达到极限应力
			C%Stress=matmul(C%Ds,C%Strain);
		else                                         !达到极限应力后应力不变
			C%Stress=C%SIG
		end if

C--------------------------------------------------------------------------------------------
		call Get_D(D, C%Et, C%ENUs, ngens)                  !设置迭代刚度矩阵  !!!ty
		D=D+1.d-6*C%Dela                             !提高收敛性
		s=C%Stress

c		write(77,*) 'Step: ',time, dtime
c		write(77,*) Beta1, cptim, timinc
c		write(77,*) C%NCrack
c		write(77,*) C%Strain(1:3)
c		write(77,*) C%Stress(1:3)
c		write(77,*)
		close(77)
		return
	end subroutine

	subroutine Simple_Concrete(D,E,DE,S,NGENS,N,NN,KC,matus,cptim, timinc)
		use typedef 
		use My_Mod
		implicit none
		integer ngens, N(2), NN, KC, matus
		real*8 D(ngens,ngens), E(ngens), DE(ngens), S(ngens)
		real*8 cptim, timinc
	
		type(typ_concrete) :: C
c
		if(cptim==0..and.timinc==0.) then                                 ! 初始化
			call Initial_ALL(cptim, timinc)
		end if

		if(cptim>time) then                                               ! 新的荷载步开始
			call New_Step (cptim,timinc)
		end if

		if(cptim==time.and. timinc< dtime.and.cptim.ne.0) then            ! 发生荷载步折半
			call Cut_Step (cptim,timinc)
		end if

		C=My_ConA(n(1),nn)
		call Get_DS(C, D, E,DE,S,NGENS,cptim, timinc)                     ! 调用核心计算程序
		My_ConB(n(1),nn)=C
		return
	end

	subroutine Cut_Step(cptim,timinc)      ! 发生荷载步折半，恢复数据
		use typedef  
		use My_Mod
		implicit none
		real*8 :: cptim, timinc
		dtime=timinc;
		My_ConB=My_ConA
		return
	end subroutine Cut_Step
	
	subroutine New_Step(cptim,timinc)      ! 新荷载步，更新变量
		use typedef 
		use My_Mod
		implicit none
		real*8 :: cptim, timinc
		time=cptim; dtime=timinc;
		My_ConA=My_ConB
		return
	end subroutine New_Step

	subroutine Initial_All(cptim, timinc)  ! 第一次计算，赋予初值
		use typedef 
		use My_Mod
		implicit none
		integer I,J
		real*8 :: cptim, timinc
		time=cptim; dtime=timinc
		return
	end subroutine Initial_ALl

c	marc 接口程序
      subroutine hypela2(d,g,e,de,s,t,dt,ngens,n,nn,kc,matus,ndi,nshear,
     2                   disp,dispt,coord,ffn,frotn,strechn,eigvn,ffn1,
     3                   frotn1,strechn1,eigvn1,ncrd,itel,ndeg,ndm,
     4                   nnode,jtype,lclass,ifr,ifu)
		IMPLICIT REAL *8 (A-H, O-Z)
      common/creeps/cptim,timinc,timinc_p,timinc_s,timincm,
     1 timinc_a,timinc_b,creept(33),icptim,icfte,icfst,
     2 icfeq,icftm,icetem,mcreep,jcreep,icpa,icftmp,icfstr,
     3 icfqcp,icfcpm,icrppr,icrcha,icpb,iicpmt,iicpa

		dimension e(ngens),de(ngens),t(*),dt(*),g(ngens),
	1		d(ngens,ngens),s(ngens)
		dimension n(2),coord(ncrd,*),disp(ndeg,*),matus(2),
     *          dispt(ndeg,*),ffn(itel,*),frotn(itel,*),
     *          strechn(itel),eigvn(itel,*),ffn1(itel,*),
     *          frotn1(itel,*),strechn1(itel),eigvn1(itel,*)

	call Simple_Concrete(D,E,DE,S,NGENS,N,NN,KC,matus(1),cptim, timinc)

	return
	end

	SUBROUTINE UBGINC(INC,INCSUB) ! 赋予总体参数
		use typedef 
		use My_Mod
		IMPLICIT REAL *8 (A-H, O-Z)

      common/creeps/cptim,timinc,timinc_p,timinc_s,timincm,
     1 timinc_a,timinc_b,creept(33),icptim,icfte,icfst,
     2 icfeq,icftm,icetem,mcreep,jcreep,icpa,icftmp,icfstr,
     3 icfqcp,icfcpm,icrppr,icrcha,icpb,iicpmt,iicpa

		type(typ_concrete) :: C
		integer :: I,J,II,NP
		if(FisrtStep.ne.123.and.cptim.eq.0.and.timinc.eq.0 ) then ! 第一次计算
			FirstStep=123
			open(66,file='debug0.txt')        !!!ty
				write(66,*)	cptim, timinc
			close(66)
c		初次计算，清零并赋值
			C%fc=30.; C%E0=30e3; C%ENU=0.18;
			C%NCrack=0; C%NCrack=0; C%EPS_Crush=-0.0033; 
			C%Beta=0;	C%Beta=0;
			C%SIG=0.; C%EPS=0.; C%Stress=0.; C%Strain=0.;
			do I=1, size(My_ConA,1)
				do J=1, size(My_ConA,2)
					My_ConA(I,J)=C; My_ConB(I,J)=C; 
				end do
			end do
		end if
c	----------------------------------------------------


	RETURN
	END

	subroutine Get_D(D, E, ENU,ngens)
c	计算弹性矩阵
		implicit none
		integer ngens, K1, K2
		real*8 D(ngens,ngens), E, ENU
		D=0.
		do K1=1, 3
			do K2=1, 3
				D(K2,K1)=ENU
			end do
			D(K1,K1)=1.-ENU
		END do
		do K1=4,6
			D(K1,K1)=(1.-2.*ENU)*0.5
		end do
		D=D*E/(1.+ENU)/(1.-2.*ENU)
		return
	end subroutine
