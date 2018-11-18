
      Program DFSFNFFS
      IMPLICIT NONE
!     Data decleration
      INTEGER           N1,N3
      PARAMETER         (N1=4)
      INTEGER           LWORK,LRWORK
      PARAMETER         (LWORK=2*N1,LRWORK=2*N1)
      INTEGER           LWORKhat
      PARAMETER         (LWORKhat=4*2*N1)
!
      INTEGER           INFO,thth,tg,jg,gj
!!
      COMPLEX*16        VL(1,1),VR(N1,N1),W(N1),WORK(LWORK)
      COMPLEX*16        WORK2(LWORKhat)
      DOUBLE PRECISION  RWORK(LRWORK)
!
      EXTERNAL ZGEEV,ZGETRI,ZGETRF
      COMPLEX*16  Imxkn1(4,4,110),Imxfkn(4,4,110),Imxkn(4,4,110),&
SgImkn(4,4,110),iSgImkn(4,4,110),Gkn(4,4,110),Idmxkn(4,4,110),&
Gnewkn(4,4,110),Mkn(4,4,110),AD(4,4),EIG(4,4)
      COMPLEX*16  ANCOMMsmspkn1(4,4,110),ANCOMMsmspkn2(4,4,110),&
ANCOMMsmspkn(4,4,110),SdgImkn(4,4,110)
      COMPLEX*16  ANCOMMsmdeckn1(4,4,110),ANCOMMsmdeckn2(4,4,110),&
ANCOMMsmdeckn(4,4,110),TriSgImkn(160),TdriSgImkn(160),&
TriSgImkn2(160),TriSgImkn3(160),TriSgImkn4(160)
      COMPLEX*16  TriSgxImkn(160),TriSgyImkn(160),TriSgzzImkn(160)
      COMPLEX*16  TriSgx4Imkn(160),Sgx4Imkn(4,4,110)
      COMPLEX*16  SgxImkn(4,4,110),SgyImkn(4,4,110),SgzzImkn(4,4,110),&
SgImkn2(4,4,110),SgImkn3(4,4,110),SgImkn4(4,4,110)
      COMPLEX*16  ANCOMMdecspkn1(4,4,110),ANCOMMdecspkn2(4,4,110),&
ANCOMMdecspkn(4,4,110),rtGkn(4,4,110),itGkn(4,4,110)
      COMPLEX*16  Dkn(4,4,110),Dknt(4,4,110),Dkn2(4,4,110),DknI(4,4,110)&
,Dkn1v2(4,4,110),Dkn1vp1(4,4,110),INEIG(4,4,110),EIGs(4,4,110),&
inEIGs(4,4),G_A(4,4,110),G_K(4,4,110)
      COMPLEX*16  Egvec(4,4,110),inEgvec(4,4,110),landa(4,4,110),&
invlanda(4,4,110),tolGkn(4,4,110),inEIGs2(4,4,110)
      COMPLEX*16  Gdecn(4,4,110),Gdcn2(4,4,110),Gdecn1(4,4),Gdecn2&
(4,4,110),nsigma(4,4,110),DknII(4,4)
      COMPLEX*16  tau3sigma0(4,4),tau3sigma1(4,4),sigmay(4,4),sigmazz(4&
,4),tau3sigma3(4,4),tau0sigma0(4,4),tau0sigma3(4,4)
      COMPLEX*16  tau3sigma2(4,4)
      DOUBLE PRECISION t(2),scale,theta3,OMEGA,EJ
      DOUBLE PRECISION TcvEthF1,TcvEthF2,TcvEthN,hvTc,Gn,lkg,lkg2
      DOUBLE PRECISION RelConductance(125),Conductance(125)
      INTEGER tt,ee,vv,i,j,yuy
      COMPLEX*16 errIn,hg
      INTEGER dd,s,q,rr,zx,f,h,u,v,y,ii,jj,ll,lll,qq,ww,yy,qqq,www,&
tmpr,nphi,kj,kjjj,k,l,hhh
      DOUBLE PRECISION phi
      DOUBLE PRECISION pi
      DOUBLE PRECISION AIn(220,160),AInsx(220,160),errGkn(4,4,110),&
AIn2(220,160),AIn3(220,160),AIn4(220,160)
	  DOUBLE PRECISION A10In(220,160),A20In(220,160),A30In(220,160),&
A40In(220,160),A50In(220,160)
	  DOUBLE PRECISION A10In2(220,160),A20In2(220,160),A30In2(220,160),&
A40In2(220,160),A50In2(220,160)
	  DOUBLE PRECISION A10In3(220,160),A20In3(220,160),A30In3(220,160),&
A40In3(220,160),A50In3(220,160)
	  DOUBLE PRECISION A10In4(220,160),A20In4(220,160),A30In4(220,160),&
A40In4(220,160),A50In4(220,160),uyu
!!C	 
	  DOUBLE PRECISION A10Insx(220,160),A20Insx(220,160),A30Insx(220,160),A40Insx(220,160),A50Insx(220,160)
	  DOUBLE PRECISION A10Insy(220,160),A20Insy(220,160),A30Insy(220,160),A40Insy(220,160),A50Insy(220,160)
	  DOUBLE PRECISION A10Insz(220,160),A20Insz(220,160),A30Insz(220,160),A40Insz(220,160),A50Insz(220,160)
	  DOUBLE PRECISION A10Insx4(220,160),A20Insx4(220,160),A30Insx4(220,160),A40Insx4(220,160),A50Insx4(220,160)
	  DOUBLE PRECISION ksA10In(220,160),ksA20In(220,160),ksA30In(220,160),ksA40In(220,160),ksA50In(220,160)
      DOUBLE PRECISION ksA10Insx(220,160),ksA20Insx(220,160),ksA30Insx(220,160),ksA40Insx(220,160),ksA50Insx(220,160)
      DOUBLE PRECISION ksA10Insy(220,160),ksA20Insy(220,160),ksA30Insy(220,160),ksA40Insy(220,160),ksA50Insy(220,160)
	  DOUBLE PRECISION ksA10Insz(220,160),ksA20Insz(220,160),ksA30Insz(220,160),ksA40Insz(220,160),ksA50Insz(220,160)
	  DOUBLE PRECISION AInsy(220,160),AInsz(220,160)
      DOUBLE PRECISION mazx10(160),mazx20(160),mazx30(160),mazx40(160),mazx50(160)
	  DOUBLE PRECISION mazxsx10(160),mazxsx20(160),mazxsx30(160),mazxsx40(160),mazxsx50(160)
	  DOUBLE PRECISION mazxsy10(160),mazxsy20(160),mazxsy30(160),mazxsy40(160),mazxsy50(160)
	  DOUBLE PRECISION mazxsz10(160),mazxsz20(160),mazxsz30(160),mazxsz40(160),mazxsz50(160)
      DOUBLE PRECISION mazx4(160)
      DOUBLE PRECISION h1vTc,h0vTc,cte,lkgncte,lvxiF1,Qxi,Qxi1,Qxi2
      DOUBLE PRECISION TH(2),WkvDelta(2),lvxiF2,lvxiN
      DOUBLE PRECISION alpha(1000000),Ikn(160),errIkn(160),tolIn(160)
      DOUBLE PRECISION In(160),ksAIn(220,160),Insx(160),ksAInsx(220,160)&
,In2(160),In3(160),In4(160)
      DOUBLE PRECISION Insx2(160),ksAInsx4(220,160),AInsx4(220,160)
      DOUBLE PRECISION Iknsx(160),errIknsx(160),tolInsx(160)
      DOUBLE PRECISION Iknsx4(160),errIknsx4(160),tolInsx4(160)
      DOUBLE PRECISION Iknsy(160),errIknsy(160),tolInsy(160)
      DOUBLE PRECISION Iknsz(160),errIknsz(160),tolInsz(160),tolIn2(160),&
tolIn3(160),tolIn4(160),Ikn2(160),Ikn3(160),Ikn4(160)
      DOUBLE PRECISION Insy(160),ksAInsy(220,160),Meanrho1(160,100)
      DOUBLE PRECISION Insz(160),ksAInsz(220,160),Insx4(160)
      INTEGER d,kk,kkk,ttt,s1,s2,s3,kjj
      INTEGER IPIV(4)
      DOUBLE PRECISION betha,gamma(160),lkgncte1,Meanrho2(160,100)
	  DOUBLE PRECISION gamma1(160),gamma2(160)
	  COMPLEX* 16 lkgncte2
      COMPLEX*16 one1
      COMPLEX*16 one(4,4)  
	
!     Open output files

      OPEN (UNIT=10,FILE='t19xyh10t2L5p1is0tstpdi.out',STATUS='unknown')
      OPEN (UNIT=81,FILE='integrall1h1025ebtedaas.dat',STATUS='unknown')
	  OPEN (UNIT=82,FILE='theta.dat',STATUS='unknown')
	  OPEN (UNIT=83,FILE='phi.dat',STATUS='unknown')
      OPEN (UNIT=11,FILE='Fthetapi4phipi.out',STATUS='unknown')
	  OPEN (UNIT=12,FILE='t19xyh10t2L5p13tstr2pdi.out',STATUS='unknown')
	  OPEN (UNIT=13,FILE='fT19hth053L1t0tstr3pdi.out',STATUS='unknown')
	  OPEN (UNIT=14,FILE='t19xyh10t2L5p1tstxpdi.out',STATUS='unknown')
      OPEN (UNIT=15,FILE='t19xyh10t2L5p1tstypdi.out',STATUS='unknown')
      OPEN (UNIT=16,FILE='t19xyh10t2L5p1tstzpdi.out',STATUS='unknown')
	  OPEN (UNIT=17,FILE='t19xyh10t2L5p1qpx4pdi.out',STATUS='unknown')
      OPEN (UNIT=21,FILE='fn10xyh10t2L5p1tstsphi.out',STATUS='unknown')
      OPEN (UNIT=31,FILE='fn10xyh10t2L5p1tstxphi.out',STATUS='unknown')
	  OPEN (UNIT=41,FILE='fn10xyh10t2L5p1tstyphi.out',STATUS='unknown')
      OPEN (UNIT=51,FILE='fn10xyh10t2L5p1tstzphi.out',STATUS='unknown')
	  OPEN (UNIT=22,FILE='fn20xyh10t2L5p1tstsphi.out',STATUS='unknown')
      OPEN (UNIT=32,FILE='fn20xyh10t2L5p1tstxphi.out',STATUS='unknown')
	  OPEN (UNIT=42,FILE='fn20xyh10t2L5p1tstyphi.out',STATUS='unknown')
	  OPEN (UNIT=52,FILE='fn20xyh10t2L5p1tstzphi.out',STATUS='unknown')
	  OPEN (UNIT=23,FILE='fn30xyh10t2L5p1tstsphi.out',STATUS='unknown')
      OPEN (UNIT=33,FILE='fn30xyh10t2L5p1tstxphi.out',STATUS='unknown')
	  OPEN (UNIT=43,FILE='fn30xyh10t2L5p1tstyphi.out',STATUS='unknown')
	  OPEN (UNIT=53,FILE='fn30xyh10t2L5p1tstzphi.out',STATUS='unknown')
	  OPEN (UNIT=24,FILE='fn40xyh10t2L5p1tstsphi.out',STATUS='unknown')
      OPEN (UNIT=34,FILE='fn40xyh10t2L5p1tstxphi.out',STATUS='unknown')
	  OPEN (UNIT=44,FILE='fn40xyh10t2L5p1tstyphi.out',STATUS='unknown')
	  OPEN (UNIT=54,FILE='fn40xyh10t2L5p1tstzphi.out',STATUS='unknown')
	  OPEN (UNIT=25,FILE='fn50xyh10t2L5p1tstsphi.out',STATUS='unknown')
      OPEN (UNIT=35,FILE='fn50xyh10t2L5p1tstxphi.out',STATUS='unknown')
	  OPEN (UNIT=45,FILE='fn50xyh10t2L5p1tstyphi.out',STATUS='unknown')
	  OPEN (UNIT=55,FILE='fn50xyh10t2L5p1tstzphi.out',STATUS='unknown')
	  OPEN (UNIT=61,FILE='fn10xyh10t2L5p1qpx4phi.out',STATUS='unknown')
	  OPEN (UNIT=62,FILE='fn20xyh10t2L5p1qpx4phi.out',STATUS='unknown')
      OPEN (UNIT=63,FILE='fn30xyh10t2L5p1qpx4phi.out',STATUS='unknown')
	  OPEN (UNIT=64,FILE='fn40xyh10t2L5p1qpx4phi.out',STATUS='unknown')
	  OPEN (UNIT=65,FILE='fn50xyh10t2L5p1qpx4phi.out',STATUS='unknown')
!
      OPEN (UNIT=20,FILE='maxfn10xyh10t2L5p1tsts.out',STATUS='unknown')
      OPEN (UNIT=30,FILE='maxfn10xyh10t2L5p1tstx.out',STATUS='unknown')
	  OPEN (UNIT=40,FILE='maxfn10xyh10t2L5p1tsty.out',STATUS='unknown')
	  OPEN (UNIT=50,FILE='maxfn10xyh10t2L5p1tstz.out',STATUS='unknown')
      OPEN (UNIT=26,FILE='maxfn20xyh10t2L5p1tsts.out',STATUS='unknown')
      OPEN (UNIT=36,FILE='maxfn20xyh10t2L5p1tstx.out',STATUS='unknown')
      OPEN (UNIT=46,FILE='maxfn20xyh10t2L5p1tsty.out',STATUS='unknown')
	  OPEN (UNIT=56,FILE='maxfn20xyh10t2L5p1tstz.out',STATUS='unknown')
      OPEN (UNIT=27,FILE='maxfn30xyh10t2L5p1tsts.out',STATUS='unknown')
      OPEN (UNIT=37,FILE='maxfn30xyh10t2L5p1tstx.out',STATUS='unknown')
	  OPEN (UNIT=47,FILE='maxfn30xyh10t2L5p1tsty.out',STATUS='unknown')
	  OPEN (UNIT=57,FILE='maxfn30xyh10t2L5p1tstz.out',STATUS='unknown')
	  OPEN (UNIT=28,FILE='maxfn40xyh10t2L5p1tsts.out',STATUS='unknown')
      OPEN (UNIT=38,FILE='maxfn40xyh10t2L5p1tstx.out',STATUS='unknown')
	  OPEN (UNIT=48,FILE='maxfn40xyh10t2L5p1tsty.out',STATUS='unknown')
	  OPEN (UNIT=58,FILE='maxfn40xyh10t2L5p1tstz.out',STATUS='unknown')
	  OPEN (UNIT=29,FILE='maxfn50xyh10t2L5p1tsts.out',STATUS='unknown')
      OPEN (UNIT=39,FILE='maxfn50xyh10t2L5p1tstx.out',STATUS='unknown')
	  OPEN (UNIT=49,FILE='maxfn50xyh10t2L5p1tsty.out',STATUS='unknown')
	  OPEN (UNIT=59,FILE='maxfn50xyh10t2L5p1tstz.out',STATUS='unknown')

      WRITE(10,*) 'TcvEth, is:'
	  pi=3.14159D0
      N3=5

      DO 1000 i=1,N1,1
	    DO 20 j=1,N1,1
	     IF (i.NE.j) THEN
	      one(i,j)=(0.0,0.0)
	     ELSE
  	      one(i,j)=(1.0,0.0)
	     ENDIF
20     CONTINUE
1000  CONTINUE

	  tau3sigma0(1,1)=(1.0,0.0)
	  tau3sigma0(1,2)=(0.0,0.0)
      tau3sigma0(1,3)=(0.0,0.0)
	  tau3sigma0(1,4)=(0.0,0.0)

      tau3sigma0(2,1)=(0.0,0.0)
 	  tau3sigma0(2,2)=(1.0,0.0)
      tau3sigma0(2,3)=(0.0,0.0)
	  tau3sigma0(2,4)=(0.0,0.0)
!
1     tau3sigma0(3,1)=(0.0,0.0)
	  tau3sigma0(3,2)=(0.0,0.0)
      tau3sigma0(3,3)=(-1.0,0.0)
      tau3sigma0(3,4)=(0.0,0.0)
!
      tau3sigma0(4,1)=(0.0,0.0)
  	  tau3sigma0(4,2)=(0.0,0.0)
      tau3sigma0(4,3)=(0.0,0.0)
	  tau3sigma0(4,4)=(-1.0,0.0)
!
	  tau0sigma0(1,1)=(1.0,0.0)
	  tau0sigma0(1,2)=(0.0,0.0)
      tau0sigma0(1,3)=(0.0,0.0)
	  tau0sigma0(1,4)=(0.0,0.0)
!
      tau0sigma0(2,1)=(0.0,0.0)
	  tau0sigma0(2,2)=(1.0,0.0)
      tau0sigma0(2,3)=(0.0,0.0)
	  tau0sigma0(2,4)=(0.0,0.0)
!
      tau0sigma0(3,1)=(0.0,0.0)
	  tau0sigma0(3,2)=(0.0,0.0)
      tau0sigma0(3,3)=(1.0,0.0)
      tau0sigma0(3,4)=(0.0,0.0)
!
      tau0sigma0(4,1)=(0.0,0.0)
 	  tau0sigma0(4,2)=(0.0,0.0)
      tau0sigma0(4,3)=(0.0,0.0)
	  tau0sigma0(4,4)=(1.0,0.0)
!
  	  tau3sigma1(1,1)=(0.0,0.0)
	  tau3sigma1(1,2)=(1.0,0.0)
      tau3sigma1(1,3)=(0.0,0.0)
	  tau3sigma1(1,4)=(0.0,0.0)
!
      tau3sigma1(2,1)=(1.0,0.0)
	  tau3sigma1(2,2)=(0.0,0.0)
      tau3sigma1(2,3)=(0.0,0.0)
	  tau3sigma1(2,4)=(0.0,0.0)
!
      tau3sigma1(3,1)=(0.0,0.0)
	  tau3sigma1(3,2)=(0.0,0.0)
      tau3sigma1(3,3)=(0.0,0.0)
	  tau3sigma1(3,4)=(-1.0,0.0)
!
      tau3sigma1(4,1)=(0.0,0.0)
	  tau3sigma1(4,2)=(0.0,0.0)
      tau3sigma1(4,3)=(-1.0,0.0)
	  tau3sigma1(4,4)=(0.0,0.0)
!
	  tau3sigma2(1,1)=(0.0,0.0)
	  tau3sigma2(1,2)=(0.0,-1.0)
      tau3sigma2(1,3)=(0.0,0.0)
	  tau3sigma2(1,4)=(0.0,0.0)
!
      tau3sigma2(2,1)=(0.0,1.0)
      tau3sigma2(2,2)=(0.0,0.0)
      tau3sigma2(2,3)=(0.0,0.0)
	  tau3sigma2(2,4)=(0.0,0.0)
!
      tau3sigma2(3,1)=(0.0,0.0)
	  tau3sigma2(3,2)=(0.0,0.0)
      tau3sigma2(3,3)=(0.0,0.0)
	  tau3sigma2(3,4)=(0.0,1.0)
!
      tau3sigma2(4,1)=(0.0,0.0)
	  tau3sigma2(4,2)=(0.0,0.0)
      tau3sigma2(4,3)=(0.0,-1.0)
	  tau3sigma2(4,4)=(0.0,0.0)
!
  	  tau3sigma3(1,1)=(1.0,0.0)
	  tau3sigma3(1,2)=(0.0,0.0)
      tau3sigma3(1,3)=(0.0,0.0)
	  tau3sigma3(1,4)=(0.0,0.0)
!
      tau3sigma3(2,1)=(0.0,0.0)
	  tau3sigma3(2,2)=(-1.0,0.0)
      tau3sigma3(2,3)=(0.0,0.0)
	  tau3sigma3(2,4)=(0.0,0.0)
!
      tau3sigma3(3,1)=(0.0,0.0)
	  tau3sigma3(3,2)=(0.0,0.0)
      tau3sigma3(3,3)=(-1.0,0.0)
	  tau3sigma3(3,4)=(0.0,0.0)
!
      tau3sigma3(4,1)=(0.0,0.0)
	  tau3sigma3(4,2)=(0.0,0.0)
      tau3sigma3(4,3)=(0.0,0.0)
	  tau3sigma3(4,4)=(1.0,0.0)
!
 	  tau0sigma3(1,1)=(1.0,0.0)
	  tau0sigma3(1,2)=(0.0,0.0)
      tau0sigma3(1,3)=(0.0,0.0)
	  tau0sigma3(1,4)=(0.0,0.0)
!
      tau0sigma3(2,1)=(0.0,0.0)
	  tau0sigma3(2,2)=(-1.0,0.0)
      tau0sigma3(2,3)=(0.0,0.0)
	  tau0sigma3(2,4)=(0.0,0.0)
!
      tau0sigma3(3,1)=(0.0,0.0)
	  tau0sigma3(3,2)=(0.0,0.0)
      tau0sigma3(3,3)=(1.0,0.0)
	  tau0sigma3(3,4)=(0.0,0.0)
!
      tau0sigma3(4,1)=(0.0,0.0)
	  tau0sigma3(4,2)=(0.0,0.0)
      tau0sigma3(4,3)=(0.0,0.0)
	  tau0sigma3(4,4)=(-1.0,0.0)
!
	  hg=(0.0,1.0)
!
	  DO 30  kjj=1,1,1
	   DO 50 dd=1,1,1
	    t(1)=0.1*dd
	    t(2)=0.1*dd
	    zx=1
		write(81,*)'Temperature is:', 0.1D0*dd
60	    IF (zx.EQ.1)THEN
	     DO 40 tt=1,41,1
          write(81,*)'phase difference counter is', tt
		  write(56,*)'phi is',1.0D0*tt*pi/4.0D0
  	      Qxi1=pi
 	      Qxi2=1.0D0*pi
	      hvTc=10.0D0
	      lvxiF1=0.1D0
	      lvxiN=0.5D0*kjj
	      lvxiF2=0.1D0
	      TcvEthF1=(lvxiF1*lvxiF1)/(2.0D0*1.76D0)
	      TcvEthN=(lvxiN*lvxiN)/(2.0D0*1.76D0)
	      TcvEthF2=(lvxiF2*lvxiF2)/(2.0D0*1.76D0)
          write(19,*)'TcvEth is',TcvEthN,' hvTc is ',hvTC
	      DO 23 i=1,N1,1
           DO 33 j=1,N1,1
            Gdecn(i,j,1)=(0.0,0.0)
	        Gdecn(i,j,N3)=(0.0,0.0)
33 	       CONTINUE
23	      CONTINUE

!         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         PART 2
!         This part will find the number of nodes which converges
!         This loop will continue till the n which satisfies the
!        		
	      DO 120 qqq=1,N1,1
           DO 130  qq=1,N1,1
            DO 140  q=1,N3,1
	         errGkn(qqq,qq,q)=0.0000001D0
140         CONTINUE
130        CONTINUE
120       CONTINUE   
	      DO 1111 q=1,N3,1
           errIkn(q)=0.00001D0
	       errIknsx(q)=0.00001D0
	       errIknsx4(q)=0.00001D0
           errIknsy(q)=0.00001D0
           errIknsz(q)=0.00001D0
1111      CONTINUE	
!
!         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         PART
!
          Gn=1.0D0
          RelConductance(1)=1.0D0*(N3-1)
          RelConductance(N3-1)=1.0D0*(N3-1)
!         Here, number of nodes in F_c is N3-1
          DO 80 i=2,N3-2,1
           RelConductance(i)=1.0D0*(N3-1)
80        CONTINUE

!         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         PART 3
!
          nphi=21	
          DO 800 kj=1,41,1
           phi=((tt-1)*pi)/(nphi-1)
           k=0
	       DO 99999 ii=1,N3,1
	        In(ii)=0.0D0
		    In2(ii)=0.0D0
		    In3(ii)=0.0D0
		    In4(ii)=0.0D0
	        Insx(ii)=0.0D0
	        Insy(ii)=0.0D0
	        Insz(ii)=0.0D0
	        Insx4(ii)=0.0D0
99999      CONTINUE
!          Initial guess of Green's functions of nodes in Normal and ferromagnetic layers
           DO 90 q=2,N3-1,1
!
            Gkn(1,1,q)=(1.0,0.0)
            Gkn(1,2,q)=(0.0,0.0)
            Gkn(1,3,q)=(0.0,0.0)
            Gkn(1,4,q)=(0.0,0.0)
!
            Gkn(2,1,q)=(0.0,0.0)
            Gkn(2,2,q)=(1.0,0.0)
            Gkn(2,3,q)=(0.0,0.0)
            Gkn(2,4,q)=(0.0,0.0)
!
            Gkn(3,1,q)=(0.0,0.0)
            Gkn(3,2,q)=(0.0,0.0)
            Gkn(3,3,q)=(-1.0,0.0)
            Gkn(3,4,q)=(0.0,0.0)
!
            Gkn(4,1,q)=(0.0,0.0)
            Gkn(4,2,q)=(0.0,0.0)
            Gkn(4,3,q)=(0.0,0.0)
            Gkn(4,4,q)=(-1.0,0.0)
!
90         CONTINUE
!
           rr=1
!          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          PART 4
!          This part will iterate till the convergence of current 
!          for specIFied phi occurs
!
100        IF (rr.EQ.1) THEN
            f=0
            k=k+1	     
            alpha(k)=((2.0D0*(k)-1.0D0)*pi)/(1.76D0)
            DO  110 h=1,2,1
             TH(h)=(tanh(1.74D0*(SQRT((1.0D0/t(h))-1.0D0))))
             WkvDelta(h)=(alpha(k)*t(h))/TH(h)
110         CONTINUE
            cte=(1.0D0/(SQRT((WkvDelta(1)*WkvDelta(1))+1.0D0)))
            Gkn(1,1,1)=cte*WkvDelta(1)
	        Gkn(1,2,1)=(0.0,0.0)
	        Gkn(1,4,1)=(0.0,0.0)
	        Gkn(1,3,1)=-hg*cte
	        Gkn(2,1,1)=(0.0,0.0)
	        Gkn(2,2,1)=cte*WkvDelta(1)
	        Gkn(2,4,1)=hg*cte
	        Gkn(2,3,1)=(0.0,0.0)
	        Gkn(3,2,1)=(0.0,0.0)
	        Gkn(3,1,1)=hg*cte
	        Gkn(3,3,1)=-cte*WkvDelta(1)
	        Gkn(3,4,1)=(0.0,0.0)
	        Gkn(4,2,1)=-hg*cte
	        Gkn(4,1,1)=(0.0,0.0)
	        Gkn(4,3,1)=(0.0,0.0)
	        Gkn(4,4,1)=-cte*WkvDelta(1)
!
            if ((tt.EQ.1).OR.(tt.EQ.21)) then
	         Gkn(1,1,N3)=cte*WkvDelta(1)
             Gkn(1,2,N3)=(0.0,0.0)
             Gkn(1,4,N3)=(0.0,0.0)
             Gkn(1,3,N3)=-hg*cte
             Gkn(2,1,N3)=(0.0,0.0)
             Gkn(2,2,N3)=cte*WkvDelta(1)
             Gkn(2,4,N3)=hg*cte
             Gkn(2,3,N3)=(0.0,0.0)
             Gkn(3,2,N3)=(0.0,0.0)
             Gkn(3,1,N3)=hg*cte
             Gkn(3,3,N3)=-cte*WkvDelta(1)
             Gkn(3,4,N3)=(0.0,0.0)
             Gkn(4,2,N3)=-hg*cte
             Gkn(4,1,N3)=(0.0,0.0)
             Gkn(4,3,N3)=(0.0,0.0)
             Gkn(4,4,N3)=-cte*WkvDelta(1)
            elseif ((tt.EQ.11)) then
             Gkn(1,1,N3)=cte*WkvDelta(1)
             Gkn(1,2,N3)=(0.0,0.0)
             Gkn(1,4,N3)=(0.0,0.0)
             Gkn(1,3,N3)=hg*cte
             Gkn(2,1,N3)=(0.0,0.0)
             Gkn(2,2,N3)=cte*WkvDelta(1)
             Gkn(2,4,N3)=-hg*cte
             Gkn(2,3,N3)=(0.0,0.0)
             Gkn(3,2,N3)=(0.0,0.0)
             Gkn(3,1,N3)=-hg*cte
             Gkn(3,3,N3)=-cte*WkvDelta(1)
             Gkn(3,4,N3)=(0.0,0.0)
             Gkn(4,2,N3)=hg*cte
             Gkn(4,1,N3)=(0.0,0.0)
             Gkn(4,3,N3)=(0.0,0.0)
             Gkn(4,4,N3)=-cte*WkvDelta(1)
            else
	         Gkn(1,1,N3)=cte*WkvDelta(1)
	         Gkn(1,2,N3)=(0.0,0.0)
	         Gkn(1,4,N3)=(0.0,0.0)
	         Gkn(1,3,N3)=-hg*cte*exp(hg*phi)
             Gkn(2,1,N3)=(0.0,0.0)
             Gkn(2,2,N3)=cte*WkvDelta(1)
             Gkn(2,4,N3)=hg*cte*exp(hg*phi)
             Gkn(2,3,N3)=(0.0,0.0)
             Gkn(3,2,N3)=(0.0,0.0)
             Gkn(3,1,N3)=hg*cte*exp(-hg*phi)
             Gkn(3,3,N3)=-cte*WkvDelta(1)
             Gkn(3,4,N3)=(0.0,0.0)
             Gkn(4,2,N3)=-hg*cte*exp(-hg*phi)
             Gkn(4,1,N3)=(0.0,0.0)
             Gkn(4,3,N3)=(0.0,0.0)
             Gkn(4,4,N3)=-cte*WkvDelta(1)
            end if
!
             lkgncte1=-((2.0D0*1.76D0)/(N3-2))*(TcvEthN)
             lkg=lkgncte1*(-alpha(k)*t(1))
             lkg2=lkgncte1*(alpha(k)*t(1))
!
             Gdecn1(1,1)=CMPLX(lkg,0.0)
             Gdecn1(2,2)=CMPLX(lkg,0.0)
             Gdecn1(3,3)=CMPLX(lkg2,0.0)
             Gdecn1(4,4)=CMPLX(lkg2,0.0)
             Gdecn1(1,2)=(0.0,0.0)
             Gdecn1(1,3)=(0.0,0.0)
             Gdecn1(1,4)=(0.0,0.0)
             Gdecn1(2,1)=(0.0,0.0)
             Gdecn1(2,3)=(0.0,0.0)
             Gdecn1(2,4)=(0.0,0.0)
             Gdecn1(3,1)=(0.0,0.0)
             Gdecn1(3,2)=(0.0,0.0)
             Gdecn1(3,4)=(0.0,0.0)
             Gdecn1(4,1)=(0.0,0.0)
             Gdecn1(4,2)=(0.0,0.0)
             Gdecn1(4,3)=(0.0,0.0)

	         DO 2100 s=2,N3-1,1
	          if (s.EQ.2) then
               gamma2(s)=0.0D0
               nsigma(1,1,s)=(0.0,0.0)
               nsigma(1,2,s)=CMPLX(sin(gamma2(s)),-cos(gamma2(s)))
               nsigma(1,3,s)=(0.0,0.0)
               nsigma(1,4,s)=(0.0,0.0)
               nsigma(2,1,s)=CMPLX(sin(gamma2(s)),cos(gamma2(s)))
               nsigma(2,2,s)=(0.0,0.0)
               nsigma(2,3,s)=(0.0,0.0)
               nsigma(2,4,s)=(0.0,0.0)
               nsigma(3,1,s)=(0.0,0.0)
               nsigma(3,2,s)=(0.0,0.0)
               nsigma(3,3,s)=(0.0,0.0)
               nsigma(3,4,s)=CMPLX(sin(gamma2(s)),-cos(gamma2(s)))
               nsigma(4,1,s)=(0.0,0.0)
               nsigma(4,2,s)=(0.0,0.0)
               nsigma(4,3,s)=CMPLX(sin(gamma2(s)),cos(gamma2(s)))
               nsigma(4,4,s)=(0.0,0.0)

              elseif (s.EQ.N3-1) then
                
               gamma2(s)=(1.0D0*(kj-1)*pi/20.0D0)
               nsigma(1,1,s)=(0.0,0.0)
               nsigma(1,2,s)=CMPLX(sin(gamma2(s)),-cos(gamma2(s)))
               nsigma(1,3,s)=(0.0,0.0)
               nsigma(1,4,s)=(0.0,0.0)
               nsigma(2,1,s)=CMPLX(sin(gamma2(s)),cos(gamma2(s)))
               nsigma(2,2,s)=(0.0,0.0)
               nsigma(2,3,s)=(0.0,0.0)
               nsigma(2,4,s)=(0.0,0.0)
               nsigma(3,1,s)=(0.0,0.0)
               nsigma(3,2,s)=(0.0,0.0)
               nsigma(3,3,s)=(0.0,0.0)
               nsigma(3,4,s)=CMPLX(sin(gamma2(s)),-cos(gamma2(s)))
               nsigma(4,1,s)=(0.0,0.0)
               nsigma(4,2,s)=(0.0,0.0)
               nsigma(4,3,s)=CMPLX(sin(gamma2(s)),cos(gamma2(s)))
               nsigma(4,4,s)=(0.0,0.0)

              else
	 
	           nsigma(1,1,s)=(0.0,0.0)
               nsigma(1,2,s)=(0.0,0.0)
               nsigma(1,3,s)=(0.0,0.0)
               nsigma(1,4,s)=(0.0,0.0)
               nsigma(2,1,s)=(0.0,0.0)
               nsigma(2,2,s)=(0.0,0.0)
               nsigma(2,3,s)=(0.0,0.0)
               nsigma(2,4,s)=(0.0,0.0)
               nsigma(3,1,s)=(0.0,0.0)
               nsigma(3,2,s)=(0.0,0.0)
               nsigma(3,3,s)=(0.0,0.0)
               nsigma(3,4,s)=(0.0,0.0)
               nsigma(4,1,s)=(0.0,0.0)
               nsigma(4,2,s)=(0.0,0.0)
               nsigma(4,3,s)=(0.0,0.0)
               nsigma(4,4,s)=(0.0,0.0)

              endif

            DO 190 i=1,N1,1
             DO 200 j=1,N1,1
			  Gdecn2(i,j,s)=lkgncte2*nsigma(i,j,s)
200		     CONTINUE
190		    CONTINUE
!C                                                             
            DO 210 i=1,N1,1
             DO 220 j=1,N1,1
 		      Gdecn(i,j,s)=Gdecn1(i,j)+Gdecn2(i,j,s)
220		     CONTINUE
210		    CONTINUE
2100   CONTINUE
  
	   DO 270 ii=1,N1,1
	    DO 280 jj=1,N1,1
	     DO 290 s=1,N3,1
		   Gdcn2(ii,jj,s)=(0.0,0.0)
		   DO 310 kkk=1,N1,1
             Gdcn2(ii,jj,s)=Gdecn(ii,kkk,s)*Gdecn(kkk,jj,s)+Gdcn2(ii,jj,s)
310        CONTINUE
290	     CONTINUE
280     CONTINUE
270    CONTINUE

       d=1
150    IF (d.EQ.1)THEN
	    DO 160 kk=1,N1,1
         DO 170 jj=1,N1,1
	 	  DO 180 s=2,N3-1,1
		   Mkn(kk,jj,s)=RelConductance(s-1)*Gkn(kk,jj,s-1)+RelConductance(s)*Gkn(kk,jj,s+1)+Gdecn(kk,jj,s)
           ANCOMMsmspkn1(kk,jj,s)=(0.0,0.0)
           ANCOMMsmspkn2(kk,jj,s)=(0.0,0.0)
           ANCOMMsmdeckn1(kk,jj,s)=(0.0,0.0)
           ANCOMMsmdeckn2(kk,jj,s)=(0.0,0.0)
           ANCOMMdecspkn1(kk,jj,s)=(0.0,0.0)
           ANCOMMdecspkn2(kk,jj,s)=(0.0,0.0)
		   DO 230 i=1,N1,1
            ANCOMMsmspkn1(kk,jj,s)=Gkn(kk,i,s+1)*Gkn(i,jj,s-1)+ANCOMMsmspkn1(kk,jj,s)
            ANCOMMsmspkn2(kk,jj,s)=Gkn(kk,i,s-1)*Gkn(i,jj,s+1)+ANCOMMsmspkn2(kk,jj,s)
            ANCOMMsmdeckn1(kk,jj,s)=Gkn(kk,i,s-1)*Gdecn(i,jj,s)+ANCOMMsmdeckn1(kk,jj,s)
            ANCOMMsmdeckn2(kk,jj,s)=Gdecn(kk,i,s)*Gkn(i,jj,s-1)+ANCOMMsmdeckn2(kk,jj,s)
            ANCOMMdecspkn1(kk,jj,s)=Gkn(kk,i,s+1)*Gdecn(i,jj,s)+ANCOMMdecspkn1(kk,jj,s)
            ANCOMMdecspkn2(kk,jj,s)=Gdecn(kk,i,s)*Gkn(i,jj,s+1)+ANCOMMdecspkn2(kk,jj,s)
230		   CONTINUE
		   ANCOMMsmdeckn(kk,jj,s)=ANCOMMsmdeckn1(kk,jj,s)+ANCOMMsmdeckn2(kk,jj,s)
		   ANCOMMdecspkn(kk,jj,s)=ANCOMMdecspkn1(kk,jj,s)+ANCOMMdecspkn2(kk,jj,s)
		   ANCOMMsmspkn(kk,jj,s)=ANCOMMsmspkn1(kk,jj,s)+ANCOMMsmspkn2(kk,jj,s)
180	      CONTINUE
170      CONTINUE
160     CONTINUE

        DO 240 ii=1,N1,1
         DO 250 jj=1,N1,1
          DO 260 s=2,N3-1,1
          Dkn(ii,jj,s)=RelConductance(s-1)*RelConductance(s-1)*one(ii,jj)+&
RelConductance(s)*RelConductance(s)*one(ii,jj)+Gdcn2(ii,jj,s)+&
RelConductance(s)*RelConductance(s-1)*ANCOMMsmspkn(ii,jj,s)+RelConductance(s-1)*&
ANCOMMsmdeckn(ii,jj,s)+RelConductance(s)*&
ANCOMMdecspkn(ii,jj,s)
260       CONTINUE
250      CONTINUE
240     CONTINUE

!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!       This part will calculate eigenvalues and eigenvector
!       P is the matrix of eigenvectors and landa is the mat
!
	    DO 350 s=2,N3-1,1
         DO 360 ii=1,N1,1
          DO 370 jj=1,N1,1
           AD(ii,jj)=Dkn(ii,jj,s)
370       CONTINUE
360      CONTINUE
         CALL ZGEEV('N','V',N1,AD,N1,W,VL,1,VR,N1,WORK,N1*2,RWORK,INFO)
         DO 390 ii=1,N1,1
          DO 400 jj=1,N1,1
           IF (ii.NE.jj) THEN
		    landa(ii,jj,s)=(0.0,0.0)
           ELSE
		    landa(ii,jj,s)=W(ii)
           ENDIF
400	      CONTINUE
390 	 CONTINUE
!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!	     This part will calculate the inverse
!        matrix of P. P is the matrix of eigenvectors
!
         DO  410 ii=1,N1,1
          DO 420 kk=1,N1,1
		   EIGs(ii,kk,s)=VR(ii,kk)
		   inEIGs(ii,kk)=VR(ii,kk)
420       CONTINUE
410	     CONTINUE
!
         CALL ZGETRF( N1, N1, inEIGs, N1, IPIV, INFO )
	     CALL ZGETRI( N1, inEIGs, N1, IPIV, WORK2, LWORK, INFO )

!          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          D^(d)=P^(-1)*landa^(d)*P###########################
!
		 DO 570 ii=1,N1,1
		  DO 5700 jj=1,N1,1
		   inEIGs2(ii,jj,s)=inEIGs(ii,jj)
!
		   IF (ii.NE.jj) THEN
 		    invlanda(ii,jj,s)=(0.0,0.0)
		   ELSE
		    invlanda(ii,jj,s)=(1.0,0.0)/(SQRT(landa(ii,jj,s)))
		   ENDIF
5700 	  CONTINUE
570	     CONTINUE
         DO 590 ii=1,N1,1
          DO 600 jj=1,N1,1
           Dkn1v2(ii,jj,s)=(0.0,0.0)
           DO 610 ll=1,N1,1
            DO 620 l=1,N1,1
             Dkn1v2(ii,jj,s)=EIGs(ii,ll,s)*(invlanda(ll,l,s))*inEIGs2(l,jj,s)+Dkn1v2(ii,jj,s)
620 		CONTINUE
610        CONTINUE
600		  CONTINUE
590      CONTINUE
!
	     DO 630 ii=1,N1,1
          DO 640 jj=1,N1,1
           Gnewkn(ii,jj,s) =(0.0,0.0)
           DO 660 kk=1,N1,1
            Gnewkn(ii,jj,s)=Mkn(ii,kk,s)*Dkn1v2(kk,jj,s)+Gnewkn(ii,jj,s)
660        CONTINUE
           tolGkn(ii,jj,s)=Gnewkn(ii,jj,s)-Gkn(ii,jj,s)
640       CONTINUE
630      CONTINUE

350     CONTINUE

!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        f=f+1
! 	    Here, improved Green's function is replaced with the
!        old one
        DO 670 ii=1,N1,1
         DO 680 jj=1,N1,1
          DO 690 s=2,N3-1,1
	       Gkn(ii,jj,s)=Gnewkn(ii,jj,s)
690		  CONTINUE
680		 CONTINUE
670	    CONTINUE
        d=0
!      This part will check IF the tolerance is smaller
!	   than errG
        DO 700 ii=1,N1,1
         DO 710 jj=1,N1,1
          DO 720 s=2,N3-1,1
		   IF  ((abs(real(tolGkn(ii,jj,s))).GT.errGkn(ii,jj,s)).OR.&
(abs(AIMAG(tolGkn(ii,jj,s))).GT.errGkn(ii,jj,s))) THEN
             d=1
           ENDIF
720		  CONTINUE
710	     CONTINUE
700	    CONTINUE
       GOTO 150
	  ENDIF

!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	  DO  830 ii=1,N1,1
       DO 840 jj=1,N1,1
        DO 850 s=1,N3-1,1
         Imxkn(ii,jj,s)=(0.0,0.0)
         DO 860 l=1,N1,1
	      Imxkn(ii,jj,s)=(Gkn(ii,l,s)*Gkn(l,jj,s+1)-Gkn(ii,l,s+1)*&
Gkn(l,jj,s))*((exp((hg*(2*k-1)*pi)/2)-exp((-hg*(2*k-1)*pi)/2))/&
(exp((hg*(2*k-1)*pi)/2)+exp((-hg*(2*k-1)*pi)/2)))+Imxkn(ii,jj,s)
860      CONTINUE
850     CONTINUE
840    CONTINUE
830   CONTINUE
	  DO  8301 ii=1,N1,1
       DO 8401 jj=1,N1,1
        DO 8501 s=1,N3,1
         Imxfkn(ii,jj,s)=(0.0,0.0)
         DO 8601 l=1,N1,1
	      Imxfkn(ii,jj,s)=(Gkn(ii,l,s)*Gdecn2(l,jj,s)-Gdecn2(ii,l,s)*&
Gkn(l,jj,s))*((exp((hg*(2*k-1)*pi)/2)-exp((-hg*(2*k-1)*pi)/2))/&
(exp((hg*(2*k-1)*pi)/2)+exp((-hg*(2*k-1)*pi)/2)))+Imxfkn(ii,jj,s)
8601     CONTINUE
8501    CONTINUE
8401   CONTINUE
8301  CONTINUE

      DO 870 ii=1,N1,1
       DO 880 jj=1,N1,1
        DO 890 s=1,N3-1,1
         SgImkn(ii,jj,s)=(0.0,0.0)
         DO 900 l=1,N1,1
          SgImkn(ii,jj,s)=hg*tau3sigma0(ii,l)*Imxkn(l,jj,s)+SgImkn(ii,jj,s)
900      CONTINUE
890     CONTINUE
880	   CONTINUE
870   CONTINUE

      DO 8170 ii=1,N1,1
       DO 8180 jj=1,N1,1
        DO 8190 s=1,N3-1,1
         SgImkn2(ii,jj,s)=(0.0,0.0)
         DO 9100 l=1,N1,1
          SgImkn2(ii,jj,s)=hg*tau3sigma1(ii,l)*Imxkn(l,jj,s)+SgImkn2(ii,jj,s)
9100     CONTINUE
8190    CONTINUE
8180   CONTINUE
8170  CONTINUE

      DO 8270 ii=1,N1,1
       DO 8280 jj=1,N1,1
        DO 8290 s=1,N3-1,1
         SgImkn3(ii,jj,s)=(0.0,0.0)
         DO 9200 l=1,N1,1
          SgImkn3(ii,jj,s)=hg*tau3sigma2(ii,l)*Imxkn(l,jj,s)+SgImkn3(ii,jj,s)
9200     CONTINUE
8290    CONTINUE
8280   CONTINUE
8270  CONTINUE

      DO 8370 ii=1,N1,1
       DO 8380 jj=1,N1,1
        DO 8390 s=1,N3-1,1
         SgImkn4(ii,jj,s)=(0.0,0.0)
         DO 9300 l=1,N1,1
          SgImkn4(ii,jj,s)=hg*tau0sigma3(ii,l)*Imxkn(l,jj,s)+SgImkn4(ii,jj,s)
9300     CONTINUE
8390    CONTINUE
8380   CONTINUE
8370  CONTINUE

      DO 8703 ii=1,N1,1
       DO 8803 jj=1,N1,1
        DO 8903 s=1,N3,1
         SgxImkn(ii,jj,s)=(0.0,0.0)
         DO 9003 l=1,N1,1
          SgxImkn(ii,jj,s)=hg*tau3sigma1(ii,l)*Imxfkn(l,jj,s+1)+SgxImkn(ii,jj,s)
9003     CONTINUE
8903    CONTINUE
8803   CONTINUE
8703  CONTINUE

      DO 8704 ii=1,N1,1
       DO 8804 jj=1,N1,1
        DO 8904 s=1,N3,1
         SgyImkn(ii,jj,s)=(0.0,0.0)
         DO 9004 l=1,N1,1
          SgyImkn(ii,jj,s)=hg*tau3sigma2(ii,l)*Imxfkn(l,jj,s)+SgyImkn(ii,jj,s)
9004     CONTINUE
8904    CONTINUE
8804   CONTINUE
8704  CONTINUE

      DO 8705 ii=1,N1,1
       DO 8805 jj=1,N1,1
        DO 8905 s=1,N3,1
         SgzzImkn(ii,jj,s)=(0.0,0.0)
         DO 9005 l=1,N1,1
          SgzzImkn(ii,jj,s)=hg*tau0sigma3(ii,l)*Imxfkn(l,jj,s)+SgzzImkn(ii,jj,s)
9005     CONTINUE
8905    CONTINUE
8805   CONTINUE
8705        CONTINUE

       DO 8706 ii=1,N1,1
        DO 8806 jj=1,N1,1
         DO 8906 s=1,N3,1
          Sgx4Imkn(ii,jj,s)=(0.0,0.0)
          DO 9006 l=1,N1,1
           Sgx4Imkn(ii,jj,s)=hg*tau3sigma3(ii,l)*Imxfkn(l,jj,s)+Sgx4Imkn(ii,jj,s)
9006      CONTINUE
8906     CONTINUE
8806    CONTINUE
8706   CONTINUE
       DO 920 s=1,N3-1,1
        TriSgImkn(s)=(0.0,0.0)
920    CONTINUE
       DO 921 s=1,N3-1,1
        TriSgImkn2(s)=(0.0,0.0)
921    CONTINUE
       DO 922 s=1,N3-1,1
        TriSgImkn3(s)=(0.0,0.0)
922    CONTINUE
       DO 923 s=1,N3-1,1
        TriSgImkn4(s)=(0.0,0.0)
923    CONTINUE
	   DO 9203 s=1,N3,1
        TriSgxImkn(s)=(0.0,0.0)
9203   CONTINUE
	   DO 9204 s=1,N3,1
        TriSgyImkn(s)=(0.0,0.0)
9204   CONTINUE
	   DO 9207 s=1,N3,1
        TriSgx4Imkn(s)=(0.0,0.0)
9207   CONTINUE
	   DO 9205 s=1,N3,1
        TriSgzzImkn(s)=(0.0,0.0)
9205   CONTINUE
       DO 930 ii=1,N1,1
        DO 940 jj=1,N1,1
         DO 950 s=1,N3-1,1
          IF (ii.EQ.jj) THEN
           TriSgImkn(s)=SgImkn(ii,jj,s)+TriSgImkn(s)
          ENDIF
950      CONTINUE
940     CONTINUE
930    CONTINUE

       DO 932 ii=1,N1,1
        DO 942 jj=1,N1,1
         DO 952 s=1,N3-1,1
          IF (ii.EQ.jj) THEN
           TriSgImkn2(s)=SgImkn2(ii,jj,s)+TriSgImkn2(s)
          ENDIF
952      CONTINUE
942     CONTINUE
932    CONTINUE

       DO 933 ii=1,N1,1
        DO 943 jj=1,N1,1
         DO 953 s=1,N3-1,1
          IF (ii.EQ.jj) THEN
           TriSgImkn3(s)=SgImkn3(ii,jj,s)+TriSgImkn3(s)
          ENDIF
953      CONTINUE
943     CONTINUE
933    CONTINUE

       DO 934 ii=1,N1,1
        DO 944 jj=1,N1,1
         DO 954 s=1,N3-1,1
          IF (ii.EQ.jj) THEN
           TriSgImkn4(s)=SgImkn4(ii,jj,s)+TriSgImkn4(s)
          ENDIF
954      CONTINUE
944     CONTINUE
934    CONTINUE

       DO 9303 ii=1,N1,1
        DO 9403 jj=1,N1,1
         DO 9503 s=1,N3,1
          IF (ii.EQ.jj) THEN
           TriSgxImkn(s)=SgxImkn(ii,jj,s)+TriSgxImkn(s)
          ENDIF
9503     CONTINUE
9403    CONTINUE
9303   CONTINUE

       DO 9304 ii=1,N1,1
        DO 9404 jj=1,N1,1
         DO 9504 s=1,N3,1
          IF (ii.EQ.jj) THEN
           TriSgyImkn(s)=SgyImkn(ii,jj,s)+TriSgyImkn(s)
          ENDIF
9504     CONTINUE
9404    CONTINUE
9304   CONTINUE

       DO 9305 ii=1,N1,1
        DO 9405 jj=1,N1,1
         DO 9505 s=1,N3,1
          IF (ii.EQ.jj) THEN
           TriSgzzImkn(s)=SgzzImkn(ii,jj,s)+TriSgzzImkn(s)
          ENDIF
9505     CONTINUE
9405    CONTINUE
9305   CONTINUE

       DO 9306 ii=1,N1,1
        DO 9406 jj=1,N1,1
         DO 9506 s=1,N3,1
          IF (ii.EQ.jj) THEN
           TriSgx4Imkn(s)=Sgx4Imkn(ii,jj,s)+TriSgx4Imkn(s)
          ENDIF
9506     CONTINUE
9406    CONTINUE
9306   CONTINUE
       DO 960 s=1,N3-1
        Ikn(s)=AIMAG(((2.0D0*t(1)*RelConductance(s)*(1.0D0/(1.76D0*4.0D0))*(TriSgImkn(s)))))
        In(s)=Ikn(s)+In(s)
        IF (In(s).NE.0.0D0) THEN
         tolIn(s)=abs(Ikn(s)/In(s))
        ENDIF
960    CONTINUE
       DO 962 s=1,N3-1
        Ikn2(s)=AIMAG(((2.0D0*t(1)*RelConductance(s)*(1.0D0/(1.76D0*4.0D0))*(TriSgImkn2(s)))))
        In2(s)=Ikn2(s)+In2(s)
        IF (In2(s).NE.0.0D0) THEN
         tolIn2(s)=abs(Ikn2(s)/In2(s))
        ENDIF
962    CONTINUE
       DO 963 s=1,N3-1
        Ikn3(s)=AIMAG(((2.0D0*t(1)*RelConductance(s)*(1.0D0/(1.76D0*4.0D0))*(TriSgImkn3(s)))))
        In3(s)=Ikn3(s)+In3(s)
        IF (In3(s).NE.0.0D0) THEN
         tolIn3(s)=abs(Ikn3(s)/In3(s))
        ENDIF
963    CONTINUE
       DO 964 s=1,N3-1
        Ikn4(s)=AIMAG(((2.0D0*t(1)*RelConductance(s)*(1.0D0/(1.76D0*4.0D0))*(TriSgImkn4(s)))))
        In4(s)=Ikn4(s)+In4(s)
        IF (In4(s).NE.0.0D0) THEN
         tolIn4(s)=abs(Ikn4(s)/In4(s))
        ENDIF
964    CONTINUE

       rr=0
       DO 970 s=1,N3-1,1
        IF (tolIn(s).GT.errIkn(s)) THEN
         rr=1
        ENDIF
        IF (tolIn2(s).GT.errIkn(s)) THEN
         rr=1
        ENDIF
        IF (tolIn3(s).GT.errIkn(s)) THEN
         rr=1
        ENDIF
        IF (tolIn4(s).GT.errIkn(s)) THEN
         rr=1
        ENDIF
970    CONTINUE
       DO 9603 s=1,N3
        Iknsx(s)=aimag(((2.0D0*t(1)*(1.0D0/(1.76D0*4.0D0))*(TriSgxImkn(s)))))
         Insx(s)=Iknsx(s)+Insx(s)
         IF (Insx(s).NE.0.0D0) THEN
          tolInsx(s)=abs(Iknsx(s)/Insx(s))
         ENDIF
9603    CONTINUE

        DO 9604 s=1,N3
         Iknsy(s)=aimag(((2.0D0*t(1)*(1.0D0/(1.76D0*4.0D0))*(TriSgyImkn(s)))))
         Insy(s)=Iknsy(s)+Insy(s)
         IF (Insy(s).NE.0.0D0) THEN
          tolInsy(s)=abs(Iknsy(s)/Insy(s))
         ENDIF
9604    CONTINUE

        DO 9605 s=1,N3
         Iknsz(s)=aimag(((2.0D0*t(1)*(1.0D0/(1.76D0*4.0D0))*(TriSgzzImkn(s)))))
         Insz(s)=Iknsz(s)+Insz(s)
         IF (Insz(s).NE.0.0D0) THEN
          tolInsz(s)=abs(Iknsz(s)/Insz(s))
         ENDIF
9605    CONTINUE

		DO 9606 s=1,N3
         Iknsx4(s)=aimag(((2*t(1)*(1/(1.76*4))*(TriSgx4Imkn(s)))))
         Insx4(s)=Iknsx4(s)+Insx4(s)
         IF (Insx4(s).NE.0.0D0) THEN
          tolInsx4(s)=abs(Iknsx4(s)/Insx4(s))
         ENDIF
9606    CONTINUE
        GOTO 100
        ENDIF
        A10In(kj,tt)=In(2)
        A10In2(kj,tt)=In2(3)
!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        A20In(kj,tt)=In(4)
        A20In2(kj,tt)=In2(4)
        A20In3(kj,tt)=In3(2)
        A20In4(kj,tt)=In4(2)
        A20Insx(kj,tt)=Insx(2)
        A20Insy(kj,tt)=Insy(2)
        A20Insz(kj,tt)=Insz(4)
        A20Insx4(kj,tt)=Insx4(2)
        A10In3(kj,tt)=In3(3)
        A10In4(kj,tt)=In4(3)
        A10Insx(kj,tt)=Insx(3)
        A10Insy(kj,tt)=Insy(3)
        A10Insz(kj,tt)=Insz(3)
        A10Insx4(kj,tt)=Insx4(3)
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        write(10,*)'current is',In(1),In(2),In(3),In(4)
        write(21,*) In(2)/93000
        write(22,*) In(3)
        write(14,*)'current is',Insx(1),Insx(2),Insx(3),Insx(4)
    	write(31,*) Insx(2)
        write(32,*) Insx(3)
        write(15,*)'current is',Insy(1),Insy(2),Insy(3),Insy(4)
    	write(41,*) Insy(2)
        write(42,*) Insy(3)
        write(16,*)'current is',Insz(1),Insz(2),Insz(3),Insz(4)
        write(51,*) Insz(2)
        write(52,*) Insz(3)
        write(17,*)'current is',Insx4(1),Insx4(2),Insx4(3),Insx4(4)
        write(61,*) Insx4(2)
        write(62,*) Insx4(3)
	 	write(13,*)'current is',In2(1),In2(2),In2(3),In2(4)
	    WRITE(12,*)'current is',In3(1),In3(2),In3(3),In3(4)
		WRITE(11,*)'current is',In4(1),In4(2),In4(3),In4(4)
    	write(17,*) Insz(2)
!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

800    CONTINUE
40    CONTINUE
      zx=0
      GOTO 60
      ENDIF
!     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     Calculation of the Free energy
	  DO 4444 jg=1,41,1
       DO 5555 tg=1,41,1
	    Meanrho1(jg,tg) = 0.0d0
	    scale=93000.0D0*24.0D0
	    Do hhh=1, tg , 2
         Meanrho1(jg,tg)=Meanrho1(jg,tg)+ A20In(1,hhh-1)/(3.0D0*scale)+&
4.0d0*A20In(1,hhh)/(3.0D0*scale)+A20In(1,hhh+1)/(3.0D0*scale)
	    End do
        Meanrho1(jg,tg) = Meanrho1(jg,tg) / 3.d0
!       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	    Meanrho2(jg,tg) = 0.0d0
	    Do hhh=1, jg , 2
	     Meanrho2(jg,tg)=Meanrho2(jg,tg)+ A20Insz(hhh-1,tg)/scale+&
4.0d0*A20Insz(hhh,tg)/scale+A20Insz(hhh+1,tg)/scale
	    End do
	    Meanrho2(jg,tg) = Meanrho2(jg,tg) / 3.d0

        write(81,*)'a(',jg,',',tg,',',kjj,')=',Meanrho1(jg,tg)+Meanrho2(jg,tg),';'

5555   CONTINUE
4444  CONTINUE

      mazx10(tt)=abs(A10In(tt,1))
      mazx20(tt)=abs(A20In(tt,1))
      thth=2
      DO 971 vv=2,21,1
       ksA10In(tt,vv)=A10In(tt,vv)
       ksA20In(tt,vv)=A20In(tt,vv)
       IF (A10In(tt,vv).LT.0.0D0) THEN
        A10In(tt,vv)=abs(A10In(tt,vv))
       ENDIF
       IF (A10In(tt,vv).GT.abs(mazx10(tt))) THEN
        mazx10(tt)=ksA10In(tt,vv)
        thth=vv
       ENDIF
       IF (A20In(tt,vv).LT.0.0D0) THEN
        A20In(tt,vv)=abs(A20In(tt,vv))
       ENDIF
971	  CONTINUE
      write(20,*)mazx10(tt)
      write(26,*)mazx20(tt)
      mazxsx10(tt)=abs(A10Insx(tt,1))
      mazxsx20(tt)=abs(A20Insx(tt,1))

      DO 9713 vv=2,21,1
       ksA10Insx(tt,vv)=A10Insx(tt,vv)
       ksA20Insx(tt,vv)=A20Insx(tt,vv)
       IF (A10Insx(tt,vv).LT.0.0D0) THEN
        A10Insx(tt,vv)=abs(A10Insx(tt,vv))
       ENDIF
       IF (A10Insx(tt,vv).GT.abs(mazxsx10(tt))) THEN
        mazxsx10(tt)=ksA10Insx(tt,vv)
       ENDIF
       IF (A20Insx(tt,vv).LT.0.0D0) THEN
        A20Insx(tt,vv)=abs(A20Insx(tt,vv))
       ENDIF
       IF (A20Insx(tt,vv).GT.abs(mazxsx20(tt))) THEN
        mazxsx20(tt)=ksA20Insx(tt,vv)
       ENDIF
9713  CONTINUE
      write(30,*)mazxsx10(tt)
      write(36,*)mazxsx20(tt)
      mazxsy10(tt)=abs(A10Insy(tt,1))
      mazxsy20(tt)=abs(A20Insy(tt,1))

      DO 9714 vv=2,21,1
       ksA10Insy(tt,vv)=A10Insy(tt,vv)
       ksA20Insy(tt,vv)=A20Insy(tt,vv)
        IF (A10Insy(tt,vv).LT.0.0D0) THEN
         A10Insy(tt,vv)=abs(A10Insy(tt,vv))
        ENDIF
        IF (A10Insy(tt,vv).GT.abs(mazxsy10(tt))) THEN
         mazxsy10(tt)=ksA10Insy(tt,vv)
        ENDIF
        IF (A20Insy(tt,vv).LT.0.0D0) THEN
         A20Insy(tt,vv)=abs(A20Insy(tt,vv))
        ENDIF
        IF (A20Insy(tt,vv).GT.abs(mazxsy20(tt))) THEN
         mazxsy20(tt)=ksA20Insy(tt,vv)
        ENDIF
9714  CONTINUE

      write(40,*)mazxsy10(tt)
      write(46,*)mazxsy20(tt)

      mazxsz10(tt)=abs(A10Insz(tt,1))
      mazxsz20(tt)=abs(A20Insz(tt,1))
	  
      DO 9715 vv=2,21,1
       ksA10Insz(tt,vv)=A10Insz(tt,vv)
       ksA20Insz(tt,vv)=A20Insz(tt,vv)
       IF (A10Insz(tt,vv).LT.0.0D0) THEN
        A10Insz(tt,vv)=abs(A10Insz(tt,vv))
       ENDIF
       IF (A10Insz(tt,vv).GT.mazxsz10(tt)) THEN
        mazxsz10(tt)=A10Insz(tt,vv)
		write(50,*) mazxsz10(tt)
       ENDIF
       IF (A20Insz(tt,vv).LT.0.0D0) THEN
        A20Insz(tt,vv)=abs(A20Insz(tt,vv))
       ENDIF
       IF (A20Insz(tt,vv).GT.abs(mazxsz20(tt))) THEN
        mazxsz20(tt)=ksA20Insz(tt,vv)
       ENDIF
9715  CONTINUE
      write(56,*)mazxsz10(tt)

50    CONTINUE
30    CONTINUE
      END PROGRAM DFSFNFFS
