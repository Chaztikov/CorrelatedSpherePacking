c------------------------------------------------------------------------
c------------------------------------------------------------------------
      
      PROGRAM RandomPack
	
  	c	WRITE (*,*) 'Start the rnpack...'

      OPEN (3, FILE = 'packin.par', STATUS = 'OLD')
      OPEN (7, FILE = 'pack.par')
      OPEN (8, FILE = 'rntmp.dat')
      OPEN (9, FILE = 'pack.out')
 
      WRITE (*,*) 'Read the parameters...'

      CALL ReadParam

      WRITE (*,*) 'Generate the random packing with overlaps...'

      CALL RndmPack

      WRITE (*,*) 'Allocate particles into compartments...'

      CALL AlloctCmpt

      WRITE (*,*) 'Seek overlap-free packing conditions...'

      CALL Sort

      WRITE (*,*) 'End the rnpack!!'

      END

c------------------------------------------------------------------------

		SUBROUTINE ReadParam

      IMPLICIT REAL (a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
      COMMON /local/ intensity, lshift, move
		COMMON /corr/ ncor, gridln, Ngrid, Nzgrid, mr
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		READ(3,1000)
		READ(3,1010) alphar, betar, cnfint, porosity
      READ(3,1020) intensity, mpBC, ncor, nc, ncz, nrndsd, sdnrm
      READ(3,1030) pi, eps, tolrnce

		IF(ncor.EQ.1) THEN
         READ(3,1040) gridln, Ngrid, Nzgrid, mr
		ENDIF

1000  FORMAT()
1010  FORMAT(/,4f9.7)
1020  FORMAT(/,6i7,e13.5)
1030  FORMAT(/,f13.7,2e13.5)
1040  FORMAT(/,f13.7,3i7)

      END

c------------------------------------------------------------------------

		SUBROUTINE RndmPack

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

      COMMON /corr/ ncor, gridln, Ngrid, Nzgrid, mr

c.....Generate the initial random packing

      IF(ncor.EQ.0) THEN
		   CALL GnrtRdii
		   CALL GnrtSite

c.....Generate the initial random packing with spatial correlations

      ELSE
			CALL GnrtPack
      ENDIF

c.....Output the random packing with overlaps

		CALL Output 

      END

c------------------------------------------------------------------------

		SUBROUTINE GnrtRdii
      
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
      COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

c.....Set the radius range, based on the spacified confidence interval

		rminlog = betar*RInvNorm(0.5*(1.0-cnfint))+alphar
		rmaxlog = betar*RInvNorm(cnfint+0.5*(1.0-cnfint))+alphar

c.....Initialize

		npart = 0
		Vpack = 0.0
		rlmt = 0.05

c.....Generate the log-normal distribution of particle radii

      rmax = EXP(rmaxlog)
      rmin = EXP(rminlog)
		cmptln = 2.0*rmax+rlmt
		size = cmptln*FLOAT(nc)
		height = cmptln*FLOAT(ncz)
		ravg = EXP(alphar)
		Vt = size*size*height
		Pd = 0.0
     
		IF(ABS(betar).LE.tolrnce) THEN
			Vsphere = 4.0*pi*ravg*ravg*ravg/3.0
			npart = INT(Vt*(1-porosity)/Vsphere)
			DO i = 1, npart
				r(i) = ravg
         ENDDO
      ELSE 
			DO WHILE(Pd.LT.(1-porosity))
			   npart = npart+1
			   sample = rmaxlog
			   DO WHILE(sample.GE.rmaxlog.OR.sample.LE.rminlog)
			   	CALL RNorm(alphar,betar,sample,sdnrm)
            ENDDO
			   r(npart) = EXP(sample)
			   Vpack = Vpack+4.0*pi*r(npart)*r(npart)*r(npart)/3.0
			   Pd = Vpack/Vt
         ENDDO
      ENDIF
		
	c	WRITE(*,*) 'The maximum radius = ', rmax
	c	WRITE(*,*) 'The average radius = ', ravg
	c	WRITE(*,*) 'The sample size = ', size
      WRITE(*,*) 'The sampple height = ', height
	c	WRITE(*,*) 'The number of spheres in the initial packing = ',
     1           npart

		bound = cmptln-2.0*rmax

		END

c------------------------------------------------------------------------

      SUBROUTINE RNorm(Mu,Sigma,Sample,Seed)
      REAL Umin,Umax,Seed,Mu,U1,U2,Sigma,Sample

      Umin = 0.00
      Umax = 1.00
      
100   U1=(MOD(24298.00*Seed + 99991.00, 199017.00))
      Seed = U1
      U1=U1*(Umax-Umin)/199017.00 + Umin
      U2=(MOD(24298.00*Seed + 99991.00,199017.00))
      Seed=U2
      U2=U2*(Umax-Umin)/199017.00+Umin
      IF(U1.EQ.0.00) GOTO 100
      Sample=SQRT(-2.00*LOG(U1))*COS(2.00*3.1415900*U2)*Sigma+Mu

		END

c------------------------------------------------------------------------
      
		FUNCTION RInvNorm(P)
      IMPLICIT REAL(A-H,O-Z)
C
C         INVERSE CUMULATIVE NORMAL DISTRIBUTION FUNCTION.
C            P = CNORM(RInvNorm(P))
C     THIS IS A COPY OF ALGORITHM AS 24,
C     APPLIED STATISTICS, VOLUME 18, NUMBER 3, 1969
C     WRITTEN BY S.W. CUNNINGHAM.
C
      COMMON /MCHCOM/ BASE,DNOT(3),FINITY
      DIMENSION A(5),CONNOR(17),HSTNGS(6)
      DATA CONNOR/8.0327350124E-17, 1.4483264644E-15, 2.4668270103E-14,
     *            3.9554295164E-13, 5.9477940136E-12, 8.3507027951E-11,
     *            1.0892221037E-9 , 1.3122532964E-8 , 1.4503852223E-7 ,
     *            1.4589169001E-6 , 1.3227513228E-5 , 1.0683760684E-4 ,
     *            7.5757575758E-4 , 4.6296296296E-3 , 2.3809523810E-2 ,
     *            0.10           , 0.333333333330  /
      DATA RTHFPI/1.25331413730/
      DATA RRT2PI/0.39894228040/
      DATA TERMIN/1.0E-11/
      DATA HSTNGS/2.5155170, 0.8028530, 0.0103280,
     *            1.4327880, 0.1892690, 0.0013080/
C
C         IF P IS OUT OF RANGE, RETURN POSITIVE OR
C     NEGATIVE INFINITY.
C
      IF (P .GT. 0.00) GOTO 91
        ORD=-(FINITY/(BASE**2))
        GOTO 10
   91 IF (P .LT. 1.00) GOTO 92
        ORD = FINITY
        GOTO 10
   92 CONTINUE
C
C         GET FIRST APPROXIMATION, X0, TO ORDINATE BY
C     HASTINGS' FORMULA.
C
      B = P
      IF (B .GT. 0.50) B=1.00-B
      F = -LOG(B)
      E = SQRT(F+F)
      X0 = -E + ((HSTNGS(3)*E+HSTNGS(2))*E+HSTNGS(1))/
     *          (((HSTNGS(6)*E+HSTNGS(5))*E+HSTNGS(4))*E+1.00)
      IF (X0 .LT. 0.00) GO TO 1
      X0 = 0.00
      P0 = 0.50
      X1 = -RTHFPI
      GO TO 7
C
C         FIND THE AREA, P0, CORRESPONDING TO X0.
C
   1  Y = X0**2
      IF (X0.LE.-1.90) GO TO 3
      Y = -0.50*Y
C
C         (1) SERIES APPROXIMATION.
C
      P0 = CONNOR(1)
      DO 2 L=2,17
        P0 = P0*Y + CONNOR(L)
   2  CONTINUE
      P0 = (P0*Y+1.00) * X0
      X1 = -(P0*RTHFPI) * EXP(-Y)
      P0 = P0*RRT2PI + 0.50
      GO TO 7
C
C         (2) CONTINUED FRACTION APPROXIMATION.
C
   3  Z = 1.00/Y
      A(2) = 1.00
      A(3) = 1.00
      A(4) = Z + 1.00
      A(5) = 1.00
      W = 2.00
   4  DO 6 L=1,3,2
        DO 5 J=1,2
          K = L+J
          KA = 7-K
          A(K) = A(KA) + A(K)*W*Z
   5    CONTINUE
        W = W+1.00
   6  CONTINUE
      APPRXU = A(2)/A(3)
      APPRXL = A(5)/A(4)
      C = APPRXU-APPRXL
      IF (C.GE.TERMIN) GO TO 4
      X1 = APPRXL/X0
      P0 = -X1*RRT2PI*EXP(-0.50*Y)
C
C         GET ACCURATE VALUE OF ORDINATE BY TAYLOR SERIES.
C     X1, X2, AND X3 ARE DERIVATIVES FOR THE TAYLOR SERIES.
C
   7  D = F + LOG(P0)
      X2 = X0*X1*X1 - X1
      X3 = X1**3 + 2.00*X0*X1*X2 - X2
      X = ((X3*D/3.00 + X2)*D/2.00 + X1)*D + X0
      ORD = X
      IF (P .GT. 0.50) ORD=-X
   10 RInvNorm = ORD
      
		END

c------------------------------------------------------------------------

      SUBROUTINE GnrtSite
		
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

      COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

      DO i = 1, npart
			x(i) = ABS(size*RAN(nrndsd))
			y(i) = ABS(size*RAN(nrndsd))
			z(i) = ABS(height*RAN(nrndsd))
      ENDDO

		END

c------------------------------------------------------------------------

      SUBROUTINE GnrtPack

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
      COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 
      COMMON /corr/ ncor, gridln, Ngrid, Nzgrid, mr

		DIMENSION cpdf(100), nradius(100), ncradius(100), 
     1          Pfun(100,100,100)

      OPEN(4, FILE = 'pfield.dat', STATUS = 'OLD')
		OPEN(12, FILE = 'cpsd.out')

c.....Generate the cumulative radius size distribution
		
		rminlog = betar*RInvNorm(0.5*(1.0-cnfint))+alphar
		rmaxlog = betar*RInvNorm(cnfint+0.5*(1.0-cnfint))+alphar
		rmax = EXP(rmaxlog)
		rmin = EXP(rminlog)
		ravg = EXP(alphar)
		delr = (rmax-rmin)/FLOAT(mr)
		rlmt = 0.05

      cmptln = 2.0*rmax+rlmt
		size = cmptln*FLOAT(nc)
		height = cmptln*FLOAT(ncz)

c.....Initialize

		DO k = 1, mr
			nradius(k) = 0
      ENDDO

		npart = 0
		Pd = 0.0
		Vpack = 0.0
		Vt = size*size*height

		DO isample = 1, 100000
			sample = rmaxlog
			DO WHILE(sample.GE.rmaxlog.OR.sample.LE.rminlog)
				CALL RNorm(alphar,betar,sample,sdnrm)
         ENDDO
			radius = EXP(sample)
			DO k = 1, mr
				IF(radius.GE.(rmin+delr*(k-1)).AND.
     1         radius.LT.(rmin+delr*k)) THEN
               nradius(k) = nradius(k)+1
            ENDIF
         ENDDO
      ENDDO

		ncradius(1) = nradius(1)
		
		DO k = 1, mr-1
			ncradius(k+1) = nradius(k+1)+ncradius(k)
      ENDDO

c.....Normalize the distribution
      
		DO k = 1, mr
			cpdf(k) = FLOAT(ncradius(k))/FLOAT(ncradius(mr))
		c	WRITE(12,2000) (rmin+delr*(k-1)), cpdf(k)
      ENDDO

c.....Read in the random P-field data

		READ(4,1000)
		DO kz = 1, Nzgrid
			DO jy = 1, Ngrid
				DO i1x = 1, Ngrid
					READ(4,1010) Pfun(i1x,jy,kz)
c    			c	WRITE(*,*) i1x,jy,kz,Pfun(i1x,jy,kz)
            ENDDO
         ENDDO
      ENDDO

c.....Create the initial random packing with spatial correlations

		DO WHILE(Pd.LT.(1-porosity))
			npart = npart+1
			
			xsample = ABS(size*RAN(nrndsd))
			DO l = 1, Ngrid
				IF(xsample.GE.gridln*(l-1).AND.
     1         xsample.LT.gridln*l) THEN
					i1x = l 
            ENDIF
         ENDDO
			ysample = ABS(size*RAN(nrndsd))
			DO l = 1, Ngrid
				IF(ysample.GE.gridln*(l-1).AND.
     1         ysample.LT.gridln*l) THEN
					jy = l 
            ENDIF
         ENDDO
			zsample = ABS(height*RAN(nrndsd))
			DO l = 1, Nzgrid
				IF(zsample.GE.gridln*(l-1).AND.
     1         zsample.LT.gridln*l) THEN
					kz = l 
            ENDIF
         ENDDO
		
		   x(npart) = xsample
		   y(npart) = ysample
		   z(npart) = zsample
 
			pf = Pfun(i1x,jy,kz)
			IF(pf.LT.cpdf(1)) THEN
				r(npart) = rmin-delr*RAN(nrndsd)
         ELSEIF(pf.GE.cpdf(mr)) THEN
				r(npart) = rmax+delr*RAN(nrndsd)
         ELSE
			   DO k = 1, mr-1
			   	IF(pf.GE.cpdf(k).AND.pf.LT.cpdf(k+1)) THEN
			   		r(npart) = rmin+delr*(k-1)+delr*RAN(nrndsd)
               ENDIF
            ENDDO
			ENDIF
			Vpack = Vpack+4.0*pi*r(npart)*r(npart)*r(npart)/3.0
			Pd = Vpack/Vt
      ENDDO

	c	WRITE(*,*) 'The maximum radius = ', rmax
	c	WRITE(*,*) 'The minimum radius = ', rmin
	c	WRITE(*,*) 'The average radius = ', ravg
	c	WRITE(*,*) 'The sample size = ', size
	c	WRITE(*,*) 'The sampple height = ', height
	c	WRITE(*,*) 'The number of spheres in the initial packing = ',
     1           npart


      bound = cmptln-2.0*rmax

	c	WRITE(*,*) 'Do you want to start rearranging process(1--Yes, 0--No)?'
		READ(*,*) lsignal
      IF(lsignal.EQ.0) stop 

1000  FORMAT(///)
1010  FORMAT(10x,f10.4)
2000  FORMAT(2(2x,e13.5))

		END

c------------------------------------------------------------------------

      SUBROUTINE Output

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

c     WRITE(8,2000) 
c	c	WRITE(8,*) (x(i), y(i), z(i), r(i), i = 1, npart)

		DO i = 1, npart
		   WRITE(8,2010) x(i), y(i), z(i), r(i) 
      ENDDO

2000  FORMAT('x(i)  y(i)  z(i)  r(i) ... ')
2010  FORMAT(4(2x,e13.5))

		END

c------------------------------------------------------------------------

      SUBROUTINE AlloctCmpt
		
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)
	
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

c.....Read in particles and allocate them into compartments

		REWIND 8
c		READ(8,*)
c		READ(8,*) (x(i), y(i), z(i), r(i), i = 1, npart)
     
		DO i = 1, npart
		   READ(8,1000) x(i), y(i), z(i), r(i)
		ENDDO

	c	WRITE(9,2000)

		mz = 0
		DO kz = 1, ncz
			CALL InitialCnt
		   DO i = 1, npart
				CALL Allocating(i,mz)
         ENDDO

c.....Sum up the total counts in the layer

			DO jy = 1, nc
				DO i1x = 1, nc
					npart_kz = npart_kz+icnt(i1x,jy)
            ENDDO
         ENDDO

c.....Write particles in file after allocation

		c	WRITE(9,2010) npart_kz
			DO jy = 1, nc
            DO i1x = 1, nc
               ncpart = icnt(i1x,jy)
				c	WRITE(9,2020) i1x, jy, ncpart
c				c	WRITE(9,*) (xc(i1x,jy,ic), yc(i1x,jy,ic), zc(i1x,jy,ic), 
c     1                     rc(i1x,jy,ic), ic = 1, ncpart) 
				   DO ic = 1, ncpart
					c	WRITE(9,2030) xc(i1x,jy,ic), yc(i1x,jy,ic),
     1                          zc(i1x,jy,ic), rc(i1x,jy,ic)
					ENDDO
				ENDDO
         ENDDO
			mz = mz+1
      ENDDO

1000  FORMAT(4(2x,e13.5))
2000  FORMAT('npart_kz',/,'i1x  jy  ncpart',/,'xc  yc  zc  rc')
2010  FORMAT(i7)
2020  FORMAT(3i7)
2030  FORMAT(4(2x,e13.5))

		END

c------------------------------------------------------------------------

		SUBROUTINE InitialCnt

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

c.....Initiallize counts in each compartment

      DO jy = 1, nc
         DO i1x = 1, nc
				icnt(i1x,jy) = 0
         ENDDO
      ENDDO

c.....Initialize count in each layer

		npart_kz = 0

		END

c------------------------------------------------------------------------

      SUBROUTINE Allocating(i,mz)
		
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 
		
		my = 0
      DO jy = 1, nc
			mx = 0
         DO i1x = 1, nc
            IF(x(i).GE.mx*cmptln.AND.x(i).LT.(mx+1)*cmptln.AND.
     1         y(i).GE.my*cmptln.AND.y(i).LT.(my+1)*cmptln.AND.
     2         z(i).GE.mz*cmptln.AND.z(i).LT.(mz+1)*cmptln) THEN
					icnt(i1x,jy) = icnt(i1x,jy)+1
               ic = icnt(i1x,jy)
					xc(i1x,jy,ic) = x(i)
					yc(i1x,jy,ic) = y(i)
					zc(i1x,jy,ic) = z(i)
					rc(i1x,jy,ic) = r(i)
            ENDIF
				mx = mx+1
         ENDDO
			my = my+1
      ENDDO

		END

c------------------------------------------------------------------------

		SUBROUTINE Sort

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
 
c.....Initialize after-sort counts

		DO lz = 1, ncz+1 
	      CALL InitialSort(lz)
      ENDDO

		niterate = 0

c.....Sort particles layer by layer

		REWIND 9
		READ(9,1000)

		DO lz = 1, ncz
		   READ(9,*) npart_kz
			DO jy = 1, nc
				DO i1x = 1, nc
			      READ(9,*) dummy, dummy, ncpart
               icnt(i1x,jy) = ncpart
               DO ic = 1, ncpart
						READ(9, 1010) xc(i1x,jy,ic), yc(i1x,jy,ic),
     1                          zc(i1x,jy,ic), rc(i1x,jy,ic)
               ENDDO
				ENDDO
         ENDDO

   	c	WRITE(*,*) 'Working on the vertical cell: ', lz

			DO ipart = 1, npart_kz
				CALL FindZmin(ic1,ncpart1)
				CALL Reorder(ic1,ncpart1)
				CALL Sorting(lz)
         ENDDO

      ENDDO

c.....Calculate the final porosity and write out the packing data
		
      CALL FinalSet


1000  FORMAT(//)
1010  FORMAT(4(2x,e13.5))

		END

c------------------------------------------------------------------------

		SUBROUTINE InitialSort(lz)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 

c.....Initialize counts in a layer of lz after sort

      DO jy = 1, nc
         DO i1x = 1, nc
            icntst(i1x,jy,lz) = 0
         ENDDO
      ENDDO

		END

c------------------------------------------------------------------------

		SUBROUTINE FindZmin(ic1,ncpart1)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

		z1 = size
      
		DO jy = 1, nc
         DO i1x = 1, nc
            ncpart = icnt(i1x,jy)
				DO ic = 1, ncpart
					IF(zc(i1x,jy,ic).LT.z1) THEN
						z1 = zc(i1x,jy,ic)
						ic1 = ic
						ncpart1 = ncpart
						i1x1 = i1x
						jy1 = jy
						
               ENDIF
            ENDDO
			ENDDO
      ENDDO

		END

c------------------------------------------------------------------------

		SUBROUTINE Reorder(ic1,ncpart1)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 

      x1 = xc(i1x1,jy1,ic1)
      y1 = yc(i1x1,jy1,ic1)
      z1 = zc(i1x1,jy1,ic1)
      r1 = rc(i1x1,jy1,ic1)

      icnt(i1x1,jy1) = icnt(i1x1,jy1)-1
		DO ic = ic1+1, ncpart1
			xc(i1x1,jy1,ic-1) = xc(i1x1,jy1,ic)
			yc(i1x1,jy1,ic-1) = yc(i1x1,jy1,ic)
			zc(i1x1,jy1,ic-1) = zc(i1x1,jy1,ic)
			rc(i1x1,jy1,ic-1) = rc(i1x1,jy1,ic)
      ENDDO

		END

c------------------------------------------------------------------------

      SUBROUTINE Sorting(lz)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /local/ intensity, lshift, move
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		kz1 = lz
		ntshift = 0

c.....Check overlap conditions

		mroll = 0
10		CALL FindXYZ
		lshift = 1
      
		IF(z1.GE.bound.AND.x1.GE.bound.AND.y1.GE.bound.AND.
     1   x1.LE.(size-bound).AND.y1.LE.(size-bound)) THEN
         CALL FindLoc(lz)
			IF(nrec.EQ.1) GOTO 20
      ENDIF

		IF(z1.GE.bound.AND.(x1.LT.bound.OR.x1.GT.(size-bound).OR.
     1   y1.LT.bound.OR.y1.GT.(size-bound))) THEN
			CALL FindLocBC(lz)
			IF(mcount.EQ.11) GOTO 10
      ENDIF

20    CALL Register

		END

c------------------------------------------------------------------------

      SUBROUTINE FindXYZ

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /tmp/ mrolltmp

      mflag = 1
		movlpcond = 0
		mrolltmp = 0
		zorig = z1

		DO WHILE(mflag.EQ.1)
			IF(i1x1.EQ.nc) i1xx = nc-1
			IF(i1x1.LT.nc) i1xx = i1x1
			IF(jy1.EQ.nc) jyy = nc-1
			IF(jy1.LT.nc) jyy = jy1
			IF(kz1.EQ.(ncz+1)) kzz = ncz
         IF(kz1.LT.(ncz+1)) kzz = kz1
			mflag = 0
            IF(kz1.EQ.1.OR.kz1.EQ.(ncz+1)) THEN
                  IF(i1x1.EQ.1.AND.jy1.EQ.1.OR.i1x1.EQ.nc.AND.jy1.EQ.1.OR.i1x1.EQ.1.AND.jy1.EQ.nc.OR.i1x1.EQ.nc.AND.jy1.EQ.nc) THEN
                        CALL Shift1(i1xx,jyy,kzz,mflag)
            ELSEIF(i1x1.EQ.1.OR.i1x1.EQ.nc.OR.jy1.EQ.1.OR.jy1.EQ.nc) THEN
					CALL Shift2(i1xx,jyy,kzz,mflag)
            ELSE
					CALL Shift3(i1xx,jyy,kzz,mflag)
            ENDIF
         ELSEIF(jy1.EQ.1.OR.jy1.EQ.nc) THEN
            IF(i1x1.EQ.1.OR.i1x1.EQ.nc) THEN
					CALL Shift2(i1xx,jyy,kzz,mflag)
            ELSE
					CALL Shift3(i1xx,jyy,kzz,mflag)
            ENDIF
         ELSEIF(i1x1.EQ.1.OR.i1x1.EQ.nc) THEN
				CALL Shift3(i1xx,jyy,kzz,mflag)
         ELSE
				CALL Shift4(i1xx,jyy,kzz,mflag)
         ENDIF
      ENDDO

		IF(mroll.EQ.0.AND.ABS(zorig-z1).GT.tolrnce) mroll = 1
      IF(mrolltmp.EQ.3.AND.movlpcond.EQ.0) mroll = mrolltmp

		END

c------------------------------------------------------------------------

      SUBROUTINE Shift1(i1xx,jyy,kzz,mflag)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzz, kzz+1
			DO jy = jyy, jyy+1
				DO i1x = i1xx, i1xx+1
					CALL ChkOvlp(mflag)
					IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
            ENDDO
         ENDDO
      ENDDO

10    END

c------------------------------------------------------------------------

      SUBROUTINE Shift2(i1xx,jyy,kzz,mflag)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      IF(i1x1.GE.2.AND.i1x1.LT.nc) THEN
			DO kz = kzz, kzz+1
				DO jy = jyy, jyy+1
					DO i1x = i1xx-1, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(jy1.GE.2.AND.jy1.LT.nc) THEN
			DO kz = kzz, kzz+1
				DO jy = jyy-1, jyy+1
					DO i1x = i1xx, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(kz1.GE.2.AND.kz1.LE.ncz) THEN
			DO kz = kzz-1, kzz+1
				DO jy = jyy, jyy+1
					DO i1x = i1xx, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ENDIF

10    END

c------------------------------------------------------------------------

      SUBROUTINE Shift3(i1xx,jyy,kzz,mflag)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      IF(i1x1.EQ.1.OR.i1x1.EQ.nc) THEN
			DO kz = kzz-1, kzz+1
				DO jy = jyy-1, jyy+1
					DO i1x = i1xx, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(jy1.EQ.1.OR.jy1.EQ.nc) THEN
			DO kz = kzz-1, kzz+1
				DO jy = jyy, jyy+1
					DO i1x = i1xx-1, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(kz1.EQ.1.OR.kz1.EQ.(ncz+1)) THEN
			DO kz = kzz, kzz+1
				DO jy = jyy-1, jyy+1
					DO i1x = i1xx-1, i1xx+1
						CALL ChkOvlp(mflag)
					   IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
               ENDDO
            ENDDO
         ENDDO
      ENDIF

10    END

c------------------------------------------------------------------------

      SUBROUTINE Shift4(i1xx,jyy,kzz,mflag)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzz-1, kzz+1
			DO jy = jyy-1, jyy+1
				DO i1x = i1xx-1, i1xx+1
					CALL ChkOvlp(mflag)
					IF(mflag.EQ.1.OR.movlpcond.EQ.1) GOTO 10
            ENDDO
         ENDDO
      ENDDO

10    END

c------------------------------------------------------------------------

      SUBROUTINE ChkOvlp(mflag)

      IMPLICIT REAL(a-h,o-z)

      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate


      IF(mroll.EQ.0) THEN
			CALL FindZmax(mflag) 
	   ELSEIF(mroll.EQ.1) THEN
			CALL Roll1Loc
			GOTO 10
      ELSEIF(mroll.EQ.2) THEN
			CALL Roll2Loc
      ENDIF

10    END

c------------------------------------------------------------------------
      
		SUBROUTINE FindZmax(mflag)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 

		ncntst = icntst(i1x,jy,kz)

		IF(ncntst.EQ.0) GOTO 20

		DO icst = 1, ncntst
			xr = xcst(i1x,jy,kz,icst)
			yr = ycst(i1x,jy,kz,icst)
			zr = zcst(i1x,jy,kz,icst)
			rr = rcst(i1x,jy,kz,icst)
         d = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr)+(z1-zr)*(z1-zr))
			sumr = r1+rr
			IF(d.LT.(sumr-eps)) THEN
10          dxy = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr))
				IF(dxy.LE.tolrnce) THEN
               x1 = x1+1.0e-4
               y1 = y1+1.0e-4
					IF(x1.GE.cmptln*i1x1.AND.i1x1.LT.nc) i1x1 = i1x1+1
					IF(y1.GE.cmptln*jy1.AND.jy1.LT.nc) jy1 = jy1+1
				   GOTO 10
            ENDIF
				cosineA = dxy/sumr
				IF(cosineA.GT.1.0) cosineA = 1.0
				sineA = SQRT(1.0-cosineA*cosineA)
				zp = sumr*sineA
				z1 = zr+zp
				x2 = xr
				y2 = yr
				z2 = zr
            r2 = rr
				i1x2 = i1x
				jy2 = jy 
				kz2 = kz
				icst2 = icst 
				IF(z1.GE.cmptln*kz.AND.kz1.LE.ncz) kz1 = kz+1
				mflag = 1
				GOTO 20
         ENDIF
      ENDDO

20    END

c------------------------------------------------------------------------

      SUBROUTINE Roll1Loc

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      ncntst = icntst(i1x,jy,kz)
      
      DO icst = 1, ncntst
			IF(.NOT.(i1x.EQ.i1x2.AND.jy.EQ.jy2.AND.kz.EQ.kz2.AND.
     1      icst.EQ.icst2)) THEN
				CALL CalR1Loc(icst)
				IF(movlpcond.EQ.1) GOTO 10
	      ENDIF
      ENDDO
		
10    END

c------------------------------------------------------------------------

		SUBROUTINE CalR1Loc(icst)
	   
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

   	xr = xcst(i1x,jy,kz,icst)
      yr = ycst(i1x,jy,kz,icst)
      zr = zcst(i1x,jy,kz,icst)
      rr = rcst(i1x,jy,kz,icst)
      d = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr)+(z1-zr)*(z1-zr))
      sumr = r1+rr
      IF(d.GE.(sumr-eps).AND.d.LE.sumr) THEN
		   x3 = xr
		   y3 = yr
		   z3 = zr
         r3 = rr
	      i1x3 = i1x
		   jy3 = jy 
		   kz3 = kz
		   icst3 = icst 
		   mroll = 2
	   ELSEIF(d.LT.(sumr-eps)) THEN
			IF(dtheta.LE.eps) THEN
	   	   x3 = xr
	   	   y3 = yr
	   	   z3 = zr
            r3 = rr
	         i1x3 = i1x
		      jy3 = jy 
		      kz3 = kz
		      icst3 = icst 
		      mroll = 2
         ELSE
            theta = theta-dtheta
            dtheta = 0.5*dtheta
			   movlpcond = 1
         ENDIF
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE Roll2Loc

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      ncntst = icntst(i1x,jy,kz)

      DO icst = 1, ncntst
			IF((.NOT.(i1x.EQ.i1x2.AND.jy.EQ.jy2.AND.kz.EQ.kz2.AND.
     1      icst.EQ.icst2)).AND.(.NOT.(i1x.EQ.i1x3.AND.
     2      jy.EQ.jy3.AND.kz.EQ.kz3.AND.icst.EQ.icst3))) THEN
			   CALL CalR2Loc(icst)
	   	   IF(movlpcond.EQ.1) GOTO 10
         ENDIF
      ENDDO

10    END

c------------------------------------------------------------------------

		SUBROUTINE CalR2Loc(icst)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /tmp/ mrolltmp

		xr = xcst(i1x,jy,kz,icst)
      yr = ycst(i1x,jy,kz,icst)
      zr = zcst(i1x,jy,kz,icst)
      rr = rcst(i1x,jy,kz,icst)
      d = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr)+(z1-zr)*(z1-zr))
      sumr = r1+rr
      IF(d.GE.(sumr-eps).AND.d.LE.sumr) THEN
		   x4 = xr
         y4 = yr
         z4 = zr
         r4 = rr
         i1x4 = i1x
         jy4 = jy 
         kz4 = kz
         icst4 = icst 
	      mrolltmp = 3 
      ELSEIF(d.LT.(sumr-eps)) THEN
			IF(dtheta.LE.eps) THEN
		      x4 = xr
            y4 = yr
            z4 = zr
            r4 = rr
            i1x4 = i1x
            jy4 = jy 
            kz4 = kz
            icst4 = icst 
	         mrolltmp = 3 
         ELSE
            theta = theta-dtheta
		      dtheta = 0.5*dtheta
		      movlpcond = 1
         ENDIF
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE Register

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      icst1 = icntst(i1x1,jy1,kz1)+1
      icntst(i1x1,jy1,kz1) = icst1
      xcst(i1x1,jy1,kz1,icst1) = x1
      ycst(i1x1,jy1,kz1,icst1) = y1
      zcst(i1x1,jy1,kz1,icst1) = z1
      rcst(i1x1,jy1,kz1,icst1) = r1

c      write(*,2010) xcst(i1x1,jy1,kz1,icst1), ycst(i1x1,jy1,kz1,icst1),
c     1              zcst(i1x1,jy1,kz1,icst1), rcst(i1x1,jy1,kz1,icst1)
c      write(*,2020) i1x1,jy1,kz1,icst1

c      write(*,2000) rcst(i1x1,jy1,kz1,icst1), i1x1,jy1,kz1,icst1

2000	FORMAT(f8.5, 1x, 4(i3,i1x))
2010	FORMAT(4(e13.5,2x))
2020  FORMAT(4(i3,i1x))

		END

c------------------------------------------------------------------------

      SUBROUTINE ShiftDown

      IMPLICIT REAL(a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
 
		ravg = EXP(alphar)
		ncdiv = 30
		delz = ravg/FLOAT(ncdiv)

		DO WHILE(z1.GE.bound)
			z1 = z1-delz
			IF(z1.LT.cmptln*(kz1-1).AND.kz1.GT.1) THEN
				kz1 = kz1-1
         ENDIF
			mroll = 0
			CALL FindXYZ
			IF(mroll.EQ.1) GOTO 10
      ENDDO

10    END

c------------------------------------------------------------------------

      SUBROUTINE FindLocBC(lz)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		mcount = 0

      nrec = 0
   	mroll = 0
		DO WHILE(nrec.NE.1)
		   IF(mroll.EQ.0) CALL ShiftDown
		   IF(z1.LT.bound) THEN
		   	nrec = 1
		   	GOTO 20
         ENDIF
		   IF(mroll.EQ.1) THEN
		   	IF((x1.LT.bound.OR.x1.GT.(size-bound)).AND.
     1         (y1.LT.bound.OR.y1.GT.(size-bound))) THEN
					IF(kz1.GT.lz) THEN
						GOTO 10
               ELSE
					   nrec = 1
			         GOTO 20
               ENDIF
				ELSE
			   	CALL RollItBC
            ENDIF
         ENDIF
			IF(mroll.GE.2) THEN
				niterate = niterate+1
				IF(kz1.GT.lz) THEN
					GOTO 10
            ELSE
				   nrec = 1
		         GOTO 20
            ENDIF
			ENDIF
10      	IF(kz1.GT.lz) CALL RelocateBC
			IF(mcount.EQ.11) GOTO 20
      ENDDO

20    END

c------------------------------------------------------------------------

      SUBROUTINE RollItBC

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /flag/ mark
		
      IF(x1.LT.bound.OR.x1.GT.(size-bound)) THEN
			mark = 1
			o1 = x1
			o2 = x2
			p1 = y1
			p2 = y2
			lp1 = jy1
      ELSE
			mark = 2
			o1 = y1
			o2 = y2
			p1 = x1
			p2 = x2
			lp1 = i1x1
      ENDIF

		IF(ABS(p1-p2).LE.tolrnce) THEN
			p1 = p1+1.0e-4
			IF(p1.GE.cmptln*lp1.AND.lp1.LT.nc) lp1 = lp1+1
			rtrack = z1-z2
			sineT = 1.0e-4/rtrack
			IF(sineT.GT.1.0) sineT = 1.0
			cosineT = SQRT(1.0-sineT*sineT)
			z1 = z2+rtrack*cosineT
      ENDIF

		d12 = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
		o12 = ABS(o1-o2)
		rtrack = SQRT(d12*d12-o12*o12)
		delp = ABS(p1-p2)
      dptort = delp/rtrack
		IF(dptort.GT.1.0) THEN
			theta = 0.5*pi
      ELSE
         theta = ASIN(dptort)
      ENDIF
		
		dtheta = 0.15
		angle = 0.5*pi

		DO WHILE(theta.LT.angle)
			p1tmp = p1
			z1tmp = z1
			lp1tmp = lp1
			kz1tmp = kz1
10       theta = theta+dtheta
         IF(p1.GT.p2) THEN 
				p1 = p2+rtrack*SIN(theta)
			ELSE
				p1 = p2-rtrack*SIN(theta)
         ENDIF

			z1 = z2+rtrack*COS(theta)
			CALL ChkBCmpt(p1,lp1)
			CALL FindXYZBC(p1,lp1)
			IF(movlpcond.EQ.1) GOTO 10

			IF(mroll.EQ.2.OR.z1.LT.bound) THEN
				IF(mark.EQ.1) THEN
					IF(dtheta.LE.eps) THEN
						y1 = p1tmp
						z1 = z1tmp
						jy1 = lp1tmp
						kz1 = kz1tmp
               ELSE
					   y1 = p1
					   jy1 = lp1
               ENDIF
				ELSE
					IF(dtheta.LE.eps) THEN
						x1 = p1tmp
						z1 = z1tmp
						i1x1 = lp1tmp
						kz1 = kz1tmp
               ELSE
					   x1 = p1
					   i1x1 = lp1
               ENDIF
            ENDIF
				nrec = 1
				GOTO 20
         ENDIF
      ENDDO

		mroll = 0
      IF(mark.EQ.1) THEN
		   y1 = p1
		   jy1 = lp1
		ELSE
		   x1 = p1
		   i1x1 = lp1
      ENDIF

20    END

c------------------------------------------------------------------------

      SUBROUTINE ChkBCmpt(p1,lp1)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 

      IF(p1.GE.(cmptln*lp1).AND.lp1.LT.nc) lp1 = lp1+1

		IF(p1.LT.(cmptln*(lp1-1)).AND.lp1.GE.2) lp1 = lp1-1
 
		IF(z1.LT.(cmptln*(kz1-1)).AND.kz1.GE.2) kz1 = kz1-1

      END

c------------------------------------------------------------------------

      SUBROUTINE FindXYZBC(p1,lp1)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      movlpcond = 0

      IF(lp1.EQ.nc) lpp = nc-1
      IF(lp1.LT.nc) lpp = lp1
      IF(kz1.EQ.(ncz+1)) kzp = ncz
      IF(kz1.LT.(ncz+1)) kzp = kz1

      IF(lp1.EQ.1.AND.kz1.EQ.1.OR.lp1.EQ.nc.AND.kz1.EQ.1.OR.
     1   lp1.EQ.1.AND.kz1.EQ.(ncz+1).OR.lp1.EQ.nc.AND.
     2   kz1.EQ.(ncz+1)) THEN
         CALL BCRoll1(p1,lp1,lpp,kzp)
      ELSEIF(lp1.EQ.1.OR.lp1.EQ.nc) THEN
         CALL BCRoll2(p1,lp1,lpp,kzp)
      ELSEIF(kz1.EQ.1.OR.kz1.EQ.(ncz+1)) THEN
         CALL BCRoll3(p1,lp1,lpp,kzp)
      ELSE
         CALL BCRoll4(p1,lp1,lpp,kzp)
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE BCRoll1(p1,lp1,lpp,kzp)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzp, kzp+1
			DO lp = lpp, lpp+1
				CALL ChkBCOvlp(p1,lp1,lp)
				IF(movlpcond.EQ.1) GOTO 10
         ENDDO
      ENDDO
		 
10    END

c------------------------------------------------------------------------

      SUBROUTINE BCRoll2(p1,lp1,lpp,kzp)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzp-1, kzp+1
			DO lp = lpp, lpp+1
				CALL ChkBCOvlp(p1,lp1,lp)
				IF(movlpcond.EQ.1) GOTO 10
         ENDDO
      ENDDO
		 
10    END

c------------------------------------------------------------------------

      SUBROUTINE BCRoll3(p1,lp1,lpp,kzp)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzp, kzp+1
			DO lp = lpp-1, lpp+1
				CALL ChkBCOvlp(p1,lp1,lp)
				IF(movlpcond.EQ.1) GOTO 10
         ENDDO
      ENDDO
		 
10    END

c------------------------------------------------------------------------

      SUBROUTINE BCRoll4(p1,lp1,lpp,kzp)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      DO kz = kzp-1, kzp+1
			DO lp = lpp-1, lpp+1
				CALL ChkBCOvlp(p1,lp1,lp)
				IF(movlpcond.EQ.1) GOTO 10
         ENDDO
      ENDDO
		 
10    END

c------------------------------------------------------------------------

      SUBROUTINE ChkBCOvlp(p1,lp1,lp)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /flag/ mark

      IF(mark.EQ.1) THEN
			y1 = p1
			jy = lp
			i1x = i1x1
      ELSE
         x1 = p1
			i1x = lp
			jy = jy1
      ENDIF

		ncntst = icntst(i1x,jy,kz)

      DO icst = 1, ncntst
			IF(.NOT.(i1x.EQ.i1x2.AND.jy.EQ.jy2.AND.kz.EQ.kz2.AND.
     1      icst.EQ.icst2)) THEN
				CALL CalRLocBC(icst)
				IF(movlpcond.EQ.1) GOTO 10
	      ENDIF
      ENDDO
		
10    END

c------------------------------------------------------------------------

		SUBROUTINE CalRLocBC(icst)
		
      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

	   xr = xcst(i1x,jy,kz,icst)
	   yr = ycst(i1x,jy,kz,icst)
	   zr = zcst(i1x,jy,kz,icst)
	   rr = rcst(i1x,jy,kz,icst)
      d = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr)+(z1-zr)*(z1-zr))
      sumr = r1+rr
      IF(d.GE.(sumr-eps).AND.d.LE.sumr) THEN
			mroll = 2
	   ELSEIF(d.LT.(sumr-eps)) THEN
			IF(dtheta.LE.eps) THEN
				mroll = 2
         ELSE 
			   theta = theta-dtheta
			   dtheta = 0.5*dtheta
			   movlpcond = 1
         ENDIF
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE RelocateBC
      
		IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /sign/ lsign
      COMMON /flag/ mark

      nchalf = INT(0.5*nc)
      IF(mcount.EQ.10) THEN
			IF(mark.EQ.1) THEN
				IF(x1.LT.bound) THEN
					x1 = x1+cmptln
					IF(x1.GE.cmptln) i1x1 = i1x1+1
            ELSE
					x1 = x1-cmptln
					IF(x1.LT.(size-cmptln)) i1x1 = i1x1-1
            ENDIF
	      ELSE				
				IF(y1.LT.bound) THEN
					y1 = y1+cmptln
					IF(y1.GE.cmptln) jy1 = jy1+1
            ELSE
					y1 = y1-cmptln
					IF(y1.LT.(size-cmptln)) jy1 = jy1-1
            ENDIF
         ENDIF
			GOTO 20
      ENDIF

		IF(ntshift.GE.60) THEN
         r1 = 0.9*r1
         mcount = 0
			ntshift = 0
      ENDIF
			
10    IF(mcount.EQ.0) THEN
			IF((x1.LT.bound.OR.x1.GT.(size-bound)).AND.
     1      (y1.LT.bound.OR.y1.GT.(size-bound))) THEN
			   IF(jy1.LE.nchalf) THEN
					lsign = 1
            ELSE
					lsign = -1
            ENDIF
			ELSEIF(x1.LT.bound.OR.x1.GT.(size-bound)) THEN
			   IF(jy1.LE.nchalf) THEN
					lsign = 1
            ELSE
					lsign = -1
            ENDIF
         ELSE
				IF(i1x1.LE.nchalf) THEN
					lsign = 1
            ELSE
					lsign = -1
            ENDIF
         ENDIF
      ENDIF

		IF(mark.EQ.1) THEN
			IF(y1.GT.size) THEN
				y1 = 2.0*size-y1
         ELSEIF(y1.LT.0.0) THEN
				y1 = -y1
         ENDIF
			y1 = y1+cmptln*lsign
			jy1 = jy1+lsign
      ELSE
			IF(x1.GT.size) THEN
				x1 = 2.0*size-x1
         ELSEIF(x1.LT.0.0) THEN
				x1 = -x1
         ENDIF
			x1 = x1+cmptln*lsign
			i1x1 = i1x1+lsign
      ENDIF

		IF(i1x1.EQ.0.OR.i1x1.EQ.(nc+1)) THEN
			x1 = x1-2.0*cmptln*FLOAT(lsign)
			i1x1 = i1x1-2*lsign
         r1 = 0.9*r1
         mcount = 0
			ntshift = 0
		   GOTO 10
      ELSEIF(jy1.EQ.0.OR.jy1.EQ.(nc+1)) THEN
			y1 = y1-2.0*cmptln*FLOAT(lsign)
			jy1 = jy1-2*lsign
         r1 = 0.9*r1
         mcount = 0
			ntshift = 0
		   GOTO 10
      ENDIF

20    mcount = mcount+1
		mroll = 0
		nrec = 0
		ntshift = ntshift+1

      END

c------------------------------------------------------------------------

      SUBROUTINE FindLoc(lz)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /local/ intensity, lshift, move
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      mcount = 0
      nrec = 0

10    IF(ABS(r1-rlmt).LE.tolrnce) THEN
			DO WHILE(nrec.NE.1) 
			   IF(mroll.EQ.0) CALL ShiftDown
				IF(z1.LT.bound) THEN
					nrec = 1
					GOTO 20
            ENDIF
				IF(nrec.EQ.2) GOTO 20
				IF(mroll.EQ.2) THEN
					CALL SettleRlmt
					niterate = niterate+1
            ENDIF
         ENDDO
      ELSE
         nloop = 0
         DO WHILE(nrec.NE.1)
				IF(mroll.EQ.0) CALL ShiftDown
				IF(z1.LT.bound) THEN
					nrec = 1
					GOTO 20
            ENDIF
				IF(mroll.EQ.1) CALL RollOn1
				IF(nrec.EQ.2) GOTO 20
				IF(mroll.EQ.2) CALL RollOn2
				IF(mroll.EQ.3) CALL ChkCond
				IF(nrec.NE.1) nloop = nloop+1
				IF(nloop.GE.30) THEN
					CALL Relocate
					nloop = 0
					lshift = 1
            ENDIF
				IF(nrec.EQ.1) THEN
					niterate = niterate+1
					IF(z1.GE.bound) THEN
				   	CALL TmpRec
               ENDIF
            ENDIF
				IF(nrec.EQ.1.AND.lshift.EQ.move.AND.kz1.GT.lz) THEN
				   CALL Relocate
				   nloop = 0
				   lshift = 1
            ENDIF
				IF(ABS(r1-rlmt).LE.tolrnce) THEN
					GOTO 10
            ENDIF
         ENDDO
      ENDIF

20    END

c------------------------------------------------------------------------

      SUBROUTINE RollOn1

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		dx = ABS(x1-x2)
		dxy = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
		dxtodxy = dx/dxy
		IF(dxtodxy.GT.1.0) THEN
			beta = 0.5*pi
      ELSE
		   beta = ASIN(dxtodxy)
      ENDIF

		d12 = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
		dxytod12 = dxy/d12
		IF(dxytod12.GT.1.0) dxytod12 = 1.0
		IF(z1.GT.z2) THEN
			theta = ASIN(dxytod12)
		ELSE
			mroll = 0
			GOTO 20
      ENDIF
		dtheta = 0.15
		angle = 0.5*pi

		DO WHILE(theta.LT.angle)
10       theta = theta+dtheta
			dxy = d12*SIN(theta)

         IF(x1.GE.x2) THEN
				x1 = x2+dxy*SIN(beta)
         ELSE
				x1 = x2-dxy*SIN(beta)
         ENDIF

         IF(y1.GE.y2) THEN
				y1 = y2+dxy*COS(beta)
         ELSE
				y1 = y2-dxy*COS(beta)
         ENDIF

         z1 = z2+d12*COS(theta)
			
			CALL ChkCmpt
			mroll = 1
			CALL FindXYZ
			IF(movlpcond.EQ.1) GOTO 10
 
			IF(z1.LT.bound) THEN
				nrec = 1
				GOTO 20
         ELSEIF(x1.LT.bound.OR.x1.GT.(size-bound).OR.
     1          y1.LT.bound.OR.y1.GT.(size-bound)) THEN
				IF(z1.LE.z2) THEN
					mroll = 0
				ENDIF
				nrec = 2
				GOTO 20
         ELSEIF(mroll.GE.2) THEN
				GOTO 20
         ENDIF
      ENDDO

		mroll = 0

20    END

c------------------------------------------------------------------------

      SUBROUTINE ChkCmpt

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 

      IF(x1.GE.(cmptln*i1x1).AND.i1x1.LT.nc) i1x1 = i1x1+1

		IF(x1.LT.(cmptln*(i1x1-1)).AND.i1x1.GE.2) i1x1 = i1x1-1

      IF(y1.GE.(cmptln*jy1).AND.jy1.LT.nc) jy1 = jy1+1

		IF(y1.LT.(cmptln*(jy1-1)).AND.jy1.GE.2) jy1 = jy1-1

      IF(z1.GE.(cmptln*kz1).AND.kz1.LE.ncz) kz1 = kz1+1

		IF(z1.LT.(cmptln*(kz1-1)).AND.kz1.GE.2) kz1 = kz1-1

      END

c------------------------------------------------------------------------

      SUBROUTINE SettleRlmt

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		CALL CalRlmtLoc(alpha,alpha2,alpha3)
		CALL ChkCmpt
      mroll = 0
		CALL FindXYZ
		angle = 0.5*pi

		IF(mroll.EQ.0.AND.(alpha+alpha3).LE.angle) THEN
			nrec = 1
      ELSEIF(mroll.EQ.0.AND.(alpha+alpha3).GT.angle) THEN
         x2 = x3
			y2 = y3
			z2 = z3
			r2 = r3
			i1x2 = i1x3
			jy2 = jy3
			kz2 = kz3
			icst2 = icst3
         mroll = 1
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE CalRlmtLoc(alpha,alpha2,alpha3)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 

		d23 = SQRT((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3))
		d12 = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
      d13 = SQRT((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3))
		dxy23 = SQRT(d23*d23-(z2-z3)*(z2-z3))
		z23tod23 = ABS(z2-z3)/d23
		IF(z23tod23.GT.1.0) z23tod23 = 1.0
		alpha = ASIN(z23tod23)
		ratio2 = 0.5*(d12*d12+d23*d23-d13*d13)/(d12*d23)
		ratio3 = 0.5*(d13*d13+d23*d23-d12*d12)/(d13*d23)
		IF(ratio2.GT.1.0) ratio2 = 1.0
		IF(ratio3.GT.1.0) ratio3 = 1.0
		alpha2 = ACOS(ratio2)
		alpha3 = ACOS(ratio3)
		ry23 = ABS(y2-y3)/dxy23
		IF(ABS(1.0-ry23).LE.eps)  THEN
			beta = 0.5*pi
      ELSE
         beta = ASIN(ry23)
      ENDIF
		
		ang = 0.5*pi
		IF(z2.GE.z3) THEN
			dz1 = d12*SIN(alpha2-alpha)
			z1 = z2+dz1
			dz13 = z1-z3
			dxy13 = SQRT(d13*d13-dz13*dz13)
			dx1 = dxy13*COS(beta)
			dy1 = dxy13*SIN(beta)
      ELSE
			dz1 = d13*SIN(alpha3-alpha)
         z1 = z3+dz1
			dz12 = z1-z2
			dxy12 = SQRT(d12*d12-dz12*dz12)
			IF((alpha+alpha2).LE.ang) THEN
			   dx1 = (dxy23-dxy12)*COS(beta)
			   dy1 = (dxy23-dxy12)*SIN(beta)
         ELSE
			   dx1 = (dxy23+dxy12)*COS(beta)
			   dy1 = (dxy23+dxy12)*SIN(beta)
         ENDIF
      ENDIF

		IF(z2.GE.z3) THEN
			IF((alpha+alpha3).GT.ang) THEN
				nsign = -1
         ELSE
				nsign = 1
         ENDIF
      ELSE 
		  	nsign = 1
      ENDIF
	
		IF(x2.GE.x3.AND.y2.GE.y3) THEN
			x1 = x3+nsign*dx1
			y1 = y3+nsign*dy1
      ENDIF

		IF(x2.GE.x3.AND.y2.LT.y3) THEN
			x1 = x3+nsign*dx1
			y1 = y3-nsign*dy1
      ENDIF

		IF(x2.LT.x3.AND.y2.GE.y3) THEN
			x1 = x3-nsign*dx1
			y1 = y3+nsign*dy1
      ENDIF

		IF(x2.LT.x3.AND.y2.LT.y3) THEN
			x1 = x3-nsign*dx1
			y1 = y3-nsign*dy1
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE RollOn2

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 

c.....Find the highest point of the rolling track by using subroutine
c     CalRlmtLoc

      xtmp = x1
      ytmp = y1
      ztmp = z1
      
		CALL CalRlmtLoc(alpha,alpha2,alpha3)
		
		x1max = x1
		y1max = y1
		z1max = z1
      x1 = xtmp
      y1 = ytmp
      z1 = ztmp
		
c.....Find the center and radius of the track

      d13 = SQRT((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3))
      d12 = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
      rtrack = d13*SIN(alpha3)
      rn = d13*COS(alpha3)
      rm = d12*COS(alpha2)
		ratio = rm/rn
		angle = 0.5*pi
      IF(ABS(alpha3-angle).LE.tolrnce) THEN
			x0 = x3
			y0 = y3
			z0 = z3
		ELSE
		   x0 = (x2+ratio*x3)/(1+ratio)
		   y0 = (y2+ratio*y3)/(1+ratio)
		   z0 = (z2+ratio*z3)/(1+ratio)
      ENDIF

c.....Start rolling

		IF(ABS(z2-z3).LE.0.01) THEN
c		IF(ABS(z2-z3).LE.eps) THEN
			CALL VertRoll(x1max,y1max,z1max)
      ELSE
			CALL CalRLoc(x1max,y1max,z1max,x0,y0,z0,rtrack)
      ENDIF

      END

c------------------------------------------------------------------------

      SUBROUTINE VertRoll(x1max,y1max,z1max)

      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

      rtrack = z1max-z2
      IF(ABS(z1max-z1).LE.tolrnce) THEN
			IF(ABS(y2-y3).LE.tolrnce) THEN
				x1 = x1+1.0e-5
				IF(x1.GE.cmptln*i1x1.AND.i1x1.LT.nc) i1x1 = i1x1+1
				beta = 0.0
				sineT = (1.0e-5/rtrack)
				IF(sineT.GT.1.0) sineT = 1.0
				theta = ASIN(sineT)
				z1 = z2+rtrack*COS(theta)
         ELSE 
			   theta = 1.0e-6
				dxy = rtrack*SIN(theta)
				beta1 = ATAN(ABS(x3-x2)/ABS(y3-y2))
				beta = 0.5*pi-beta1
				dx = dxy*SIN(beta)
         ENDIF
      ELSE
	   	dx = ABS(x1-x1max)
	   	dxy = SQRT((x1-x1max)*(x1-x1max)+(y1-y1max)*(y1-y1max))
			dxtodxy = dx/dxy 
			IF(dxtodxy.GT.1.0) dxtodxy = 1.0
         beta = ASIN(dxtodxy)
			dxytort = dxy/rtrack
			IF(dxytort.GT.1.0) dxytort = 1.0
		   theta = ASIN(dxy/rtrack)
      ENDIF
		dtheta = 0.15
		angle = 0.5*pi

		DO WHILE(theta.LT.angle)
10       theta = theta+dtheta
			dxy = rtrack*SIN(theta)
         
			IF(x1.GE.x1max) THEN
				x1 = x1max+dxy*SIN(beta)
         ELSE
            x1 = x1max-dxy*SIN(beta)
         ENDIF

			IF(y1.GE.y1max) THEN
				y1 = y1max+dxy*COS(beta)
         ELSE
            y1 = y1max-dxy*COS(beta)
         ENDIF

			z1 = z2+rtrack*COS(theta)
			CALL ChkCmpt
			mroll = 2
			CALL FindXYZ
			IF(movlpcond.EQ.1) GOTO 10

			IF(z1.LT.bound) THEN
				nrec = 1
				GOTO 20
         ELSEIF(x1.LT.bound.OR.x1.GT.(size-bound).OR.
     1          y1.LT.bound.OR.y1.GT.(size-bound)) THEN
				nrec = 1
				GOTO 20
         ELSEIF(mroll.EQ.3) THEN
				GOTO 20
			ENDIF
      ENDDO

		mroll = 0

20    END

c------------------------------------------------------------------------

      SUBROUTINE CalRLoc(x1max,y1max,z1max,x0,y0,z0,rtrack)

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /data1/ x11, x12, x13, x11max, y11, y12, y13, y11max, 
     1               z11, z12, z11max, rxy, 
     2               sineA, cosineA, sineB, cosineB
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      
c.....Transform axes
 
		CALL TransFwd(x1max,y1max,z1max,x0,y0,z0,rtrack)

c.....Roll and check overlaps

		x13tort = x13/rtrack
		IF(x13tort.GT.1.0) x13tort = 1.0
		IF(x13tort.LT.(-1.0)) x13tort = -1.0
      theta = ACOS(x13/rtrack)
      IF(y13.GE.0.0) THEN
			msign = 1
      ELSE 
			msign = -1
      ENDIF

		dtheta = 0.15
		angle = 0.5*pi

		IF(theta.GT.angle) THEN
			mroll = 0
			GOTO 20
      ENDIF
	  
		DO WHILE(theta.LT.angle)
			x1tmp = x1
			y1tmp = y1
			z1tmp = z1
			i1x1tmp = i1x1
			jy1tmp = jy1
			kz1tmp = kz1
10       theta = theta+dtheta
			x13 = rtrack*COS(theta)
			y13 = msign*rtrack*SIN(theta)
			CALL TransBwd(x0,y0,z0)
			CALL ChkCmpt
			mroll = 2
			CALL FindXYZ
			IF(movlpcond.EQ.1) GOTO 10

			IF(z1.LT.bound) THEN
            nrec = 1
            GOTO 20
         ELSEIF(x1.LT.bound.OR.x1.GT.(size-bound).OR.
     1          y1.LT.bound.OR.y1.GT.(size-bound)) THEN
				nrec = 1 
            GOTO 20
         ELSEIF(mroll.EQ.3) THEN
				IF(dtheta.LE.eps) THEN
					x1 = x1tmp
					y1 = y1tmp
					z1 = z1tmp
					i1x1 = i1x1tmp
					jy1 = jy1tmp
					kz1 = kz1tmp
            ENDIF
				GOTO 20
         ENDIF
      ENDDO

		mroll = 1

		IF(z3.LT.z2) THEN
			x2 = x3
			y2 = y3
			z2 = z3
			r2 = r3
			i1x2 = i1x3
			jy2 = jy3
			kz2 = kz3
			icst2 = icst3
      ENDIF

20    END

c------------------------------------------------------------------------

      SUBROUTINE TransFwd(x1max,y1max,z1max,x0,y0,z0,rtrack)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /data1/ x11, x12, x13, x11max, y11, y12, y13, y11max, 
     1               z11, z12, z11max, rxy, 
     2               sineA, cosineA, sineB, cosineB
      
c.....Parallel translation

		x11max = x1max-x0
		y11max = y1max-y0
		z11max = z1max-z0

      x11 = x1-x0
      y11 = y1-y0
      z11 = z1-z0

c.....Rotation of x and y axes by angle A

		rxy = SQRT(rtrack*rtrack-z11max*z11max)
		sineA = y11max/rxy
		IF(sineA.GT.1.0) sineA = 1.0
		IF(sineA.LT.(-1.0)) sineA = -1.0
		IF(x11max.LT.0.0) THEN
			cosineA = -SQRT(1.0-sineA*sineA)
      ELSE
			cosineA = SQRT(1.0-sineA*sineA)
      ENDIF

		x12 = x11*cosineA+y11*sineA
		y12 = -x11*sineA+y11*cosineA
		z12 = z11

c.....Rotation of x and z axes by angle B

		sineB = z11max/rtrack
		IF(sineB.GT.1.0) sineB = 1.0
		cosineB = SQRT(1.0-sineB*sineB)

		x13 = x12*cosineB+z12*sineB
		y13 = y12

      END

c------------------------------------------------------------------------

      SUBROUTINE TransBwd(x0,y0,z0)

      IMPLICIT REAL(a-h,o-z)

		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /data1/ x11, x12, x13, x11max, y11, y12, y13, y11max, 
     1               z11, z12, z11max, rxy, 
     2               sineA, cosineA, sineB, cosineB

c.....Rotation of x and z axes by angle B

		x12 = x13*cosineB
		y12 = y13
		z12 = x13*sineB

c.....Rotation of x and y axes by angle A

		x11 = x12*cosineA-y12*sineA
		y11 = x12*sineA+y12*cosineA
		z11 = z12

c.....Parallel translation

		x1 = x11+x0
		y1 = y11+y0
		z1 = z11+z0

      END

c------------------------------------------------------------------------

      SUBROUTINE ChkCond

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		d12 = SQRT((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
		d13 = SQRT((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3))
		d14 = SQRT((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4))
		d23 = SQRT((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3))
		d24 = SQRT((x2-x4)*(x2-x4)+(y2-y4)*(y2-y4))
		d34 = SQRT((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4))

      ratio1 = 0.5*(d23*d23+d12*d12-d13*d13)/(d23*d12)
		ratio2 = 0.5*(d24*d24+d12*d12-d14*d14)/(d24*d12)
		ratio3 = 0.5*(d23*d23+d13*d13-d12*d12)/(d23*d13)
		ratio4 = 0.5*(d34*d34+d13*d13-d14*d14)/(d34*d13)
		ratio12 = 0.5*(d23*d23+d24*d24-d34*d34)/(d23*d24)
		ratio34 = 0.5*(d23*d23+d34*d34-d24*d24)/(d23*d34)
		
		IF( ratio1.GT.1.0) ratio1 = 1.0
		IF( ratio2.GT.1.0) ratio2 = 1.0
		IF( ratio3.GT.1.0) ratio3 = 1.0
		IF( ratio4.GT.1.0) ratio4 = 1.0
		IF( ratio12.GT.1.0) ratio12 = 1.0
		IF( ratio34.GT.1.0) ratio34 = 1.0
		
		IF(d14.LE.tolrnce) THEN
			nrec = 1
      ELSE
		   alpha1 = ACOS(ratio1)
         alpha2 = ACOS(ratio2)
         alpha3 = ACOS(ratio3)
         alpha4 = ACOS(ratio4)
         alpha12 = ACOS(ratio12)
         alpha34 = ACOS(ratio34)

			case1 = alpha2-alpha12
			case2 = alpha1-alpha12
			case3 = alpha3-alpha34

    		IF(case1.GT.tolrnce) THEN
			   x1 = x1+1.0e-4
				IF(x1.GE.cmptln*i1x1.AND.i1x1.LT.nc) i1x1 = i1x1+1
			   y1 = y1+1.0e-4
				IF(y1.GE.cmptln*jy1.AND.jy1.LT.nc) jy1 = jy1+1
				mroll = 0
    		ELSEIF(case2.GT.tolrnce) THEN
         	x3 = x4
	   		y3 = y4
            z3 = z4
            r3 = r4
            i1x3 = i1x4
            jy3 = jy4
            kz3 = kz4
            icst3 = icst4
            mroll = 2
         ELSEIF(case3.GT.tolrnce) THEN
            x2 = x4
            y2 = y4
            z2 = z4
            r2 = r4
            i1x2 = i1x4
            jy2 = jy4
            kz2 = kz4
            icst2 = icst4
            mroll = 2
         ELSE
				nrec = 1
         ENDIF
      ENDIF

		END

c------------------------------------------------------------------------

      SUBROUTINE TmpRec

		IMPLICIT REAL(a-h,o-z)

      COMMON /local/ intensity, lshift, move

		move = 10

		IF(intensity.EQ.0) GOTO 10
		IF(lshift.LE.move) CALL LocalShift
		IF(lshift.EQ.move) THEN
			CALL ChooseSite 
			GOTO 10
      ENDIF
		lshift = lshift+1

10		END

c------------------------------------------------------------------------

      SUBROUTINE LocalShift

		IMPLICIT REAL(a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /local/ intensity, lshift, move
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /shift/ xshift(30), yshift(30), zshift(30), i1xshift(30),
     1               jyshift(30), kzshift(30), indx(30)
			  
      xshift(lshift) = x1
      yshift(lshift) = y1
      zshift(lshift) = z1
      i1xshift(lshift) = i1x1
      jyshift(lshift) = jy1
      kzshift(lshift) = kz1

		IF(lshift.LT.move) THEN
			sample = size
			DO WHILE(sample.GE.(size-bound).OR.sample.LE.bound)
				sample = x1-rmax+ABS(2.0*rmax*RAN(nrndsd))
         ENDDO
			x1 = sample
			
			sample = size
			DO WHILE(sample.GE.(size-bound).OR.sample.LE.bound)
				sample = y1-rmax+ABS(2.0*rmax*RAN(nrndsd))
         ENDDO
			y1 = sample
		   
			CALL ChkCmpt
		   mroll = 0
		   nrec = 0
      ENDIF

		END 

c------------------------------------------------------------------------

      SUBROUTINE ChooseSite

      IMPLICIT REAL (a-h,o-z)

      COMMON /local/ intensity, lshift, move
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
		COMMON /shift/ xshift(30), yshift(30), zshift(30), i1xshift(30),
     1               jyshift(30), kzshift(30), indx(30)
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

		DO j = 1, move
			indx(j) = j
      ENDDO

		l = move/2+1
		ir = move

10    CONTINUE
		IF(l.GT.1) THEN
			l = l-1
			indxt = indx(l)
			q = zshift(indxt)
      ELSE
			indxt = indx(ir)
			q = zshift(indxt)
			indx(ir) = indx(1)
			ir = ir-1
			IF(ir.EQ.1) THEN
				indx(1) = indxt
				GOTO 30
         ENDIF
      ENDIF
		i = l
		j = l+l

20    IF(j.LE.ir) THEN
	      IF(j.LT.ir) THEN
		 	   IF(zshift(indx(j)).LT.zshift(indx(j+1))) j = j+1
         ENDIF
		   IF(q.LT.zshift(indx(j))) THEN
		      indx(i) = indx(j)
		      i = j
			   j = j+j
         ELSE
		 	   j = ir+1
         ENDIF
		   GOTO 20
      ENDIF
		indx(i) = indxt
		GOTO 10

30		i = move+1-intensity
		x1 = xshift(indx(i))
		y1 = yshift(indx(i))
		z1 = zshift(indx(i))
		i1x1 = i1xshift(indx(i))
		jy1 = jyshift(indx(i))
		kz1 = kzshift(indx(i))
      nrec = 1

      END

c------------------------------------------------------------------------

      SUBROUTINE Relocate
      
		IMPLICIT REAL(a-h,o-z)

		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
		COMMON /data/ icst2, icst3, icst4, i1x, i1x1, i1x2, i1x3, i1x4,
     1              jy, jy1, jy2, jy3, jy4, kz, kz1, kz2, kz3, kz4,
     2              r1, r2, r3, r4, x1, x2, x3, x4, 
     3              y1, y2, y3, y4, z1, z2, z3, z4 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate
      COMMON /sign/ lsign

      nchalf = INT(0.5*nc)
      IF(mcount.EQ.40.OR.ntshift.GE.60) THEN
         r1 = 0.9*r1
         mcount = 0
			ntshift = 0
      ENDIF
			
10    IF(mcount.EQ.0) THEN
			IF(i1x1.LE.nchalf) THEN
				lsign = 1
         ELSE
				lsign = -1
         ENDIF
      ENDIF

		IF(jy1.LT.nc) THEN
			y1 = y1+cmptln
			jy1 = jy1+1
      ELSE
			x1 = x1+cmptln*lsign
			i1x1 = i1x1+lsign
			y1 = y1+cmptln-size
			jy1 = 1
      ENDIF

		IF(i1x1.EQ.0.OR.i1x1.EQ.(nc+1)) THEN
			x1 = x1-2.0*cmptln*FLOAT(lsign)
			i1x1 = i1x1-2*lsign
         r1 = 0.9*r1
         mcount = 0
			ntshift = 0
		   GOTO 10
      ENDIF

		mcount = mcount+1
		mroll = 0
		nrec = 0
		ntshift = ntshift+1

      END

c------------------------------------------------------------------------

		SUBROUTINE FinalSet

      INCLUDE "paramt"
      IMPLICIT REAL(a-h,o-z)

      COMMON /dist/ alphar, betar, cnfint, rmin, rmax, sdnrm
		COMMON /param/ bound, cmptln, dtheta, eps, height, nc, npart, 
     1               npart_kz, pi, porosity, rlmt, nrndsd, size, 
     2               theta, tolrnce, ncz
      COMMON /local/ intensity, lshift, move
      COMMON /origdim/ r(nmax), rc(ncmax,ncmax,icmax),
     1                 x(nmax), xc(ncmax,ncmax,icmax), 
     2                 y(nmax), yc(ncmax,ncmax,icmax), 
     3                 z(nmax), zc(ncmax,ncmax,icmax), 
     4                 icnt(ncmax,ncmax) 
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /index/ mcount, movlpcond, mpBC, mroll, nrec, ntshift,
     1               niterate

c.....Calculate the final porosity, number of spheres in the whole
c     domain
		
		Vpack = 0.0
	   npart = 0
		ravg = EXP(alphar)
		rmin = ravg
		rmax = ravg

		DO kz = 1, ncz-2
			DO jy = 1, nc
				DO i1x = 1, nc
					ncntst = icntst(i1x,jy,kz)
					DO icst = 1, ncntst
	               r1 = rcst(i1x,jy,kz,icst)
					   Vpack = Vpack+4.0*pi*r1*r1*r1/3.0
						IF(r1.LT.rmin) rmin = r1
						IF(r1.GT.rmax) rmax = r1
						npart = npart+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO

		avgVc = Vpack/(ncz-2.0)

		Vctop = 0.0

		DO kz = ncz-1, ncz
         DO jy = 1, nc
            DO i1x = 1, nc
			   	ncntst = icntst(i1x,jy,kz)
			   	DO icst = 1, ncntst
                  r1 = rcst(i1x,jy,kz,icst)
		   		   Vctop = Vctop+4.0*pi*r1*r1*r1/3.0
						IF(r1.LT.rmin) rmin = r1
						IF(r1.GT.rmax) rmax = r1
						npart = npart+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
	
		Hctop = cmptln*Vctop/avgVc
		IF(Hctop.GE.2*cmptln) Hctop = 2*cmptln
		trueH = (ncz-2)*cmptln+Hctop
      porosity = 1.0-avgVc/(size*size*cmptln)

      WRITE(7,2000) porosity, size, size, trueH, cmptln,
     1              npart, nc, ncz, mpBC
	c	WRITE(7,2010) alphar,betar, rmax, rmin, ravg, rlmt
	c	WRITE(7,2015) niterate, intensity

c.....Write out the final extracted packing data!

		OPEN (5, FILE = 'subdim.par',STATUS = 'OLD')
		READ(5,1000) ncl, nch, nczl, nczh
		REWIND 9
		npart = 0
		
		DO kz = nczl, nczh
			DO jy = ncl, nch
				DO i1x = ncl, nch
					ncntst = icntst(i1x,jy,kz)
					ic = 0
					DO icst = 1, ncntst
                  x1 = xcst(i1x,jy,kz,icst)
                  y1 = ycst(i1x,jy,kz,icst)
                  z1 = zcst(i1x,jy,kz,icst)
                  r1 = rcst(i1x,jy,kz,icst)
					   IF(x1.GE.tolrnce.AND.x1.LE.(size+tolrnce).AND.
     1               y1.GE.tolrnce.AND.y1.LE.(size+tolrnce).AND.
     2               z1.GE.tolrnce.AND.z1.LE.(height+tolrnce)) THEN
							npart = npart+1
							ic = ic+1
							x(ic) = x1
							y(ic) = y1
							z(ic) = z1
							r(ic) = r1
                  ENDIF
               ENDDO
				c	WRITE(9,2020) i1x, jy, kz, ic
					DO iic = 1, ic
					c	WRITE(9,2030) x(iic), y(iic), z(iic), r(iic)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

	c	WRITE(7,2016) npart, ncl, nch, nczl, nczh

1000  FORMAT(/,4i7)

2000  FORMAT('                     Porosity of the sample: ',e13.5,/,
     1       '                              Sample length: ',e13.5,/,
     2       '                               Sample width: ',e13.5,/,
     3       '                               Sample depth: ',e13.5,/,
     4       '                                Cell length: ',e13.5,/,
     5       '                        Number of particles: ',i13,/,
     6       '       Number of cells in each x and y axis: ',i13,/,
     7       '                  Number of cells in z axis: ',i13,/,
     8       '                Periodical BC (Yes=1, No=0): ',i13)
2010  FORMAT('                           Average of ln(r): ',e13.5,/,
     1       '                             S.D. for ln(r): ',e13.5,/,
     2       '                        Maximum samplable r: ',e13.5,/,
     3       '                        Minimum samplable r: ',e13.5,/,
     4       '                        Average samplable r: ',e13.5,/,
     5       '                Small particle radius limit: ',e13.5)
2015  FORMAT('                        Number of iteration: ',i13,/,
     1       '                             Porosity index: ',i13)
2016  FORMAT('              Number of particles extracted: ',i13,/,
     1       '        Starting cell (in x or y direction): ',i13,/, 
     2       '          Ending cell (in x or y direction): ',i13,/,
     3       '             Starting cell (in z direction): ',i13,/, 
     4       '               Ending cell (in z direction): ',i13)
2020  FORMAT(4i7)
2030  FORMAT(4(e13.5,2x))

		END

c------------------------------------------------------------------------
c------------------------------------------------------------------------

         



		







