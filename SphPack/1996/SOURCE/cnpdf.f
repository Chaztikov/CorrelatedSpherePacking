c------------------------------------------------------------------------
c------------------------------------------------------------------------
      
      PROGRAM CNdistribution

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max

		CALL DataIn

		CALL ComputeCN

		CALL CNPDF

		IF(nsignal.EQ.1) THEN
			CALL SizeCntribute
      ENDIF

		END

c------------------------------------------------------------------------

      SUBROUTINE DataIn

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 

      OPEN (3, FILE = 'cnparam.in', STATUS = 'OLD')
      OPEN (4, FILE = 'cndata.in', STATUS = 'OLD')
      OPEN (7, FILE = 'cnpdf.out')

c.....Read in  parameters

		READ(3,1000)
		READ(3,1010) nc, ncz, ncnz, eps
		READ(3,1015) nsignal, ncontact_l, ncontact_h

c.....Read in packing data

		DO kz = 1, ncz
			DO jy = 1, nc
            DO ix = 1, nc
					READ(4,1020) ncntst
					DO icst = 1, ncntst
						READ(4,1030) x1, y1, z1, r1
						xcst(ix,jy,kz,icst) = x1
						ycst(ix,jy,kz,icst) = y1
						zcst(ix,jy,kz,icst) = z1
						rcst(ix,jy,kz,icst) = r1
               ENDDO
					icntst(ix,jy,kz) = ncntst
            ENDDO
         ENDDO
      ENDDO

1000  FORMAT()
1010  FORMAT(/,3i7,e13.5)
1015  FORMAT(/,3i11)
1020  FORMAT(21x,i7)
1030  FORMAT(4(e13.5,2x))

		END
 
c------------------------------------------------------------------------

      SUBROUTINE ComputeCN

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /cnarray/ r(nmax), ncn(nmax), intcnt(50), pd(50)

		npart = 0

		DO kz = 2, ncnz-1
			DO jy = 2, nc-1
				DO ix = 2, nc-1
					ncntst = icntst(ix,jy,kz)
					DO icst = 1, ncntst
						x1 = xcst(ix,jy,kz,icst)
						y1 = ycst(ix,jy,kz,icst)
						z1 = zcst(ix,jy,kz,icst)
						r1 = rcst(ix,jy,kz,icst)
                  npart = npart+1
						r(npart) = r1
						CALL CalcCN
               ENDDO
            ENDDO
         ENDDO
      ENDDO

		END

c------------------------------------------------------------------------

      SUBROUTINE CalcCN

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max
		COMMON /sortdim/ rcst(ncmax,ncmax,ncmax,icmax),
     1                 xcst(ncmax,ncmax,ncmax,icmax), 
     2                 ycst(ncmax,ncmax,ncmax,icmax), 
     3                 zcst(ncmax,ncmax,ncmax,icmax), 
     4                 icntst(ncmax,ncmax,ncmax) 
      COMMON /cnarray/ r(nmax), ncn(nmax), intcnt(50), pd(50)

		ncontact = 0
 
		DO k = kz-1, kz+1
			DO j = jy-1, jy+1
				DO i = ix-1, ix+1
					num = icntst(i,j,k)
					DO ic = 1, num
						IF(.NOT.(i.EQ.ix.AND.j.EQ.jy.AND.k.EQ.kz.AND.
     1                     ic.EQ.icst)) THEN
							xr = xcst(i,j,k,ic)
							yr = ycst(i,j,k,ic)
							zr = zcst(i,j,k,ic)
							rr = rcst(i,j,k,ic)
                     d = SQRT((x1-xr)*(x1-xr)+(y1-yr)*(y1-yr)
     1                        +(z1-zr)*(z1-zr))
			            sumr = r1+rr
			            IF((d-sumr).LE.eps) THEN
								ncontact = ncontact+1
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

		ncn(npart) = ncontact

		END

c------------------------------------------------------------------------

      SUBROUTINE CNPDF

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max
      COMMON /cnarray/ r(nmax), ncn(nmax), intcnt(50), pd(50)

c.....Initialize the counting

		DO i = ncontact_l, ncontact_h
			intcnt(i) = 0
      ENDDO

c.....Count in each division

		DO i = ncontact_l, ncontact_h
			DO ipart = 1, npart
				IF(ncn(ipart).EQ.i) THEN
					intcnt(i) = intcnt(i)+1
            ENDIF
         ENDDO
      ENDDO

c.....Compute the mean CN

		sum1 = 0.0
		sum2 = 0.0
		DO i = ncontact_l, ncontact_h
		   sum1 = sum1+FLOAT(i)*FLOAT(intcnt(i))
			sum2 = sum2+FLOAT(intcnt(i))
      ENDDO

		cnmean = sum1/sum2

c.....Normalize the PDF

		DO i = ncontact_l, ncontact_h
			pd(i) = FLOAT(intcnt(i))/FLOAT(npart)
      ENDDO
      
		
		WRITE(7,2000) cnmean 
		IF(nsignal.EQ.0) THEN
			WRITE(7,2010)
		   DO i = ncontact_l, ncontact_h
		   	WRITE(7,2020) i, pd(i)
         ENDDO
      ENDIF

2000  FORMAT('            CN mean = ',f7.3)
2010  FORMAT('  CN ',2x,'  PDF  ',/,
     1       '-----',2x,'-------')
2020  FORMAT(i5,2x,f7.3)

		END

c------------------------------------------------------------------------

      SUBROUTINE SizeCntribute

      INCLUDE "cnparam"
      IMPLICIT REAL(a-h,o-z)

		COMMON /param/ ix, jy, kz, icst, x1, y1,z1,r1, 
     1               cnmean, nc, ncz, ncnz, npart, eps,
     2               nsignal, ncontact_l, ncontact_h, 
     3               rmin, rmax, intcnt_max
      COMMON /cnarray/ r(nmax), ncn(nmax), intcnt(50), pd(50)
      
		DIMENSION mcnfr(50,10), pdfr(50,10)

c.....Find the rmin and rmax

		rmin = 0.05
		rsubmax = 0.55
      rmax = 1.55

		nseg = 5
		delr = (rsubmax-rmin)/FLOAT(nseg)

c.....Initialize the counting array

		DO i = ncontact_l, ncontact_h
			DO k = 1, nseg+1
				mcnfr(i,k) = 0
         ENDDO
      ENDDO

		DO k = 1, nseg
			DO ipart = 1, npart
				IF(r(ipart).GE.(rmin+delr*(k-1)).AND.
     1         r(ipart).LT.(rmin+delr*k)) THEN
               DO i = ncontact_l, ncontact_h
						IF(ncn(ipart).EQ.i) THEN
                     mcnfr(i,k) = mcnfr(i,k)+1
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

		DO ipart = 1, npart
    		IF(r(ipart).GE.rsubmax) THEN
            DO i = ncontact_l, ncontact_h
	            IF(ncn(ipart).EQ.i) THEN
				   	mcnfr(i,nseg+1) = mcnfr(i,nseg+1)+1
               ENDIF
            ENDDO
         ENDIF
      ENDDO

c.....Normalize the CN distributions

		DO k = 1, nseg+1
			DO i = ncontact_l, ncontact_h
				pdfr(i,k) = FLOAT(mcnfr(i,k))/FLOAT(npart)
         ENDDO
      ENDDO

c.....Write out the CN distributions

      WRITE(7,2000)
		DO k = 1, nseg
			r_low = rmin+delr*(k-1)
			r_high = rmin+delr*k
			WRITE(7,2010) k, r_low, r_high
      ENDDO

      WRITE(7,2010) nseg+1, rsubmax, rmax


		WRITE(7,2020)
		DO i = ncontact_l, ncontact_h
			WRITE(7,2030) i,pdfr(i,1), pdfr(i,2), pdfr(i,3),
     1                   pdfr(i,4), pdfr(i,5), pdfr(i,6), pd(i)
      ENDDO

2000  FORMAT(11x,'k',3x,'r_low',4x,'r_high',/,
     1       10x,'---',x,'-------',2x,'-------')
2010  FORMAT(10x,i2,2x,2(f7.3,2x))
2020  FORMAT(5x,'CN',4x,'pdf(1)',3x,'pdf(2)',3x,'pdf(3)',
     1       3x,'pdf(4)',3x,'pdf(5)',3x,'pdf(6)',4x,'pdf',/,
     2       3x,'------',7(x,'--------'))
2030  FORMAT(5x,i2,2x,7(x,f8.3))

		END

c------------------------------------------------------------------------
c------------------------------------------------------------------------













