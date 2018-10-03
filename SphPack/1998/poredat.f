c------------------------------------------------------------------------
c------------------------------------------------------------------------

		PROGRAM PoreConstruct

		OPEN (3, FILE = 'pack.par', STATUS = 'OLD')
		OPEN (4, FILE = 'pack.out', STATUS = 'OLD')
		OPEN (7, FILE = 'pore.par')
		OPEN (8, FILE = 'pore.out')
		OPEN (9, FILE = 'pore2d.out')

		CALL DataIn

		CALL BuildPore

		CALL ParamOut

		END

c------------------------------------------------------------------------

		SUBROUTINE DataIn

		INCLUDE "pore.inc"
		IMPLICIT REAL(a-h,o-z)

      COMMON /param/ cmptln, eps, rpix, ncell, nc, ncl, nch, 
     1               ncz, nczl, nczh, 
     2               npix_x, npix_z, npore, porosity
		COMMON /data/ ix, jy, kz, mx1, mx2, my1, my2, mz1, mz2
		COMMON /array/rcst(ncmax,ncmax,ncmax,icmax),
     1              xcst(ncmax,ncmax,ncmax,icmax),
     2              ycst(ncmax,ncmax,ncmax,icmax),
     3              zcst(ncmax,ncmax,ncmax,icmax),
     4              icntst(ncmax,ncmax,ncmax)

c.....Read in parameters

		READ(3,1000) cmptln, betar, rmax, rmin, r0, npart, 
     1             ncl, nch, nczl, nczh

c.....Read in packing data

10    READ(4,1010) ix, jy, kz, icpart
      DO icst = 1, icpart
			READ(4,1020) x1, y1, z1, r1
			xcst(ix,jy,kz,icst) = x1
			ycst(ix,jy,kz,icst) = y1
			zcst(ix,jy,kz,icst) = z1
			rcst(ix,jy,kz,icst) = r1
      ENDDO
      icntst(ix,jy,kz) = icpart
		IF(ix.EQ.nch.AND.jy.EQ.nch.AND.kz.EQ.nczh) GOTO 20
		GOTO 10

20    eps = 1.0e-6
  
		nc = nch-ncl+1
		ncz = nczh-nczl+1

		size = cmptln*FLOAT(nc)
		depth = cmptln*FLOAT(ncz)
		temp = rmax+r0

c.....Pixel resolution is inputed here

		ncell = INT(24.5*cmptln/temp)

		pixsize = cmptln/ncell
		rpix = 0.5*pixsize
		npix_x = ncell*nc
		npix_z = ncell*ncz

c.....Write out packing domain, pore info

      WRITE(*,*) 'Number of pixels in each packing cell: ', ncell

		WRITE(7,2000) size, depth, nc, ncz, npart,
     1              rmax, rmin, r0, betar
		WRITE(7,2010) pixsize, ncell, npix_x, npix_z

1000  FORMAT(////,45x,e13.5,
     1       //////,45x,e13.5,
     2       /,45x,e13.5,
     3       /,45x,e13.5,
     4       /,45x,e13.5,
     5       ////,45x,i13,
     6       /,45x,i13,
     7       /,45x,i13,
     8       /,45x,i13,
     9       /,45x,i13)
1010  FORMAT(4i7)
1020  FORMAT(4(e13.5,2x))

2000  FORMAT('                                Domain Size: ',e13.5,/,
     1       '                               Domain depth: ',e13.5,/,
     2       '             Number of cells in x or y axis: ',i13,/
     3       '                  Number of cells in z axis: ',i13,/
     4       '                        Number of particles: ',i13,/
     5       '                        Maximum samplable r: ',e13.5,/,
     6       '                        Minimum samplable r: ',e13.5,/,
     7       '                        Average samplable r: ',e13.5,/,
     8       '                             S.D. for ln(r): ',e13.5)
2010  FORMAT('                               Pixel length: ',e13.5,/,
     1       '    Number of pixels in each axis in a cell: ',i13,/
     2       '            Number of pixels in x or y axis: ',i13,/
     3       '                 Number of pixels in z axis: ',i13)

		END

c------------------------------------------------------------------------

      SUBROUTINE BuildPore

		INCLUDE "pore.inc"
		IMPLICIT REAL(a-h,o-z)
		INTEGER*1 ipix

      COMMON /param/ cmptln, eps, rpix, ncell, nc, ncl, nch, 
     1               ncz, nczl, nczh, 
     2               npix_x, npix_z, npore, porosity
		COMMON /data/ ix, jy, kz, mx1, mx2, my1, my2, mz1, mz2
		COMMON /array/rcst(ncmax,ncmax,ncmax,icmax),
     1              xcst(ncmax,ncmax,ncmax,icmax),
     2              ycst(ncmax,ncmax,ncmax,icmax),
     3              zcst(ncmax,ncmax,ncmax,icmax),
     4              icntst(ncmax,ncmax,ncmax)
      COMMON /index/ ipix(maxpix,maxpix,maxcell)

		COMMON /subpore/nclow, nchigh, layer, nthpix

c.....Open file to get the size (cube!) data for pore output

      OPEN (5, FILE = 'domain.par',STATUS = 'OLD')

c.....nclow and nchigh are cell bounds for x, y and z directions
c     nflag=2 means output 2D data, layer specifies kz value for 
c     2D image, nthpix specifies the nth pixel along z dir in layer 
c     cell

		READ(5,1000) nclow, nchigh
		READ(5,1010) nflag, layer, nthpix

		npore = 0

   	npixel = ((nchigh-nclow)+1)*ncell

		DO kz = nczl, nczh
         IF(kz.EQ.nczl) THEN
				mz1 = nczl
				mz2 = nczl+1
         ELSEIF(kz.EQ.nczh) THEN
				mz1 = nczh-1
				mz2 = nczh
         ELSE
				mz1 = kz-1
				mz2 = kz+1
         ENDIF
 
			CALL LayerProcess

         IF(kz.GE.nclow.AND.kz.LE.nchigh) THEN
            npix_low = (nclow-1)*ncell+1
            npix_high = nchigh*ncell
			   DO numpz = 1, ncell
					indxpz= numpz+(kz-1)*ncell
			   	DO numpy = npix_low, npix_high
cj			   		DO numpx = npix_low, npix_high
cj						WRITE(8,*)  
cj    2                          ipix(numpx,numpy,numpz)
                  WRITE(8,1006) (ipix(numpx,numpy,numpz), 
     1                       numpx = npix_low, npix_high) 
c					   ENDDO
               ENDDO
            ENDDO
         ENDIF
1006   format(128i2)
			IF(nflag.EQ.2.AND.kz.EQ.layer) THEN
				numpz = nthpix 
				DO numpy = npix_low, npix_high
					DO numpx = npix_low, npix_high
						WRITE(9,2020) (numpx-npix_low+1), (numpy-npix_low+1),
     1                          ipix(numpx,numpy,numpz)
					ENDDO
            ENDDO
         ENDIF

      ENDDO

1000  FORMAT(/,2(i7))
1010  FORMAT(/,3(i7))

2000  FORMAT('         Number of pixels in each direction: ',i13)
2010  FORMAT(4(i7))
2020  FORMAT(3(i7))

      END

c------------------------------------------------------------------------

      SUBROUTINE LayerProcess

		IMPLICIT REAL(a-h,o-z)

      COMMON /param/ cmptln, eps, rpix, ncell, nc, ncl, nch, 
     1               ncz, nczl, nczh, 
     2               npix_x, npix_z, npore, porosity
		COMMON /data/ ix, jy, kz, mx1, mx2, my1, my2, mz1, mz2

		DO jy = ncl, nch
			DO ix = ncl, nch
				IF(ix.EQ.ncl) THEN
					mx1 = ncl
					mx2 = ncl+1
            ELSEIF(ix.EQ.nch) THEN
					mx1 = nch-1
					mx2 = nch
            ELSE
					mx1 = ix-1
					mx2 = ix+1
            ENDIF
				IF(jy.EQ.ncl) THEN
					my1 = ncl
					my2 = ncl+1
            ELSEIF(jy.EQ.nch) THEN
					my1 = nch-1
					my2 = nch
            ELSE
					my1 = jy-1
					my2 = jy+1
            ENDIF

				CALL CellCount

         ENDDO
      ENDDO

		END

c------------------------------------------------------------------------

      SUBROUTINE CellCount

		INCLUDE "pore.inc"
		IMPLICIT REAL(a-h,o-z)
		INTEGER*1 ipix

      COMMON /param/ cmptln, eps, rpix, ncell, nc, ncl, nch, 
     1               ncz, nczl, nczh, 
     2               npix_x, npix_z, npore, porosity
		COMMON /data/ ix, jy, kz, mx1, mx2, my1, my2, mz1, mz2
		COMMON /array/rcst(ncmax,ncmax,ncmax,icmax),
     1              xcst(ncmax,ncmax,ncmax,icmax),
     2              ycst(ncmax,ncmax,ncmax,icmax),
     3              zcst(ncmax,ncmax,ncmax,icmax),
     4              icntst(ncmax,ncmax,ncmax)
      COMMON /index/ ipix(maxpix,maxpix,maxcell)

      DO numpz = 1, ncell
			DO numpy = 1, ncell
				DO numpx = 1, ncell
               xpix = 2.0*rpix*numpx-rpix+(ix-1)*cmptln
               ypix = 2.0*rpix*numpy-rpix+(jy-1)*cmptln
               zpix = 2.0*rpix*numpz-rpix+(kz-1)*cmptln
					indxpx = numpx+(ix-1)*ncell
					indxpy = numpy+(jy-1)*ncell
					indxpz = numpz
					sumovlp = 0.0
               DO icz = mz1, mz2
						DO icy = my1, my2
							DO icx = mx1, mx2
								icpart = icntst(icx,icy,icz)
								DO ic = 1, icpart
									x1 = xcst(icx,icy,icz,ic)
									y1 = ycst(icx,icy,icz,ic)
									z1 = zcst(icx,icy,icz,ic)
									r1 = rcst(icx,icy,icz,ic)
									dpix1 = SQRT((xpix-x1)*(xpix-x1)
     1                                +(ypix-y1)*(ypix-y1)
     2                                +(zpix-z1)*(zpix-z1))
									sumpix1 = r1+rpix
									dovlp = sumpix1-dpix1
									IF(dovlp.GT.eps) THEN
										sumovlp = sumovlp+dovlp
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

c.....Pixel value must be defined for visualization, ipix = 255 red

					IF(sumovlp.GT.(1.05*rpix)) THEN
c.....solid
						ipix(indxpx,indxpy,indxpz) = 0
               ELSE
c.....pore
						ipix(indxpx,indxpy,indxpz) = 1
						npore = npore+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

		END

c------------------------------------------------------------------------

      SUBROUTINE ParamOut

		IMPLICIT REAL(a-h,o-z)

      COMMON /param/ cmptln, eps, rpix, ncell, nc, ncl, nch, 
     1               ncz, nczl, nczh, 
     2               npix_x, npix_z, npore, porosity
		COMMON /subpore/nclow, nchigh, layer, nthpix

      porosity = FLOAT(npore)/FLOAT(npix_x*npix_x*npix_z)
      npixel = ((nchigh-nclow)+1)*ncell

		WRITE(7,2000) porosity, nclow, nchigh, npixel,
     1              layer, nthpix

2000  FORMAT('                     Porosity of the sample: ',e13.5,//,
     1       '              Pore data extracted from cell: ',i13,/,
     2       '                Pore data extracted to cell: ',i13,/,
     3       '         Number of pixels in each direction: ',i13,/,
     4       '            2-D pore data extracted at cell: ',i13,/,
     5       '      2-D pore data extracted at x_th pixel: ',i13)

		END

c------------------------------------------------------------------------
c------------------------------------------------------------------------









