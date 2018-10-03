c------------------------------------------------------------------------
c------------------------------------------------------------------------

		PROGRAM ViewPore

      PARAMETER(max=200)
      IMPLICIT REAL(a-h,o-z)
		INTEGER*1 ipix

		DIMENSION ipix(max,max,max)

		OPEN (3, FILE = 'domain.par', STATUS = 'OLD')
      OPEN (4, FILE = 'pore.out', STATUS = 'OLD')
      OPEN (1, FILE = 'avsdim.par', STATUS = 'OLD')
    	OPEN (9, FILE = 'avspore.out')

c.....Read in number of pixels along each direction in pore.out

		READ(3,1000) nclow, nchigh

c.....Read the bounds for avs view

      READ(1,1010) nxl, nxh, nyl, nyh, nzl, nzh

c.....Read in pore data

		READ(4,1015) np

		DO nz = 1, np
			WRITE(*,*) 'Start reading nz = ', nz
			DO ny = 1, np
				DO nx = 1, np
   				READ(4,1020) ipix(nx,ny,nz)
c   				write(*,*) ipix(nx,ny,nz)
            ENDDO
         ENDDO
      ENDDO

c.....Output extracted pore data
		
		DO nz = nzl, nzh
			DO ny = nyl, nyh
				DO nx = nxl, nxh
					ipixel = ipix(nx,ny,nz)

c.....Assign pixel color for visualization

					IF(ipixel.EQ.1) ipixel = 255
					WRITE(9,2000) nx, ny, nz, ipixel
            ENDDO
         ENDDO
      ENDDO

1000  FORMAT(/,2(i7))
1010  FORMAT(/,6(i6))
1015  FORMAT(45x,i13)
1020  FORMAT(21x,i7)

2000  FORMAT(4(i7))

		END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	  
