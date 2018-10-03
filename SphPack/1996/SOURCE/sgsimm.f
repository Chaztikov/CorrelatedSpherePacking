      program sgsimm
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1992 Stanford Center for Reservoir Forecasting.  All   %
C rights reserved.  Distributed with: C.V. Deutsch and A.G. Journel.   %
C ``GSLIB: Geostatistical Software Library and User's Guide,'' Oxford  %
C University Press, New York, 1992.                                    %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
c
c           Conditional Simulation of a 3-D Rectangular Grid
c           ************************************************
c
c This is a template driver program for GSLIB's "sgsim" subroutine.
c 3-D realizations of a Gaussian process with a given autocovariance
c model and conditional to input Gaussian data are created.  The
c conditional simulation is achieved by the sequential approach.  A
c normal score transform is performed if necessary.
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example sgsim.par)
c and should contain the following information:
c
c   -  Name of the data file (GEOEAS format)
c   -  column numbers for x, y, z, wt, and variable (if no wt set to 0)
c   -  Minimum acceptable value (used to flag missing values)
c   -  Are the data Gaussian (1=yes,0=no)?
c   -  A file for the transformation table (used if transformation
c        is necessary).
c   -  An output file (may be overwritten)
c   -  The debugging level (integer code - larger means more)
c   -  A file for the debugging output
c   -  Random Number Seed
c   -  Kriging Type: 0=simple kriging, 1=ordinary kriging
c   -  The number of simulations
c   -  X grid definition (number, minimum, size): nx,xmn,xsiz
c   -  Y grid definition (number, minimum, size): ny,ymn,ysiz
c   -  Z grid definition (number, minimum, size): nz,zmn,zsiz
c   -  The data can be assimilated to grid nodes which means that only
c        a spiral search is necessary and a covariance table look-up
c        has to be used.  The default is a two part search.
c   -  minimum and maximum number of data values to use in simulation.
c   -  maximum data per octant (0 -> no octant search)
c   -  maximum search radius
c   -  parameters for anisotropic search distance calculation:
c        a) azimuth principal direction (measured clockwise from Y).
c        b) dip of principal direction (measured negative down from X).
c        c) a third rotation of the two minor directions about the
c             principal direction. This angle acts counterclocwise
c             when looking in the principal direction.
c        d) search radius in minor direction at 90 degrees from the
c             principal direction divided by the principal radius.
c        e) search radius in minor direction at 90 degrees vertical from
c             the principal direction divided by the principal radius.
c   -  maximum number of previously simuulated grid nodes to use
c   -  Variogram Definition:  number of structures(nst), and nugget
c   -  the next "nst*2" lines require:
c      First line:
c        a) an integer code for variogram (1=sph,2=exp,3=gaus,4=pow)
c        b) "a" parameter (range except for power model)
c        b) "c" parameter (contribution except for power model).
c      Second line:
c        a) azimuth principal direction (measured clockwise from Y).
c        b) dip of principal direction (measured negative down from X).
c        c) a third rotation of the two minor directions about the
c             principal direction. This angle acts counterclocwise
c             when looking in the principal direction.
c        Two anisotropy factors are required to complete the definition
c        of the geometric anisotropy of each nested structure:
c        d) radius in minor direction at 90 degrees from the principal
c             direction divided by the principal radius.
c        e) radius in minor direction at 90 degrees vertical from the
c             principal direction divided by the principal radius.
c
c
c
c The output file will be a GEOEAS file containing the simulated values
c The file is ordered by x,y,z, and then simulation (i.e., x cycles
c fastest, then y, then z, then simulation number).  The values will be
c backtransformed to the original data values if a normal scores
c transform was performed (i.e., you don't even have to know that the
c data were transformed).
c
c
c
c Original: C.V. Deutsch                               Date: August 1990
c-----------------------------------------------------------------------
      include  'sgsim.inc'
      dimension vrg(MAXDAT),tmp(MAXDAT)
      character transfl*40,datafl*40,outfl*40,dbgfl*40
c
c Input/Output units used:
c
      lin  = 1
      lout = 2
      ldbg = 3
      ltr  = 4
      lint = 7
c
c Read the Parameter File and the Data:
c
      call readparm(tmin,tmax,iwt,itrans,transfl,datafl,outfl,dbgfl)
c
c If required: create normal scores and write the transformation table
c out for the backtransformation and for future use:
c
      if(itrans.eq.0) then
            call nscore(nd,vr,tmin,tmax,iwt,wt,tmp,vrg,ierror)
            do 1 i=1,nd
                  tmp(i) = vr(i)
                  vr(i)  = vrg(i)
 1          continue
c
c Sort the transformation table and write it out for future reference:
c
            call sortem(1,nd,vrg,1,tmp,c,d,e,f,g,h)
            open(ltr,file=transfl,status='UNKNOWN')
            do 2 i=1,nd
                  write(ltr,'(2f8.4)') tmp(i),vrg(i)
 2          continue
            close(ltr)
      endif
c
c Call sgsim for the simulation:
c
      call sgsim
c
c If required: read the transformation table
c
      if(itrans.eq.0) then
            open(ltr,file=transfl,status='OLD')
            do 3 i=1,nd
                  read(ltr,'(2f8.4)') vr(i),vrg(i)
 3          continue
            close(ltr)
      endif
c
c Now read everything from the interim file, back transform the data
c if required and write it to the permanent output file.  The first five
c elements of the "tmp" array will be used to keep track of the mean
c and variance of the normal scores and backtransformed points:
c
      rewind lint
 4    read(lint,*,end=6) simval
c
c If this node has been simulated then accumulate statistics and back
c transform if necessary:
c
      if(simval.ne.UNEST) then
            if(itrans.eq.0) then
                  bac = backtr(simval,nd,vr,vrg,zmin,zmax,ltail,
     +                         ltpar,utail,utpar)
                  simval = bac
            endif
      endif
      write(lout,'(f8.4)') simval
      go to 4
 6    close(lint)
      close(lout)
c
c Finished:
c
      write(ldbg,101) outfl,dbgfl,transfl
      write(*   ,101) outfl,dbgfl,transfl
 101  format(/' Finished SGSIM: simulated results in ',a40,/,
     +        '                 debugging output in  ',a40,/,
     +        '                 transformation table ',a40,//)
      stop
      end



      subroutine readparm(tmin,tmax,iwt,itrans,transfl,datafl,outfl,
     +                    dbgfl)
c-----------------------------------------------------------------------
c
c                  Initialization and Read Parameters
c                  **********************************
c
c The input parameters and data are read in from their files. Some quick
c error checking is performed and the statistics of all the variables
c being considered are written to standard output.
c
c
c
c Original: C.V. Deutsch                                 Date: July 1990
c-----------------------------------------------------------------------
      include  'sgsim.inc'
      parameter(MV=20)
      real      var(MV)
      character transfl*40,datafl*40,outfl*40,dbgfl*40,str*40,title*80
      logical   testfl
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' SGSIM Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      write(*,*) 'Which parameter file do you want to use?'
      read (*,'(a40)') str
      if(str(1:1).eq.' ')str='sgsim.par                               '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a40)',err=98) datafl
      write(*,*) ' data file:           ',datafl
      read(lin,*,err=98)       ixl,iyl,izl,ivrl,iwt
      write(*,*) ' input columns:       ',ixl,iyl,izl,ivrl,iwt
      read(lin,*,err=98)       tmin,tmax
      write(*,*) ' trimming limits:     ',tmin,tmax
      read(lin,*,err=98)       itrans
      write(*,*) ' transformation flag: ',itrans
      read(lin,'(a40)',err=98) transfl
      write(*,*) ' transformation file: ',transfl
      read(lin,*,err=98)       zmin,zmax
      write(*,*) ' data limits (tails): ',zmin,zmax
      read(lin,*,err=98)       ltail,ltpar
      write(*,*) ' lower tail:          ',ltail,ltpar
      read(lin,*,err=98)       utail,utpar
      write(*,*) ' upper tail:          ',utail,utpar
      read(lin,'(a40)',err=98) outfl
      write(*,*) ' output file:         ',outfl
      read(lin,*,err=98)       idbg
      write(*,*) ' debugging level:     ',idbg
      read(lin,'(a40)',err=98) dbgfl
      write(*,*) ' debugging file:      ',dbgfl
      read(lin,*,err=98)       seed
      write(*,*) ' random number seed:  ',seed
      read(lin,*,err=98)       ktype
      write(*,*) ' kriging type:        ',ktype
      read(lin,*,err=98)       nsim
      write(*,*) ' number of simulations',nsim
      read(lin,*,err=98)       nx,xmn,xsiz
      write(*,*) ' X grid specification:',nx,xmn,xsiz
      read(lin,*,err=98)       ny,ymn,ysiz
      write(*,*) ' Y grid specification:',ny,ymn,ysiz
      read(lin,*,err=98)       nz,zmn,zsiz
      write(*,*) ' Z grid specification:',nz,zmn,zsiz
      read(lin,*,err=98)       sstrat
      write(*,*) ' search strategy flag:',sstrat
      read(lin,*,err=98)       noct
      write(*,*) ' number of octants:   ',noct
      read(lin,*,err=98)       radius
      write(*,*) ' search radius:       ',radius
      radsqd = radius * radius
      read(lin,*,err=98)       sang1,sang2,sang3,sanis1,sanis2
      write(*,*) ' search anisotropy: ',sang1,sang2,sang3,sanis1,sanis2
      read(lin,*,err=98)       ndmin,ndmax
      write(*,*) ' min and max data:    ',ndmin,ndmax
      read(lin,*,err=98)       nodmax
      write(*,*) ' maximum previous node',nodmax
      read(lin,*,err=98)       nst,c0
      sill = c0
      write(*,*) ' nst, c0:             ',nst,c0
      if(nst.le.0) then
            write(*,9997) nst
 9997       format(' nst must be at least 1, it has been set to ',i4,/,
     +             ' The c or a values can be set to zero')
            stop
      endif
      do 2 i=1,nst
            read(lin,*,err=98) it(i),aa(i),cc(i)
            write(*,*) ' it,aa,cc:            ',it(i),aa(i),cc(i)
            sill = sill + cc(i)
            if(it(i).eq.4) then
                  write(*,*) ' A power model is NOT allowed '
                  write(*,*) ' Choose a different model and re start '
                  stop
            endif
            read(lin,*,err=98) ang1(i),ang2(i),ang3(i),anis1(i),anis2(i)
            write(*,*) ' ang123,anis12:       ',
     +                   ang1(i),ang2(i),ang3(i),anis1(i),anis2(i)
 2    continue
      write(*,*)
      close(lin)
c
c Warn the user if the sill is different than 1.0:
c
      if(sill.gt.(1.0+EPSLON).or.sill.lt.(1.0-EPSLON)) then
            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
            write(*,*) '        the sill = ',sill
      end if
c
c Make sure the flag for no weights is set properly for "nscore":
c
      if(iwt.le.0) iwt = 0
c
c Perform some quick error checking:
c
      if(nx.gt.MAXX) stop 'nx is too big - modify .inc file'
      if(ny.gt.MAXY) stop 'ny is too big - modify .inc file'
      if(nz.gt.MAXZ) stop 'nz is too big - modify .inc file'
      if(ltail.ne.1.and.ltail.ne.2) then
            write(*,*) 'ERROR invalid lower tail option ',ltail
            write(*,*) '      only allow 1 or 2 - see manual '
            stop
      endif
      if(utail.ne.1.and.utail.ne.2.and.utail.ne.4) then
            write(*,*) 'ERROR invalid upper tail option ',ltail
            write(*,*) '      only allow 1,2 or 4 - see manual '
            stop
      endif
      if(utail.eq.4.and.utpar.lt.1.0) then
            write(*,*) 'ERROR invalid power for hyperbolic tail',utpar
            write(*,*) '      must be greater than 1.0!'
            stop
      endif
      if(ltail.eq.2.and.ltpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',ltpar
            write(*,*) '      must be greater than 0.0!'
            stop
      endif
      if(utail.eq.2.and.utpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',utpar
            write(*,*) '      must be greater than 0.0!'
            stop
      endif
      if(ixl.le.0) nx = 1
      if(iyl.le.0) ny = 1
      if(izl.le.0) nz = 1
c
c Check debugging level:
c
      if(idbg.ge.3.and.(nsim*nx*ny*nz).gt.100) then
            write(*,100) idbg,(nsim*nx*ny*nz)
 100        format(/' The debugging level is set high: ',i3,' and there'
     +            ,/' are many points to simulate:     ',i9
     +           ,//' Do you want to proceed? (y/n)')
            read (*,'(a40)') str
            if(str(1:1).ne.'y'.and.str(1:1).ne.'Y') stop
      end if
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      title = 'SGSIM SIMULATIONS:                      '//
     +        '                                        '
      nd = 0
      av = 0.0
      ss = 0.0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'WARNING data file ',datafl,' does not exist!'
            write(*,*) '   - Hope your intention was to create an ',
     +                       'unconditional simulation'
            write(*,*) '   - Resetting ndmin,  ndmax to 0 '
            write(*,*) '   - Resetting itrans, sstrat to 1 '
            ndmin  = 0
            ndmax  = 0
            itrans = 1
            sstrat = 1
      else
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
            open(lin,file=datafl,status='OLD')
            read(lin,'(a60)',err=99) title(21:80)
            read(lin,*,err=99)       nvari
            do 3 i=1,nvari
 3          read(lin,'(a40)',err=99) str
c
c Read all the data until the end of the file:
c
            nd = 0
 4          read(lin,*,end=8,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 4
            nd = nd + 1
            if(nd.gt.MAXDAT) then
                  write(*,*) ' ERROR exceeded MAXDAT - check inc file'
                  stop
            end if
            if(ixl.le.0) then
                  x(nd) = xmn
            else
                  x(nd) = var(ixl)
            endif
            if(iyl.le.0) then
                  y(nd) = ymn
            else
                  y(nd) = var(iyl)
            endif
            if(izl.le.0) then
                  z(nd) = zmn
            else
                  z(nd) = var(izl)
            endif
            vr(nd) = var(ivrl)
            if(iwt.gt.0) wt(nd) = var(iwt)
            av     = av + var(ivrl)
            ss     = ss + var(ivrl)*var(ivrl)
            go to 4
 8          close(lin)
c
c Note if there are less than 10 data:
c
            if(nd.le.10.and.itrans.eq.0) write(*,101) 
 101        format(/'WARNING too few data for a decent transform - ',/,
     +        '         continuing anyway.')
c
c Compute the averages and variances as an error check for the user:
c
            av = av / amax1(real(nd),1.0)
            ss =(ss / amax1(real(nd),1.0)) - av * av
      endif
c
c Open the debugging and output files:
c
      open(ldbg,file=dbgfl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
      write(ldbg,*) 'Data for SGSIM: Variable number ',ivrl
      write(ldbg,*) '  Number of acceptable data  = ',nd
      write(ldbg,*) '  Equal Weighted Average     = ',av
      write(ldbg,*) '  Equal Weighted Variance    = ',ss
c
c Write a header on the output file and return:
c
      write(lout,102) title
 102  format(a80,/,'1',/,'Simulated Value')
      return
c
c Error in an Input File Somewhere:
c
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end
