      program backtrm
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
c                     Gaussian Back Transformation
c                     ****************************
c
c INPUT/OUTPUT Parameters:
c
c   datafl                 the input data
c   ivr                    normal variable of interest
c   tmin,tmax              trimming limits: acceptable tmin<=z<tmax
c   outfl                  the output file for results
c   transfl                the input file with transformation table
c   zmin,zmax              the minimum and maximum data value
c   ltail,ltpar            lower tail option and parameter
c   utail,utpar            upper tail option and parameter
c
c
c
c PROGRAM NOTES:
c
c  1. ltail, utail options: 1=linear interpolation, 2=power model
c     interpolation, and 4=hyperbolic model interpolation (only for
c     upper tail)
c
c
c
c EXTERNAL REFERENCES:
c
c   gcum     Inverse of Gaussian cumulative distribution function
c   locate   locate a position in an array
c   powint   power law interpolation
c
c
c
c Original: C.V. Deutsch                               Date: August 1990
c-----------------------------------------------------------------------
      parameter(MAXDAT=500, EPSLON=0.00001, MV=20, VERSION=1.200)

      real      vr(MAXDAT),vrg(MAXDAT),var(MV),ltpar,utpar
      character datafl*40,outfl*40,transfl*40,str*40
      integer   ltail,utail
      logical   testfl
      data      lin/1/,lout/2/
c
c Note VERSION number before anything else:
c
      write(*,9999) VERSION
 9999 format(/' BACKTR Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      write(*,*) 'Which parameter file do you want to use?'
      read (*,'(a40)') str
      if(str(1:1).eq.' ')str='backtr.par                              '
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
 1    read(lin,'(a40)',end=98) str
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a40)',err=98) datafl
      write(*,*) ' data file:            ',datafl
      read(lin,*,err=98)       ivr
      write(*,*) ' column of normal:     ',ivr
      read(lin,*,err=98)       tmin,tmax
      write(*,*) ' trimming limits:      ',tmin,tmax
      read(lin,'(a40)',err=98) outfl
      write(*,*) ' output file:          ',outfl
      read(lin,'(a40)',err=98) transfl
      write(*,*) ' transformation file:  ',transfl
      read(lin,*,err=98)       zmin,zmax  
      write(*,*) ' data limits:          ',zmin,zmax
      read(lin,*,err=98)       ltail,ltpar
      write(*,*) ' lower tail option:    ',ltail,ltpar
      read(lin,*,err=98)       utail,utpar
      write(*,*) ' upper tail option:    ',utail,utpar
      close(lin)
c
c Check for error situation:
c
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
c
c Read in the transformation table:
c
      inquire(file=transfl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR transformation file does not exist'
            stop
      else
            open(lin,file=transfl,status='OLD')
            nt = 0
 2          read(lin,*,end=3) (var(i),i=1,2)
            nt = nt + 1
            if(nt.ge.MAXDAT) then
                  write(*,*) 'ERROR too many data in transformation '
                  write(*,*) '      have room for ',MAXDAT
                  stop
            endif
            vr(nt)  = var(1)
            vrg(nt) = var(2)
            if(nt.gt.1) then
            if(vr(nt) .lt.vr(nt-1) .or.
     +         vrg(nt).lt.vrg(nt-1)) then
                  write(*,*) 'ERROR transformation table must be '
                  write(*,*) '      monotonic increasing! '
                  write(*,*) '      Inconsistency at line ',nt
                  stop
            endif
            endif
            go to 2
 3          close(lin)
      endif
c
c Check for error situation:
c
      if(utail.eq.4.and.vr(nt).le.0.0) then
            write(*,*) 'ERROR can not use hyperbolic tail with '
            write(*,*) '      negative values! - see manual '
            stop
      endif
      if(zmin.gt.vr(1)) then
            write(*,*) 'ERROR zmin should be no larger than the first'
            write(*,*) '      entry in the transformation table '
            write(*,*) '      zmin = ',zmin,' vr1 ',vr(1)
            stop
      endif
      if(zmax.lt.vr(nt)) then
            write(*,*) 'ERROR zmax should be no less than the last'
            write(*,*) '      entry in the transformation table '
            write(*,*) '      zmax = ',zmax,' vrnt ',vr(nt)
            stop
      endif
c
c Now read through the data:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
c
c The data file exists so open the file and read in the header and
c write a header on the output file:
c
      open(lin, file=datafl,status='OLD')
      open(lout,file=outfl, status='UNKNOWN')
      read(lin,'(a40)',err=99) str
      write(lout,'(a40)')      str
      read(lin,*,err=99)       nvari
      write(lout,'(i2)')       nvari+1
      do 6 i=1,nvari
            read(lin,'(a40)',err=99) str
            write(lout,'(a40)')      str
 6    continue
      str='Back Transform                          '
      write(lout,'(a40)') str
c
c Read and write all the data until the end of the file:
c
 7    read(lin,*,end=8,err=99) (var(i),i=1,nvari)
      bac = backtr(var(ivr),nt,vr,vrg,zmin,zmax,ltail,ltpar,utail,utpar)
      write(lout,'(10(f8.4,1x))') (var(i),i=1,nvari),bac
      go to 7
 8    close(lin)
      close(lout)
      stop 'Succesful finish'
c
c Error in an Input File Somewhere:
c
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      real function backtr(vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar,
     +                     utail,utpar)
c-----------------------------------------------------------------------
c
c           Back Transform Univariate Data from Normal Scores
c           *************************************************
c
c This subroutine backtransforms a standard normal deviate from a
c specified back transform table and option for the tails of the
c distribution.  Call once with "first" set to true then set to false
c unless one of the options for the tail changes.
c
c
c
c INPUT VARIABLES:
c
c   vrgs             normal score value to be back transformed
c   nt               number of values in the back transform tbale
c   vr(nd)           original data values that were transformed
c   vrg(nd)          the corresponding transformed values
c   zmin,zmax        limits possibly used for linear or power model
c   ltail            option to handle values less than vrg(1):
c   ltpar            parameter required for option ltail
c   utail            option to handle values greater than vrg(nt):
c   utpar            parameter required for option utail
c   first            set to true the first time called
c
c
c
c Original: C.V. Deutsch                                     August 1990
c-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)
      dimension vr(nt),vrg(nt)
      real      ltpar,utpar,lambda
      integer   ltail,utail
c
c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
c
      if(vrgs.le.vrg(1)) then
            backtr = vr(1)
            cdflo  = gcum(vrg(1))
            cdfbt  = gcum(vrgs)
            if(ltail.eq.1) then
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,1.0)
            else if(ltail.eq.2) then
                  cpow   = 1.0 / ltpar
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,cpow)
            endif
c
c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
c
      else if(vrgs.ge.vrg(nt)) then
            backtr = vr(nt)
            cdfhi  = gcum(vrg(nt))
            cdfbt  = gcum(vrgs)
            if(utail.eq.1) then
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,1.0)
            else if(utail.eq.2) then
                  cpow   = 1.0 / utpar
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,cpow)
            else if(utail.eq.4) then
                  lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                  backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
            endif
      else
c
c Value within the transformation table:
c
            call locate(vrg,nt,vrgs,j)
            j = max0(min0((nt-1),j),1)
            backtr = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,1.0)
      endif
      return
      end



      real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
c-----------------------------------------------------------------------
c
c Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
c                 for a value of x and a power pow.
c
c-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)
      powint = ylow +        (yhigh-ylow)*
     +        (((xval-xlow)/amax1(EPSLON,(xhigh-xlow)))**pow)
      return
      end



      subroutine locate(xx,n,x,j)
c-----------------------------------------------------------------------
c
c Given an array "xx" of length "n", and given a value "x", this routine
c returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
c must be monotonic, either increasing or decreasing.  j=0 or j=n is
c returned to indicate that x is out of range.
c
c Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
c-----------------------------------------------------------------------
      dimension xx(n)
c
c Initialize lower and upper methods:
c
      jl = 0
      ju = n
c
c If we are not done then compute a midpoint:
c
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
c
c Replace the lower or upper limit with the midpoint:
c
            if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
c
c Return with the array index:
c
      j = jl
      return
      end



      real function gcum(x)
c-----------------------------------------------------------------------
c
c Evaluate the standard normal cdf given a normal deviate x.  gcum is
c the area under a unit normal curve to the left of x.  The results are
c accurate only to about 5 decimal places.
c
c
c-----------------------------------------------------------------------
      z = x
      if(z.lt.0.) z = -z
      t    = 1./(1.+ 0.2316419*z)
      gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 +
     +       t*(-1.821255978 + t*1.330274429))))
      e2   = 0.
c
c  6 standard deviations out gets treated as infinity:
c
      if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
      gcum = 1.0- e2 * gcum
      if(x.ge.0.) return
      gcum = 1.0 - gcum
      return
      end
