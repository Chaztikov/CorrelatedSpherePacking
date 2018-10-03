      program histplt
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
c                       Histogram of Some Data
c                       **********************
c
c This program generates a PostScript file with a histogram and summary
c statistics.
c
c
c
c PROGRAM NOTES:
c
c 1. The program is executed with no command line arguments.  The user
c    will be prompted for the name of a parameter file.  The parameter 
c    file is described in the documentation (see the example hist.par).
c
c 2. If data are <tmin , >=tmax , or the weight is <EPSLON the data will
c    not be considered.
c
c
c
c INPUT/OUTPUT Parameters:
c
c   datafl      the input data
c   ivr,iwt     columns for the variable and the weight (0 if non)
c   outfl       output file with histogram
c   tmin,tmax   trimming limts (acceptable: >= tmin and < tmax)
c   hmin,hmax   plotting limts (will choose automatically if hmax<hmin)
c   ncl         the number of classes
c   ilog        1=log scale, 0=arithmetic
c   title       title for output PostScript file
c
c
c
c Original: C.V. Deutsch                                 Date: June 1990
c-----------------------------------------------------------------------
      parameter (MAXDAT=130000, MV=20, EPSLON=1.0e-20, VERSION=1.201)

      character  datafl*40,outfl*40,title*40,str*80
      real       ar1(MAXDAT),ar2(MAXDAT),var(20)
      logical    testfl
c
c Common variables for Postscript Plotting:
c
      common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                xmax,ymin,ymax
c
c Hardwire many of the plot parameters:
c
      data       lin/1/,lpsout/1/,pscl/0.24/,pxmin/0.0/,pxmax/288.0/
     +           pymin/0.0/,pymax/216.0/,xmin/-10.0/,xmax/60.0/,
     +           ymin/-10.0/,ymax/60.0/,hpxmin/1.0/,hpxmax/59.5/,
     +           hpymin/0.0/,hpymax/58.0/
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' HISTPLT Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      write(*,*) 'Which parameter file do you want to use?'
      read (*,'(a20)') str(1:20)
      if(str(1:1).eq.' ') str(1:20) = 'histplt.par         '
      inquire(file=str(1:20),exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            stop
      endif
      open(lin,file=str(1:20),status='old')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a40)',err=98) datafl
      write(*,*) ' data file:        ',datafl
      read(lin,*,err=98)       ivr,iwt
      write(*,*) ' columns:          ',ivr,iwt
      read(lin,'(a40)',err=98) outfl
      write(*,*) ' output file:      ',outfl
      read(lin,*,err=98)       tmin,tmax
      write(*,*) ' trimming  limits: ',tmin,tmax
      read(lin,*,err=98)       hmin,hmax
      write(*,*) ' histogram limits: ',hmin,hmax
      read(lin,*,err=98)       ncl
      write(*,*) ' number of classes:',ncl
      read(lin,*,err=98)       ilog
      write(*,*) ' log scale option: ',ilog
      read(lin,'(a40)',err=98) title
      write(*,*) ' title:            ',title
      close(lin)
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the data file does not exist,'
            write(*,*) '        check for the file and try again  '
            stop
      endif
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
      open(lin,file=datafl)
      read(lin,'(a)',err=99) str
      read(lin,*,err=99)     nvari
      do 6 i=1,nvari
 6    read(lin,'()',err=99)
      if(ivr.gt.nvari) then
            write(*,*) ' ivr is greater than the number in data file!'
            stop
      end if
c
c Read as much data as possible:
c
      vrmin = 1.0e21
      vrmax =-1.0e21
      nd = 0
      nt = 0
 7    read(lin,*,end=8,err=99) (var(j),j=1,nvari)
c
c Trim this data?
c
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax) then
            nt = nt + 1
            go to 7
      endif
      if(iwt.ge.1) then
            if(var(iwt).le.EPSLON) then
                  nt = nt + 1
                  go to 7
            endif
      endif
c
c Accept this data:
c
      nd = nd + 1
      if(nd.gt.MAXDAT) then
            write(*,*) 'ERROR: exceeded available storage'
            write(*,*) '       have ',MAXDAT,' available'
            stop
      endif
      ar1(nd) = var(ivr)
      vrmin = amin1(ar1(nd),vrmin)
      vrmax = amax1(ar1(nd),vrmax)
      if(iwt.ge.1) then
            ar2(nd) = var(iwt)
      else
            ar2(nd) = 1.0
      endif
c
c Go back for another data:
c
      go to 7
 8    close(lin)
c
c Assign the defaults if necessary:
c
      if(hmax.le.hmin) then
            hmin = vrmin
            hmax = vrmax
      endif
      if(iwt.ge.1) then
            iwt = 1
      else
            iwt = 0
      endif
c
c Plot the histogram of the first and second variable:
c
      open(lpsout,file=outfl)
      call plthist(ar1,nd,nt,ar2,ilog,hmin,hmax,ncl,
     +             hpxmin,hpxmax,hpymin,hpymax,title)
      close(lpsout)
c
c Finished:
c
      stop
 98   stop 'ERROR in parameters somewhere'
 99   stop 'ERROR in data somewhere'
      end



        subroutine plthist(ar1,nd,nt,ar2,ilog,hmin,hmax,ncl,
     +                     hpxmin,hpxmax,hpymin,hpymax,title)
c-----------------------------------------------------------------------
c
c                      Histogram and Statistics
c                      ************************
c
c This program computes univariate stats of a data set and displays the
c  histogram and statistics on postscript output.
c
c PARAMETERS:
c      ar1     The variable for the histogram
c      nd      The number of data
c      nt      The number of data trimmed
c      ar2     A weight variable for each data
c      ilog    =1, then take base 10 logarithms, otherwise leave alone
c      hmin    the minimum value for the histogram (<0 - then default)
c      hmax    the maximum value for the histogram (<0 - then default)
c      ncl     the number of classes (<0 - then default)
c      hpxmin  minimum X to use for histogram plotting
c      hpxmax  maximum X to use for histogram plotting
c      hpymin  minimum Y to use for histogram plotting
c      hpymax  maximum Y to use for histogram plotting
c
c
c
c AUTHOR: Clayton Deutsch                               DATE: April 1990
c-----------------------------------------------------------------------
      parameter(BIGNUM=99999999.,EPSLON=0.00001,ZERO=0.0,MAXCLS=102)
      dimension ar1(*),ar2(*),ar3(MAXCLS)
      real      xh(5),yh(5)
      character str*80,title*40
c
c Common variables for Postscript Plotting:
c
      common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,wxmin,
     +                wxmax,wymin,wymax
c
c Add a header:
c
      write(lpsout,998)
 998  format('%!PS',/,'/m {moveto} def',/,
     +         '/l {lineto} def',/,'/r {rlineto} def',/,
     +         '/s {stroke} def',/,'/n {newpath} def',/,
     +         '/c {closepath} def',/,'2 setlinejoin ',/,
     +         '/rtext{ dup stringwidth pop -1 div 0 ',
     +         ' rmoveto show } def',/,'/ltext{show} def',/,
     +         '/ctext{ dup stringwidth pop -2 div 0 ',
     +         ' rmoveto show } def',/,
     +         '/bullet{ 6 0 360 arc c fill } def',/,
     +         '0.240000 0.240000 scale')
c
c Stop right away if there are too many classes:
c
      if((ncl+2).gt.MAXCLS) then
            write(*,*) 'ERROR: exceeded available number of classes'
            write(*,*) '       have ',MAXCLS-2,' available'
            stop
      endif
      if(nd.le.1) then
            write(*,*) 'ERROR: there is less than one datum ',nd
            write(*,*) '       check your data file'
            stop
      endif
c
c Get mean and total weight:
c
      xtwt = 0.0
      xmen = 0.0
      do 2 i=1,nd
            xmen = xmen + ar1(i)*ar2(i)
            xtwt = xtwt + ar2(i)
 2    continue
      xtwti = 1.0  / xtwt
      xmen  = xmen * xtwti
c
c Get the variance:
c
      xvar = 0.0
      do 3 i=1,nd
            xvar = xvar + (ar1(i)-xmen) * (ar1(i)-xmen) * ar2(i)
 3    continue
      xvar  = xvar * xtwti
c
c Infer some of the histogram parameters:
c
      if(ncl.lt.0) then
            ncl             = 100
            if(nd.le.500) ncl= 50
            if(nd.lt.200) ncl= 20
      endif
      ncl = min0((MAXCLS-2),(ncl+2))
      dcl = (hmax-hmin)/real(ncl-2)
      if(ilog.eq.1) then
            if(hmin.lt.0.0) hmin = xmin
            if(hmax.lt.0.0) hmax = xmax
            hmin = alog10(amax1(hmin,xmin,EPSLON))
            hmax = alog10(amax1(hmax,EPSLON))
            dcl  = (hmax-hmin)/real(ncl-2)
      endif
c
c Determine the histogram class structure:
c
      do 4 i=1,ncl
 4    ar3(i) = 0
      do 5 i=1,nd
            if(ilog.eq.1) then
                  art = alog10(amax1(ar1(i),EPSLON))
            else
                  art = ar1(i)
            endif
            wt1 = ar2(i) * xtwti
            j = ((art-hmin)/dcl)+2
            if(j.lt.1)   j = 1
            if(j.gt.ncl) j = ncl
            ar3(j) =ar3(j) + wt1
 5    continue
c
c Sort the Data in Ascending Order:
c
      call sortem(1,1,nd,ar1,1,ar2,c,d,e,f,g,h)
c
c Turn the weights into a cdf:
c
      oldcp = 0.0
      cp    = 0.0
      do 6 i=1,nd
            cp     = cp + ar2(i) * xtwti
            ar2(i) =(cp + oldcp) * 0.5
            oldcp  = cp
 6    continue
c
c Obtain the quantiles:
c
      call locate(ar2,nd,1,nd,0.25,i)
      if(i.eq.0) then
            xlqt = ar1(1)
      else if(i.eq.nd) then
            xlqt = ar1(nd)
      else
            xlqt = ar1(i) +      (ar1(i+1)-ar1(i)) *
     +                      (0.25-ar2(i))/(ar2(i+1)-ar2(i))
      endif
      call locate(ar2,nd,1,nd,0.50,i)
      if(i.eq.0) then
            xmed = ar1(1)
      else if(i.eq.nd) then
            xmed = ar1(nd)
      else
            xmed = ar1(i) +      (ar1(i+1)-ar1(i)) *
     +                      (0.50-ar2(i))/(ar2(i+1)-ar2(i))
      endif
      call locate(ar2,nd,1,nd,0.75,i)
      if(i.eq.0) then
            xuqt = ar1(1)
      else if(i.eq.nd) then
            xuqt = ar1(nd)
      else
            xuqt = ar1(i) +      (ar1(i+1)-ar1(i)) *
     +                      (0.75-ar2(i))/(ar2(i+1)-ar2(i))
      endif
c
c Calculate some statistics of the data:
c
      xmin = ar1(1)
      xmax = ar1(nd)
      if(xmin.lt.0.0.or.xmen.le.EPSLON) then
            xcvr = -1.0
      else
            xcvr = sqrt(xvar)/xmen
      endif
c
c Find the Maximum Class Frequency:
c
      xfrmx  = ar3(1)
      do 10 i =2,ncl
 10   xfrmx  = amax1(xfrmx,ar3(i))
c
c Set some scaling parameters:
c
      xhsmin = hmin
      xhsmax = hmin + real(ncl+1)*dcl
      yhsmin = 0.0
      yhsmax = 1.03*xfrmx
      xrange = hpxmax - hpxmin
      yrange = hpymax - hpymin
c
c Write the title and the boxes:
c
      ts   = 7.5
      xloc = hpxmin - 0.15*xrange
      yloc = hpymin + 0.50*yrange
      call pstext(xloc,yloc,9,'Frequency',ts,1,90.0,1)
      xloc = hpxmin + 0.50*xrange
      yloc = hpymin - 0.15*yrange
      call pstext(xloc,yloc,8,'Variable',ts,1,0.0,1)
c
c Scale and Draw the histogram:
c
      call scal(xhsmin,xhsmax,yhsmin,yhsmax,hpxmin,hpxmax,hpymin,
     +          hpymax)
      x1 = hmin-dcl
      x2 = hmin
      nh = 5
      do 11 i=1,ncl
            xh(1) = resc(xhsmin,xhsmax,hpxmin,hpxmax,x1)
            xh(2) = xh(1)
            xh(3) = resc(xhsmin,xhsmax,hpxmin,hpxmax,x2)
            xh(4) = xh(3)
            xh(5) = xh(1)
            yh(1) = hpymin
            yh(2) = resc(yhsmin,yhsmax,hpymin,hpymax,ar3(i))
            yh(3) = yh(2)
            yh(4) = yh(1)
            yh(5) = yh(1)
            call psfill(nh,xh,yh,0.3,0.9)
            x1=x1+dcl
            x2=x2+dcl
 11   continue
c
c Write the Statistics to the Output Device:
c
      x1 = hpxmin + 0.74*xrange
      x2 = hpxmin + 0.76*xrange
      yd = 0.044*yrange
      yd2= 0.060*yrange
      yy = hpymax
      ts =  7.0
      ls = 40
      yy = hpymax + 0.01*(hpxmax-hpxmin)
      call pstext(hpxmin,yy,40,title,8.0,3,0.0,0)
      str = 'Number of Data'
      call histtext(ls,str,x1,yy,ts,0.0,2)
      write(str,'(i8)') nd
      call histtext(ls,str,x2,yy,ts,yd2,0)
      if(nt.gt.0) then
            str = 'number trimmed'
            call histtext(ls,str,x1,yy,ts,0.0,2)
            write(str,'(i8)') nt
            call histtext(ls,str,x2,yy,ts,yd2,0)
      endif
      str = 'mean'
      call histtext(ls,str,x1,yy,ts,0.0,2)
      write(str,'(f14.4)') xmen
      call histtext(ls,str,x2,yy,ts,yd,0)
      
      str = 'std. dev.'
      call histtext(ls,str,x1,yy,ts,0.0,2)
      write(str,'(f14.4)') sqrt(xvar)
      call histtext(ls,str,x2,yy,ts,yd,0)
      str = 'coef. of var'
      call histtext(ls,str,x1,yy,ts,0.0,2)
      if(xcvr.lt.0.0) then
            str = ' undefined          '
      else
            write(str,'(f14.4)') xcvr
      endif
      call histtext(ls,str,x2,yy,ts,yd2,0)
      str = 'maximum'
      call histtext(ls,str,x1,yy,ts,0.,2)
      write(str,'(f14.4)') xmax
      call histtext(ls,str,x2,yy,ts,yd,0)
      str = 'upper quartile'
      call histtext(ls,str,x1,yy,ts,0.,2)
      write(str,'(f14.4)') xuqt
      call histtext(ls,str,x2,yy,ts,yd,0)
      str = 'median'
      call histtext(ls,str,x1,yy,ts,0.,2)
      write(str,'(f14.4)') xmed
      call histtext(ls,str,x2,yy,ts,yd,0)
      str = 'lower quartile'
      call histtext(ls,str,x1,yy,ts,0.,2)
      write(str,'(f14.4)') xlqt
       call histtext(ls,str,x2,yy,ts,yd,0)
      str = 'minimum'
      call histtext(ls,str,x1,yy,ts,0.,2)
      write(str,'(f14.4)') xmin
      call histtext(ls,str,x2,yy,ts,yd2,0)
      if(ilog.eq.1) then
            str = 'Histogram has Logarithmic Scale'
            call histtext(ls,str,x1,yy,ts,yd,2)
      endif
c
c Add a footer to the Postscript plot file:
c
      write(lpsout,999)
 999  format('%END OF HISTOGRAM',/,'4.166667 4.166667 scale',/,
     +       '%showpage')
c
c Return to Calling Program:
c
      return
      end



      subroutine histtext(lostr,str,x,y,ts,yd,iadj)
c-----------------------------------------------------------------------
c
c A space saving routine that either writes the text to the screen or
c writes the text to a postscript file and then decrements the y 
c position for the next time.
c
c INPUT:
c      lostr - length of the string
c      str   - the character string
c      x,y   - the location to position the text
c      ts    - text size
c      yd    - amount top decrement the y position
c      iadj  - adjustment: 0 left, 1 center, 2 right
c
c-----------------------------------------------------------------------
      character str*80
c
c Common variables for Postscript Plotting:
c
      common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                xmax,ymin,ymax
c
c Plot the text on the requested device:
c
      it = 1
      rt = 0.0
      call pstext(x,y,lostr,str,ts,it,rt,iadj)
      y = y - yd
      return
      end



      subroutine sortem(iflag,ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c
c This is a subroutine for sorting a real array in ascending or
c  descending order. This routine is a translation into fortran iv
c  of algorithm 271, quickersort, by r.s. scowen in collected
c  algorithms of the acm. consult this for further references.
c
c The method is used is that of continually splitting the array
c  into parts such that all elements of one part are less than
c  all elements of the other, with a third part in the middle
c  consisting of one element. An element with value t is chosen 
c  arbitrarily (here we choose the  middle element). i and j give
c  the lower and upper limits of the segment being split. after 
c  the split a value q will have been found such that a(q)=t and
c  a(l)<=t<=a(m) for all i<=l<q<m<=j. the program then performs
c  operations on the two segments (i,q-1) and (q+1,j) as follows.
c  the smaller segment is split and the position of the larger
c  segment is stored in the lt and ut arrays. if the segment to
c  be split contains two or fewer elements, it is sorted and
c  another segment is obtained from the lt and ut arrays. when
c  no more segments remain, the array is completely sorted.
c
c
c INPUT PARAMETERS:
c
c    iflag  =1 ascending order,=0 decending order
c
c    a      =the array, a portion of which has to be sorted.
c
c    ib,ie  =the positions in array a of the first and last elements of 
c            the portion to be sorted.
c
c    iperm  =0 no other array is permuted.
c           =1 array b is permuted according to array a
c           =2 arrays b,c are permuted.
c           =3 arrays b,c,d are permuted.
c           =4 arrays b,c,d,e are permuted.
c           =5 arrays b,c,d,e,f are permuted.
c           =6 arrays b,c,d,e,f,g are permuted.
c           =7 arrays b,c,d,e,f,g,h are permuted.
c           >7 no other array is permuted.
c
c    b,c,d,e,f,g,h  =arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c Taken from G. Verly's MGLIB                         Revised: June 1989
c-----------------------------------------------------------------------
      dimension a(1),b(1),c(1),d(1),e(1),f(1),g(1),h(1)
      real ta,tb,tc,td,te,tf,tg,th,xa,xb,xc,xd,xe,xf,xg,xh
      integer   lt(25),ut(25),i,j,k,m,p,q
c
c       the dimensions for lt and ut have to be at least log  n.
c                                                           2
c       17 was chosen to handle n<131,073.
c
      j=ie
      m=1
      i=ib
      iring=iperm+1
      if (iperm.gt.7) iring=1
c
c if this segment has more than two elements  we split it
c
 10      if (j-i-1) 100,90,15
c
c p is the position of an arbitrary element in the segment
c we choose the middle element. under certain circumstances
c it may be advantageous to choose p at random.
c
 15      p=(j+i)/2
      ta=a(p)
      a(p)=a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163      th=h(p)
      h(p)=h(i)
 162      tg=g(p)
      g(p)=g(i)
 161      tf=f(p)
      f(p)=f(i)
 16      te=e(p)
      e(p)=e(i)
 17      td=d(p)
      d(p)=d(i)
 18      tc=c(p)
      c(p)=c(i)
 19      tb=b(p)
      b(p)=b(i)
 21      continue
c
c Starting at the beginning of the segment, search for k
c such that a(k)>t
c
      q=j
      k=i
 20      k=k+1
c
c Modification so that we can sort ascending or decending:
c
      if (k.gt.q) go to 60
      if(iflag.eq.1) then
            if (a(k).le.ta) go to 20
      else
            if (a(k).ge.ta) go to 20
      endif
c
c Such an element has now been found now search for a q such that
c a(q)<t starting at the end of the segment.
c
 30      continue
      if(iflag.eq.1) then
            if (a(q).lt.ta) go to 40
      else
            if (a(q).gt.ta) go to 40
      endif
      q=q-1
      if (q.gt.k) go to 30
      go to 50
c
c a(q) has now been found. we interchange a(q) and a(k)
c
 40      xa=a(k)
      a(k)=a(q)
      a(q)=xa
      go to (45,44,43,42,41,411,412,413),iring
 413      xh=h(k)
      h(k)=h(q)
      h(q)=xh
 412      xg=g(k)
      g(k)=g(q)
      g(q)=xg
 411      xf=f(k)
      f(k)=f(q)
      f(q)=xf
 41      xe=e(k)
      e(k)=e(q)
      e(q)=xe
 42      xd=d(k)
      d(k)=d(q)
      d(q)=xd
 43      xc=c(k)
      c(k)=c(q)
      c(q)=xc
 44      xb=b(k)
      b(k)=b(q)
      b(q)=xb
 45      continue
c
c Update q and search for another pair to interchange:
c
      q=q-1
      go to 20
 50      q=k-1
 60      continue
c
c The upwards search has now met the downwards search:
c
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613      h(i)=h(q)
      h(q)=th
 612      g(i)=g(q)
      g(q)=tg
 611      f(i)=f(q)
      f(q)=tf
 61      e(i)=e(q)
      e(q)=te
 62      d(i)=d(q)
      d(q)=td
 63      c(i)=c(q)
      c(q)=tc
 64      b(i)=b(q)
      b(q)=tb
 65      continue
c
c The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
c store the position of the largest segment in lt and ut
c
      if (2*q.le.i+j) go to 70
        lt(m)=i
        ut(m)=q-1
        i=q+1
      go to 80
 70     lt(m)=q+1
        ut(m)=j
        j=q-1
c
c       update m and split the new smaller segment
 80   m=m+1
      go to 10
c
c       we arrive here if the segment has  two elements
c       we test to see if the segment is properly ordered
c       if not, we perform an interchange
c
c Modification:
c
 90      continue
      if(iflag.eq.1) then
            if (a(i).le.a(j)) go to 100
      else
            if (a(i).ge.a(j)) go to 100
      endif
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913  xh=h(i)
      h(i)=h(j)
      h(j)=xh
 912  xg=g(i)
      g(i)=g(j)
      g(j)=xg
 911  xf=f(i)
      f(i)=f(j)
      f(j)=xf
   91 xe=e(i)
      e(i)=e(j)
      e(j)=xe
   92 xd=d(i)
      d(i)=d(j)
      d(j)=xd
   93 xc=c(i)
      c(i)=c(j)
      c(j)=xc
   94 xb=b(i)
      b(i)=b(j)
      b(j)=xb
   95 continue
c
c       if lt and ut contain more segments to be sorted repeat process
 100  m=m-1
      if (m.le.0) goto 110
      i=lt(m)
      j=ut(m)
      go to 10
  110 continue
      return
      end



      subroutine pstext(xs,ys,lostr,str,tsiz,ifont,rot,iadj)
c-----------------------------------------------------------------------
c
c              Write Postscript Text commands to a file
c              ****************************************
c
c
c CALLING ARGUMENTS:
c
c  xs            starting value of x in the range xmin to xmax
c  ys            starting value of y in the range ymin to ymax
c  lostr      number of characters in str to print
c  str            the character string
c  tsiz            Text size in 1/72 of an inch
c  ifont      Font Number: See font number below
c  rot             Rotation Angle to post the text (default to 0.0)
c  iadj            Adjustment: 0=left adjusted, 1=centre, 2=right
c
c
c AUTHOR: C. Deutsch                                 DATE: November 1988
c-----------------------------------------------------------------------
      character str*80,fnnt(10)*32,line*132,size*4,part1*1,part2*7
c
c Common Block for Postscript Output Unit and Scaling:
c
        common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                  xmax,ymin,ymax
      save      fnnt,ifold,tsold,izip
c
c Preset 10 different fonts:
c
      data fnnt/'/Helvetica             findfont ',
     +            '/Helvetica-Bold        findfont ',
     +            '/Helvetica-BoldOblique findfont ',
     +            '/Times-Roman           findfont ',
     +            '/Times-Bold            findfont ',
     +            '/Times-Italic          findfont ',
     +            '/Times-BoldItalic      findfont ',
     +            '/Courier               findfont ',
     +            '/Courier-Bold          findfont ',
     +            '/Courier-BoldOblique   findfont '/
      data ifold/0/,tsold/0.0/,izip/0/
      part1 = '('
      part2 = ')  text'
c
c Remove leading and trailing blanks:
c
      lost = lostr
      do 1 i=1,lostr
            if(str(1:1).eq.' ') then
                  lost = lost - 1
                  do 2 j=1,lost
                        k = j + 1
                        str(j:j) = str(k:k)
 2                continue
            else
                  go to 3
            endif
 1    continue
 3    continue
      k = lost
      do 4 i=1,k
            ix = k - i + 1
            if(str(ix:ix).ne.' ') go to 5
            lost = lost - 1
 4    continue
 5    if(lost.le.0) return
c
c Create line to set the text size and type:
c
      if(ifont.ne.ifold.or.tsiz.ne.tsold) then
            isiz=int(tsiz/pscl)
            write(size,'(i4)') isiz
            line=fnnt(ifont)//size//' scalefont setfont'      
            write(lpsout,'(a)')line(1:54)
            ifold = ifont
            tsold = tsiz
      endif
c
c Set the correct adjustment:
c
      part2(3:3) = 'l'
      if(iadj.eq.1) part2(3:3) = 'c'
      if(iadj.eq.2) part2(3:3) = 'r'
c
c Write the lines and position to the Postscript file:
c                  
      ix = int((resc(xmin,xmax,pxmin,pxmax,xs))/pscl)
      iy = int((resc(ymin,ymax,pymin,pymax,ys))/pscl)
c
c Rotate if Necessary:
c
      line = part1//str(1:lost)//part2
      if(rot.ne.0.0) then
            irot = int(rot)
            write(lpsout,102) ix,iy
            write(lpsout,103) irot 
            write(lpsout,100) izip,izip
            write(lpsout,'(a)')line(1:lost+8)
            ix   = -1.0 * ix
            iy   = -1.0 * iy
            irot = -1.0 * irot
            write(lpsout,103) irot 
            write(lpsout,102) ix,iy
      else
c
c Just write out the text if no rotation:
c
            write(lpsout,100)ix,iy
            write(lpsout,'(a)')line(1:lost+8)
      endif
 100  format(i4.4,1x,i4.4,1x,'m')
 102  format(i5,1x,i5,1x,'translate')
 103  format(i5,1x,'rotate')
c
c Finished - Return to calling program:
c
      return
      end



      subroutine psfill(np,x,y,lwidt,gray)
c-----------------------------------------------------------------------
c
c
c CALLING ARGUMENTS:
c
c  x            X location of the center of the box
c  y            Y location of the center of the box
c  np           number of points
c  lwidt        The width of the line (1.0 = dark, 0.5 = light)
c  gray         the grayness of the fill area
c
c NOTES:
c
c  1. The pxmin,pxmax,.. variables are in the standard 1/72 inch 
c     resolution of the postscript page. If a different scale is 
c     going to be used in the printing set pscl to the scale.
c
c
c
c AUTHOR: C. Deutsch                                     DATE: June 1989
c-----------------------------------------------------------------------
      parameter(EPSLON=0.0001)
      real lwidt,x(1),y(1)
c
c Common Block for Postscript Output Unit and Scaling:
c
        common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                  xmax,ymin,ymax
c
c Change the line width:
c
      if(pscl.lt.0.01) pscl = 1.0
      width = lwidt/pscl
      write(lpsout,103) width
c
c Start a new path and loop through the points:
c
      write(lpsout,100)
      do 1 i=1,np
            ix = int(resc(xmin,xmax,pxmin,pxmax,x(i))/pscl)
            iy = int(resc(ymin,ymax,pymin,pymax,y(i))/pscl)
            if(i.eq.1) then
                  write(lpsout,101) ix,iy
            else
                  write(lpsout,102) ix,iy
            endif
 1    continue         
      if(lwidt.le.EPSLON) then
            write(lpsout,104) gray
      else
            write(lpsout,105) gray
      endif
 100  format('n')
 101  format(i4.4,1x,i4.4,1x,'m')
 102  format(i4.4,1x,i4.4,1x,'l')
 103  format(f6.3,' setlinewidth')
 104  format('c',/,f4.2,' setgray',/,'fill',/,'0.0 setgray')
 105  format('c gsave ',/,f4.2,' setgray',/,'fill',/,'grestore s',
     +         /,'0.00 setgray')
      return
      end



      subroutine psline(np,x,y,lwidt,idsh)
c-----------------------------------------------------------------------
c
c              Write Postscript line commands to a file
c              ****************************************
c
c
c CALLING ARGUMENTS:
c
c  np           the number of points in the x and y array to join
c  x()          array of x values in the range xmin to xmax
c  y()          array of y values in the range ymin to ymax
c  lwidt        the width of the line (1.0 = dark, 0.5 = light)
c  idsh         Dashing Index
c
c NOTES:
c
c  1. The pxmin,pxmax,.. variables are in the standard 1/72 inch 
c     resolution of the postscript page. If a different scale is 
c     going to be used in the printing set pscl to the scale.
c
c  2. If "idsh" is zero then no dashing is perfomed
c
c
c AUTHOR: C. Deutsch                                 DATE: November 1988
c-----------------------------------------------------------------------
      real      x(*),y(*),lwidt,lwold
      character dash(10)*24
c
c Common Block for Postscript Output Unit and Scaling:
c
      common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                xmax,ymin,ymax
      save   lwold
c
c Dash Patterns:
c
        data dash/'[40 20] 0 setdash       ',
     +            '[13 14 13 20] 0 setdash ',
     +            '[12 21 4 21] 0 setdash  ',
     +            '[10 10] 0 setdash       ',
     +            '[20 20] 0 setdash       ',
     +            '[30 30] 0 setdash       ',
     +            '[40 40] 0 setdash       ',
     +            '[50 50] 0 setdash       ',
     +            '[50 50] 0 setdash       ',
     +            '[50 50] 0 setdash       '/
c
c Change the line width if necessary:
c
      if(pscl.lt.0.01) pscl = 1.0
      if(idsh.gt.10)   idsh = 10 
      if(lwidt.ne.lwold) then
            width = lwidt/pscl
            write(lpsout,100) width
 100        format(f6.3,' setlinewidth')
            lwold = lwidt
      endif
c
c Start a new path and loop through the points:
c
      if(idsh.gt.0) write(lpsout,'(a24)') dash(idsh)
      write(lpsout,101)
 101  format('n')
      do 1 i=1,np
            ix = int(resc(xmin,xmax,pxmin,pxmax,x(i))/pscl)
            iy = int(resc(ymin,ymax,pymin,pymax,y(i))/pscl)
            if(i.eq.1) then
                  write(lpsout,102) ix,iy
 102              format(i4.4,1x,i4.4,' m')
            else
                  write(lpsout,103) ix,iy
 103              format(i4.4,1x,i4.4,' l')
            endif
 1    continue      
      write(lpsout,104)
 104  format('s')
      if(idsh.gt.0) write(lpsout,105)
 105  format('[] 0 setdash')
c
c Finished - Return to calling program:
c
      return
      end



        subroutine scal(xmin,xmax,ymin,ymax,xaxmin,xaxmax,yaxmin,
     +                  yaxmax)
c-----------------------------------------------------------------------
c
c Determines an appropriate labelling and tic mark interval for a set of
c axis and writes the results to either an interactive graphics device
c or a postscript file.
c
c INPUT VARIABLES:
c       xmin   - the minimum of the x axis (labeled on axis)
c       xmax   - the maximum of the x axis (labeled on axis)
c       ymin   - the minimum of the y axis (labeled on axis)
c       ymax   - the maximum of the y axis (labeled on axis)
c       xaxmin - the minimum of the x axis (on GKS window)
c       xaxmax - the maximum of the x axis (on GKS window)
c       yaxmin - the minimum of the y axis (on GKS window)
c       yaxmax - the maximum of the y axis (on GKS window)
c   
c
c EXTERNAL REFERENCES: The cdps library of subroutines
c
c AUTHOR: Clayton Deutsch                      Last Revision: April 1990
c-----------------------------------------------------------------------
      real       xloc(5),yloc(5)
      character  label*8,lfmt*8
c
c Common Block for Postscript Output Unit and Scaling:
c
      common /psdata/ lpsout,psscl,pxmin,pxmax,pymin,pymax,wxmin,
     +                  wxmax,wymin,wymax
c
c Check to make sure that the scale can be plotted:
c
      if((xmax-xmin).le.0.0001.or.(ymax-ymin).le.0.0001) return
      if(((xmin-.05*(xmax-xmin))*(xmin+.05*(xmax-xmin))).le.0.0)
     +  xmin = 0.0
      if(((ymin-.05*(ymax-ymin))*(ymin+.05*(ymax-ymin))).le.0.0)
     +  ymin = 0.0
c
c set up some of the parameters:
c
      tlng = 0.013 * ( (xaxmax-xaxmin) + (yaxmax-yaxmin) )
      tsht = 0.007 * ( (xaxmax-xaxmin) + (yaxmax-yaxmin) )
      psz  = 4.0 + 0.06  * (xaxmax-xaxmin)
      pl1  = 0.6 + 0.005 * (xaxmax-xaxmin)
      pl2  = 0.3 + 0.003 * (xaxmax-xaxmin)
c
c Write in the axis:
c
      xloc(1) = xaxmin
      yloc(1) = yaxmax
      xloc(2) = xaxmin
      yloc(2) = yaxmin
      xloc(3) = xaxmax
      yloc(3) = yaxmin
      num     = 3
      call psline(num,xloc,yloc,pl1,0)
c
c CONSTRUCT THE X AXIS:
c
c
c      Determine a labelling and tic mark increment:
c
            do 1 i=1,20
                  test = (xmax-xmin)/(10.0**(6-i))
                  if(test.gt.0.9) go to 2
 1            continue
 2            if(test.gt.3.0)                 zval = 1.0
            if(test.le.3.0.and.test.gt.2.0) zval = 0.5
            if(test.le.2.0.and.test.gt.1.2) zval = 0.4
            if(test.le.1.2)                 zval = 0.2
            nval = 5
            if(zval.eq.0.4.or.zval.eq.0.2)  nval = 4
            zval = zval * 10.0**(6-i)
            tval = zval / real(nval)
            if(i.ge.12) lfmt = '(f8.8)'
            if(i.eq.11) lfmt = '(f8.7)'
            if(i.eq.10) lfmt = '(f8.6)'
            if(i.eq.9)  lfmt = '(f8.5)'
            if(i.eq.8)  lfmt = '(f8.4)'
            if(i.eq.7)  lfmt = '(f8.3)'
            if(i.eq.6)  lfmt = '(f8.2)'
            if(i.eq.5)  lfmt = '(f8.1)'
            if(i.le.4)  lfmt = '(f8.0)'
c
c      Loop along the axis drawing the tic marks and labels:
c
            yloc(1) = yaxmin
            pos     = xmin
            num     = 2
            do 3 i=1,100
                  yloc(2) = yaxmin - tlng
                  xloc(1) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                  xloc(2) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                  call psline(num,xloc,yloc,pl2,0)
                  yloc(2) = yloc(1) - 2.5*tlng
                  write(label,lfmt) pos                        
                  call pstext(xloc(1),yloc(2),8,label,psz,1,0.0,1)
                  yloc(2) = yaxmin - tsht
                  do 4 j=1,nval-1
                       pos     = pos + tval
                       if(pos.gt.xmax) go to 5
                       xloc(1) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                       xloc(2) = resc(xmin,xmax,xaxmin,xaxmax,pos)
                       call psline(num,xloc,yloc,pl2,0)
 4                  continue
                  pos = pos + tval
                  if(pos.gt.xmax) go to 5
 3            continue
 5            continue
c
c CONSTRUCT THE Y AXIS:
c
c
c      Determine a labelling and tic mark increment:
c
            do 11 i=1,20
                  test = (ymax-ymin)/(10.0**(6-i))
                  if(test.gt.0.9) go to 12
 11            continue
 12             if(test.ge.3.0)                 zval = 1.0
            if(test.le.3.0.and.test.gt.2.0) zval = 0.5
            if(test.le.2.0.and.test.gt.1.2) zval = 0.4
            if(test.le.1.2)                 zval = 0.2
            nval = 5
            if(zval.eq.0.4.or.zval.eq.0.2)  nval = 4
            zval = zval * 10.0**(6-i)
            tval = zval / real(nval)
            if(i.ge.12) lfmt = '(f8.8)'
            if(i.eq.11) lfmt = '(f8.7)'
            if(i.eq.10) lfmt = '(f8.6)'
            if(i.eq.9)  lfmt = '(f8.5)'
            if(i.eq.8)  lfmt = '(f8.4)'
            if(i.eq.7)  lfmt = '(f8.3)'
            if(i.eq.6)  lfmt = '(f8.2)'
            if(i.eq.5)  lfmt = '(f8.1)'
            if(i.le.4)  lfmt = '(f8.0)'
c
c      Loop along the axis drawing the tic marks and labels:
c
            xloc(1) = xaxmin
            pos     = ymin
            num     = 2
            do 13 i=1,100
                  xloc(2) = xaxmin - tlng
                  yloc(1) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                  yloc(2) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                  call psline(num,xloc,yloc,pl2,0)
                  xloc(2) = xloc(2) - 0.2*tlng
                  write(label,lfmt) pos                        
                  call pstext(xloc(2),yloc(2),8,label,psz,1,0.0,2)
                  xloc(2) = xaxmin - tsht
                  do 14 j=1,nval-1
                       pos     = pos + tval
                       if(pos.gt.ymax) go to 15
                       yloc(1) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                       yloc(2) = resc(ymin,ymax,yaxmin,yaxmax,pos)
                       call psline(num,xloc,yloc,pl2,0)
 14                  continue
                  pos = pos + tval
                  if(pos.gt.ymax) go to 15
 13            continue
 15            continue
c
c Return to calling program:
c
      return
      end
        real function resc(xmin1,xmax1,xmin2,xmax2,x111)
        rsc=(xmax2-xmin2)/(xmax1-xmin1)
        resc=xmin2 + (x111 - xmin1 ) * rsc
        return
        end



      subroutine locate(xx,n,is,ie,x,j)
c-----------------------------------------------------------------------
c
c Given an array "xx" of length "n", and given a value "x", this routine
c returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
c must be monotonic, either increasing or decreasing.  j=0 or j=n is
c returned to indicate that x is out of range.
c
c Modified to set the start and end points by "is" and "ie" 
c
c Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
c-----------------------------------------------------------------------
      dimension xx(n)
c
c Initialize lower and upper methods:
c
      jl = is
      ju = ie
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
