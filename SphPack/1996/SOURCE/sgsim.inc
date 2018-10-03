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
c The following Parameters control static dimensioning for SGSIM:
c
c   MAXX      maximum nodes in X
c   MAXY      maximum nodes in Y
c   MAXZ      maximum nodes in Z
c   MAXCTX    maximum X points in covariance table (ODD number)
c   MAXCTY    maximum Y points in covariance table (ODD number)
c   MAXCTZ    maximum Z points in covariance table (ODD number)
c   MAXSBX    maximum super block nodes in X direction
c   MAXSBY    maximum super block nodes in Y direction
c   MAXSBZ    maximum super block nodes in Z direction
c   MAXDAT    maximum number of data
c   MAXNOD    maximum number of previously simulated nodes to use
c   MAXSAM    maximum number of data for one simulation
c   MAXNST    maximum number of nested structures for variogram
c
c-----------------------------------------------------------------------
c
c User Adjustable Parameters:
c
      parameter(MAXX   =    80, MAXY   =  80, MAXZ   = 80,
     +          MAXCTX =    61, MAXCTY =  61, MAXCTZ = 61,
     +          MAXSBX =     1, MAXSBY =   1, MAXSBZ =  1,
     +          MAXDAT =   200, MAXNOD =  24, MAXSAM = 24,
     +          MAXNST =     4)
c
c Fixed Parameters:
c
      parameter(MAXKR1=MAXNOD+MAXSAM+1,MAXKR2=MAXKR1*MAXKR1,
     +          MAXROT=MAXNST+1,MAXCXY=MAXCTX*MAXCTY,
     +          MAXXYZ=MAXCTX*MAXCTY*MAXCTZ,MAXSB=MAXSBX*MAXSBY*MAXSBZ,
     +          MXSXY=4*MAXSBX*MAXSBY,MXSX=2*MAXSBX,UNEST=-99.0,
     +          EPSLON=1.0e-20,VERSION=1.201)
c
c Variable Declaration:
c
      real      ltpar,x(MAXDAT),y(MAXDAT),z(MAXDAT),vr(MAXDAT),
     +          wt(MAXDAT),close(MAXDAT),sim(MAXX,MAXY,MAXZ),
     +          cc(MAXNST),aa(MAXNST),ang1(MAXNST),ang2(MAXNST),
     +          ang3(MAXNST),anis1(MAXNST),anis2(MAXNST),
     +          rotmat(MAXROT,3,3),covtab(MAXCTX,MAXCTY,MAXCTZ),
     +          cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     +          cnodev(MAXNOD),vra(MAXKR1)
      real*8    r(MAXKR1),rr(MAXKR1),s(MAXKR1),a(MAXKR2)
      integer*2 ixnode(MAXXYZ),iynode(MAXXYZ),iznode(MAXXYZ)
      integer   seed,ltail,utail,sstrat,nisb(MAXSB),sbtosr(8*MAXSB),
     +          inoct(8),it(MAXNST),icnode(MAXNOD)
      logical   atnode(MAXDAT)
c
c Common Blocks:
c
      common /grid3d/ nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn

      common /data3d/ nd,x,y,z,vr,wt,ltail,ltpar,utail,utpar,zmin,zmax

      common /datalg/ atnode

      common /simula/ nsim,seed,ktype,sim,idbg,ldbg,lint,lin,lout

      common /superb/ nsbx,nsby,nsbz,xsbsiz,ysbsiz,zsbsiz,nisb,sbtosr,
     +                nsbtosr

      common /srch3d/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,
     +                noct,icros,nclose,close,ndmin,ndmax,sstrat,
     +                infoct,inoct

      common /cova3d/ nst,c0,it,cmax,cc,aa,ang1,ang2,ang3,anis1,anis2,
     +                isrot,rotmat

      common /clooku/ nctx,ncty,nctz,covtab,nlooku,ncnode,icnode,
     +                cnodex,cnodey,cnodez,cnodev,nodmax,ixnode,iynode,
     +                iznode

      common /krigev/ r,rr,s,a,vra,cbb
