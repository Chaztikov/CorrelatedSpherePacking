      subroutine sgsim
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
c This subroutine generates 3-D realizations of a Gaussian process with
c a given autocovariance model, and conditional to input Gaussian data.
c The conditional simulation is achieved by sequential simulation of all
c the nodes visited by a random path.
c
c
c
c PROGRAM NOTES:
c
c  1. The three dimensional anisotropy parameters, i.e., of the search
c     ellipse and variogram ranges are described in section 2.3 of the
c     manual.   The variogram parameters are described in the same place
c
c  2. The original data and previously simulated grid nodes can be
c     searched separately.  There can be a different maximum number of
c     each and a minimum number of original data can be specified
c     to restrict simulation beyond the limits of the data.  The
c     closeness of previously simulated grid nodes is measured according
c     to the variogram structural distance.
c
c
c
c INPUT VARIABLES:
c
c   nd               Number of data (no missing values)
c   x,y,z(nd)        coordinates of the data
c   vr(nd)           gaussian data (normal scores)
c
c   nx,ny,nz         Number of blocks in X,Y, and Z
c   xmn,ymn,zmn      Coordinate at the center of the first Block
c   xsiz,ysiz,zsiz   spacing of the grid nodes (block size)
c
c   nsim             number of simulations
c   seed             starting random number seed
c   ktype            =1, ordinary kriging; =0, simple kriging
c   sim(nx,ny,nz)    one simulation
c   idbg             integer debugging level (0=none,2=normal,4=serious)
c   ldbg             unit number for the debugging output
c   lout             unit number for the output
c   lint             interim storage of the simulated normal scores
c
c   radius           Maximum search radius
c   sang1            Azimuth angle of the principal search direction
c   sang2            Dip angle of the principal search direction
c   sang3            Third rotation angle of the search ellipse
c   sanis1           Anisotropy for the dip angle
c   sanis2           Anisotropy for the plunge angle
c   ndmin            Minimum number of data required before sim
c   ndmax            Maximum number of samples for simulation
c   noct             Maximum number per octant if an octant search is
c                      desired (if <= 0, then no octant search)
c
c   nodmax           Maximum number of previously simulated grid nodes
c                      to consider in the simulation.  The structural
c                      variogram distance is used to identify close ones
c
c   c0               Nugget constant (isotropic).
c   cc(nst)          Multiplicative factor of each nested structure.
c   aa(nst)          Parameter "a" of each nested structure.
c   it(nst)          Type of nested structures (1=sph,2=exp,3=gau,4=pow)
c   ang1(nst)        Azimuth angle for the principal direction
c   ang2(nst)        Dip angle for the principal direction
c   ang3(nst)        Third rotation angle to rotate the two minor
c                      directions around the principal direction
c   anis1(nst)       Anisotropy (radius in minor direction at 90
c                      degrees from "ang1" divided by the principal
c                      radius in direction "ang1")
c   anis2(nst)       Anisotropy (radius in minor direction at 90 degrees
c                      vertical from "ang1" divided by the principal
c                      radius in direction "ang1")
c
c
c OUTPUT VARIABLES:  Simulated Values are written to "lint"
c
c
c
c EXTERNAL REFERENCES:
c
c   super            Sets up the super block search of original data
c   search           Search for nearby data values
c   ctable           Builds a covariance table and "spiral" search
c   srchnd           Search for nearby simulated grid nodes
c   sqdist           computes anisotropic squared distance
c   sortem           sorts multiple arrays in ascending order (separate)
c   cova3            Calculates the covariance given a variogram model
c   krige            Sets up and solves either the SK or OK system
c   ksol             Linear system solver using Gaussian elimination
c
c
c
c Concepts taken from F. Alabert and E. Isaaks
c
c Original:  C.V. Deutsch                                      Aug. 1990
c-----------------------------------------------------------------------
      include  'sgsim.inc'
      real      randnu(1)
      real*8    p
c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search.
c
      do 50 is=1,nst
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
 50   continue
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
c
c Set up the super block search, the covariance table, and the spiral
c search for previously simulated grid nodes:
c
      if(sstrat.ne.1) then
            call super
            call picksup
      else
            ndmax = 0
      endif
      call ctable
      call rand(seed,1,randnu)
      seed = 0
c
c The random path generation procedure of Srivastava and Gomez has been
c adopted in this subroutine.  A linear congruential generator of the
c form: r(i) = 5*r(i-1)+1 mod(2**n) has a cycle length of 2**n.  By
c choosing the smallest power of 2 that is still larger than the total
c number of points to be simulated, the method ensures that all indices
c will be generated once and only once:
c
      nxy  = nx * ny
      nxyz = nx * ny * nz
      modlus = 1
 1    modlus = modlus * 2
      if(modlus.le.nxyz) go to 1
      write(ldbg,*)
      write(ldbg,*) 'The random number generator to generate the path:'
      write(ldbg,*) '        r(i) = 5*r(i-1)+1 mod ',modlus
c
c Open a temporary file for the untransformed realization:
c
      lint = 7
      open(lint,status='SCRATCH')
c
c MAIN LOOP OVER ALL THE SIMULAUTIONS:
c
      do 2 isim=1,nsim
c
c Initialize the simulation:
c
            do 3 iz=1,nz
            do 3 iy=1,ny
            do 3 ix=1,nx
                  sim(ix,iy,iz) = UNEST
 3          continue
            write(*,101) isim
 101        format(/' Working on simulation number: ',i4)
c
c Make sure that data which fall on grid nodes are assigned to
c the nodes before proceeding (eliminates singular matrices).  Also,
c assign them anyway if a two part search strategy is chosen:
c
            TINY = 0.0001
            do 30 id=1,nd
                  atnode(id) = .false.
                  ix = min0(max0(int((x(id)-xmn)/xsiz+1.5),0),nx)
                  iy = min0(max0(int((y(id)-ymn)/ysiz+1.5),0),ny)
                  iz = min0(max0(int((z(id)-zmn)/zsiz+1.5),0),nz)
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
                  test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
                  if(sstrat.eq.1.or.test.le.TINY) then
c
c This data should be assigned to the nearest node location (unless
c there is a closer data):
c
                        if(sim(ix,iy,iz).ge.0.0) then
                              id2 = int(sim(ix,iy,iz)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2)) 
     +                                               + abs(zz-z(id2))
                              if(test.le.test2) sim(ix,iy,iz) = real(id)
                              write(ldbg,102) id,id2
                        else
                              sim(ix,iy,iz) = real(id)
                        end if
                        atnode(id) = .true.
                  end if
 30         continue
 102        format(' WARNING data values ',2i5,' are both assigned to ',
     +           /,'         the same node - taking the closest')
c
c Now, enter data values into the simulated grid:
c
            do 31 iz=1,nz
            do 31 iy=1,ny
            do 31 ix=1,nx
                  id = int(sim(ix,iy,iz)+0.5)
                  if(id.gt.0) sim(ix,iy,iz) = vr(id)
 31         continue
c
c Locate a random starting point:
c
            call rand(seed,1,randnu)
            index = 1 + int(randnu(1)*(nx*ny*nz-1))
            irepo = max0(1,min0((nxyz/10),10000))
c
c MAIN LOOP OVER ALL THE NODES:
c
            do 4 in=1,nxyz
                  if((int(in/irepo)*irepo).eq.in) write(*,103) in
 103              format('   currently on node ',i9)
c
c Figure out the location of this point and make sure it has
c not been assigned a value already:
c
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  if(sim(ix,iy,iz).ne.UNEST) go to 5
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
c
c Now, we'll simulate the point ix,iy,iz.  First, get the close data
c and make sure that there are enough to actually simulate a value,
c we'll only keep the closest "ndmax" data, and look for previously
c simulated grid nodes:
c
                  if(sstrat.ne.1) then
                        call search(xx,yy,zz)
                        if(nclose.lt.ndmin) go to 5
                        if(nclose.gt.ndmax) nclose = ndmax
                  endif
                  call srchnd(ix,iy,iz)
c
c Calculate the conditional mean and standard deviation.  This will be
c done with kriging if there are data, otherwise, the global mean and
c standard deviation will be used:
c
                  if((nclose+ncnode).lt.1) then
                        cmean  = 0.0
                        cstdev = 1.0
                  else
                        call krige(ix,iy,iz,xx,yy,zz,cmean,cstdev)
                  endif
c
c Draw a random number and assign a value to this node:
c
                  call rand(seed,1,randnu)
                  p = dble(randnu(1))
                  call gauinv(p,xp,ierr)
                  sim(ix,iy,iz) = xp * cstdev + cmean
c
c Quick check for far out results:
c
                  if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or.
     +               abs(sim(ix,iy,iz)).gt.6.0) then
                  write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(ix,iy,iz)
  104             format('WARNING: grid node location: ',3i5,/,
     +                   '         conditional mean:   ',f12.5,/,
     +                   '         conditional stdev:  ',f12.5,/,
     +                   '         simulated value:    ',f12.5)
                  endif
c
c Get a new node to simulate:
c
 5                index = mod(5*index+1,modlus)
                  if(index.gt.nxyz.or.index.lt.1) go to 5
c
c END MAIN LOOP OVER NODES:
c
 4          continue
c
c Write this simulation to the output file:
c
            do 10 iz=1,nz
            do 10 iy=1,ny
            do 10 ix=1,nx
                  write(lint,'(f8.4)') sim(ix,iy,iz)
 10         continue
c
c END MAIN LOOP OVER SIMULATIONS:
c
 2    continue
c
c Return to the main program:
c
      return
      end



      subroutine super
c-----------------------------------------------------------------------
c
c           Establish Super Block Search Limits and Sort Data
c           *************************************************
c
c This subroutine sets up a 3-D "super block" model and orders the data
c by super block number.  The limits of the super block is set to the
c minimum and maximum limits of the grid; data outside are assigned to
c the nearest edge block.
c
c The idea is to establish a 3-D block network that contains all the
c relevant data.  The data are then sorted by their index location in
c the search network, i.e., the index location is given after knowing
c the block index in each coordinate direction (ix,iy,iz):
c          ii = (iz-1)*nsbx*nsby + (iy-1)*nsbx + ix
c An array, the same size as the number of super blocks, is constructed
c that contains the cumulative number of data in the model.  With this
c array it is easy to quickly check what data are located near any given
c location.
c
c
c
c INPUT VARIABLES:
c
c   nx,xmn,xsiz      Definition of the X grid being considered
c   ny,ymn,ysiz      Definition of the Y grid being considered
c   nz,zmn,zsiz      Definition of the Z grid being considered
c   nsbx             Number of "super blocks" in X direction
c   nsby             Number of "super blocks" in Y direction
c   nsbz             Number of "super blocks" in Z direction
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   vr(nd)           Variable at each location. WARNING: if cokriging or
c                      cosimulation is being considered then all arrays
c                      must be sorted!
c   tmp(nd)          Temporary storage to keep track of the super block
c                      index associated to each data (uses the same
c                      storage already allocated for the simulation)
c
c
c
c OUTPUT VARIABLES:
c
c   nisb()           An array that tells how many data are in each super
c                      super block.
c   xsbsiz           X size of each super block
c   ysbsiz           Y size of each super block
c   zsbsiz           Z size of each super block
c
c
c
c EXTERNAL REFERENCES:
c
c   sortem           Sorting routine to sort the data
c   sqdist           Computes anisotropic squared distance
c
c
c
c Author: C. Deutsch                                     Date: July 1990
c-----------------------------------------------------------------------
      include    'sgsim.inc'
      dimension   tmp(MAXDAT)
c
c Establish the number and size of the super blocks:
c
      nsbx = min0(nx,MAXSBX)
      nsby = min0(ny,MAXSBY)
      nsbz = min0(nz,MAXSBZ)
      xsbsiz = real(nx)*xsiz/real(nsbx)
      ysbsiz = real(ny)*ysiz/real(nsby)
      zsbsiz = real(nz)*zsiz/real(nsbz)
      if(idbg.ge.3) then
            write(ldbg,*) 'SUPER BLOCK SETUP: X ',nsbx,xsbsiz
            write(ldbg,*) '                   Y ',nsby,ysbsiz
            write(ldbg,*) '                   Z ',nsbz,zsbsiz
      endif
c
c Initialize the extra super block array to zeros:
c
      do 1 i=1,(nsbx*nsby*nsbz)
 1    nisb(i) = 0
c
c Loop over all the data assigning the data to a super block and
c accumulating how many data are in each super block:
c
      do 2 i=1,nd
            ix = int((x(i)-xmn)/xsbsiz) + 1
            iy = int((y(i)-ymn)/ysbsiz) + 1
            iz = int((z(i)-zmn)/zsbsiz) + 1
c
c Make sure that each coordinate index is acceptable.  If a point is
c outside thenn it will be assigned to an edge block:
c
            ix = min0((max0(ix,1)),nsbx)
            iy = min0((max0(iy,1)),nsby)
            iz = min0((max0(iz,1)),nsbz)
c
c Assign the super block index:
c
            ii = (iz-1)*nsbx*nsby + (iy-1)*nsbx + ix
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
 2    continue
c
c Sort the holes by ascending super block number:
c
      call sortem(1,nd,tmp,4,x,y,z,vr,f,g,h)
c
c Set up array nisb with the starting address of the block data:
c
      do 3 i=1,(nsbx*nsby*nsbz-1)
 3    nisb(i+1) = nisb(i) + nisb(i+1)
c
c Set up is complete - Return to calling program:
c
      return
      end



      subroutine picksup
c-----------------------------------------------------------------------
c
c             Establish Which Super Blocks to Search
c             **************************************
c
c This subroutine establishes which super blocks must be searched given
c that a point being estimated/simulated falls within a super block
c centered at 0,0,0:
c
c
c Author: C. Deutsch                                    Date: April 1992
c-----------------------------------------------------------------------
      include 'sgsim.inc'
c
c Set up the limits of the possible super blocks:
c
      zero    = 0.0
      nsbtosr = 0
c
c MAIN Loop over all possible super blocks:
c
      do 1 i=-(nsbx-1),(nsbx-1)
      do 1 j=-(nsby-1),(nsby-1)
      do 1 k=-(nsbz-1),(nsbz-1)
c
c Find the closest distance between the corners of the super blocks:
c
            shortest = 1.0e21
            do 2 i1=0,1
            do 2 j1=0,1
            do 2 k1=0,1
                  do 3 i2=i,i+1
                  do 3 j2=j,j+1
                  do 3 k2=k,k+1
                        xdis = real(i2-i1)*xsbsiz
                        ydis = real(j2-j1)*ysbsiz
                        zdis = real(k2-k1)*zsbsiz
                        hsqd = sqdist(zero,zero,zero,xdis,ydis,zdis,
     +                                isrot,MAXROT,rotmat)
                        if(hsqd.lt.shortest) shortest = hsqd
 3                continue
 2          continue
c
c Keep this super block if it is close enoutgh:
c
            if(shortest.le.1.2*radsqd) then
                  ic = nsbx + i
                  jc = nsby + j
                  kc = nsbz + k
                  lc = (kc-1)*MXSXY + (jc-1)*MXSX + ic
                  nsbtosr = nsbtosr + 1
                  sbtosr(nsbtosr) = lc
            endif
 1    continue
c
c Set up is complete - Return to calling program:
c
      return
      end



      subroutine search(xloc,yloc,zloc)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radius           search radius
c   sang1,2,3        Define orientation of search ellipse
c   sanis1,2         Define anisotropy of search ellipse
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   vr(nd)           Variable at each location. WARNING: if cokriging or
c                      cosimulation is being considered then all arrays
c                      must be sorted!
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   xmn,ymn,zmn      Origin of the super block grid system
c   nsbx,nsby,nsbz   Number of "super blocks" in X,Y, and Z directions
c   xsbsiz           X size of each super block
c   ysbsiz           Y size of each super block
c   zsbsiz           Z size of each super block
c   nisb()           The array specifying the data in each super block
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c AUTHOR: C. Deutsch                                     DATE: July 1990
c-----------------------------------------------------------------------
      include  'sgsim.inc'
      dimension tmp(MAXDAT)
c
c Determine Super Block Location of Point being Estimated:
c
      ix = int((xloc-xmn)/xsbsiz) + 1
      iy = int((yloc-ymn)/ysbsiz) + 1
      iz = int((zloc-zmn)/zsbsiz) + 1
      ix = min0((max0(ix,1)),nsbx)
      iy = min0((max0(iy,1)),nsby)
      iz = min0((max0(iz,1)),nsbz)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isb=1,nsbtosr
c
c Is this super block within the grid system:
c
            lc = sbtosr(isb)
            kc = int(lc/MXSXY) + 1
            jc = int((lc-(kc-1)*MXSXY)/MXSX) + 1
            ic = lc - (kc-1)*MXSXY - (jc-1)*MXSX
            ixsb = ix + (ic-nsbx)
            iysb = iy + (jc-nsby)
            izsb = iz + (kc-nsbz)
            if(ixsb.le.0.or.ixsb.gt.nsbx.or.
     +         iysb.le.0.or.iysb.gt.nsby.or.
     +         izsb.le.0.or.izsb.gt.nsbz) go to 1
            ii = (izsb-1)*nsbx*nsby + (iysb-1)*nsbx + ixsb
c
c Figure out how many samples in this super block:
c
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check for a data coinciding with a node location:
c
                  if(atnode(i)) go to 2
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),isrot,
     +                          MAXROT,rotmat)
                  if(hsqd.gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = hsqd
 2          continue
 1    continue
c
c Sort the nearby samples by distance to point being estimated:
c
      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do 3 i=1,8
 3    inoct(i) = 0
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do 4 j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
 4    continue
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do 8 i=1,8
 8    if(inoct(i).gt.0) infoct = infoct + 1
      return
      end



      subroutine ctable
c-----------------------------------------------------------------------
c
c               Establish the Covariance Look up Table
c               **************************************
c
c The idea is to establish a 3-D network that contains the covariance
c value for a range of grid node offsets that should be at as large
c as twice the search radius in each direction.  The reason it has to
c be twice as large as the search radius is because we want to use it
c to compute the data covariance matrix as well as the data-block
c covariance matrix.
c
c Secondly, we want to establish a search for nearby nodes that
c in order of closeness as defined by the variogram.
c
c
c
c INPUT VARIABLES:
c
c   xsiz,ysiz,zsiz  Definition of the grid being considered
c   MAXCTX,Y,Z      Number of blocks in covariance table
c
c   covariance table parameters
c
c
c
c OUTPUT VARIABLES:  covtab()         Covariance table
c
c EXTERNAL REFERENCES:
c
c   sqdist          Computes 3-D anisotropic squared distance
c   sortem          Sorts multiple arrays in ascending order
c   cova3           Computes the covariance according to a 3-D model
c
c
c
c Author: C. Deutsch                                     Date: July 1990
c-----------------------------------------------------------------------
      parameter(TINY=1.0e-10)
      include  'sgsim.inc'
      logical   first
      dimension tmp1(MAXXYZ),tmp2(MAXXYZ)
c
c Size of the look-up table:
c
      nctx = (MAXCTX-1)/2
      ncty = (MAXCTY-1)/2
      nctz = min0(((MAXCTZ-1)/2),(nz-1))
c
c Debugging output:
c
      write(ldbg,*)
      write(ldbg,*) 'Covariance Look up table and search for previously'
      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
      write(ldbg,*) 'coordinate direction for covariance look up is:'
      write(ldbg,*) '          X direction: ',nctx*xsiz
      write(ldbg,*) '          Y direction: ',ncty*ysiz
      write(ldbg,*) '          Z direction: ',nctz*zsiz
      write(ldbg,*) 'Node Values are not searched beyond this distance!'
      write(ldbg,*)
c
c NOTE: If dynamically allocating memory, and if there is no shortage
c       it would a good idea to go at least as far as the radius and
c       twice that far if you wanted to be sure that all covariances
c       in the left hand covariance matrix are within the table look-up.
c
c Initialize the covariance subroutine and cbb at the same time:
c
      first  = .true.
      cbb    = cova3(0.0,0.0,0.0,0.0,0.0,0.0,first)
      first  = .false.
c
c Now, set up the table and keep track of the node offsets that are
c within the search radius:
c
      nlooku = 0
      do 1 i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do 1 j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do 1 k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            covtab(ic,jc,kc) = cova3(0.0,0.0,0.0,xx,yy,zz,first)
            hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,
     +                          MAXROT,rotmat)
            if(hsqd.le.radsqd) then
                  nlooku         = nlooku + 1
c
c We subtract the covariance from a large value so that the ascending
c sort subroutine will accomplish the sort we want.  Furthermore, a
c fraction of the distance is also taken off so that we search by
c anisotropic distance once we are beyond the range:
c
                  tmp1(nlooku) = - (covtab(ic,jc,kc)-TINY*hsqd)
                  tmp2(nlooku) = real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
            endif
 1    continue
c
c Finished setting up the look-up table, now order the nodes such
c that the closest ones, according to variogram distance, are searched
c first. Note: the "loc" array is used because I didn't want to make
c special allowance for 2 byte integers in the sorting subroutine:
c
      call sortem(1,nlooku,tmp1,1,tmp2,c,d,e,f,g,h)
      do 2 il=1,nlooku
            loc = int(tmp2(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = int(iz)
            iynode(il) = int(iy)
            ixnode(il) = int(ix)
 2    continue
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
c
c Debugging output if requested:
c
      if(idbg.lt.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.14) return
      do 3 i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
 3    continue
 100  format('Point ',i3,' at ',3f12.4)
c
c All finished:
c
      return
      end



      subroutine srchnd(ix,iy,iz)
c-----------------------------------------------------------------------
c
c               Search for nearby Simulated Grid nodes
c               **************************************
c
c The idea is to spiral away from the node being simulated and note all
c the nearby nodes that have been simulated.
c
c
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   sim(nx,ny,nz)   the simulation so far
c   nodmax          the maximum number of nodes that we want
c   nlooku          the number of nodes in the look up table
c   i[x,y,z]node    the relative indices of those nodes.
c   [x,y,z]mn       the origin of the global grid netwrok
c   [x,y,z]siz      the spacing of the grid nodes.
c
c
c
c OUTPUT VARIABLES:
c
c   ncnode          the number of close nodes
c   icnode()        the number in the look up table
c   cnode[x,y,z]()  the location of the nodes
c   cnodev()        the values at the nodes
c
c
c
c Author: C. Deutsch                                     Date: July 1990
c-----------------------------------------------------------------------
      include  'sgsim.inc'
      integer   ninoct(8)
c
c Consider all the nearby nodes until enough have been found:
c
      ncnode = 0
      if(noct.gt.0) then
            do 1 i=1,8
                  ninoct(i) = 0
 1          continue
      end if
      do 2 il=2,nlooku
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            if(sim(i,j,k).ne.UNEST) then
c
c Check the of data already taken from this octant:
c
                  if(noct.gt.0) then                
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                              iq = 4
                              if(idx.le.0 .and. idy.gt.0) iq = 1
                              if(idx.gt.0 .and. idy.ge.0) iq = 2
                              if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                              iq = 8
                              if(idx.le.0 .and. idy.gt.0) iq = 5
                              if(idx.gt.0 .and. idy.ge.0) iq = 6
                              if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(i,j,k)
                  if(ncnode.eq.nodmax) return
            endif
 2    continue
c
c Return to calling program:
c
      return
      end



      subroutine krige(ix,iy,iz,xx,yy,zz,cmean,cstdev)
c-----------------------------------------------------------------------
c
c            Builds and Solves the SK or OK Kriging System
c            *********************************************
c
c INPUT VARIABLES:
c
c   ix,iy,iz        index of the point currently being simulated
c   xx,yy,zz        location of the point currently being simulated
c
c
c
c OUTPUT VARIABLES:
c
c   cmean           kriged estimate
c   cstdev          kriged standard deviation
c
c
c
c EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
c
c
c
c ORIGINAL: C.V. Deutsch                               DATE: August 1990
c-----------------------------------------------------------------------
      include 'sgsim.inc'
      logical first
c
c Size of the kriging system:
c
      first = .false.
      na    = nclose + ncnode
      neq   = na + ktype
c
c Set up kriging matrices:
c
      in=0
      do 1 j=1,na
c
c Sort out the actual location of point "j"
c
            if(j.le.nclose) then
                  index  = int(close(j))
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  vra(j) = vr(index)
            else
c
c It is a previously simulated node (keep index for table look-up):
c
                  index  = j-nclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
                  vra(j) = cnodev(index)
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
            endif
            do 2 i=1,j
c
c Sort out the actual location of point "i"
c
                  if(i.le.nclose) then
                        index  = int(close(i))
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
c
c It is a previously simulated node (keep index for table look-up):
c
                        index  = i-nclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
c
c Now, get the covariance value:
c
                  in = in + 1
c
c Decide whether or not to use the covariance look-up table:
c
                  if(j.le.nclose.or.i.le.nclose) then
                        cov   = cova3(x1,y1,z1,x2,y2,z2,first)
                        a(in) = dble(cov)
                  else
c
c Try to use the covariance look-up (if the distance is in range):
c
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +                     jj.lt.1.or.jj.gt.MAXCTY.or.
     +                     kk.lt.1.or.kk.gt.MAXCTZ) then
                              cov = cova3(x1,y1,z1,x2,y2,z2,first)
                        else
                              cov = covtab(ii,jj,kk)
                        endif
                        a(in) = dble(cov)
                  endif
 2          continue
c
c Get the RHS value (possibly with covariance look-up table):
c
            if(j.le.nclose) then
                  cov  = cova3(xx,yy,zz,x1,y1,z1,first)
                  r(j) = dble(cov)
            else
c
c Try to use the covariance look-up (if the distance is in range):
c
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +               jj.lt.1.or.jj.gt.MAXCTY.or.
     +               kk.lt.1.or.kk.gt.MAXCTZ) then
                        cov = cova3(xx,yy,zz,x1,y1,z1,first)
                  else
                        cov = covtab(ii,jj,kk)
                  endif
                  r(j) = dble(cov)
            endif
            rr(j) = r(j)
 1    continue
c
c Addition of OK constraint:
c
      if(ktype.eq.1) then
            do 3 i=1,na
                  in    = in + 1
                  a(in) = 1.0
 3          continue
            in      = in + 1
            a(in)   = 0.0
            r(neq)  = 1.0
            rr(neq) = 1.0
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
      if(idbg.ge.3) then
            write(ldbg,100) ix,iy,iz
            is = 1
            do 4 i=1,neq
                  ie = is + i - 1
                  write(ldbg,101) i,r(i),(a(j),j=is,ie)
                  is = is + i
 4          continue
 100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
c
c Solve the Kriging System:
c
      if(neq.eq.1.and.ktype.eq.0) then
            s(1)  = r(1) / a(1)
            ising = 0
      else
            call ksol(1,neq,1,a,r,s,ising)
      endif
c
c Write a warning if the matrix is singular:
c
      if(ising.ne.0) then
            if(idbg.ge.1) then
                  write(ldbg,*) 'WARNING SGSIM: singular matrix'
                  write(ldbg,*) '               for node',ix,iy,iz
            endif
            cmean  = 0.0
            cstdev = 1.0
            return
      endif
c
c Write out the kriging Matrix if Seriously Debugging:
c
      if(idbg.ge.3) then
            do 40 i=1,na
                  write(ldbg,140) i,s(i)
 40         continue
 140        format(' Kriging weight for data: ',i4,' = ',f8.4)
      endif
c
c Compute the estimate, the kriging variance, and return:
c
      cmean  = 0.0
      cstdev = cbb
      do 6 i=1,na
            cmean  = cmean  + real(s(i))*vra(i)
            cstdev = cstdev - real(s(i))*rr(i)
 6    continue
      if(ktype.eq.1) cstdev = cstdev - real(s(na+1))
      if(cstdev.lt.0.0) then
            write(ldbg,*) 'NEGATIVE VARIANCE: ',cstdev
            cstdev = 0.0
      endif
      cstdev = sqrt(cstdev)
      return
      end



      real function cova3(x1,y1,z1,x2,y2,z2,first)
c-----------------------------------------------------------------------
c
c              Covariance Between Two Points (3-D Version)
c              *******************************************
c
c This function returns the covariance associated with a variogram model
c that is specified by a nugget effect and possibly four different
c nested varigoram structures.  The anisotropy definition can be
c different for each of the nested structures (spherical, exponential,
c gaussian, or power).
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   nst              Number of nested structures (max. 4).
c   c0               Nugget constant (isotropic).
c   cmax             Maximum variogram value needed for kriging when
c                      using power model.  A unique value of cmax is
c                      used for all nested structures which use the
c                      power model.  therefore, cmax should be chosen
c                      large enough to account for the largest single
c                      structure which uses the power model.
c   cc(nst)          Multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(nst)          Parameter "a" of each nested structure.
c   it(nst)          Type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c   ang1             Azimuth angle for the principal direction of
c                    continuity (measured clockwise in degrees from Y)
c   ang2             Dip angle for the principal direction of continuity
c                    (measured in negative degrees down from horizontal)
c   ang3             Third rotation angle to rotate the two minor
c                    directions around the principal direction defined
c                    by ang1 and ang2.  A positive angle acts clockwise
c                    while looking in the principal direction.
c   anis1            Anisotropy (radius in minor direction at 90
c                    degrees from "ang1" divided by the principal radius
c                    in direction "ang1")
c   anis2            Anisotropy (radius in minor direction at 90 degrees
c                    vertical from "ang1" divided by the principal
c                    radius in direction "ang1")
c   first            A logical variable which is set to true if the
c                      direction specifications have changed - causes
c                      the rotation matrices to be recomputed.
c
c
c
c OUTPUT VARIABLES: returns "cova3" the covariance obtained from the
c                   variogram model.
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,DTOR=PI/180.0,PMX=9999.)
      include  'sgsim.inc'
      logical   first
c
c The first time around calculate the maximum covariance value:
c
      if(first) then
            cmax = c0
            do 1 is=1,nst
                  if(it(is).eq.5) then
                        cmax = cmax + 1.0
                  else if(it(is).eq.4) then
                        cmax = cmax + PMX
                  else
                        cmax = cmax + cc(is)
                  endif
 1          continue
      endif
c
c Check for very small distance:
c
      is   = 1
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,is,MAXROT,rotmat)
      if(hsqd.lt.EPSLON) then
            cova3 = cmax
            return
      endif
c
c Non-zero distance, loop over all the structures:
c
      cova3 = 0.0
      do 2 is=1,nst
c
c Compute the appropriate structural distance:
c
            if(is.ne.1) hsqd=sqdist(x1,y1,z1,x2,y2,z2,is,MAXROT,rotmat)
            h = sqrt(hsqd)
            if(it(is).eq.1) then
c
c Spherical model:
c
                  hr = h/aa(is)
                  if(hr.ge.1.0) go to 2
                  cova3 = cova3 + cc(is)*(1.-hr*(1.5-.5*hr*hr))
            else if(it(is).eq.2) then
c
c Exponential model:
c
                  cova3 = cova3 +cc(is)*exp(-h/aa(is))
            else if(it(is).eq.3) then
c
c Gaussian model:
c
                  hh=-(h*h)/(aa(is)*aa(is))
                  cova3 = cova3 +cc(is)*exp(hh)
            else if(it(is).eq.4) then
c
c Power model:
c
                  cov1  = cmax - cc(is)*(h**aa(is))
                  cova3 = cova3 + cov1
            else if(it(is).eq.5) then
c
c 2-D bombing model:
c
                  hr = h/aa(is)
                  if(hr.ge.1.0) go to 2
                  si = cc(is)*(1.0-cc(is))
                  circ = 2/(PI*aa(is)*aa(is)) * (h*sqrt(aa(is)*aa(is) -
     +                   h*h) + aa(is)*aa(is)*asin(h/aa(is)))
                  cova3 = ( si - cc(is)*(1.0-cc(is)**circ) ) / si
            endif
 2    continue
      return
      end



      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              The matrix indicator to initialize
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      parameter(DEG2RAD=3.14159265/180.0,EPSLON=1.0e-20)
      real      rotmat(MAXROT,3,3)
c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
c
c Get the required sines and cosines:
c
      sina  = sin(alpha)
      sinb  = sin(beta)
      sint  = sin(theta)
      cosa  = cos(alpha)
      cosb  = cos(beta)
      cost  = cos(theta)
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / amax1(anis1,EPSLON)
      afac2 = 1.0 / amax1(anis2,EPSLON)
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end



      real function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The matrix indicator to initialize
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      real rotmat(MAXROT,3,3)
c
c Compute component distance vectors and the squared distance:
c
      dx = x1 - x2
      dy = y1 - y2
      dz = z1 - z2
      sqdist = 0.0
      do 1 i=1,3
            cont   = rotmat(ind,i,1) * dx
     +                 + rotmat(ind,i,2) * dy
     +                 + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
 1      continue
      return
      end
 
 
 
      subroutine ksol(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for OKB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for okb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do 70 k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do 60 iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do 50 i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do 30 j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
 30                     continue
                        do 40 llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
 40                     continue
 50               continue
 60         continue
 70   continue
c
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
      do 140 iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do 100 llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
 100        continue
            i=neq
            kk=ij
            do 130 ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do 120 llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do 110 j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
 110                    continue
                        s(ll1)=ap*piv
 120              continue
 130        continue
 140  continue
c
c Finished solving back, return:
c
      return
      end



      subroutine rand(seed,n,vector)
c---------------------------------------------------------------------
c
c This random number generator generates random numbers in ]0,1[
c Note that if the seed value is zero on the first call, a default
c value of 1369 will be used  in a linear congruential generator to
c generate 55 odd integers for the array 'itab()'. These values are
c preserved by a common statement, so that they may be used in sub-
c sequent calls by setting the seed to zero.If the value of 'seed'
c is greater than zero in a call to the subroutine, then the array
c 'itab' will be initialized and a new seed value will be returned
c by the subroutine. Best results are obtained by making the initial
c call with a seed of your choice and then setting the seed to '0'
c for all subsequent calls.
c
c---------------------------------------------------------------------
      dimension vector(*)
      common /unusual/itab(55),n1,n2,nseed
      integer rn1,seed
c
c Test to see if 55 odd integers must be generated.
c
      if((seed.gt.0).or.(nseed.lt.1)) then
            nseed = seed
            if(seed.le.0) nseed  = 7931
            do 10 i=1,55
                  rn1=mod(nseed*9069,32768)
                  if(mod(rn1,2).eq.0) rn1 = rn1-1
                  itab(i) = rn1
                  nseed = rn1
 10         continue
            n1 = 0
            n2 = 24
      endif
c
c generate "n" random components for the vector "VECTOR"
c
      do 30 i=1,n
            itab(55-n1) = mod(itab(55-n2)*itab(55-n1),32768)
            vector(i)   = abs(float(itab(55-n1))/float(32768))
            n1 = mod(n1+1,55)
            n2 = mod(n2+1,55)
 30   continue
      if(seed.gt.0) seed=nseed
      return
      end



      subroutine nscore(nd,vr,tmin,tmax,iwt,wt,tmp,vrg,ierror)
c-----------------------------------------------------------------------
c
c              Transform Univariate Data to Normal Scores
c              ******************************************
c
c This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
c by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
c as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
c so that the data can be returned in the same order (i.e., just in
c case there are associated arrays like the coordinate location).
c
c
c
c INPUT VARIABLES:
c
c   nd               Number of data (no missing values)
c   vr(nd)           Data values to be transformed
c   tmin,tmax        data trimming limits
c   iwt              =0, equal weighted; =1, then apply weight
c   wt(nd)           Weight for each data (don't have to sum to 1.0)
c   tmp(nd)          Temporary storage space for sorting
c
c
c
c OUTPUT VARIABLES:
c
c   vrg(nd)          Normal Scores
c   ierror           Error Flag: 0 - Error Free
c                                1 - No data to transform
c                                2 - Problem with weights (not used)
c                                3 - Unacceptable cdf value
c                                4 - Despiking is necessary
c
c
c
c EXTERNAL REFERENCES:
c
c   gauinv           Calculates the inverse of a Gaussian cdf
c   sortem           sorts a number of arrays according to a key array
c
c
c
c Original:  R.M. Srivastava                                   Aug. 1987
c Revisions: C.V. Deutsch                                      July 1990
c-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)
      dimension vr(*),wt(*),vrg(*),tmp(*)
      real*8    pd
c
c Compute the factor so that the weights sum to 1.0:
c
      ierror = 1
      wtfac  = 0.0
      do 1 i=1,nd
            if(vr(i).le.tmin.or.vr(i).gt.tmax) go to 1
            if(iwt.ge.1) then
                  wtfac = wtfac + wt(i)
            else
                  wtfac = wtfac + 1
            endif
 1    continue
c
c The sum of the weights must be greater than epsilon:
c
      if(wtfac.lt.EPSLON) return
      ierror = 0
      wtfac  = 1.0 / wtfac
c
c Prepare tmp array (so that we can get back where we started):
c
      do 3 i=1,nd
 3    tmp(i) =i
c
c Sort values in ascending order:
c
      ib = 1
      ip = 2
      call sortem(ib,nd,vr,ip,wt,tmp,d,e,f,g,h)
c
c Calculate cumulative probabilities and corresponding gaussian
c quantiles:
c
      oldcp = 0.0
      cp    = 0.0
      do 4 i=1,nd
            if(vr(i).gt.tmin.and.vr(i).le.tmax) then
                  w  = wtfac
                  if(iwt.ge.1) w = w*wt(i)
                  cp = cp + w
                  p  =(cp + oldcp)/2.0
                  if(p.lt.0.0.or.p.gt.1.0) then
                        ierror = 3
                        return
                  endif
                  pd = dble(p)
                  call gauinv(pd,vrg(i),ierr)
                  oldcp = cp
                  if(i.gt.1.and.vr(i).eq.vr(i-1)) ierror=4
            endif
 4    continue
c
c Get the arrays back in original order:
c
      ip = 3
      call sortem(ib,nd,tmp,ip,wt,vr,vrg,e,f,g,h)
c
c Now, return to calling routine:
c
      return
      end



      subroutine gauinv(p,xp,ierr)
c-----------------------------------------------------------------------
c
c Computes the inverse of the standard normal cumulative distribution
c function with a numerical approximation from : Statistical Computing,
c by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
c
c
c
c INPUT/OUTPUT:
c
c   p    = double precision cumulative probability value: dble(psingle)
c   xp   = G^-1 (p) in single precision
c   ierr = 1 - then error situation (p out of range), 0 - OK
c
c
c-----------------------------------------------------------------------
      real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p
      save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
c
c Coefficients of approximation:
c
      data lim/1.0e-20/
      data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/,
     +     p3/-0.0204231210245/,p4/-0.0000453642210148/
      data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/,
     +     q3/0.103537752850/,q4/0.0038560700634/
c
c Initialize:
c
      ierr = 1
      xp   = 0.0
      pp   = dble(p)
      if(p.gt.0.5) pp = 1 - pp
      if(dble(p).lt.lim) return
      ierr = 0      
      if(p.eq.0.5) return
c
c Approximate the function:
c
      y  = dsqrt(dlog(1.0/(pp*pp)))
      xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) /
     +               ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
      if(real(p).eq.real(pp)) xp = -xp
c
c Return with G^-1(p):
c
      return
      end



      subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c
c This is a subroutine for sorting a real array in ascending order. This
c routine is a translation into Fortran of algorithm 271, quickersort,
c by R.S. Scowen in collected algorithms of the acm. consult this for
c further references.
c
c The method is used is that of continually splitting the array into
c parts such that all elements of one part are less than all elements
c of the other, with a third part in the middle consisting of one
c element.  An element with value t is chosen arbitrarily (here we
c choose the middle element). i and j give the lower and upper limits
c of the segment being split.  After the split a value q will have
c been found such that a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.
c The program then performs operations on the two segments (i,q-1) and
c (q+1,j) as follows.  The smaller segment is split and the position
c of the larger segment is stored in the lt and ut arrays.  If the
c segment to be split contains two or fewer elements, it is sorted and
c another segment is obtained from the lt and ut arrays.  When no more
c segments remain, the array is completely sorted.
c
c
c INPUT PARAMETERS:
c
c   ib,ie        start and end index of the array to be sorteda
c   a            array, a portion of which has to be sorted.
c   iperm        0 no other array is permuted.
c                1 array b is permuted according to array a
c                2 arrays b,c are permuted.
c                3 arrays b,c,d are permuted.
c                4 arrays b,c,d,e are permuted.
c                5 arrays b,c,d,e,f are permuted.
c                6 arrays b,c,d,e,f,g are permuted.
c                7 arrays b,c,d,e,f,g,h are permuted.
c               >7 no other array is permuted.
c
c   b,c,d,e,f,g,h  arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c-----------------------------------------------------------------------
      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
      integer   lt(17),ut(17),i,j,k,m,p,q
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
      if (a(k).le.ta) go to 20
c
c Such an element has now been found now search for a q such that
c a(q)<t starting at the end of the segment.
c
 30      continue
      if (a(q).lt.ta) go to 40
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
      if (a(i).le.a(j)) go to 100
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
