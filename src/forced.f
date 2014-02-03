      subroutine forced(beta,npts,dom,ntar,zz,VLP)
c  ****** DESCRIPTION:
c  Do volume integral to solve forced problem  
c

c  ****** INPUT:
c  beta     = parameter in the argument of K_0
c  ff       = forcing function
c  zz       = coordinates of target locations as complex *16
c  Npt      = number of target locations
c  dom      = solid walls as complex *16
c  un       = solid walls as complex *16
c  Nbd      = number of boundary points

c  ****** OUTPUT:
c  VLP      = volume potential at the target points


c      implicit none

      integer maxpts,maxbd

      integer npts,ntar
      real *8 beta
      complex *16 zz(ntar)
      complex *16 dom(npts)

      real *8 VLP(npts)


      integer *4 i
      integer *4 stat,ierror
      integer *4 maxboxes,maxlevel,maxwrk
      parameter (maxboxes=2000)
      parameter (maxlevel=10)
      parameter (maxiwrk=9*maxboxes+1)
      parameter (maxrwrk=100000)
      integer *4 nboxes,nlev
      integer *4 levelbox(maxboxes)
      integer *4 iparentbox(maxboxes)
      integer *4 ichildbox(4,maxboxes)
      integer *4 icolbox(maxboxes)
      integer *4 irowbox(maxboxes)
      integer *4 iboxlev(maxboxes)
      integer *4 icolleagbox(9,maxboxes)
      integer *4 nblevel(0:maxlevel)
      integer *4 istartlev(0:maxlevel)

c      real *8 xf(16,maxboxes),yf(16,maxboxes)
      real *8 fGrid(16,maxboxes),pot(16,maxboxes)
      integer *4 iwork(maxiwrk)
      real *8 work(maxrwrk)




      call formTree(levelbox,icolbox,irowbox,nboxes,nlev,
     1      iparentbox,ichildbox,nblevel,iboxlev,istartlev,
     2      maxboxes,maxlevel)
c     Form tree that is subdivided if error between a 16x16
c     grid and 64x64 grid differ by a tolerance

      call balanceTree(levelbox,iparentbox,ichildbox,
     1    icolbox,irowbox,icolleagbox,nboxes,nlev,
     2    nblevel,iboxlev,istartlev,maxboxes)
c     Balance tree 2:1 so that neighbouring boxes differ by
c     no more than one leve

      call writeTree2Matlab(nboxes,ichildbox,
     1      irowbox,icolbox,levelbox,
     2      '../output/buildtreeForced.m')
c     write the tree structure to a matlab file

      call formRHS(fGrid,levelbox,nboxes,icolbox,irowbox,
     1      ichildbox,nlev,nblevel,iboxlev,istartlev) 
c     Get x and y coordinates corresponding to the balanced tree

      call startVolumeFMM(beta,work,maxrwrk,iwork,maxiwrk,
     1    nlev,levelbox,iparentbox,ichildbox,
     2    icolbox,irowbox,nboxes,nblevel,iboxlev,
     3    istartlev,fGrid,pot)




      do i = 1,ntar
        VLP(i) = 1.d0
      enddo


      end 

c***********************************************************************
      subroutine startVolumeFMM(beta,work,lenw,iwork,ilen,
     1    nlev,levelbox,iparentbox,ichildbox,
     2    icolbox,irowbox,nboxes,nblevel,iboxlev,
     3    istartlev,fGrid,pot)
      implicit none

      real *8 beta
      integer *4 lenw,ilen
      integer *4 nboxes,nlev
      real *8 work(lenw)
      integer *4 iwork(ilen)
      integer *4 levelbox(nboxes),iparentbox(nboxes)
      integer *4 ichildbox(4,nboxes)
      integer *4 icolbox(nboxes),irowbox(nboxes)
      integer *4 nblevel(0:nlev),iboxlev(nboxes)
      integer *4 istartlev(0:nlev)
      real *8 fGrid(16,nboxes),pot(16,nboxes)

      integer *4 nnodesmax
      integer *4 nterms,iprec
      integer *4 ninbox,ndeg
      integer *4 lcolleag,mcolleag,itot
      integer *4 lmpole,llocxp,lexp,lcoeff
      integer *4 lnodes,lscale,lcomp,lchoose,lzs
      integer *4 iscale,ibetascal
      integer *4 mpole,locexp
      integer *4 iexpn,iexps,iexpe,iexpw
      integer *4 mcomp,mcomp2,mchoose
      integer *4 mzs,mxnodes,mcoeff

      ndeg = 3
      nterms = 32
      iprec = 12
      nnodesmax = 50

      lcolleag = 9*nboxes
      mcolleag = 1
      itot = mcolleag + lcolleag

      if (itot .ge. ilen) then
        write(6,*) 'THE INTEGER WORKSPACE REQUIRED', itot
        write(6,*) 'THE TOTAL WORKSPACD NEEDED EXCEEDS THE'
        write(6,*) 'AMOUNT ALOTTED.  PLEASE READJUST MAXIWORK'
        stop
      endif


      lmpole = 2*(nterms+1)*nboxes
      llocxp = 2*(nterms+1)*nboxes
      lexp = 2*nnodesmax*nboxes
      lcoeff = 16*nboxes
      lnodes = nnodesmax
      lscale = nlev + 1
      lcomp = nnodesmax*(2*nterms+1)
      lchoose = 2*(nterms+1) * (2*nterms+1)
      lzs = 2*7*7*nnodesmax


      iscale = 1
      ibetascal = iscale + lscale
      mpole = ibetascal + lscale
      locexp = mpole + lmpole
      iexpn = locexp + llocxp
      iexps = iexpn + lexp
      iexpe = iexps + lexp
      iexpw = iexpe + lexp
      mcomp = iexpw + lcomp
      mcomp2 = mcomp + lcomp
      mchoose = mcomp2 + lcomp
      mzs = mchoose + lchoose
      mxnodes = mzs + lzs
      mcoeff = mxnodes + lnodes
      itot = mcoeff + lcoeff

      if (itot .ge. lenw) then
        write(6,*) 'THE REAL WORKSPACE REQUIRED', itot
        write(6,*) 'THE TOTAL WORKSPACD NEEDED EXCEEDS THE'
        write(6,*) 'AMOUNT ALOTTED.  PLEASE READJUST MAXRWORK'
        stop
      endif

      call evalVolumeFMM(beta,work(ibetascal),work(iscale),
     1    iprec,nlev,ndeg,work(mxnodes),work(mzs),work(mcomp),
     2    work(mcomp2),nterms,pot,
     3    work(mpole),work(locexp),work(iexpn),work(iexps),
     4    work(iexpe),work(iexpw),work(mchoose),work(mcoeff),
     5    levelbox,iparentbox,ichildbox,
     6    icolbox,irowbox,iwork(mcolleag),nboxes,
     7    nblevel,iboxlev,istartlev,fGrid)



      return
      end




c***********************************************************************
      subroutine evalVolumeFMM(beta,betascal,scale,
     1    iprec,nlev,ndeg,xnodes,zs,comp,
     2    comp2,nterms,pot,
     3    mpole,locexp,expn,exps,
     4    expe,expw,c,coeffs,
     5    levelbox,iparentbox,ichildbox,
     6    icolbox,irowbox,icolleagbox,nboxes,
     7    nblevel,iboxlev,istartlev,fGrid)
c***********************************************************************
      implicit none

      integer *4 nnodesmax
      parameter (nnodesmax = 50)

      integer iprec,nlev,ndeg,nterms,nboxes
      real *8 beta
      integer *4 icolbox(1),irowbox(1)
      integer *4 nblevel(0:1),iboxlev(1)
      integer *4 istartlev(0:1)
      integer *4 levelbox(1),iparentbox(1)
      integer *4 ichildbox(4,1),icolleagbox(9,1)
      real *8 betascal(0:nlev),scale(0:nlev)
      real *8 pot(16,1),fGrid(16,1)
      real *8 xnodes(1)
      real *8 coeffs(0:ndeg,0:ndeg,1)
      real *8 comp(nnodesmax,-nterms:nterms)
      real *8 comp2(nnodesmax,-nterms:nterms)
      complex *16 c(0:nterms,-nterms:nterms)
      complex *16 mpole(0:nterms,1),locexp(0:nterms,1)
      complex *16 expn(nnodesmax,1),exps(nnodesmax,1)
      complex *16 expe(nnodesmax,1),expw(nnodesmax,1)
      complex *16 zs(-3:3,-3:3,nnodesmax)


c      complex *16 wint(nterms/4+1,10)
      integer *4 i,j,info
      real *8 A(10,16)
      real *8 xp(16),yp(16),thetaxy(16),rouxy(16)


      betascal(0) = beta
      scale(0) = 1.d0
      do i=1,nlev
        betascal(i) = betascal(i-1)/2.d0
        scale(i) = scale(i)/2.d0
      enddo

      do i = 1,nboxes
        do j = 1,16
          pot(i,j) = 0.d0
        enddo
      enddo

      call interpMatrix(A)
c     Get interpolation matrix that takes 16 function values
c     to 10 polynomial coefficients

      call getCoeffs(coeffs,fGrid,A,nlev,
     1    ichildbox,nblevel,iboxlev,istartlev)
c     Get the polynomial coefficients for each of the
c     childless boxes


      call mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev,istartlev)
c     make colleagues of the final 2:1 balanced tree


      do i=1,4
        xp((i-1)*4+1) = -3.75e-1
        xp((i-1)*4+2) = -1.25e-1
        xp((i-1)*4+3) = 1.25e-1
        xp((i-1)*4+4) = 3.75e-1
        yp(i) = -3.75e-1
        yp(4+i) = -1.25e-1
        yp(8+i) = 1.25e-1
        yp(12+i) = 3.75e-1
      enddo

      call coorTransform(1,0.d0,0.d0,16,xp,yp,rouxy,thetaxy,info)



c     start of upward pass

      do i = 1,nboxes
        do j = 0,nterms
          mpole(j,i) = 0.d0
          locexp(j,i) = 0.d0
        enddo
        do j = 1,nnodesmax
          expe(j,i) = 0.d0
          expw(j,i) = 0.d0
          expn(j,i) = 0.d0
          exps(j,i) = 0.d0
        enddo
      enddo



      return
      end

c***********************************************************************
      subroutine coorTransform(iflag,x0,y0,N,x,y,rou,theta,info)
      implicit none

      real *8 x0,y0
      integer iflag,N
      real *8 x(N),y(N),rou(N),theta(N)
      integer *4 info

      integer *4 i
      real *8 pi
      real *8 xx,yy

      pi = 4.d0*datan(1.d0)
      if (iflag .eq. 1) then
        do i=1,N
          xx = x(i) - x0
          yy = y(i) - y0
          rou(i) = dsqrt(xx*xx + yy*yy)
          if (xx .ge. 0.d0 .and. yy .ge. 0.d0) then
            if (yy .eq. 0.d0) then
              theta(i) = 0.d0
            else
              theta(i) = datan(yy/xx)
            endif
          elseif (xx .ge. 0.d0 .and. yy .le. 0.d0) then
            if (yy .eq. 0.d0) then
              theta(i) = 0.d0
            else
              theta(i) = 2.0d0*pi - datan(-yy/xx)
            endif
          elseif (xx .le. 0.d0 .and. yy .ge. 0.d0) then
            if (yy .eq. 0.d0) then
              theta(i) = pi 
            else
              theta(i) = pi - datan(-yy/xx)
            endif
          else
            if (yy .eq. 0.0d0) then
              theta(i) = pi
            else
              theta(i) = pi + datan(yy/xx)
            endif
          endif
        enddo
      elseif (iflag .eq. 2) then
        do i = 1,N
          x(i) = x0 + rou(i)*dcos(theta(i))
          y(i) = y0 + rou(i)*dsin(theta(i))
        enddo
      endif


      return
      end







c***********************************************************************
      subroutine getCoeffs(coeffs,fGrid,A,nlev,
     1    ichildbox,nblevel,iboxlev,istartlev)
      implicit none

      integer *4 nlev
      integer *4 ichildbox(4,1)
      integer *4 nblevel(0:1),iboxlev(1),istartlev(0:1)
      real *8 coeffs(0:3,0:3,1),fGrid(16,1),A(10,16)


      integer *4 i,j,ibox
      real *8 xlength

      xlength = 1.d0

      do j = 0,nlev
        do i = istartlev(j),istartlev(j) + nblevel(j) - 1
          ibox = iboxlev(i)

          if (ichildbox(1,ibox) .lt. 0) then
            call interpolate(A,fGrid(1,ibox),
     1              coeffs(0,0,ibox),xlength)
          endif
        enddo
        xlength = xlength/2.d0
      enddo



      return
      end


c***********************************************************************
      subroutine formTree(levelbox,icolbox,irowbox,nboxes,nlev,
     1           iparentbox,ichildbox,nblevel,iboxlev,istartlev,
     3           maxboxes,maxlevel)
c***********************************************************************
c   description:

c   input:
c     maxboxes : the maximum number of boxes allowed.
c     maxlevel : the deepest level allowed.
c     eps : the error bound that determines when refinement take place
c     h : the real function that is the right hand side of the yukawa
c         equation
c
c   output:
c     istartlev : the pointer to where each level begins in the
c                 iboxlev array
c     levelbox : an array determining the level of each box. (maxboxes)
c     nboxes : the total number of boxes
c     nlev : the finest level
c     icolbox : the column of each box. (maxboxes)
c     irowbox : the row of each box. (size : maxboxes)
c     iparentbox : the parent of each box. (size : maxboxes)
c     ichildbox : the four children of each box. (4, maxboxes)
c     nblevel : the total number of boxes per level. (0:maxlevel)
c     iboxlev : the array in which the boxes are arranged. (maxboxes).
c
c   called from : forced
c
c   subroutine called :
c     intermatrix4(a) : get the precomputed matrix which maps the
c          function value to coefficients.
c     mkgrid() : find the two grids for the box. 16x16 and 64x64.
c     geterror() : the error of the 16x16 and 64x64. 64x64 are considered
c          exact.
c     subdivide1() : when the error is larger than threshhold, subdivide.
c
c***************************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4 maxboxes
      integer *4 levelbox(maxboxes)
      integer *4 nlev, nboxes,  maxlevel
      integer *4 icolbox(maxboxes),irowbox(maxboxes)
      integer *4 iparentbox(maxboxes),ichildbox(4,maxboxes)
      integer *4 nblevel(0:maxlevel),iboxlev(maxboxes)
      integer *4 istartlev(0:maxlevel)
c
c-----local variables
c
      integer *4 i, ibox
      integer *4 j, istart, iend
      integer *4 levflag, nbound
      real *8 eps
      real *8 error, hh
      real *8 xf(16), yf(16)
      real *8 xf2(64), yf2(64)
      real *8 a(10,16)
c
      
      eps = 1.d-4
      call interpMatrix(a)
c
      do i = 0, maxlevel
        nblevel(i) = 0
        istartlev(i) = 0
      end do
c
      do i = 1, maxboxes
        iboxlev(i) = 0
      end do
c
c-----first set the big parent box to the appropriate settings:
c     (initially, there is just one box and it is the big parent box
c      at level 0)
c
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
c
      nboxes = 1
      nlev = 0
c
c-----we also need to initialize the adaptive 'ladder'
c     structures to the correct initial values:
c
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1
c
      do i = 0, maxlevel - 1
        levflag = 0
        istart = istartlev(i)
        iend = istart + nblevel(i) - 1
        nbound = nboxes
        do j = istart, iend
          ibox = iboxlev(j)
c
c---------now let's find the two grids for this box:
c         xf and yf are the 16 point grids that the
c         function is actually defined on.  xf2 and yf2
c         represent the 64 point grid that would be this
c         boxes children
c
          call mkgrid(xf,yf,xf2,yf2,icolbox(ibox),
     1            irowbox(ibox),levelbox(ibox))

c
c---------now let's compute the error in the polynomial
c         approximation for a given box:
c         (the error from the approximate polynomial from
c         xf, yf is compared to the exact functin values on
c         xf2, yf2.  the l2 error between the two values
c         is set to error)
c
          call getError(xf,yf,xf2,yf2,levelbox(ibox),a,error)
c
c---------now let's test to see if the error is above or below
c         the required threshold.  note that we must scale the
c         error to account for the difference in levels:
c         now let's readjust epsilon to account for the
c         fact that we are on different levels:
c

          hh = dble(4**levelbox(ibox))
c
          if (error .ge. eps*hh) then
c
c-----------call subdivide1
c
            if (levflag .eq. 0)then
              nlev = nlev + 1
              levflag = 1
            endif
            call subdivide1(ibox,iparentbox,
     1           ichildbox,nboxes,irowbox,icolbox,levelbox,nlev,
     2           istartlev, nblevel, iboxlev,
     3           maxboxes,maxlevel)

          endif
        end do
c
      end do
c
      return
      end




c***************************************************************************
      subroutine balanceTree(levelbox,iparentbox,ichildbox,
     1    icolbox,irowbox,icolleagbox,nboxes,nlev,
     2    nblevel,iboxlev,istartlev,maxboxes)
c***************************************************************************
c   description:
c
c   input:
c
c     levelbox : an array determining the level of each box
c     iparentbox : the parent of each box
c     ichildbox : the four children of each box
c     icolbox : the column of each box
c     irowbox : the row of each box
c     icolleagbox : the colleagues of a given box
c     nboxes : the total number of boxes
c     nlev : the finest level
c     nblevel : the total number of boxes per level
c     iboxlev : the array in which the boxes are arranged
c     istartlev : the pointer to where each level begins in the
c               iboxlev array
c
c   output:
c     icolbox, irowbox, icolleagbox, nboxes, and all other
c     arrays describing the tree may be change on output
c
c
c   called from : main.
c
c   subroutine called : sortBoxes().
c                       mkcolls().
c                       subdivide().
c*************************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4 maxboxes
      integer *4 levelbox(maxboxes), icolleagbox(9,maxboxes)
      integer *4 iparentbox(maxboxes), ichildbox(4,maxboxes)
      integer *4 icolbox(maxboxes), irowbox(maxboxes)
      integer *4 nboxes, nlev
      integer *4 nblevel(0:nlev),iboxlev(maxboxes),istartlev(0:nlev)
      integer *4 itemparray(maxboxes)
c
c-----local variables
c
      integer *4 iflag(maxboxes)
      integer *4 ichild(4),icoll(9), ibox
      integer *4 i, ipar, itest, j, nb
      integer *4 itemp, ntemp, jcntr, icntr
      integer *4 start, istop

c
c-----first let's call a subroutine that will
c     generate all of the colleagues for each
c     box.  the colleagues are generated in the
c     correct order so there is no need to 'shuffle'
c     them later on.
c
      call mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev,istartlev)
c
c-----let's initialize all of the flags to zero.
c
      do i = 1, nboxes
        iflag(i) = 0
      end do

c
c-----find all of the boxes that need to be
c     flagged.  a flagged box will be denoted by
c     setting iflag(box) = 1.
c     this refers to any box that is directly touching
c     a box that is more than one level smaller than
c     it.  it is found by performing an upward pass
c     and lokking a box's parents parents and seeing
c     if they are childless and contact the given box.
c     note that we only need to get up to level two, as
c     we will not find a violation at a coarser level
c     than that.
c
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          ipar  = iparentbox(ibox)
          itest = iparentbox(ipar)
c
          icoll(1) = icolleagbox(1,itest)
          icoll(2) = icolleagbox(2,itest)
          icoll(3) = icolleagbox(3,itest)
          icoll(4) = icolleagbox(4,itest)
          icoll(5) = icolleagbox(5,itest)
          icoll(6) = icolleagbox(6,itest)
          icoll(7) = icolleagbox(7,itest)
          icoll(8) = icolleagbox(8,itest)
          icoll(9) = icolleagbox(9,itest)
c
          ichild(1) = ichildbox(1,itest)
          ichild(2) = ichildbox(2,itest)
          ichild(3) = ichildbox(3,itest)
          ichild(4) = ichildbox(4,itest)
c
          do nb = 1, 9
            itemp = icoll(nb)
            if(ichildbox(1,itemp) .lt. 0)then
c
c-------------the neighboring box is not divided
c             we could have problems.
c
              if (nb .eq. 1)then
                if(ipar .eq. ichild(4))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 2)then
                if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 3)then
                if(ipar .eq. ichild(3))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 4)then
                if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 6)then
                if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 7)then
                if(ipar .eq. ichild(1))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 8)then
                if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                  iflag(itemp) = 1
                end if
              elseif (nb .eq. 9)then
                if(ipar .eq. ichild(2))then
                  iflag(itemp) = 1
                end if
              endif
            endif
          end do
        end do
      end do

c
c-----find all of the boxes that need to be
c     given a flag+.  a flag+ box will be denoted by
c     setting iflag(box) = 2.
c     this refers to any box that is not already flagged
c     and is bigger than and is contacting a flagged box
c     or another box that has already been given a flag+.
c     it is found by performing an upward pass
c     and looking at a flagged box's parents colleagues
c     and a flag+ box's parents colleagues and seeing if
c     they are childless and present the case where a
c     bigger box is contacting a flagged or a flag+ box.
c
      do i = nlev, 2, -1
        do j = istartlev(i), istartlev(i) + nblevel(i) - 1
          ibox = iboxlev(j)
          if(iflag(ibox) .eq. 1 .or. iflag(ibox) .eq. 2)then
c
            ipar  = iparentbox(ibox)
c
            icoll(1) = icolleagbox(1,ipar)
            icoll(2) = icolleagbox(2,ipar)
            icoll(3) = icolleagbox(3,ipar)
            icoll(4) = icolleagbox(4,ipar)
            icoll(5) = icolleagbox(5,ipar)
            icoll(6) = icolleagbox(6,ipar)
            icoll(7) = icolleagbox(7,ipar)
            icoll(8) = icolleagbox(8,ipar)
            icoll(9) = icolleagbox(9,ipar)
c
            ichild(1) = ichildbox(1,ipar)
            ichild(2) = ichildbox(2,ipar)
            ichild(3) = ichildbox(3,ipar)
            ichild(4) = ichildbox(4,ipar)
c
            do nb = 1, 9
              itemp = icoll(nb)
c
c-------------let's check using the same criteria as above, but noting that
c             a flag will take precedence over a flag+.
c
              if(ichildbox(1,itemp) .lt. 0
     1            .and. iflag(itemp) .ne. 1)then
c
c---------------the neighboring box is not divided
c               we could have problems.
c
                if (nb .eq. 1)then
                  if(ibox .eq. ichild(4))then
                      iflag(itemp) = 2
                  end if
                elseif (nb .eq. 2)then
                  if(ibox .eq. ichild(3) .or. ibox .eq. ichild(4))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 3)then
                  if(ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 4)then
                  if(ibox .eq. ichild(4) .or. ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 6)then
                  if(ibox .eq. ichild(2) .or. ibox .eq. ichild(3))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 7)then
                  if(ibox .eq. ichild(1))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 8)then
                  if(ibox .eq. ichild(1) .or. ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                  end if
                elseif (nb .eq. 9)then
                  if(ibox .eq. ichild(2))then
                    iflag(itemp) = 2
                  end if
                endif
              endif
            end do
          endif
        end do
      end do
c
c-----now let's divide the boxes that need to be immediately
c     divided up.  all of the flagged and flag+ boxes need to
c     be divided one time.  the distinction lies in the fact
c     that the children of a flag+ box will never need to be
c     divided but the children of a flagged box may need to
c     be divided further.
c     below, all flagged and flag+ boxes are divided once.  the
c     children of a flag+ box are left unflagged while those of
c     the flagged boxes are given a flag++ (denoted by setting
c     iflag(box) = 3) which will be needed in the downward pass.
c

      ntemp = nboxes
      do i = 1, ntemp
c
c-------divide flagged boxes:
c
        if (iflag(i) .eq. 1)then
c
          if(ichildbox(1,i) .lt. 0)then
            call subdivide(i,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,
     3         itemparray)
          endif
c
c---------give flag++ to children of flagged boxes.
c
          itemp = ichildbox(1,i)
          iflag(itemp) = 3
c
          itemp = ichildbox(2,i)
          iflag(itemp) = 3
c
          itemp = ichildbox(3,i)
          iflag(itemp) = 3
c
          itemp = ichildbox(4,i)
          iflag(itemp) = 3
c
c-------divide flag+ boxes.
c
        elseif (iflag(i) .eq. 2)then
c
          if(ichildbox(1,i) .lt. 0)then
c            call subdivide(i,iparentbox,ichildbox,icolleagbox,
c     1         nboxes,irowbox,icolbox,levelbox,nlev,
c     2         istartlev, nblevel, iboxlev,
c     3         itemparray)
          endif
c
        endif
      end do
c
c-----now we need to do a downward pass.
c     we will concern oursleves only with the children of
c     flagged boxes and their children.  at each level,
c     for each flag++ box, test colleagues children and see
c     if they have children that are contacting you.  if so,
c     divide and flag++ all children that are created.
c
      do i = 0, nlev
        ntemp = nboxes
        start = istartlev(i)
        istop  = istartlev(i) + nblevel(i) - 1
        do 500 j = start, istop
          ibox = iboxlev(j)
c
c---------only be concerned with boxes on this level and
c         boxes that are given a flag++:
c
          if(iflag(ibox) .ne. 3)goto 500

          icoll(1) = icolleagbox(1,ibox)
          icoll(2) = icolleagbox(2,ibox)
          icoll(3) = icolleagbox(3,ibox)
          icoll(4) = icolleagbox(4,ibox)
          icoll(5) = icolleagbox(5,ibox)
          icoll(6) = icolleagbox(6,ibox)
          icoll(7) = icolleagbox(7,ibox)
          icoll(8) = icolleagbox(8,ibox)
          icoll(9) = icolleagbox(9,ibox)
c
c---------scan colleagues.
c
          do 400 jcntr = 1, 9
            if(icoll(jcntr) .lt. 0)goto 400
            if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 400
c
            ichild(1) = ichildbox(1,icoll(jcntr))
            ichild(2) = ichildbox(2,icoll(jcntr))
            ichild(3) = ichildbox(3,icoll(jcntr))
            ichild(4) = ichildbox(4,icoll(jcntr))
c
c-----------scan colleague's children.
c
            do 300 icntr = 1, 4
              if (ichildbox(1,ichild(icntr)) .lt. 0)goto 300
c
              if(jcntr .eq. 1 .and. icntr .eq. 2)then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 2 .and.
     1          (icntr .eq. 1 .or. icntr .eq. 2))then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 3 .and. icntr .eq. 1)then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev, 
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 4 .and.
     1          (icntr .eq. 2 .or. icntr .eq. 3))then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 6 .and.
     1          (icntr .eq. 1 .or. icntr .eq. 4))then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 7 .and. icntr .eq. 3)then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 8 .and.
     1          (icntr .eq. 3 .or. icntr .eq. 4))then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
c
              elseif(jcntr .eq. 9 .and. icntr .eq. 4)then
c
c---------------call subdivide
c
                if(ichildbox(1,ibox) .lt. 0)then
                  call subdivide(ibox,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
                endif
c
c---------------flag++ all children created
c
                itemp = ichildbox(1,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(2,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(3,ibox)
                iflag(itemp) = 3
                itemp = ichildbox(4,ibox)
                iflag(itemp) = 3
              endif
300         continue
400       continue
500     continue
      end do

      return
      end



c***************************************************************************
      subroutine mkgrid(xf,yf,xf2,yf2,icol,irow,level)
c***************************************************************************
c     this subroutine is only called within the formTree subroutine.
c     this subroutine generates two grids: (xf,yf) and (xf2,yf2).
c     they represent the 16 point grid that the polynomial
c     coefficients are obtained from and the 64  point grid (the children's
c     grids) that is used to test the accuracy of the polynomial
c     approximation.  (both grids are face centered.)
c
c   input:
c     icol : the column of the given box
c     irow : the row of the given box
c     level : the level of the given box
c
c   output:
c     xf : the x values of 16 cell centered points in the box
c     yf : the y values of 16 cell centered points in the box
c     xf2 : the x values of 64 cell centered points in the box
c     yf2 : the y values of 64 cell centered points in the box
c
c   called from : formTree().
c
c   subroutine called : none.
c
c***************************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4  level
      integer *4  icol, irow
      real *8  xf(16), yf(16)
      real *8  xf2(64), yf2(64)
c
c-----local variables
c
      integer *4  ilev, j, l
      integer *4  itemp1, itemp2
      real *8  xincrem, xstart
      real *8  xincrem2, xstart2
      real *8  x, xshift
      real *8  y, yshift
c
      ilev = level
      itemp1 = 4 * 2**ilev
      itemp2 = 2 * itemp1
c
      xstart   =  1.0d0 / dble(itemp2)   - .50d0
      xstart2  =  1.0d0 / dble(2*itemp2) - .50d0
      xincrem  =  1.0d0 / dble(itemp1)
      xincrem2 =  1.0d0 / dble(2*itemp1)
c
      xshift  =  dble(icol - 1) / dble(2**ilev)
      yshift  =  dble(irow - 1) / dble(2**ilev)
c
      x = xstart + xshift
      y = xstart + yshift
c
      do j = 1, 4
        do l = 1, 4
          xf(4*(l-1)+j) = x
          yf(4*(j-1)+l) = y
        end do
        x = x + xincrem
        y = y + xincrem
      end do
c
      x = xstart2 + xshift
      y = xstart2 + yshift
c
      do j = 1, 8
        do l = 1, 8
          xf2(8*(l-1)+j) = x
          yf2(8*(j-1)+l) = y
        end do
        x = x + xincrem2
        y = y + xincrem2
      end do
c
      return
      end

c
c*************************************************************************
      subroutine getError(xf,yf,xf2,yf2,level,a,error)
c*************************************************************************
c     the following subroutine calculates the l2 error between two vectors
c     that is needed in order to generate the tree structure.  this is used
c     to determine how accurate the approximating polynomial is.
c
c   input:
c     xf : the x values of 16 cell centered points in the box
c     yf : the y values of 16 cell centered points in the box
c     xf2 : the x values of 64 cell centered points in the box
c     yf2 : the y values of 64 cell centered points in the box
c     level : the level of the given box
c     a : the matrix that maps 16 function values
c       onto the 10 basis function coefficients
c
c   output:
c     error : the l2 error between the exact solution evaluated
c     at the 64 points and the approximation given by the polynomial
c     (coefficients found from the 16 points)
c
c   called from : formTree().
c
c   subroutine called :
c     forceFn() : the right hand side.
c     interpolate() : generates the polynomial coefficients.
c
c***************************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4  level
      real *8  a(10,16), error
      real *8  xf(16), yf(16)
      real *8  xf2(64), yf2(64)
c
c-----local variables
c
      integer *4  j, k, l
      real *8  coeffs(0:3,0:3)
      real *8  f(16), xlength
      real *8  vec1(64), vec2(64), sum
      real *8  xstart, xincrem
      real *8  x, y
      real *8  xx, yy
      real *8  xf1(64), yf1(64)
c
c-----external functions
c
      real *8 forceFn
      external forceFn
c
c-----now lets set the function values on the coarser grid:
c

      do j = 1, 16
        f(j) = forceFn(xf(j),yf(j))
      end do
c
      xstart  =  -0.4375d0
      xincrem =  0.125d0
c
      x = xstart
      y = xstart
c
      do j = 1, 8
        do l = 1, 8
          xf1(8*(l-1)+j) = x
          yf1(8*(j-1)+l) = y
        end do
        x = x + xincrem
        y = y + xincrem
      end do
c
c-----now call a routine giving us the coefficients that
c     match on the above grid:
c
      xlength = 1.0d0 / dble(2**level)
      call interpolate(a,f,coeffs,xlength)
c
c-----now let's calculate vec1 on the 64 point grid,
c     this vector corresponds to the function evaluated
c     at those 64 points:
c     (this is the 'correct' result.)
c
      do j = 1, 64
        vec1(j) = forceFn(xf2(j),yf2(j))
      end do
c
c-----now let's calculate vec2 on the 64 point grid,
c     this vector corresponds to the above polynomial
c     evaluated on the 64 points:
c     (this is the approximate result.)
c
      do j = 1, 64
        sum = 0.0d0
        yy = 1.0d0
        do l = 0, 3
          xx = yy
          do k = 0, 3 - l
            sum = sum + coeffs(l,k) * xf1(j)**l * yf1(j)**k * xx
            xx = xx * xlength
          end do
          yy = yy * xlength
        end do
        vec2(j) = sum
      end do
c
c-----now let's compute the error between the two vectors:
c
      error = 0.0d0
      do j = 1, 64
        error = error + (vec1(j) - vec2(j))**2
      end do
      error = dsqrt(error)
c
      return
      end




c************************************************************************
      subroutine interpolate(a, f, coeff, xlength)
c************************************************************************
c     the following subroutine is designed to generate the constants
c     in a polynomial approximating a given function.
c     on input, f is an array that is either 16 or 36 function values.
c     n is either 16 or 36 (depending on whether you are in the fourth
c     or sixth order case).  m is the degree of the polynomial that
c     you are approximateing the function with (3 or 5).  k is the
c     smaller dimension in the polynomial coefficient matrix (the total
c     number of coefficiets needed, 10 in the fourth order case,
c     21 in the sixth order case).  u, v, and sing are just the various
c     parts of the svd of the polynomial matrix and are precomputed.
c     xlength is the dimension of the box the function values are on.
c     on output the coefficents in the approximate polynomial are
c     contained in the array coeff.  (coeff(i,j) is the coefficient of
c     x^i*y^j.)
c
c     the following subroutine is designed to generate the constants
c     in a polynomial approximating a given function.
c
c     input:
c
c     f is the array of sixteen cell centered function values
c
c     xlength is the length of a side of the box
c     a is a precomputed matrix that represents a map from
c        16 cell centered function values to 10 basis functions
c
c     output:
c
c     coeffs is the array of coefficients for the basis functions
c
c************************************************************************
c
      implicit none
c
c-----global variables
c
      real *8 f(16), coeff(0:3,0:3)
      real *8 xlength, a(10,16)
c
c-----local variables
c
      integer i, j, l, ncntr
      real *8 x, y, xlen2
      real *8 sum(10)
c
c-----the matrix a is precomputed, it represents a map from the
c     16 points where the function is defined to a vector of
c     ten polynomial basis coefficients:
c     (a is just a 10 x 16 matrix.  it is derived by taking
c     an svd of the standard least sqaures problem that is
c     used to approximate the coefficients.)
c
      do i = 1, 10
        sum(i) = 0.0d0
        do j = 1, 16
          sum(i) = sum(i) + a(i,j) * f(j)
        end do
      end do
c
c-----now place the coefficients in a 2d array
c     indexed by corresponding powers of x and y
c     and scale them appropriately:
c
      ncntr = 1
      xlen2 = xlength/2.0d0
      x = 1.0d0
      do j = 0, 3
        y = x
        do l = 0, 3-j
          coeff(j,l) = sum(ncntr)/y
          ncntr = ncntr + 1
          y = y*xlen2
        end do
        x = x*xlen2
      end do
c
      return
      end




c*****************************************************************
      function forceFn(x,y)
      real *8 x,y,forceFN

      forceFN = exp(-1.d3*(x**2.d0 + y**2.d0)) + 
     1      exp(-1.d4*((x-3.d-1)**2.d0 + (y-3.d-1)**2.d0))


      return
      end



c*****************************************************************
      subroutine formRHS(rhs,levelbox,nboxes,icolbox,irowbox,
     1      ichildbox,nlev,nblevel,iboxlev,istartlev) 
      implicit none 
      
      integer *4 nboxes,nlev
      integer *4 icolbox(1),irowbox(1)
      integer *4 ichildbox(4,1),levelbox(1)
      integer *4 nblevel(0:1),iboxlev(1),istartlev(0:1)
      real *8 xf(16),yf(16)
      real *8 rhs(16,1)

      integer *4 i,j,k,l
      integer *4 itemp1,itemp2,itemp3
      integer *4 ibox
      real *8 x,y
      real *8 xstart,xincrem,xshift,yshift

      real *8 forceFn
      external forceFn
      
      itemp3=1
      do k=0,nlev
         itemp1=4*itemp3
         itemp2=2*itemp1
         do 100 i=istartlev(k),istartlev(k)+nblevel(k)-1
            ibox=iboxlev(i)
            if (ichildbox(1,ibox) .gt. 0) goto 100

            xstart=1.0d0/dble(itemp2)-0.5d0
            xincrem=1.0d0/dble(itemp1)
            xshift=dble(icolbox(ibox)-1)/dble(itemp3)
            yshift=dble(irowbox(ibox)-1)/dble(itemp3)

            x=xstart+xshift
            y=xstart+yshift
            do j=1,4
               do l=1,4
                  xf(4*(l-1)+j)=x
                  yf(4*(j-1)+l)=y
               enddo
               x=x+xincrem
               y=y+xincrem
            enddo

          do j = 1,16
            rhs(j,ibox) = forceFn(xf(j),yf(j))
          enddo
 100     continue
         itemp3=2*itemp3
      enddo



      return
      end






c*****************************************************************
      subroutine interpMatrix(a)
      implicit none
      real*8 a(10,16)

      a( 1, 1) =-0.9375000000000d-01
      a( 1, 2) = 0.6250000000000d-01
      a( 1, 3) = 0.6250000000000d-01
      a( 1, 4) =-0.9375000000000d-01
      a( 1, 5) = 0.6250000000000d-01
      a( 1, 6) = 0.2187500000000d+00
      a( 1, 7) = 0.2187500000000d+00
      a( 1, 8) = 0.6250000000000d-01
      a( 1, 9) = 0.6250000000000d-01
      a( 1,10) = 0.2187500000000d+00
      a( 1,11) = 0.2187500000000d+00
      a( 1,12) = 0.6250000000000d-01
      a( 1,13) =-0.9375000000000d-01
      a( 1,14) = 0.6250000000000d-01
      a( 1,15) = 0.6250000000000d-01
      a( 1,16) =-0.9375000000000d-01
      a( 2, 1) = 0.2083333333333d+00
      a( 2, 2) =-0.1666666666667d+00
      a( 2, 3) =-0.1666666666667d+00
      a( 2, 4) = 0.2083333333333d+00
      a( 2, 5) =-0.5000000000000d+00
      a( 2, 6) =-0.6250000000000d+00
      a( 2, 7) =-0.6250000000000d+00
      a( 2, 8) =-0.5000000000000d+00
      a( 2, 9) = 0.5000000000000d+00
      a( 2,10) = 0.6250000000000d+00
      a( 2,11) = 0.6250000000000d+00
      a( 2,12) = 0.5000000000000d+00
      a( 2,13) =-0.2083333333333d+00
      a( 2,14) = 0.1666666666667d+00
      a( 2,15) = 0.1666666666667d+00
      a( 2,16) =-0.2083333333333d+00
      a( 3, 1) = 0.2500000000000d+00
      a( 3, 2) = 0.2500000000000d+00
      a( 3, 3) = 0.2500000000000d+00
      a( 3, 4) = 0.2500000000000d+00
      a( 3, 5) =-0.2500000000000d+00
      a( 3, 6) =-0.2500000000000d+00
      a( 3, 7) =-0.2500000000000d+00
      a( 3, 8) =-0.2500000000000d+00
      a( 3, 9) =-0.2500000000000d+00
      a( 3,10) =-0.2500000000000d+00
      a( 3,11) =-0.2500000000000d+00
      a( 3,12) =-0.2500000000000d+00
      a( 3,13) = 0.2500000000000d+00
      a( 3,14) = 0.2500000000000d+00
      a( 3,15) = 0.2500000000000d+00
      a( 3,16) = 0.2500000000000d+00
      a( 4, 1) =-0.3333333333333d+00
      a( 4, 2) =-0.3333333333333d+00
      a( 4, 3) =-0.3333333333333d+00
      a( 4, 4) =-0.3333333333333d+00
      a( 4, 5) = 0.1000000000000d+01
      a( 4, 6) = 0.1000000000000d+01
      a( 4, 7) = 0.1000000000000d+01
      a( 4, 8) = 0.1000000000000d+01
      a( 4, 9) =-0.1000000000000d+01
      a( 4,10) =-0.1000000000000d+01
      a( 4,11) =-0.1000000000000d+01
      a( 4,12) =-0.1000000000000d+01
      a( 4,13) = 0.3333333333333d+00
      a( 4,14) = 0.3333333333333d+00
      a( 4,15) = 0.3333333333333d+00
      a( 4,16) = 0.3333333333333d+00
      a( 5, 1) = 0.2083333333333d+00
      a( 5, 2) =-0.5000000000000d+00
      a( 5, 3) = 0.5000000000000d+00
      a( 5, 4) =-0.2083333333333d+00
      a( 5, 5) =-0.1666666666667d+00
      a( 5, 6) =-0.6250000000000d+00
      a( 5, 7) = 0.6250000000000d+00
      a( 5, 8) = 0.1666666666667d+00
      a( 5, 9) =-0.1666666666667d+00
      a( 5,10) =-0.6250000000000d+00
      a( 5,11) = 0.6250000000000d+00
      a( 5,12) = 0.1666666666667d+00
      a( 5,13) = 0.2083333333333d+00
      a( 5,14) =-0.5000000000000d+00
      a( 5,15) = 0.5000000000000d+00
      a( 5,16) =-0.2083333333333d+00
      a( 6, 1) = 0.3600000000000d+00
      a( 6, 2) = 0.1200000000000d+00
      a( 6, 3) =-0.1200000000000d+00
      a( 6, 4) =-0.3600000000000d+00
      a( 6, 5) = 0.1200000000000d+00
      a( 6, 6) = 0.4000000000000d-01
      a( 6, 7) =-0.4000000000000d-01
      a( 6, 8) =-0.1200000000000d+00
      a( 6, 9) =-0.1200000000000d+00
      a( 6,10) =-0.4000000000000d-01
      a( 6,11) = 0.4000000000000d-01
      a( 6,12) = 0.1200000000000d+00
      a( 6,13) =-0.3600000000000d+00
      a( 6,14) =-0.1200000000000d+00
      a( 6,15) = 0.1200000000000d+00
      a( 6,16) = 0.3600000000000d+00
      a( 7, 1) =-0.6000000000000d+00
      a( 7, 2) =-0.2000000000000d+00
      a( 7, 3) = 0.2000000000000d+00
      a( 7, 4) = 0.6000000000000d+00
      a( 7, 5) = 0.6000000000000d+00
      a( 7, 6) = 0.2000000000000d+00
      a( 7, 7) =-0.2000000000000d+00
      a( 7, 8) =-0.6000000000000d+00
      a( 7, 9) = 0.6000000000000d+00
      a( 7,10) = 0.2000000000000d+00
      a( 7,11) =-0.2000000000000d+00
      a( 7,12) =-0.6000000000000d+00
      a( 7,13) =-0.6000000000000d+00
      a( 7,14) =-0.2000000000000d+00
      a( 7,15) = 0.2000000000000d+00
      a( 7,16) = 0.6000000000000d+00
      a( 8, 1) = 0.2500000000000d+00
      a( 8, 2) =-0.2500000000000d+00
      a( 8, 3) =-0.2500000000000d+00
      a( 8, 4) = 0.2500000000000d+00
      a( 8, 5) = 0.2500000000000d+00
      a( 8, 6) =-0.2500000000000d+00
      a( 8, 7) =-0.2500000000000d+00
      a( 8, 8) = 0.2500000000000d+00
      a( 8, 9) = 0.2500000000000d+00
      a( 8,10) =-0.2500000000000d+00
      a( 8,11) =-0.2500000000000d+00
      a( 8,12) = 0.2500000000000d+00
      a( 8,13) = 0.2500000000000d+00
      a( 8,14) =-0.2500000000000d+00
      a( 8,15) =-0.2500000000000d+00
      a( 8,16) = 0.2500000000000d+00
      a( 9, 1) =-0.6000000000000d+00
      a( 9, 2) = 0.6000000000000d+00
      a( 9, 3) = 0.6000000000000d+00
      a( 9, 4) =-0.6000000000000d+00
      a( 9, 5) =-0.2000000000000d+00
      a( 9, 6) = 0.2000000000000d+00
      a( 9, 7) = 0.2000000000000d+00
      a( 9, 8) =-0.2000000000000d+00
      a( 9, 9) = 0.2000000000000d+00
      a( 9,10) =-0.2000000000000d+00
      a( 9,11) =-0.2000000000000d+00
      a( 9,12) = 0.2000000000000d+00
      a( 9,13) = 0.6000000000000d+00
      a( 9,14) =-0.6000000000000d+00
      a( 9,15) =-0.6000000000000d+00
      a( 9,16) = 0.6000000000000d+00
      a(10, 1) =-0.3333333333333d+00
      a(10, 2) = 0.1000000000000d+01
      a(10, 3) =-0.1000000000000d+01
      a(10, 4) = 0.3333333333333d+00
      a(10, 5) =-0.3333333333333d+00
      a(10, 6) = 0.1000000000000d+01
      a(10, 7) =-0.1000000000000d+01
      a(10, 8) = 0.3333333333333d+00
      a(10, 9) =-0.3333333333333d+00
      a(10,10) = 0.1000000000000d+01
      a(10,11) =-0.1000000000000d+01
      a(10,12) = 0.3333333333333d+00
      a(10,13) =-0.3333333333333d+00
      a(10,14) = 0.1000000000000d+01
      a(10,15) =-0.1000000000000d+01
      a(10,16) = 0.3333333333333d+00
c
      return
      end

c*****************************************************************
      subroutine writeTree2Matlab(nboxes,ichildbox,
     1    irowbox,icolbox,levelbox,fileName)
      implicit none

      integer *4 nboxes
      integer *4 ichildbox(4,nboxes)
      integer *4 icolbox(nboxes)
      integer *4 irowbox(nboxes)
      integer *4 levelbox(nboxes)
      character(len=*) :: fileName


      integer *4 i
      real *8 radius
      complex *16 cen
      complex *16 eye

      data eye/(0.d0,1.d0)/ 


      open(unit=1,file=fileName)
      write(1,*) 'clf'
      write(1,*) 'hold on'
      do 1900 i=1,nboxes
         if (ichildbox(1,i) .ge. 0) goto 1900
c     if the box has children, no need to plot it
         cen=dble(2*icolbox(i)-1)+eye*dble(2*irowbox(i)-1)
         cen=cen*(0.5d0)**(levelbox(i)+1)-(0.5d0,0.5d0)
         radius=(0.5d0)**(levelbox(i)+1)
         write(1,*) 'z1=[',dreal(cen)-radius, dreal(cen)+radius,'];'
         write(1,*) 'z2=[',dimag(cen)+radius, dimag(cen)+radius,'];'
         write(1,*) 'plot(z1,z2,''k'')'
         write(1,*) 'z1=[',dreal(cen)-radius, dreal(cen)+radius,'];'
         write(1,*) 'z2=[',dimag(cen)-radius, dimag(cen)-radius,'];'
         write(1,*) 'plot(z1,z2,''k'')'
         write(1,*) 'z1=[',dreal(cen)-radius, dreal(cen)-radius,'];'
         write(1,*) 'z2=[',dimag(cen)-radius, dimag(cen)+radius,'];'
         write(1,*) 'plot(z1,z2,''k'')'
         write(1,*) 'z1=[',dreal(cen)+radius, dreal(cen)+radius,'];'
         write(1,*) 'z2=[',dimag(cen)-radius, dimag(cen)+radius,'];'
         write(1,*) 'plot(z1,z2,''k'')'
         write(1,*) ''
 1900 continue
      write(1,*) 'axis equal'
      write(1,*) 'axis([-1 1 -1 1]/1.9)'
      close(unit=1)


      return
      end
