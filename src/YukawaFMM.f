      subroutine particleFMM(npts,maxptsperbox,layers,dom,density,beta,
     $    pot,potx,poty,potxx,potxy,potyy)
c     ****** DESCRPTION:
c     Compute the Yukawa potential and its gradient at all locations
c     Yukawa potential at x is \phi(x) = \sum q_{j}K_{0}(\beta|x-x_j|)
c     Note that our previous defined alpha is 1/beta

c     ****** INPUT PARAMETERS:
c     npts - Number of target/source points.  Sources = Targets
c     maxptsperbox - Maximum number of points in a box set by the FMM
c        There will be boxes with too many points if it requires a tree
c        of greater depth than maxlevel
c     layers - number of layer potentials to compute
c       layers = 1 => only compute single-layer potential
c       layers = 2 => compute single- and double-layer potential
c       layers = 3 => compute single-, double-, and triple-layer
c       potentials
c     dom - location of the targets/sources stored as complex variables
c     density - charges of the source points stores as doubles
c     beta - single parameter in the Yukawa potential.  For time
c     stepping applications, beta^2 ~ 1/dt

c     ******* OUTPUT PARAMETERS:  
c     pot - the Yukawa potential
c     potx - the first component of the gradient of the Yukawa 
c        potential
c     poty - the second component of the gradient of the Yukawa
c        potential
c     potxx - the first component of the second derivative of the      
c        Yukawa potential
c     potxy - the mixed component of the second derivative of the      
c        Yukawa potential
c     potyy - the second component of the second derivative of the      
c        Yukawa potential
      implicit none

      integer npts,layers
      complex *16 dom(npts)
      real *8 density(npts)
      real *8 beta
      complex *16 pot(npts)
      complex *16 potx(npts),poty(npts)
      complex *16 potxx(npts),potxy(npts),potyy(npts)

      integer maxpts,maxlevel,maxboxes,maxwrk
      parameter (maxpts = 2**18)
c     Maximum number of targets/sources/charges
      parameter (maxlevel = 15)
c     Maximum depth of the tree
      parameter (maxboxes = 100000)
c     Maximum number of boxes in tree
      parameter (maxwrk = 100000000)

      integer *4 nlev,nboxes
      integer *4 maxptsperbox
      integer *4 levelbox(maxboxes)
      integer *4 iparentbox(maxboxes)
      integer *4 ichildbox(4,maxboxes)
      integer *4 icolbox(maxboxes),irowbox(maxboxes)
      integer *4 icolleagbox(9,maxboxes)
      integer *4 nblevel(0:maxlevel)
      integer *4 iboxlev(maxboxes)
      integer *4 istartlev(0:maxlevel)
      integer *4 iflag(maxboxes)
      integer *4 itemparray(maxboxes)
      integer *4 iwork(maxwrk)
      real *8 work(maxwrk)
c     variables for the tree structure

      integer *4 fb(maxpts)
      integer *4 numptsbox(maxboxes)
      integer *4 permute(maxpts)
      integer *4 fpt(maxboxes)
c     variables relating the points to the tree structure
      
      integer *4 list1(50,maxboxes)
      integer *4 list3(50,maxboxes)
      integer *4 nlist1(maxboxes),nlist3(maxboxes)
c     Need to store what Greengard et al. called list1 and list3.
c     Huang et al. didn't need these since they only exist if the
c     tree has neighbouring boxes that are more than one level apart
      
      integer i,k,ic_sum
      complex *16 eye,cen
      real*8 radius
c     Some local variables

      real *4 ETIME,timeep(2)
      real *4 t0_tree,t1_tree
c     For doing some timings

      eye = (0.d0,1.d0)

c     Initialize some tables for doing FMM translations
      call inittable

c     Prameters into create_tree
c     levelbox - Level of a box in the tree
c     iparentbox - Parent of each box
c     ichildbox - 4 children of each box
c     icolbox,irowbox - column and row of box
c     icolleagbox - 9 colleagues (neighbours) of a box
c     nboxes - total number of boxes
c     numptsbox - number of points in a box
c     nlev - total number of levels
c     nblevel - number of boxes at each level
c     fpt - index of the first point in the box
c     iboxlev - 
c     istartlev - first box at each level
c     itemparray - 
c     maxlevel - maximum allowable level
c     maxpts - maximum allowable number of points
c     maxboxes - maxiumum total number of boxes
c     dom - location of target points in complex form
c         - need to be contained in [-1/2 1/2 -1/2 1/2]
c     npts - number of target points
c     maxptsperbox - maxiumum number of points in a box
c     permute - used to find points inside a box
c     fb - finest box containing a given point
      t0_tree = ETIME(timeep)
      call create_tree(levelbox,iparentbox,ichildbox,icolbox,
     $   irowbox,icolleagbox,nboxes,numptsbox,nlev,nblevel,fpt,
     $   iboxlev,istartlev,iflag,itemparray,maxlevel,maxpts,maxboxes,
     $   dom,npts,maxptsperbox,permute,fb)
c     Create tree.  No restriction on neighbouring boxes.  They can
c     be as many levels apart from one another as required.

      t1_tree = ETIME(timeep)
      write(6,1020) 'Time to create tree is ',
     $    t1_tree-t0_tree,' seconds'

      write(6,*) nboxes, 'boxes before the tree is compressed'
      write(6,*) nlev, 'levels in tree'

      if (nboxes .gt. maxboxes) then
         write(6,*) 'There are more boxes then are allowed'
         write(6,*) 'for in the memory space.'
         write(6,*) 'Readjust the memory.'
         write(6,*) 'I am stopping'
         stop
      endif
c     Check to make sure we have allocated enough space for the tree

      t0_tree = ETIME(timeep)
      call compress_tree(levelbox,iparentbox,ichildbox,icolbox,
     $     irowbox,icolleagbox,iboxlev,nboxes,numptsbox,
     $     nlev,nblevel,fpt,istartlev,fb,npts,permute)
c     compress the tree by eliminating all boxes with no points inside 
c     of them.
c     readjust all the indicies appropriatly
      t1_tree = ETIME(timeep)
      write(6,1020) 'Time to compress tree is ',
     $    t1_tree-t0_tree,' seconds'


      write(6,*) nboxes, 'boxes after the tree is compressed'
      


      open(unit=1,file='../output/buildtree.m')
      write(1,*) 'clf'
      write(1,*) 'hold on'
      write(1,*) 'axis([-1 1 -1 1]/1.9)'
      do i=1,nboxes
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
      enddo
      write(1,*) 'axis equal'
      write(1,*) 'xlim([-.53 .53])'
      write(1,*) 'ylim([-.53 .53])'
      write(1,*) 'dom_x=['
      do i=1,npts
         write(1,*) dreal(dom(i))
      enddo
      write(1,*) '];'
      write(1,*) 'dom_y=['
      do i=1,npts
         write(1,*) dimag(dom(i))
      enddo
      write(1,*) '];'
      write(1,*) ''
c      write(1,*) 'plot(dom_x(1:',
c     $     npts,'),dom_y(1:',
c     $     npts,'),''r.'');'
      write(1,*) 'set(gca,''xtick'',[]);'
      write(1,*) 'set(gca,''ytick'',[]);'
      write(1,*) 'set(gca,''xcolor'',''w'');'
      write(1,*) 'set(gca,''ycolor'',''w'');'
      close(unit=1)
c     Write a matlab script that can be used to visualize the tree

      t0_tree = ETIME(timeep)
      do i=1,nboxes
         ic_sum=ichildbox(1,i)+ichildbox(2,i)+
     $        ichildbox(3,i)+ichildbox(4,i)
         if (ic_sum .eq. -4) then
            call mklists13(i,list1,nlist1,list3,nlist3,
     $           ichildbox,icolleagbox)
            if (nlist1(i) .gt. 50 .or. nlist3(i) .ge. 50) then
              print*,'TOO MANY BOXES IN LIST1 OR LIST3'
              print*,'I AM STOPPING'
              stop
            endif
         endif
      enddo
c     Build list1 and list3 which have to do with descendents of
c     neighbouring boxes that are not well-seperated, etc...

      t1_tree = ETIME(timeep)
      write(6,1020) 'Time to form List1 and List 3 is ',
     $    t1_tree-t0_tree,' seconds'


c     Start allocating memory for FMM
      call startParticleFMM(
     1    beta,work,maxwrk,iwork,maxwrk,
     2    nlev,levelbox,iparentbox,ichildbox,
     3    icolbox,irowbox,nboxes,nblevel,iboxlev,
     4    istartlev,fb,numptsbox,permute,fpt,
     5    list1,nlist1,list3,nlist3,
     6    npts,layers,dom,density,pot,potx,poty,
     7    potxx,potxy,potyy)




 1020 format(A,ES8.2,A)

      return
      end




c**********************************************************************
      subroutine startParticleFMM(
     1     beta,work,lenw,iwork,ilen,
     2     nlev,levelbox,iparentbox,ichildbox,
     3     icolbox,irowbox,nboxes,nblevel,iboxlev,
     4     istartlev,fb,numptsbox,permute,fpt,
     5     list1,nlist1,list3,nlist3,
     6     npts,layers,dom,q,pot,potx,poty,
     7     potxx,potxy,potyy)
      implicit none

      integer *4 levelbox(nboxes),iparentbox(nboxes)
      integer *4 ichildbox(4,nboxes)
      integer *4 irowbox(nboxes),icolbox(nboxes)
      integer *4 nblevel(0:nlev),iboxlev(nboxes)
      integer *4 istartlev(0:nlev)
      integer *4 fb(npts)
      integer *4 numptsbox(nboxes)
      integer *4 permute(npts)
      integer *4 fpt(nboxes)
      integer *4 iwork(ilen)
      real *8 work(lenw)
      real *8 q(npts)
      complex *16 pot(npts)
      complex *16 potx(npts),poty(npts)
      complex *16 potxx(npts),potxy(npts),potyy(npts)
      complex *16 dom(npts)
      integer *4 list1(50,nboxes),list3(50,nboxes)
      integer *4 nlist1(nboxes),nlist3(nboxes)
      real *8 beta
      integer *4 lenw,ilen
      integer *4 nlev,nboxes,npts
      integer layers
c     Number of layer potentials to compute

     

c     Local variables
      integer *4 nterms,iprec
      integer *4 nnodes
      integer *4 nnodesmax 
      parameter (nnodesmax=50)
      integer *4 mpole,mxnodes,mcolleag,mchoose,mcomp
      integer *4 lzs,mzs,locexp,llocxp,lmpole,lexp
      integer *4 mcomp2,lscale,lcomp,lcolleag,lnodes
      integer *4 itot,iscale,lchoose
      integer *4 iexpn,iexpe,iexps,iexpw
      integer *4 ibetascal

      nterms=41
      iprec=12
c     Result if iprec=3 from Huang's code...avoids the 'if' statement

c     carve up the integer workspace
      lcolleag=9*nboxes
      mcolleag=1
      itot=mcolleag+lcolleag

      if (itot .ge. lenw) then
         write(6,*) 'The integer workspace required', itot
         write(6,*) 'The total workspace needed exeeds the'
         write(6,*) 'amount alloted.  Please readjust maxwork'
         write(6,*) 'in the beginning of the code and try again.'
         stop
      endif   
      
      lmpole=2*(nterms+1)*nboxes
      llocxp=2*(nterms+1)*nboxes
      lexp=2*nnodesmax*nboxes
      lnodes=nnodesmax
      lscale=nlev+1
      lcomp=nnodesmax*(2*nterms+1)*2
      lchoose=2*(nterms+1)*(2*nterms+1)
      lzs=2*7*7*nnodesmax

c     carve up the real workspace
      iscale=1
      ibetascal=iscale+lscale
      mpole=ibetascal+lscale
      locexp=mpole+lmpole
      iexpn=locexp+llocxp
      iexps=iexpn+lexp
      iexpe=iexps+lexp
      iexpw=iexpe+lexp
      mcomp=iexpw+lexp
      mcomp2=mcomp+lcomp
      mchoose=mcomp2+lcomp
      mzs=mchoose+lchoose
      mxnodes=mzs+lzs

      itot=mxnodes+lnodes

      if (itot .ge. lenw) then
         write(6,*) 'The real workspace required', itot
         write(6,*) 'The total workspace needed exeeds the'
         write(6,*) 'amount alloted.  Please readjust maxwork'
         write(6,*) 'in the beginning of the code and try again.'
         stop
      endif

c     Call the code that does the bulk of the FMM code
      call evalParticleFMM(
     1     beta,work(ibetascal),work(iscale),
     2     12,nlev,work(mxnodes),work(mzs),work(mcomp),
     3     work(mcomp2),41,nnodes,npts,layers,dom,q,
     4     pot,potx,poty,potxx,potxy,potyy,
     5     work(mpole),work(locexp),work(iexpn),work(iexps),
     6     work(iexpe),work(iexpw),work(mchoose),
     7     levelbox,iparentbox,ichildbox,
     8     icolbox,irowbox,iwork(mcolleag),nboxes,
     9     nblevel,iboxlev,istartlev,fb,numptsbox,permute,fpt,
     1     list1,nlist1,list3,nlist3)



      return
      end




c**********************************************************************
      subroutine evalParticleFMM(
     1    beta,betascal,scale,
     2    iprec,nlev,xnodes,zs,comp,
     3    comp2,nterms,nnodes,npts,layers,dom,q,
     4    pot,potx,poty,potxx,potxy,potyy,
     5    mpole,locexp,expn,exps,
     6    expe,expw,c,
     7    levelbox,iparentbox,ichildbox,
     8    icolbox,irowbox,icolleagbox,nboxes,
     9    nblevel,iboxlev,istartlev,fb,numptsbox,permute,fpt,
     1    list1,nlist1,list3,nlist3)
      implicit none

      integer *4 nnodesmax
      parameter (nnodesmax=50)

c     global variables
      integer *4 icolbox(nboxes),irowbox(nboxes)
      integer *4 nblevel(0:nlev),iboxlev(nboxes)
      integer *4 istartlev(0:nlev)
      integer *4 levelbox(nboxes),iparentbox(nboxes)
      integer *4 ichildbox(4,nboxes),icolleagbox(9,nboxes)
      integer *4 fb(npts),numptsbox(nboxes)
      integer *4 permute(npts),fpt(nboxes)
      integer *4 list1(50,nboxes),list3(50,nboxes)
      integer *4 nlist1(nboxes),nlist3(nboxes)

      real *8 q(npts)
      complex *16 pot(npts)
      complex *16 potx(npts),poty(npts)
      complex *16 potxx(npts),potxy(npts),potyy(npts)
      complex *16 dom(npts)

      integer *4 iprec,nlev,npts,nboxes,nterms,layers

      real *8 beta
      real*8 betascal(0:nlev),scale(0:nlev)
      real*8 xnodes(1)
      real*8 comp(nnodesmax,-nterms:nterms)
      real*8 comp2(nnodesmax,-nterms:nterms)
      complex *16 c(0:nterms,-nterms:nterms)
      complex *16 mpole(0:nterms,1)
      complex *16 locexp(0:nterms,1)
      complex *16 expn(nnodesmax,1),exps(nnodesmax,1)
      complex *16 expe(nnodesmax,1),expw(nnodesmax,1)
      complex *16 zs(-3:3,-3:3,nnodesmax)

c     local variables
      integer inall(4),iynall(4)
      integer isall(4),iysall(4)
      integer ieall(4),iyeall(4)
      integer iwall(4),iywall(4)
      integer in12(4),iy12(4)
      integer  is34(4),iy34(4)
      integer  iw24(2),iy24(2)
      integer  ie13(2),iy13(2)
      integer  iw4(2),iy4(2)
      integer  iw2(2),iy2(2)
      integer  ie1(2),iy1(2)
      integer  ie3(2),iy3(2)
      integer  in123(1),iy123(1)
      integer  in124(1),iy124(1)
      integer  is134(1),iy134(1)
      integer  is234(1),iy234(1)
      complex *16 eye,czero,z,unor
      complex *16  mexnall(1:50),mexn12(1:50),mexn123(50),mexn124(50)
      complex *16  mexsall(1:50),mexs34(1:50),mexs134(50),mexs234(50)
      complex *16  mexeall(1:50),mexe13(1:50),mexe1(1:50),mexe3(1:50)
      complex *16  mexwall(1:50),mexw24(1:50),mexw2(1:50),mexw4(1:50)

      real*8 pi2
      integer *4 i,ii,j,k,ibox
      integer *4 ic1,ic2,ic3,ic4
      integer *4 ichildsum
      integer *4 nnodes
      integer *4 istart,iend,istart1,iend1,istart2,iend2
      integer *4 nsall,ns134,ns234,nn124,nn12,nn123,nnall
      integer *4 nw4,nwall,nw24,ns34,nw2,ne3,neall,ne1,ne13
      integer *4 nb,lk,kk,jj,iout
      integer *4 iflag1,iflag2,iflag3,iflag4
      real*8 betah

c     A bunch of variables so that we can do some timings of 
c     several of the FMM components
      real *4 timeep(2),ETIME
      real *4 t1,t0
      real *4 t_L2L,t1_L2L,t0_L2L
      real *4 t_M2P,t1_M2P,t0_M2P
      real *4 t_P2P,t1_P2P,t0_P2P
      real *4 t_P2L,t1_P2L,t0_P2L
      real *4 t_direct,t1_direct,t0_direct
      real *4 t_evalmp,t1_evalmp,t0_evalmp
      real *4 t_updateloc,t1_updateloc,t0_updateloc
      real *4 t_evalloc,t1_evalloc,t0_evalloc

      data czero/(0.0d0,0.0d0)/
      data eye/(0.0d0,1.0d0)/
      pi2 = 2.0d0*dacos(-1.0d0)

      t_L2L = 0.d0
      t_M2P = 0.d0
      t_P2P = 0.d0
      t_P2L = 0.d0
      t_direct = 0.d0
      t_evalmp = 0.d0
      t_updateloc = 0.d0
      t_evalloc = 0.d0
c     Initalize all timings to zero

      betascal(0)=beta
      scale(0)=1.0d0
      do i=1,nlev
         betascal(i)=betascal(i-1)/2.0d0
         scale(i)=scale(i-1)/2.0d0
      enddo
c     needed for scaling besselk and besseli

      do i=1,npts
         pot(i)=czero
         potx(i) = czero
         poty(i) = czero
         potxx(i) = czero
         potxy(i) = czero
         potyy(i) = czero
      enddo
c     initalize the potential to zero


      t0 = ETIME(timeep)
      call mkcolls(icolbox,
     $     irowbox,icolleagbox,nboxes,nlev,
     $     iparentbox,ichildbox,nblevel,
     $     iboxlev,istartlev)
c     make the colleagues of each box.  We did it early
c     to make list1 and list3 but don't have it here.
c     TODO: This should be fixed.  Only want to call once
      t1 = ETIME(timeep)
c      write(6,1020) 'Time recreating the colleagues is ',
c     $    t1-t0,' seconds'


      do i=1,nboxes
         do j=0,nterms
            mpole(j,i)=czero
            locexp(j,i)=czero
         enddo
         do j=1,nnodesmax
            expe(j,i)=czero
            expw(j,i)=czero
            expn(j,i)=czero
            exps(j,i)=czero
         enddo
      enddo
c     initalize everything to zero

      t0 = ETIME(timeep)
      do i=nlev,0,-1
         if (betascal(0)*scale(i) .le. 36.0d0) then
            istart=istartlev(i)
            iend=istart+nblevel(i)-1
            if (i .lt. nlev) then
               call ympmp_coef(1,nterms,scale,betascal,i+1,c)
c     terms needed for mutlipole to multipole translation
            endif

            do ii=istart,iend
               j=iboxlev(ii)
               ichildsum=ichildbox(1,j)+ichildbox(2,j)+
     $              ichildbox(3,j)+ichildbox(4,j)
               if (ichildsum .eq. -4) then
                  istart2=fpt(j)
                  iend2=istart2+numptsbox(j)-1
                  do k=istart2,iend2
                     call ymultipole(nterms,mpole(0,j),
     $                    dom(permute(k)),
     $                    q(permute(k)),levelbox(j),
     $                    icolbox(j),irowbox(j),betascal)
c     update the multipole coefficients of the current childless box
                  enddo
               else
                  ic1=ichildbox(2,j)
                  ic2=ichildbox(3,j)
                  ic3=ichildbox(4,j)
                  ic4=ichildbox(1,j)

                  call ychildpar(mpole(0,j),mpole(0,ic1),
     $                 mpole(0,ic2),mpole(0,ic3),
     $                 mpole(0,ic4),nterms,c,ic1,ic2,ic3,ic4)
c     shift the (upto) four multipole expansions of the children to
c     the current box
               endif
            enddo
         endif
      enddo
c     At this point, we have the multipole coefficients for each 
c     of the boxes
      t1 = ETIME(timeep)
c      write(6,1020) 'Time forming multipole coefficients is ',
c     $    t1-t0,' seconds'



      do i=0,nlev
         betah=scale(i)*betascal(0)
         call ytata_coef(1,nterms,scale,betascal,i,c)
         if (i .lt. nlev) then
            call ymppw_coef(1,nterms,iprec,nnodes,scale,betascal,i+1,
     $           xnodes,comp)
            call ypwpw_coef(1,nnodes,scale,betascal,i+1,xnodes,zs)
c     coefficients for doing multipole to plane wave and plane wave to
c     plane wave
         endif
         istart=istartlev(i)
         iend=istart+nblevel(i)-1
         do 100 ii=istart,iend
            j=iboxlev(ii)
            ichildsum=ichildbox(1,j)+ichildbox(2,j)+
     $           ichildbox(3,j)+ichildbox(4,j)
            if (ichildsum .ne. -4) then
               ic1=ichildbox(2,j)
               iflag1=1
               ic2=ichildbox(3,j)
               iflag2=1
               ic3=ichildbox(4,j)
               iflag3=1
               ic4=ichildbox(1,j)
               iflag4=1
               if (ic1 .eq. -1) then
                  iflag1=0
               endif
               if (ic2 .eq. -1) then
                  iflag2=0
               endif
               if (ic3 .eq. -1) then
                  iflag3=0
               endif
               if (ic4 .eq. -1) then
                  iflag4=0
               endif
c
               t0_L2L = ETIME(timeep)
               call yparentchild(locexp(0,j),locexp(0,ic1),
     $              locexp(0,ic2),locexp(0,ic3),locexp(0,ic4),
     $              iflag1,iflag2,iflag3,iflag4,nterms,c)
c     shift the local expansion to the (upto) four children
               t1_L2L = ETIME(timeep)
               t_L2L = t_L2L + t1_L2L - t0_L2L

               t0_M2P = ETIME(timeep)
               call ymkexp2d(j,nterms,mpole,
     $              nnodes,mexnall,mexn12,mexn123,mexn124,
     $              mexsall,mexs34,mexs134,mexs234,mexeall,
     $              mexe13,mexe1,mexe3,mexwall,mexw24,mexw2,mexw4,
     $              zs,comp,ichildbox,irowbox,icolbox)
c     convert the multipole coefficients for the current box to 
c     plane wave expansions that converge in the different directions
               t1_M2P = ETIME(timeep)
               t_M2P = t_M2P + t1_M2P - t0_M2P

               call ymklists(j,inall,nnall,iynall,
     1              in12,nn12,iy12,in123,nn123,iy123,in124,nn124,iy124,
     2              isall,nsall,iysall,
     3              is34,ns34,iy34,is134,ns134,iy134,is234,ns234,iy234,
     4              ieall,neall,iyeall,ie13,ne13,iy13,
     5              iwall,nwall,iywall,iw24,nw24,iy24,
     6              iw2,iy2,nw2,iw4,iy4,nw4,
     7              ie1,iy1,ne1,ie3,iy3,ne3,
     8              icolleagbox,ichildbox,
     9              icolbox, irowbox,levelbox(j))
c     form all of the elements in the interaction list that live at the
c     same level as box currently being considered.

               t0_P2P = ETIME(timeep)
               call yprocessno(expn,inall,nnall,iynall,
     1              in12,nn12,iy12,in123,nn123,iy123,in124,nn124,iy124,
     2              mexnall,mexn12,mexn123,mexn124,zs,nnodes)
               call yprocessso(exps,isall,nsall,iysall,
     1              is34,ns34,iy34,is134,ns134,iy134,is234,ns234,iy234,
     2              mexsall,mexs34,mexs134,mexs234,zs,nnodes)
               call yprocessea(expe,ieall,neall,iyeall,
     1              ie13,ne13,iy13,ie1,ne1,iy1,ie3,ne3,iy3,
     2              mexeall,mexe13,mexe1,mexe3,zs,nnodes)
               call yprocesswe(expw,iwall,nwall,iywall,
     1              iw24,nw24,iy24,iw2,nw2,iy2,iw4,nw4,iy4,
     2              mexwall,mexw24,mexw2,mexw4,zs,nnodes)
c     update plane-wave expansions at all the boxes in the interaction 
c     list due to current box
               t1_P2P = ETIME(timeep)
               t_P2P = t_P2P + t1_P2P - t0_P2P

            else
               do 250 nb=1,9
                  iout=icolleagbox(nb,j)

                  if (iout .le. 0) goto 250
c     if the colleague doesn't exist, there is nothing to do
                  ichildsum=ichildbox(1,iout)+ichildbox(2,iout)+
     $                 ichildbox(3,iout)+ichildbox(4,iout)
                  if (ichildsum .eq. -4) then
                    t0_direct = ETIME(timeep)
                    call direct_eval(npts,layers,beta,dom,q,permute,
     $                  fpt(j),fpt(iout),numptsbox(j),numptsbox(iout),
     $                  pot,potx,poty,potxx,potxy,potyy)
                    t1_direct = ETIME(timeep)
                    t_direct = t_direct + t1_direct - t0_direct
c     if the colleague has no children, evaluate potential directly
                  endif
 250           continue

               do k=1,nlist1(j)
                  lk=list1(k,j)
                    t0_direct = ETIME(timeep)
                    call direct_eval(npts,layers,beta,dom,q,permute,
     $                  fpt(j),fpt(lk),numptsbox(j),numptsbox(lk),
     $                  pot,potx,poty,potxx,potxy,potyy)
                    t1_direct = ETIME(timeep)
                    t_direct = t_direct + t1_direct - t0_direct

                    t0_direct = ETIME(timeep)
                    call direct_eval(npts,layers,beta,dom,q,permute,
     $                  fpt(lk),fpt(j),numptsbox(lk),numptsbox(j),
     $                  pot,potx,poty,potxx,potxy,potyy)
                    t1_direct = ETIME(timeep)
                    t_direct = t_direct + t1_direct - t0_direct
               enddo
c     potential due to list1 is done directly.  Need to do both
c     directions as list1 only contains boxes that are finer 
c     than box j
               do k=1,nlist3(j)
                  lk=list3(k,j)
                  istart1=fpt(j)
                  iend1=istart1+numptsbox(j)-1
c     How list 3 affects current box
                  do kk=istart1,iend1
                     z=dom(permute(kk))
                     t0_evalmp = ETIME(timeep)
                     call mpole_eval(layers,z,pot(permute(kk)),
     $                    potx(permute(kk)),poty(permute(kk)),
     $                    potxx(permute(kk)),potxy(permute(kk)),
     $                    potyy(permute(kk)),mpole(0,lk),
     $                    icolbox(lk),irowbox(lk),levelbox(lk),
     $                    beta,betascal(levelbox(lk)),nterms)
                     t1_evalmp = ETIME(timeep)
                     t_evalmp = t_evalmp + t1_evalmp - t0_evalmp
                  enddo
c   Potential due to points in list 3 are computed by evaluating
c   the multipole expansion
                  
                  do kk=istart1,iend1
                     z=dom(permute(kk))

                    t0_updateloc = ETIME(timeep)
                     call update_local_coe(z,q(permute(kk)),
     $                 icolbox(lk),irowbox(lk),levelbox(lk),betascal,
     $                 nterms,locexp(0,lk))
c     How current box affects list 3
                    t1_updateloc = ETIME(timeep)
                    t_updateloc = t_updateloc + 
     $                 t1_updateloc - t0_updateloc
                  enddo
               enddo

               
            endif
 100     continue



         call ypwta_coef(nterms,nnodes,comp,comp2)
c     coefficients needed to convert plane-wave expansions to local expansions

         do 150 jj=istart,iend
            j=iboxlev(jj)
            ichildsum=ichildbox(1,j)+ichildbox(2,j)+
     $           ichildbox(3,j)+ichildbox(4,j)
            if (ichildsum .ne. -4) then
               do ii=1,4
                  ic1=ichildbox(ii,j)
                  if (ic1 .ne. -1) then
                    t0_P2L = ETIME(timeep)
                     call yexp4local(locexp(0,ic1),nterms,nnodes,
     $                    expw(1,ic1),exps(1,ic1),expe(1,ic1),
     $                    expn(1,ic1),comp2)
                    t1_P2L = ETIME(timeep)
                    t_P2L = t_P2L + t1_P2L - t0_P2L
                  endif
c     convert the plane wave expansions of the (upto) four children
c     to a local expanions
               enddo
            endif
 150     continue
      enddo


      do i=1,npts
         ibox=fb(i)
         t0_evalloc = ETIME(timeep) 
         call localexpansion(layers,locexp(0,ibox),icolbox(ibox),
     $        irowbox(ibox),levelbox(ibox),dom(i),pot(i),
     $        potx(i),poty(i),potxx(i),potxy(i),potyy(i),
     $        betascal(levelbox(ibox)),beta,nterms)
         t1_evalloc = ETIME(timeep) 
         t_evalloc = t_evalloc + t1_evalloc - t0_evalloc
      enddo
c     Evaluate the local exapnsion due to the finest box

c      write(6,1020) 'Time doing local to local translations is ',
c     $      t_L2L,' seconds'
c      write(6,1020) 'Time doing multipole to wave translations is ',
c     $      t_M2P,' seconds'
c      write(6,1020) 'Time doing wave to wave translations is ',
c     $      t_P2P,' seconds'
      write(6,1020) 'Time doing direct evaluations is ',
     $      t_direct,' seconds'
c      write(6,1020) 'Time evaluating multipole expansions is ',
c     $      t_evalmp,' seconds'
c      write(6,1020) 'Time updating local expansions is ',
c     $      t_updateloc,' seconds'
c      write(6,1020) 'Time doing wave to local translations is ',
c     $      t_P2L,' seconds'
c      write(6,1020) 'Time evaluating local translations is ',
c     $      t_evalloc,' seconds'
c     Report timings for different components of the FMM



 1020 format(A,ES8.2,A)

      return
      end





c**********************************************************************
      subroutine create_tree(levelbox,iparentbox,ichildbox,icolbox,
     $     irowbox,icolleagbox,nboxes,numptsbox,nlev,nblevel,fpt,
     $     iboxlev,istartlev,iflag,itemparray,maxlevel,maxpts,maxboxes,
     $     dom,npts,maxptsperbox,permute,fb)
      implicit none

c     variables that are passed
      integer *4 levelbox(maxboxes)
      integer *4 iparentbox(maxboxes)
      integer *4 ichildbox(4,maxboxes)
      integer *4 icolbox(maxboxes),irowbox(maxboxes)
      integer *4 icolleagbox(9,maxboxes)
      integer *4 nboxes,nlev
      integer *4 numptsbox(maxboxes)
      integer *4 nblevel(0:maxlevel)
      integer *4 fpt(maxboxes)
      integer *4 iboxlev(maxboxes)
      integer *4 istartlev(0:maxlevel)
      integer *4 iflag(maxboxes)
      integer *4 itemparray(maxboxes)
      integer *4 maxlevel,maxpts,maxboxes
      integer *4 npts,maxptsperbox
      integer *4 permute(maxpts)
      integer fb(maxpts)
      complex(8) dom(maxpts)
      

c     variables that are local
      complex *16 cen1,cen2,cen3,cen4
      complex *16 eye
      integer permute_old(2**20)
      integer *4 istart,iend,istart2,iend2
      integer *4 ibox,levflag
      integer *4 nbox1,nbox2,nbox3,nbox4
      integer *4 i,j,k,n,nn
      integer *4 ic1,ic2,ic3,ic4
      real*8 eps,d1,d2,d3,d4
c     For doing some timings for different components of
c     computing the tree
      real *4 t0,t1,ETIME,timeep(2)
      real *4 t_infnorm,t_subdivide

      t_infnorm = 0.d0
      t_subdivide = 0.d0

      eps=1.0d-15
      eye=(0.d0,1.d0)

      do i=1,maxlevel
         nblevel(i)=0
         istartlev(i)=0
      enddo

      ibox=1
      levelbox(ibox)=0
      icolbox(ibox)=1
      irowbox(ibox)=1
      iparentbox(ibox)=-1
      ichildbox(1,ibox)=-1
      ichildbox(2,ibox)=-1
      ichildbox(3,ibox)=-1
      ichildbox(4,ibox)=-1
c     initializations for box one which is the entire
c     computational domain
      nboxes=1

      nblevel(0)=1
      istartlev(0)=1
      iboxlev(1)=1

      t0 = ETIME(timeep)
      call subdivide1(1,iparentbox,ichildbox,nboxes,irowbox,
     $     icolbox,levelbox,maxlevel,istartlev,nblevel,iboxlev,
     $     maxboxes,maxlevel)
      t1 = ETIME(timeep)
      t_subdivide = t_subdivide + t1 - t0

      nlev=1

      cen1=(-0.25d0,-0.25d0)
      cen2=(0.25d0,-0.25d0)
      cen3=(-0.25d0,0.25d0)
      cen4=(0.25d0,0.25d0)

      numptsbox(1)=npts
c     box 1 is the entire computational domain

      do i=1,npts
        t0 = ETIME(timeep)
         call inf_norm(dom(i)-cen1,d1)
         call inf_norm(dom(i)-cen2,d2)
         call inf_norm(dom(i)-cen3,d3)
         call inf_norm(dom(i)-cen4,d4)
        t1 = ETIME(timeep)
        t_infnorm = t_infnorm + t1 - t0
         if (d1 .le. 0.25d0) then
            fb(i)=2
         elseif (d2 .le. 0.25d0) then
            fb(i)=3
         elseif (d3 .le. 0.25d0) then
            fb(i)=4
         elseif (d4 .le. 0.25d0) then
            fb(i)=5
         else
            print*,'BIG PROBLEM ... SOME POINT WAS NOT BINNED'
         endif
      enddo
c     All points need to be contained in [-1/2 1/2]^{2}.  This
c     can always be done by scaling in space and adjusting beta

      fpt(1)=1
      numptsbox(2)=0
      numptsbox(3)=0
      numptsbox(4)=0
      numptsbox(5)=0
      n=1
      do k=2,5
         do i=1,npts
            if (fb(i) .eq. k) then
               numptsbox(k)=numptsbox(k)+1
               if (numptsbox(k) .eq. 1) then
                  fpt(k)=n
               endif
               permute(n)=i
               permute_old(n)=i
               n=n+1
            endif
         enddo
      enddo
c     Compute everything at level 1 (four children box of the 
c     entire computational domain)

      do i=1,maxlevel-1
         levflag=0
         istart=istartlev(i)
         iend=istart+nblevel(i)-1
         do j=istart,iend
            ibox=iboxlev(j)
            if (numptsbox(ibox) .gt. maxptsperbox) then
               if (levflag .eq. 0) then
                  nlev=nlev+1
                  levflag=1
               endif

               t0 = ETIME(timeep)
               call subdivide1(ibox,iparentbox,
     $              ichildbox,nboxes,irowbox,icolbox,levelbox,maxlevel,
     $              istartlev,nblevel,iboxlev,maxboxes,maxlevel)
               t1 = ETIME(timeep)
               t_subdivide = t_subdivide + t1 - t0

c     tree is subdivided...now need to update permute
               nbox1=ichildbox(1,ibox)
               cen1=dble(2*icolbox(nbox1)-1)+
     $              eye*dble(2*irowbox(nbox1)-1)
               cen1=cen1*(0.5d0)**(levelbox(nbox1)+1)-(0.5d0,0.5d0)

               nbox2=ichildbox(2,ibox)
               cen2=dble(2*icolbox(nbox2)-1)+
     $              eye*dble(2*irowbox(nbox2)-1)
               cen2=cen2*(0.5d0)**(levelbox(nbox2)+1)-(0.5d0,0.5d0)
            
               nbox3=ichildbox(3,ibox)
               cen3=dble(2*icolbox(nbox3)-1)+
     $              eye*dble(2*irowbox(nbox3)-1)
               cen3=cen3*(0.5d0)**(levelbox(nbox3)+1)-(0.5d0,0.5d0)

               nbox4=ichildbox(4,ibox)
               cen4=dble(2*icolbox(nbox4)-1)+
     $              eye*dble(2*irowbox(nbox4)-1)
               cen4=cen4*(0.5d0)**(levelbox(nbox4)+1)-(0.5d0,0.5d0)

               numptsbox(nbox1)=0
               numptsbox(nbox2)=0
               numptsbox(nbox3)=0
               numptsbox(nbox4)=0

               do k=fpt(ibox),fpt(ibox)+numptsbox(ibox)-1
                  t0 = ETIME(timeep)
                  call inf_norm(dom(permute(k))-cen1,d1)
                  t1 = ETIME(timeep)
                  t_infnorm = t_infnorm + t1 - t0
                  if (d1 .le. 0.5d0**(levelbox(nbox1)+1)+eps) then
                     fb(permute(k))=nbox1
                  endif
                  t0 = ETIME(timeep)
                  call inf_norm(dom(permute(k))-cen2,d2)
                  t1 = ETIME(timeep)
                  t_infnorm = t_infnorm + t1 - t0
                  if (d2 .le. 0.5d0**(levelbox(nbox2)+1)+eps) then
                     fb(permute(k))=nbox2
                  endif
                  t0 = ETIME(timeep)
                  call inf_norm(dom(permute(k))-cen3,d3)
                  t1 = ETIME(timeep)
                  t_infnorm = t_infnorm + t1 - t0
                  if (d3 .le. 0.5d0**(levelbox(nbox3)+1)+eps) then
                     fb(permute(k))=nbox3
                  endif
                  t0 = ETIME(timeep)
                  call inf_norm(dom(permute(k))-cen4,d4)
                  t1 = ETIME(timeep)
                  t_infnorm = t_infnorm + t1 - t0
                  if (d4 .le. 0.5d0**(levelbox(nbox4)+1)+eps) then
                     fb(permute(k))=nbox4
                  endif
               enddo

            endif
         enddo


         do 150 j=istartlev(i),istartlev(i)+nblevel(i)-1
            nn=0
            ic1=ichildbox(4,j)
            ic2=ichildbox(3,j)
            ic3=ichildbox(1,j)
            ic4=ichildbox(2,j)
            if (ic1+ic2+ic3+ic4 .eq. -4) goto 150
            istart2=fpt(j)
            iend2=istart2+numptsbox(j)-1

            n=0
            do k=istart2,iend2
               if (fb(permute_old(k)) .eq. ic1) then
                  nn=nn+1
                  n=n+1
                  numptsbox(ic1)=numptsbox(ic1)+1
                  if (numptsbox(ic1) .eq. 1) then
                     fpt(ic1)=fpt(j)+nn-1
                  endif
                  permute(n+fpt(ic1)-1)=permute_old(k)
               endif
            enddo

            n=0
            do k=istart2,iend2
               if (fb(permute_old(k)) .eq. ic2) then
                  nn=nn+1
                  n=n+1
                  numptsbox(ic2)=numptsbox(ic2)+1
                  if (numptsbox(ic2) .eq. 1) then
                     fpt(ic2)=fpt(j)+nn-1
                  endif
                  permute(n+fpt(ic2)-1)=permute_old(k)
               endif
            enddo

            n=0
            do k=istart2,iend2
               if (fb(permute_old(k)) .eq. ic3) then
                  nn=nn+1
                  n=n+1
                  numptsbox(ic3)=numptsbox(ic3)+1
                  if (numptsbox(ic3) .eq. 1) then
                     fpt(ic3)=fpt(j)+nn-1
                  endif
                  permute(n+fpt(ic3)-1)=permute_old(k)
               endif
            enddo
      
            n=0
            do k=istart2,iend2
               if (fb(permute_old(k)) .eq. ic4) then
                  nn=nn+1
                  n=n+1
                  numptsbox(ic4)=numptsbox(ic4)+1
                  if (numptsbox(ic4) .eq. 1) then
                     fpt(ic4)=fpt(j)+nn-1
                  endif
                  permute(n+fpt(ic4)-1)=permute_old(k)
               endif
            enddo
c         Update permute which is used to decide which points are
c         in which childless boxes


 150     continue

         do j=1,npts
            permute_old(j)=permute(j)
         enddo
      enddo

      return
      end


c**********************************************************************
      subroutine subdivide(iparbox,iparentbox,ichildbox,icolleagbox,
     1         nboxes,irowbox,icolbox,levelbox,nlev,
     2         istartlev, nblevel, iboxlev,itemparray)
c**********************************************************************
c     the following subroutine is designed to divide up a childless
c     box into four children.
c     the children are placed in correct order (clockwise starting
c     from the upper left corner) so there is no need to 'shuffle'
c     the child order later on. in the periodic case, the colleagues
c     must be obtained by looking at the potential colleague numbers
c     and their row and column and seeing if they lie outside of
c     the domain. if they do it must be readjusted to account for the
c     periodicity.
c
c   input:
c     iparbox : the box being divided
c     iparentbox : the parent of each box
c     ichildbox : the four children of each box
c     icolleagbox : the colleagues of a given box
c     nboxes : the total number of boxes
c     irowbox : the row of each box
c     icolbox : the column of each box
c     levelbox : an array determining the level of each box
c     nlev : the finest level
c     istartlev : the pointer to where each level begins in the
c               iboxlev array
c     nblevel : the total number of boxes per level
c     iboxlev : the array in which the boxes are arranged
c     itemparray : just a dummy array
c
c   output:
c     nboxes and ichildbox : are altered to
c            reflect the addition of new boxes
c
c**********************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4 nlev, nboxes
      integer *4 iparentbox(nboxes+4), ichildbox(4,nboxes+4)
      integer *4 icolbox(nboxes+4), irowbox(nboxes+4)
      integer *4 levelbox(nboxes+4)
      integer *4 iparbox
      integer *4 nblevel(0:nlev),iboxlev(nboxes+4),istartlev(0:nlev)
      integer *4 icolleagbox(9,nboxes+4)
c
c-----local variables
c
      integer *4  level, ibox, itemparray(1)
      integer *4  icolumn, irow, i, j,l
      integer *4  icntr, jcntr, isister
      integer *4  icoltemp, icoltest
      integer *4  irowtemp, irowtest
      integer *4  partemp, colleague, nside, ilev, itest
c
c-----let's initialize the array itemparray to zero:
c
      do i = 1, nboxes + 4
        itemparray(i) = 0
      end do
c
c-----level, icolumn, and irow refer to the level, column,
c     and row of the parent box, respectively.
c
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)
c
c-----here are the new boxes placed in the
c     correct positions.  they are all childless.
c     there columns and rows are determined from
c     the parents columns and rows.  the level is
c     obviously one level finer than the parent.
c
      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2
c
      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2
c
      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1
c
c-----set up a temporary array to store the old one in:
c
      do i = 1, nboxes
        itemparray(i) = iboxlev(i)
      end do
c
c-----now let's rearrange the ladder structure:
c
      nblevel(level + 1) = nblevel(level + 1) + 4
c
      iboxlev(istartlev(level+2))   = nboxes + 1
      iboxlev(istartlev(level+2)+1) = nboxes + 2
      iboxlev(istartlev(level+2)+2) = nboxes + 3
      iboxlev(istartlev(level+2)+3) = nboxes + 4
c
      do i = istartlev(level+2) + 4, nboxes + 4
        iboxlev(i) = itemparray(i - 4)
      end do
c
      do i = level + 2, nlev
        istartlev(i) = istartlev(i) + 4
      end do
      nboxes = nboxes + 4
c
c------now let's go through the process of reforming any
c      neccessary colleagues.  for each of the child boxes
c      that we just formed, all we need to do is scan through
c      the boxes that are children of the above parent boxes
c      colleagues and test the column and row numbers.  we can
c      also take advantage of the fact that for every one of
c      the newly formed boxes colleagues, that box will list
c      the newly formed box as one of its colleagues.
c      the colleague numbers can be found easily if we think
c      of a 'reflection.'  colleague 1 and 9 are opposites,
c      3 and 7 etc.
c      first do the free space case:
c
       do 200 i = 1, 4
         if(ichildbox(1,iparbox) .lt. 0)goto 200
         ibox = ichildbox(i,iparbox)
         icolleagbox(5,ibox) = ibox
         do j = 1, 4
           icolleagbox(j,ibox) = -1
         end do
         do j = 6, 9
           icolleagbox(j,ibox) = -1
         end do
c
         partemp = iparentbox(ibox)
c
c--------irowtemp and icoltemp denote the
c        row and column of the test box.
c
         irowtemp = irowbox(ibox)
         icoltemp = icolbox(ibox)
c
         do 100 jcntr = 1, 9
c
c----------colleague denotes the colleague of the parent box.
c
           colleague = icolleagbox(jcntr,partemp)
c
c----------if the colleague doesn't exist
c          or is childless, skip it:
c
           if (colleague .lt. 0)goto 100
           if (ichildbox(1,colleague) .lt. 0)goto 100
c
c----------otherwise scan the four children:
c
           do icntr = 1, 4
             j = ichildbox(icntr,colleague)
c
c------------irowtest and icoltest denote the row and column of
c            the box being compared to the test box.
c
             irowtest = irowbox(j)
             icoltest = icolbox(j)
c
             if(irowtemp .eq. irowtest+1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(1,ibox) = j
                 icolleagbox(9,j) = ibox
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(2,ibox) = j
                 icolleagbox(8,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(3,ibox) = j
                 icolleagbox(7,j) = ibox
               endif
             elseif(irowtemp .eq. irowtest)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(4,ibox) = j
                 icolleagbox(6,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(6,ibox) = j
                 icolleagbox(4,j) = ibox
               endif
             elseif(irowtemp .eq. irowtest-1)then
               if(icoltemp .eq. icoltest+1)then
                 icolleagbox(7,ibox) = j
                 icolleagbox(3,j) = ibox
               elseif(icoltemp .eq. icoltest)then
                 icolleagbox(8,ibox) = j
                 icolleagbox(2,j) = ibox
               elseif(icoltemp .eq. icoltest-1)then
                 icolleagbox(9,ibox) = j
                 icolleagbox(1,j) = ibox
               endif
             endif
           end do
100      continue
200   continue



      return
      end



c**********************************************************************
      subroutine subdivide1(iparbox,iparentbox,
     $     ichildbox, nboxes,irowbox,icolbox,levelbox,nlev,
     $     istartlev, nblevel, iboxlev, 
     $     maxboxes,maxlevel)
c**********************************************************************
c     this routine is identical to the subdivide routine except that
c     it does not concern itself at all with generating colleagues.
c     this routine is used only within the mktree4 routine.
c
c   input:
c     iparbox : the box being divided
c     iparentbox : the parent of each box
c     ichildbox : the four children of each box
c     nboxes : the total number of boxes
c     irowbox : the row of each box
c     icolbox : the column of each box
c     levelbox : an array determining the level of each box
c     nlev : the finest level
c     istartlev : the pointer to where each level begins in the
c               iboxlev array
c     nblevel : the total number of boxes per level
c     iboxlev : the array in which the boxes are arranged
c
c   output:
c     nboxes and ichildbox : altered to reflect the addition
c            of new boxes
c
c   called from : mktree4().
c
c   subroutine called : none.
c
c********************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4 iparbox,nboxes,nlev,maxboxes,maxlevel
      integer *4 iparentbox(maxboxes), ichildbox(4,maxboxes)
      integer *4 icolbox(maxboxes), irowbox(maxboxes)
      integer *4 levelbox(maxboxes)
      integer *4 nblevel(0:maxlevel), iboxlev(maxboxes)
      integer *4 istartlev(0:maxlevel)
c
c-----local variables
c
      integer *4 i,j
      integer *4 ncntr,level
      integer *4 irow,icolumn,ibox

c
c-----level, icolumn, and irow refer to the level, column,
c     and row of the parent box, respectively.
c
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)
c
c-----here are the new boxes placed in the
c     correct positions.  they are all childless.
c     there columns and rows are determined from
c     the parents columns and rows.  the level is
c     obviously one level finer than the parent.
c


      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2
c
      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2
c
      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1
c
      nboxes = nboxes + 4
c     store relationships between the current box and its children


c      ncntr = 1
cc      do i = 0, nlev
cc      No need going through all levels of the tree at each subdivision
c      do i = 0,level+1
c        nblevel(i) = 0
c        istartlev(i) = ncntr
c        do j = 1, nboxes
c          if(levelbox(j) .eq. i)then
c            iboxlev(ncntr) = j
c            ncntr = ncntr + 1
c            nblevel(i) = nblevel(i) + 1
c          endif
c        end do
c      end do

      if (nblevel(level+1) == 0) then
        istartlev(level+1) = nboxes - 3
      endif
      nblevel(level+1) = nblevel(level+1) + 4
      do i=nboxes-3,nboxes
        iboxlev(i) = i
      enddo
c     Sept 17, 2012 - This does the exact same thing as the
c     above piece of code, but orders of magnitued faster
c     Don't forsee it causing any problems

      return
      end



c**********************************************************************
      subroutine compress_tree(levelbox,iparentbox,ichildbox,icolbox,
     $     irowbox,icolleagbox,iboxlev,nboxes,numptsbox,
     $     nlev,nblevel,fpt,istartlev,fb,npts,permute)
c     Get rid of empty boxes in the tree so that we don't compute
c     on them at all.  All the variables are identical to those
c     in create_tree

      implicit none

      integer *4 levelbox(nboxes)
      integer *4 iparentbox(nboxes)
      integer *4 ichildbox(4,nboxes)
      integer *4 icolbox(nboxes),irowbox(nboxes)
      integer *4 icolleagbox(9,nboxes)
      integer *4 iboxlev(nboxes)
      integer *4 nboxes
      integer *4 numptsbox(nboxes)
      integer *4 nlev
      integer *4 nblevel(0:nlev)
      integer *4 fpt(nboxes)
      integer *4 fb(npts)
      integer *4 istartlev(0:nlev)
      integer *4 npts
      integer *4 permute(npts)


      integer *4 i,j,k
      integer *4 ncount
      integer *4 irow,icol,irow2,icol2
      integer *4 istart,iend,istart2,iend2
      

      do i=0,nlev
         nblevel(i)=0
      enddo


      ncount=1
      do i=1,nboxes
         if (numptsbox(i) .ne. 0) then
            icolbox(ncount)=icolbox(i)
            irowbox(ncount)=irowbox(i)
            levelbox(ncount)=levelbox(i)
            iboxlev(ncount)=ncount
            numptsbox(ncount)=numptsbox(i)
            fpt(ncount)=fpt(i)
            do j=fpt(ncount),fpt(ncount)+numptsbox(ncount)-1
               fb(permute(j))=ncount
            enddo
            nblevel(levelbox(ncount))=nblevel(levelbox(ncount))+1
            if (nblevel(levelbox(ncount)) .eq. 1) then
               istartlev(levelbox(ncount))=ncount
            endif
            ncount=ncount+1
         endif
      enddo

      nboxes=ncount-1

      iparentbox(1)=-1
      do i=1,nboxes
         do k=1,4
            ichildbox(k,i)=-1
         enddo
      enddo
c     set the parent of the coarsest box and all the children to -1

      do i=1,nlev
         istart=istartlev(i)
         iend=istart+nblevel(i)-1
         istart2=istartlev(i-1)
         iend2=istart2+nblevel(i-1)-1
         do j=istart,iend
            icol=icolbox(j)
            irow=irowbox(j)
            do k=istart2,iend2
               icol2=2*icolbox(k)
               irow2=2*irowbox(k)
               if (irow .eq. irow2 .and. icol .eq. icol2) then
                  ichildbox(2,k)=j
                  iparentbox(j)=k
               endif
               if (irow .eq. irow2 .and. icol .eq. icol2-1) then
                  ichildbox(1,k)=j
                  iparentbox(j)=k
               endif
               if (irow .eq. irow2-1 .and. icol .eq. icol2) then
                  ichildbox(3,k)=j
                  iparentbox(j)=k
               endif
               if (irow .eq. irow2-1 .and. icol .eq. icol2-1) then
                  ichildbox(4,k)=j
                  iparentbox(j)=k
               endif
            enddo
         enddo
      enddo

      call mkcolls(icolbox,irowbox,icolleagbox,nboxes,nlev,
     $     iparentbox,ichildbox,nblevel,iboxlev,istartlev)
c     Need colleagues next for forming list1 and list3.

      return
      end




c**********************************************************************
      subroutine mkcolls(icolbox,
     1      irowbox,icolleagbox,nboxes,nlev,
     2      iparentbox,ichildbox,nblevel,
     3      iboxlev, istartlev)
c**********************************************************************
c     the following subroutine is used to generate the colleagues
c     for all of the boxes in the tree structure.  if a colleague
c     doesn't exist it is set to -1.  each box has nine colleagues
c     and they are ordered as follows:
c
c                        7     8     9
c
c                        4     5     6
c
c                        1     2     3
c
c     you are your own colleague number 5.
c     the algorithm used here is recursive and takes advantage of
c     the fact that your colleagues can only be the children of
c     your parents colleagues.  there is no need to scan all of the
c     boxes.  
c
c   input:
c     levelbox is an array determining the level of each box
c     icolbox denotes the column of each box
c     irowbox denotes the row of each box
c     nboxes is the total number of boxes
c     nlev is the finest level
c     iparentbox denotes the parent of each box
c     ichildbox denotes the four children of each box
c     nblevel is the total number of boxes per level
c     iboxlev is the array in which the boxes are arranged
c     istartlev is the pointer to where each level begins in the
c               iboxlev array
c
c   output:
c     icolleagbox denotes the colleagues of a given box
c
c**********************************************************************
c
      implicit none
c
c-----global variables
c
      integer *4 iparentbox(nboxes)
      integer *4 icolleagbox(9,nboxes)
      integer *4 icolbox(nboxes), irowbox(nboxes)
      integer *4 ichildbox(4,nboxes)
      integer *4 nblevel(0:nlev),iboxlev(nboxes),istartlev(0:nlev)
      integer *4 nboxes,nlev
c
c-----local variables
c
      integer *4 colleague, partemp
      integer *4 jcntr, ibox, itest
      integer *4 icntr, ilev, j, l, nside
      integer *4 irowtemp, icoltemp
      integer *4 irowtest, icoltest
      integer *4 ichildsum
c
c-----initialize colleague number 5 to
c     yourself and all other colleagues to
c     -1.  -1 is the flag for the case when
c     the colleagues don't exist.  it can
c     be overwritten below.
c
      do ibox = 1, nboxes
        icolleagbox(5,ibox) = ibox
        do j = 1, 4
          icolleagbox(j,ibox) = -1
        end do
        do j = 6, 9
          icolleagbox(j,ibox) = -1
        end do
      end do
c
c-----scan through all of the levels except the coarsest level.
c     the one box at this level cannot have any colleagues.
c     do the uniform case first:
c
      do ilev = 1, nlev
c
c---------scan through all of the boxes on each level.  for each test
c         box, scan the parent's colleagues and test to see if
c         their children are in contact with the box being tested.
c         each colleague is placed in the correct order, so there is
c         no need to 'shuffle' them later on.
c
         do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
            ibox    = iboxlev(l)
            partemp = iparentbox(ibox)
c
c-----------irowtemp and icoltemp denote the row and column of
c           the test box.
c
            irowtemp = irowbox(ibox)
            icoltemp = icolbox(ibox)
c
            do 100 jcntr = 1, 9
c
c-------------colleague denotes the colleague of the parent box.
c
               colleague = icolleagbox(jcntr,partemp)
c
c-------------if the colleague doesn't exist
c             or is childless, skip it:
c
               if (colleague .lt. 0)goto 100
               ichildsum=ichildbox(1,colleague)+ichildbox(2,colleague)+
     $              ichildbox(3,colleague)+ichildbox(4,colleague)
               if (ichildsum .eq. -4) goto 100
               do icntr = 1, 4
                  j = ichildbox(icntr,colleague)
c
c---------------irowtest and icoltest denote the row and column of
c               the box being compared to the test box.
c
                  irowtest = irowbox(j)
                  icoltest = icolbox(j)
c
                  if(irowtemp .eq. irowtest+1)then
                     if(icoltemp .eq. icoltest+1)then
                        icolleagbox(1,ibox) = j
                     elseif(icoltemp .eq. icoltest)then
                        icolleagbox(2,ibox) = j
                     elseif(icoltemp .eq. icoltest-1)then
                        icolleagbox(3,ibox) = j
                     endif
                  elseif(irowtemp .eq. irowtest)then
                     if(icoltemp .eq. icoltest+1)then
                        icolleagbox(4,ibox) = j
                     elseif(icoltemp .eq. icoltest-1)then
                        icolleagbox(6,ibox) = j
                     endif
                  elseif(irowtemp .eq. irowtest-1)then
                     if(icoltemp .eq. icoltest+1)then
                        icolleagbox(7,ibox) = j
                     elseif(icoltemp .eq. icoltest)then
                        icolleagbox(8,ibox) = j
                     elseif(icoltemp .eq. icoltest-1)then
                        icolleagbox(9,ibox) = j
                     endif
                  endif
               end do
 100        continue
         end do
      end do

c
      return
      end


c**********************************************************************
      subroutine mklists13(iout,list1,nlist1,list3,nlist3,
     $     ichildbox,icolleagbox)
c     Form list1 and list3 of the tree.  They are defined as follows
c     List1 - only defined for childless boxes.  The set of all
c       childless adjacent boxes
c     List3 - only defined for childless boxes.  The descendents
c       neighbouring boxes that are not adjacent to the box of 
c       interest, but whose parent is adjacent to the box of
c       interest.
      implicit none

      integer *4 iout
      integer *4 list1(50,1),list3(50,1)
      integer *4 nlist1(1),nlist3(1)
      integer *4 ichildbox(4,1)
      integer *4 icolleagbox(9,1)

c     Local variables
      integer *4 ichildren1(50),ichildren2(50)
      integer *4 ichildren3(10),ichildren4(10)
      integer *4 i,ii
      integer *4 ncount1_lev,ncount2_lev
      integer *4 ncount3_lev,ncount4_lev
      integer *4 nchildren1,nchildren2
      integer *4 nchildren3,nchildren4
      integer *4 nchildren1_temp
      integer *4 nchildren2_temp
      integer *4 nchildren3_temp
      integer *4 nchildren4_temp
      integer *4 nb,iflag
      integer *4 ic1,ic2,ic3,ic4

      nlist1(iout)=1
      nlist3(iout)=1


      if (icolleagbox(1,iout) .ne. -1) then
        nb = icolleagbox(1,iout)
        iflag=0
        nchildren2=1
        ncount2_lev=0

        ic2=ichildbox(2,nb)

        if (ic2 .ne. -1) then
          ichildren2(nchildren2+ncount2_lev)=ic2
          ncount2_lev=ncount2_lev+1
          iflag=1
        endif
c     put the top right box in the ichildren2 list if this box exists
        nchildren2_temp=nchildren2
        nchildren2=nchildren2+ncount2_lev
          
        ic1=ichildbox(1,nb)
        ic3=ichildbox(3,nb)
        ic4=ichildbox(4,nb)
          
        if (ic1 .ne. -1) then
          list3(nlist3(iout),iout)=ic1
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic3 .ne. -1) then
          list3(nlist3(iout),iout)=ic3
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic4 .ne. -1) then
          list3(nlist3(iout),iout)=ic4
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the first, third or fourth children exist, they are in list 3

        do while (iflag .eq. 1)
          iflag=0
          ncount2_lev=0
c     loop over new boxes that have been put in ichildren2
          do i=nchildren2_temp,nchildren2-1
            ic1=ichildbox(1,ichildren2(i))
            ic2=ichildbox(2,ichildren2(i))
            ic3=ichildbox(3,ichildren2(i))
            ic4=ichildbox(4,ichildren2(i))
            if (ic2 .ne. -1) then
              ichildren2(nchildren2+ncount2_lev)=ic2
              ncount2_lev=ncount2_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic3 .ne. -1) then
              list3(nlist3(iout),iout)=ic3
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren2(i)
              nlist1(iout)=nlist1(iout)+1
            endif
          enddo
          nchildren2_temp=nchildren2
          nchildren2=nchildren2+ncount2_lev
        enddo
      endif



c     If it is the second colleague, list 3 comes from third 
c     and fourth children and list 1 comes from first and 
c     second children
      if (icolleagbox(2,iout) .ne. -1) then
        nb = icolleagbox(2,iout)
        iflag=0
        nchildren1=1
        nchildren2=1
        ncount1_lev=0
        ncount2_lev=0
          
        ic1=ichildbox(1,nb)
        ic2=ichildbox(2,nb)
          
        if (ic1 .ne. -1) then
          ichildren1(nchildren1+ncount1_lev)=ic1
          ncount1_lev=ncount1_lev+1
          iflag=1
        endif
        if (ic2 .ne. -1) then
          ichildren2(nchildren2+ncount2_lev)=ic2
          ncount2_lev=ncount2_lev+1
          iflag=1
        endif
c     put the top two boxes in the lists ichildren1 and ichildren2
c     if these boxes exist
        nchildren1_temp=nchildren1
        nchildren2_temp=nchildren2
        nchildren1=nchildren1+ncount1_lev
        nchildren2=nchildren2+ncount2_lev
          
        ic3=ichildbox(3,nb)
        ic4=ichildbox(4,nb)
        if (ic3 .ne. -1) then
          list3(nlist3(iout),iout)=ic3
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic4 .ne. -1) then
          list3(nlist3(iout),iout)=ic4
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the third and fourth children exist, they are in list 3
       
        do while (iflag .eq. 1)
          iflag=0
          ncount1_lev=0
          ncount2_lev=0
c     loop over new boxes that have been put in ichildren1
          do i=nchildren1_temp,nchildren1-1
            ic1=ichildbox(1,ichildren1(i))
            ic2=ichildbox(2,ichildren1(i))
            ic3=ichildbox(3,ichildren1(i))
            ic4=ichildbox(4,ichildren1(i))
            if (ic1 .ne. -1) then
              ichildren1(nchildren1+ncount1_lev)=ic1
              ncount1_lev=ncount1_lev+1
              iflag=1
            endif
            if (ic2 .ne. -1) then
              ichildren2(nchildren2+ncount2_lev)=ic2
              ncount2_lev=ncount2_lev+1
              iflag=1
            endif
            if (ic3 .ne. -1) then
              list3(nlist3(iout),iout)=ic3
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren1(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list1
          enddo
c     loop over new boxes that have been put in ichildren2
          do i=nchildren2_temp,nchildren2-1
            ic1=ichildbox(1,ichildren2(i))
            ic2=ichildbox(2,ichildren2(i))
            ic3=ichildbox(3,ichildren2(i))
            ic4=ichildbox(4,ichildren2(i))
            if (ic1 .ne. -1) then
              ichildren1(nchildren1+ncount1_lev)=ic1
              ncount1_lev=ncount1_lev+1
              iflag=1
            endif
            if (ic2 .ne. -1) then
              ichildren2(nchildren2+ncount2_lev)=ic2
              ncount2_lev=ncount2_lev+1
              iflag=1
            endif
            if (ic3 .ne. -1) then
             list3(nlist3(iout),iout)=ic3
             nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren2(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list 1
          enddo
          nchildren1_temp=nchildren1
          nchildren2_temp=nchildren2
          nchildren1=nchildren1+ncount1_lev
          nchildren2=nchildren2+ncount2_lev
c     update the counters
        enddo
      endif



c     If it is the third colleague, list 3 comes from second, third and 
c     fourth children and list 1 comes from first child
      if (icolleagbox(3,iout) .ne. -1) then
        nb = icolleagbox(3,iout)
        iflag=0
          
        nchildren1=1
        ncount1_lev=0

        ic1=ichildbox(1,nb)
          
        if (ic1 .ne. -1) then
          ichildren1(nchildren1+ncount1_lev)=ic1
          ncount1_lev=ncount1_lev+1
          iflag=1
        endif
c     put the top right box in the ichildren2 list if this box exists
        nchildren1_temp=nchildren1
        nchildren1=nchildren1+ncount1_lev
          
        ic2=ichildbox(2,nb)
        ic3=ichildbox(3,nb)
        ic4=ichildbox(4,nb)
          
        if (ic2 .ne. -1) then
          list3(nlist3(iout),iout)=ic2
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic3 .ne. -1) then
          list3(nlist3(iout),iout)=ic3
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic4 .ne. -1) then
          list3(nlist3(iout),iout)=ic4
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the first, third or fourth children exist, they are in list 3
          
        do while (iflag .eq. 1)
          iflag=0
          ncount1_lev=0
c     loop over new boxes that have been put in ichildren2
           do i=nchildren1_temp,nchildren1-1
             ic1=ichildbox(1,ichildren1(i))
             ic2=ichildbox(2,ichildren1(i))
             ic3=ichildbox(3,ichildren1(i))
             ic4=ichildbox(4,ichildren1(i))
             if (ic1 .ne. -1) then
               ichildren1(nchildren1+ncount1_lev)=ic1
               ncount1_lev=ncount1_lev+1
               iflag=1
             endif
             if (ic2 .ne. -1) then
               list3(nlist3(iout),iout)=ic2
               nlist3(iout)=nlist3(iout)+1
             endif
             if (ic3 .ne. -1) then
               list3(nlist3(iout),iout)=ic3
               nlist3(iout)=nlist3(iout)+1
             endif
             if (ic4 .ne. -1) then
               list3(nlist3(iout),iout)=ic4
               nlist3(iout)=nlist3(iout)+1
             endif
c     put their children in appropriate bins (if they exist)
             if (ic1+ic2+ic3+ic4 .eq. -4) then
               list1(nlist1(iout),iout)=ichildren1(i)
               nlist1(iout)=nlist1(iout)+1
             endif
           enddo
           nchildren1_temp=nchildren1
           nchildren1=nchildren1+ncount1_lev
         enddo
      endif
       

c     If it is the fourth colleague, list 3 comes from first and fourth children
c     and list 1 comes from second and third children
      if (icolleagbox(4,iout) .ne. -1) then
        nb = icolleagbox(4,iout)
        iflag=0
        nchildren2=1
        nchildren3=1
        ncount2_lev=0
        ncount3_lev=0
          
        ic2=ichildbox(2,nb)
        ic3=ichildbox(3,nb)
          
        if (ic2 .ne. -1) then
          ichildren2(nchildren2+ncount2_lev)=ic2
          ncount2_lev=ncount2_lev+1
          iflag=1
        endif
        if (ic3 .ne. -1) then
          ichildren3(nchildren3+ncount3_lev)=ic3
          ncount3_lev=ncount3_lev+1
          iflag=1
        endif
c     put the top two boxes in the lists ichildren1 and ichildren2
c     if these boxes exist
        nchildren2_temp=nchildren2
        nchildren3_temp=nchildren3
        nchildren2=nchildren2+ncount2_lev
        nchildren3=nchildren3+ncount3_lev
        
        ic1=ichildbox(1,nb)
        ic4=ichildbox(4,nb)
        if (ic1 .ne. -1) then
          list3(nlist3(iout),iout)=ic1
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic4 .ne. -1) then
          list3(nlist3(iout),iout)=ic4
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the third and fourth children exist, they are in list 3
       
        do while (iflag .eq. 1)
          iflag=0
          ncount2_lev=0
          ncount3_lev=0
c     loop over new boxes that have been put in ichildren1
          do i=nchildren2_temp,nchildren2-1
            ic1=ichildbox(1,ichildren2(i))
            ic2=ichildbox(2,ichildren2(i))
            ic3=ichildbox(3,ichildren2(i))
            ic4=ichildbox(4,ichildren2(i))
            if (ic2 .ne. -1) then
              ichildren2(nchildren2+ncount2_lev)=ic2
              ncount2_lev=ncount2_lev+1
              iflag=1
            endif
            if (ic3 .ne. -1) then
              ichildren3(nchildren3+ncount3_lev)=ic3
              ncount3_lev=ncount3_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren2(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list1
          enddo
c     loop over new boxes that have been put in ichildren2
          do i=nchildren3_temp,nchildren3-1
            ic1=ichildbox(1,ichildren3(i))
            ic2=ichildbox(2,ichildren3(i))
            ic3=ichildbox(3,ichildren3(i))
            ic4=ichildbox(4,ichildren3(i))
            if (ic2 .ne. -1) then
              ichildren2(nchildren2+ncount2_lev)=ic2
              ncount2_lev=ncount2_lev+1
              iflag=1
            endif
            if (ic3 .ne. -1) then
              ichildren3(nchildren3+ncount3_lev)=ic3
              ncount3_lev=ncount3_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren3(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list 1
          enddo
          nchildren2_temp=nchildren2
          nchildren3_temp=nchildren3
          nchildren2=nchildren2+ncount2_lev
          nchildren3=nchildren3+ncount3_lev
c     update the counters
        enddo
      endif



c     If it is the sixth colleague, list 3 comes from second and third children
c     and list 1 comes from first and fourth children
      if (icolleagbox(6,iout) .ne. -1) then
        nb = icolleagbox(6,iout)
        iflag=0
        nchildren1=1
        nchildren4=1
        ncount1_lev=0
        ncount4_lev=0
          
        ic1=ichildbox(1,nb)
        ic4=ichildbox(4,nb)

        if (ic1 .ne. -1) then
          ichildren1(nchildren1+ncount1_lev)=ic1
          ncount1_lev=ncount1_lev+1
          iflag=1
        endif
        if (ic4 .ne. -1) then
          ichildren4(nchildren4+ncount4_lev)=ic4
          ncount4_lev=ncount4_lev+1
          iflag=1
        endif
c     put the top two boxes in the lists ichildren1 and ichildren2
c     if these boxes exist
        nchildren1_temp=nchildren1
        nchildren4_temp=nchildren4
        nchildren1=nchildren1+ncount1_lev
        nchildren4=nchildren4+ncount4_lev
          
        ic2=ichildbox(2,nb)
        ic3=ichildbox(3,nb)
        if (ic2 .ne. -1) then
          list3(nlist3(iout),iout)=ic2
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic3 .ne. -1) then
          list3(nlist3(iout),iout)=ic3
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the third and fourth children exist, they are in list 3
       
        do while (iflag .eq. 1)
          iflag=0
          ncount1_lev=0
          ncount4_lev=0
c     loop over new boxes that have been put in ichildren1
          do i=nchildren1_temp,nchildren1-1
            ic1=ichildbox(1,ichildren1(i))
            ic2=ichildbox(2,ichildren1(i))
            ic3=ichildbox(3,ichildren1(i))
            ic4=ichildbox(4,ichildren1(i))
            if (ic1 .ne. -1) then
              ichildren1(nchildren1+ncount1_lev)=ic1
              ncount1_lev=ncount1_lev+1
              iflag=1
            endif
            if (ic4 .ne. -1) then
              ichildren4(nchildren4+ncount4_lev)=ic4
              ncount4_lev=ncount4_lev+1
              iflag=1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic3 .ne. -1) then
              list3(nlist3(iout),iout)=ic3
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren1(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list1
          enddo
c     loop over new boxes that have been put in ichildren2
          do i=nchildren4_temp,nchildren4-1
            ic1=ichildbox(1,ichildren4(i))
            ic2=ichildbox(2,ichildren4(i))
            ic3=ichildbox(3,ichildren4(i))
            ic4=ichildbox(4,ichildren4(i))
            if (ic1 .ne. -1) then
              ichildren1(nchildren1+ncount1_lev)=ic1
              ncount1_lev=ncount1_lev+1
              iflag=1
            endif
            if (ic4 .ne. -1) then
              ichildren4(nchildren4+ncount4_lev)=ic4
              ncount4_lev=ncount4_lev+1
              iflag=1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic3 .ne. -1) then
              list3(nlist3(iout),iout)=ic3
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren4(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list 1
          enddo
          nchildren1_temp=nchildren1
          nchildren4_temp=nchildren4
          nchildren1=nchildren1+ncount1_lev
          nchildren4=nchildren4+ncount4_lev
c     update the counters
        enddo
      endif


c     If it is the seventh colleague, list 3 comes from first, second, and 
c     fourth children and list 1 comes from third child
      if (icolleagbox(7,iout) .ne. -1) then
        nb = icolleagbox(7,iout)
        iflag=0
        nchildren3=1
        ncount3_lev=0
        
        ic3=ichildbox(3,nb)
          
        if (ic3 .ne. -1) then
          ichildren3(nchildren3+ncount3_lev)=ic3
          ncount3_lev=ncount3_lev+1
          iflag=1
        endif
c     put the top right box in the ichildren2 list if this box exists
        nchildren3_temp=nchildren3
        nchildren3=nchildren3+ncount3_lev
        
        ic1=ichildbox(1,nb)
        ic2=ichildbox(2,nb)
        ic4=ichildbox(4,nb)
        
        if (ic1 .ne. -1) then
          list3(nlist3(iout),iout)=ic1
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic2 .ne. -1) then
          list3(nlist3(iout),iout)=ic2
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic4 .ne. -1) then
          list3(nlist3(iout),iout)=ic4
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the first, third or fourth children exist, they are in list 3


        do while (iflag .eq. 1)
          iflag=0
          ncount3_lev=0
c     loop over new boxes that have been put in ichildren2
          do i=nchildren3_temp,nchildren3-1
            ic1=ichildbox(1,ichildren3(i))
            ic2=ichildbox(2,ichildren3(i))
            ic3=ichildbox(3,ichildren3(i))
            ic4=ichildbox(4,ichildren3(i))
            if (ic3 .ne. -1) then
              ichildren3(nchildren3+ncount3_lev)=ic3
              ncount3_lev=ncount3_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic4 .ne. -1) then
              list3(nlist3(iout),iout)=ic4
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren3(i)
              nlist1(iout)=nlist1(iout)+1
            endif
          enddo
          nchildren3_temp=nchildren3
          nchildren3=nchildren3+ncount3_lev
        enddo
      endif



c     If it is the eighth colleague, list 3 comes from first and second children
c     and list 1 comes from third and fourth children
      if (icolleagbox(8,iout) .ne. -1) then
        nb = icolleagbox(8,iout)
        iflag=0
        nchildren3=1
        nchildren4=1
        ncount3_lev=0
        ncount4_lev=0
          
        ic3=ichildbox(3,nb)
        ic4=ichildbox(4,nb)
          
        if (ic3 .ne. -1) then
          ichildren3(nchildren3+ncount3_lev)=ic3
          ncount3_lev=ncount3_lev+1
          iflag=1
        endif
        if (ic4 .ne. -1) then
          ichildren4(nchildren4+ncount4_lev)=ic4
          ncount4_lev=ncount4_lev+1
          iflag=1
        endif
c     put the top two boxes in the lists ichildren1 and ichildren2
c     if these boxes exist
        nchildren3_temp=nchildren3
        nchildren4_temp=nchildren4
        nchildren3=nchildren3+ncount3_lev
        nchildren4=nchildren4+ncount4_lev
          
        ic1=ichildbox(1,nb)
        ic2=ichildbox(2,nb)
        if (ic1 .ne. -1) then
          list3(nlist3(iout),iout)=ic1
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic2 .ne. -1) then
          list3(nlist3(iout),iout)=ic2
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the third and fourth children exist, they are in list 3
       
        do while (iflag .eq. 1)
          iflag=0
          ncount3_lev=0
          ncount4_lev=0
c     loop over new boxes that have been put in ichildren1
          do i=nchildren3_temp,nchildren3-1
            ic1=ichildbox(1,ichildren3(i))
            ic2=ichildbox(2,ichildren3(i))
            ic3=ichildbox(3,ichildren3(i))
            ic4=ichildbox(4,ichildren3(i))
            if (ic3 .ne. -1) then
              ichildren3(nchildren3+ncount3_lev)=ic3
              ncount3_lev=ncount3_lev+1
              iflag=1
            endif
            if (ic4 .ne. -1) then
              ichildren4(nchildren4+ncount4_lev)=ic4
              ncount4_lev=ncount4_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren3(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list1
          enddo
c     loop over new boxes that have been put in ichildren2
          do i=nchildren4_temp,nchildren4-1
            ic1=ichildbox(1,ichildren4(i))
            ic2=ichildbox(2,ichildren4(i))
            ic3=ichildbox(3,ichildren4(i))
            ic4=ichildbox(4,ichildren4(i))
            if (ic3 .ne. -1) then
              ichildren3(nchildren3+ncount3_lev)=ic3
              ncount3_lev=ncount3_lev+1
              iflag=1
            endif
            if (ic4 .ne. -1) then
              ichildren4(nchildren4+ncount4_lev)=ic4
              ncount4_lev=ncount4_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren4(i)
              nlist1(iout)=nlist1(iout)+1
            endif
c     if the box is childless, it is in list 1
          enddo
          nchildren3_temp=nchildren3
          nchildren4_temp=nchildren4
          nchildren3=nchildren3+ncount3_lev
          nchildren4=nchildren4+ncount4_lev
c     update the counters
        enddo
      endif


c     If it is the ninth colleague, list 3 comes from first, second, and 
c     third children and list 1 comes from fourth child
      if (icolleagbox(9,iout) .ne. -1) then
        nb = icolleagbox(9,iout)
        iflag=0
        nchildren4=1
        ncount4_lev=0

        ic4=ichildbox(4,nb)
        
        if (ic4 .ne. -1) then
          ichildren4(nchildren4+ncount4_lev)=ic4
          ncount4_lev=ncount4_lev+1
          iflag=1
        endif
c     put the top right box in the ichildren2 list if this box exists
        nchildren4_temp=nchildren4
        nchildren4=nchildren4+ncount4_lev
          
        ic1=ichildbox(1,nb)
        ic2=ichildbox(2,nb)
        ic3=ichildbox(3,nb)
          
        if (ic1 .ne. -1) then
          list3(nlist3(iout),iout)=ic1
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic2 .ne. -1) then
          list3(nlist3(iout),iout)=ic2
          nlist3(iout)=nlist3(iout)+1
        endif
        if (ic3 .ne. -1) then
          list3(nlist3(iout),iout)=ic3
          nlist3(iout)=nlist3(iout)+1
        endif
c     if the first, third or fourth children exist, they are in list 3


        do while (iflag .eq. 1)
          iflag=0
          ncount4_lev=0
c     loop over new boxes that have been put in ichildren2
          do i=nchildren4_temp,nchildren4-1
            ic1=ichildbox(1,ichildren4(i))
            ic2=ichildbox(2,ichildren4(i))
            ic3=ichildbox(3,ichildren4(i))
            ic4=ichildbox(4,ichildren4(i))
            if (ic4 .ne. -1) then
              ichildren4(nchildren4+ncount4_lev)=ic4
              ncount4_lev=ncount4_lev+1
              iflag=1
            endif
            if (ic1 .ne. -1) then
              list3(nlist3(iout),iout)=ic1
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic2 .ne. -1) then
              list3(nlist3(iout),iout)=ic2
              nlist3(iout)=nlist3(iout)+1
            endif
            if (ic3 .ne. -1) then
              list3(nlist3(iout),iout)=ic3
              nlist3(iout)=nlist3(iout)+1
            endif
c     put their children in appropriate bins (if they exist)
            if (ic1+ic2+ic3+ic4 .eq. -4) then
              list1(nlist1(iout),iout)=ichildren4(i)
              nlist1(iout)=nlist1(iout)+1
            endif
          enddo
          nchildren4_temp=nchildren4
          nchildren4=nchildren4+ncount4_lev
        enddo
      endif


      nlist1(iout)=nlist1(iout)-1
      nlist3(iout)=nlist3(iout)-1


      return
      end


c**********************************************************************
      subroutine ymultipole(nterms,mpole,x,charge,level,icol,irow,
     $     betascal)
c     Update the multipole coefficients due to a single point x
c     This is only called for childess boxes

c     Input:
c     nterms - number of multipole coefficients we are keeping
c     mpole - multipole coefficients thus far
c     x - location of source point in complex format
c     charge - charge at source point in double format
c     level - level of box that contains x
c     icol - column of box that contains x
c     irow - row of box that contains 
c     betascal - beta divided by powers of two

c     Output:
c     mpole - update multipole coefficients due to (x,charge)
      implicit none

      integer *4 nterms,level
      integer *4 icol,irow
      complex *16 mpole(0:nterms)
      complex *16 x
      real *8 charge
      real *8 betascal(0:1)
      real *8 bi(0:50)

      complex *16 cen
      complex *16 eye
      integer p
      real *8 r,theta
      real *8 beta

      data eye/(0.0d0,1.0d0)/
      beta=betascal(0)

      cen=dble(2*icol-1)+eye*dble(2*irow-1)
      cen=cen*(0.5d0)**(level+1)-(0.5d0,0.5d0)

      call polar_convert(x-cen,r,theta)
      call besseli(nterms+1,r*beta,bi,betascal(level))

      do p=0,nterms
        mpole(p) = mpole(p)+charge*bi(p)*exp(-eye*p*theta)
      enddo
c     Updated multipole coefficients


      return
      end
c**********************************************************************
      subroutine ychildpar(mpolepar,mpole1,mpole2,mpole3,mpole4,
     1                    nterms,c,ic1,ic2,ic3,ic4)
c**********************************************************************
c     this subroutine is designed to merge the multipole expansions
c     of four child boxes together to form the multipole expansion
c     of the parent.
c     on input,
c       mpole1, mpole2, mpole3, and mpole4 : represent the
c         four child multipole expansions.  the boxes are numbered
c         in a clockwise manner starting in the upper left corner.
c       nterms is the length of the multipole expansions.
c       c : is an array of precomputed coefficients c(i,j), and it is
c         designed for the first box.
c
c     on output,
c       mpolepar : is the multipole expansion of the parent box.
c
c     the algorithm works by passing through the loop mod 4 and
c     rescaling the new coefficients so that they are on the parent's
c     level. note that all the scaling are included in the precomputed
c     coefficients.
c
c     subroutine called : dconjg()
c
c     the ordering of the children boxes:
c
c            4  1
c            3  2
c
c     note : the computation should be further simplified to
c            be efficient.
c
c     subroutine called : none.
c
c     called from : adapfmm4()
c                   boundfmm4().
c**********************************************************************
c
      implicit none
c
c-----global variables
c
      integer   nterms
      complex *16 c(0:nterms, -nterms:nterms)
      complex *16 mpole1(0:nterms), mpole2(0:nterms)
      complex *16 mpole3(0:nterms), mpole4(0:nterms)
      complex *16 mpolepar(0:nterms)
c
c-----local variables
c
      integer ic1,ic2,ic3,ic4

      integer   m, k, ind,ind2
      complex *16 temp(0:50, 0:3)
      complex *16 temp1, temp2, temp3, temp4, imag
      data imag/(0.0d0,1.0d0)/
c
c-----first compute temp(), so we can use less multiplication.
c

      if (ic1 .eq. -1) then
         do m=0,nterms
            mpole1(m)=(0.0d0,0.0d0)
         enddo
      endif
      if (ic2 .eq. -1) then
         do m=0,nterms
            mpole2(m)=(0.0d0,0.0d0)
         enddo
      endif
      if (ic3 .eq. -1) then
         do m=0,nterms
            mpole3(m)=(0.0d0,0.0d0)
         enddo
      endif
      if (ic4 .eq. -1) then
         do m=0,nterms
            mpole4(m)=(0.0d0,0.0d0)
         enddo
      endif

      do m=0, nterms
        temp1 = mpole1(m) + mpole3(m)
        temp2 = mpole1(m) - mpole3(m)
        temp3 = mpole2(m) + mpole4(m)
        temp4 = imag*(mpole2(m) - mpole4(m))
        temp(m,0)=temp1+temp3
        temp(m,1)=temp2+temp4
        temp(m,2)=temp1-temp3
        temp(m,3)=temp2-temp4
      enddo
c
c-----first the zeroth term.
c
      do 1 m=0, nterms
        ind=mod(m,4)
        mpolepar(m)=c(m,0)*temp(0,ind)
1     continue
c
c-----now all the other terms of the children multipole coefficients.
c
      do 2 m=0, nterms, 4
        do k=1, nterms
          ind=mod(3*k,4)
          mpolepar(m)=mpolepar(m)+c(m,k)*temp(k,ind)
     1                  +c(m,-k)*dconjg(temp(k,ind))
        enddo
2     continue
c
      do 3 m=1, nterms,4
        do k=1, nterms
          ind=mod(1+3*k,4)
          ind2=mod(ind+2,4)
          mpolepar(m)=mpolepar(m)+c(m,k)*temp(k,ind)
     1                  +c(m,-k)*dconjg(temp(k,ind2))
        enddo
3     continue
c
      do 4 m=2, nterms,4
        do k=1, nterms
          ind=mod(2+3*k,4)
          mpolepar(m)=mpolepar(m)+c(m,k)*temp(k,ind)
     1                  +c(m,-k)*dconjg(temp(k,ind))
        enddo
4     continue
c
      do 5 m=3, nterms, 4
        do k=1, nterms
          ind=mod(3+3*k,4)
          ind2=mod(ind+2,4)
          mpolepar(m)=mpolepar(m)+c(m,k)*temp(k,ind)
     1                  +c(m,-k)*dconjg(temp(k,ind2))
        enddo
5     continue
c


      return
      end

c**********************************************************************
       subroutine yparentchild(betahatpar, beta1hat, beta2hat, beta3hat,
     $     beta4hat, iflag1,iflag2,iflag3,iflag4,
     $     nterms, c)
c**********************************************************************
c     this subroutine is designed to shift the local expansions of
c     a parent box to its four children.  this is used in the downward
c     pass.
c     on input,
c       betahatpar is the local expansion for the parent box.
c       nterms is the number of terms in the local expansions, and
c       c is an array of precomputed coefficients for the shifting, it is
c         computed in the subroutine tatacoeff().
c
c     on output,
c       beta1hat, beta2hat, beta3hat, and beta4hat are the
c         local expansions of the four child boxes.  they are ordered
c         as in the atandard convention.  clockwise starting from the
c         upper left hand corner.
c
c     the ordering of the children boxes :
c
c         4  1
c         3  2
c
c     subroutine called : none.
c
c     called from :
c
c       adapfmm4().
c       boundfmm4().
c
c**********************************************************************
      implicit none
c
c-----global variables
c
      integer  nterms
      complex *16 c(0:nterms, -nterms:nterms)
      complex *16 beta1hat(0:nterms), beta2hat(0:nterms)
      complex *16 beta3hat(0:nterms), beta4hat(0:nterms)
      complex *16 betahatpar(0:nterms)
      integer iflag1,iflag2,iflag3,iflag4
c
c-----local variables
c
      integer   m, k, ind1,ind2
      real *8 pn1(0:3)
      complex *16 dconjg,tempc1,tempc2
      complex *16 imag
      complex *16 pimag(0:3),pnimag(0:3)
      data imag/(0.0d0,1.0d0)/
      data pimag/(1.0d0, 0.0d0),(0.0d0,1.0d0),
     1           (-1.0d0,0.0d0),(0.0d0,-1.0d0)/
      data pnimag/(1.0d0, 0.0d0),(0.0d0,-1.0d0),
     1           (-1.0d0,0.0d0),(0.0d0,+1.0d0)/
      data pn1/1.0d0,-1.0d0,1.0d0,-1.0d0/
      save pimag, pnimag, pn1
c
c-----first the zeroth term.
c
      if (iflag1 .eq. 1) then
         do m=0,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta1hat(m)=beta1hat(m)+tempc1
         enddo
         do m=1,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta1hat(m)=beta1hat(m)+tempc1
         enddo
         do m=2,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta1hat(m)=beta1hat(m)+tempc1
         enddo
         do m=3,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta1hat(m)=beta1hat(m)+tempc1
         enddo
         
         do m=0,nterms
            do k=1,nterms
               tempc1=betahatpar(k)*c(m,k)
               tempc2=dconjg(betahatpar(k))*c(m,-k)
               beta1hat(m)=beta1hat(m)+tempc1+tempc2
            enddo
         enddo
      endif

      if (iflag2 .eq. 1) then
         do m=0,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta2hat(m)=beta2hat(m)+tempc1
         enddo
         do m=1,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            tempc2=tempc1*imag
            beta2hat(m)=beta2hat(m)+tempc2
         enddo
         do m=2,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta2hat(m)=beta2hat(m)-tempc1
         enddo
         do m=3,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            tempc2=tempc1*imag
            beta2hat(m)=beta2hat(m)-tempc2
         enddo
         
         do m=0,nterms
            do k=1,nterms
               tempc1=betahatpar(k)*c(m,k)
               tempc2=dconjg(betahatpar(k))*c(m,-k)
               ind1=mod(m+3*k,4)
               ind2=mod(m+k,4)
               beta2hat(m)=beta2hat(m)+tempc1*pimag(ind1)
     $              +tempc2*pimag(ind2)
            enddo
         enddo
      endif

      if (iflag3 .eq. 1) then
         do m=0,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta3hat(m)=beta3hat(m)+tempc1
         enddo
         do m=1,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta3hat(m)=beta3hat(m)-tempc1
         enddo
         do m=2,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta3hat(m)=beta3hat(m)+tempc1
         enddo
         do m=3,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta3hat(m)=beta3hat(m)-tempc1
         enddo
         
         do m=0,nterms
            do k=1,nterms
               tempc1=betahatpar(k)*c(m,k)
               tempc2=dconjg(betahatpar(k))*c(m,-k)
               ind1=mod(m+3*k,4)
               ind2=mod(m+k,4)
               beta3hat(m)=beta3hat(m)+tempc1*pn1(ind1)
     $              +tempc2*pn1(ind2)
            enddo
         enddo
      endif

      if (iflag4 .eq. 1) then
         do m=0,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta4hat(m)=beta4hat(m)+tempc1
         enddo
         do m=1,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            tempc2=tempc1*imag
            beta4hat(m)=beta4hat(m)-tempc2
         enddo
         do m=2,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            beta4hat(m)=beta4hat(m)-tempc1
         enddo
         do m=3,nterms,4
            tempc1=betahatpar(0)*c(m,0)
            tempc2=tempc1*imag
            beta4hat(m)=beta4hat(m)+tempc2
         enddo
         
         do m=0,nterms
            do k=1,nterms
               tempc1=betahatpar(k)*c(m,k)
               tempc2=dconjg(betahatpar(k))*c(m,-k)
               ind1=mod(m+3*k,4)
               ind2=mod(m+k,4)
               beta4hat(m)=beta4hat(m)+tempc1*pnimag(ind1)
     $              +tempc2*pnimag(ind2)
            enddo
         enddo
      endif


      return
      end


c***********************************************************************
      subroutine ymkexp2d(ibox,nterms,mpole,
     1           nnodes,mexnall,mexn12,mexn123,mexn124,
     2           mexsall,mexs34,mexs134,mexs234,
     3           mexeall,mexe13,mexe1,mexe3,mexwall,mexw24,mexw2,mexw4,
     4           zs,comp,ichildbox,irowbox,icolbox)
c***********************************************************************
c
c     this subroutine creates the north (+y)  and south (-y) exponential
c     expansions for a parent box due to all four children.
c
c     some intelligence is used in the order of summation.
c
c     on input :
c       ibox : the index of current box.
c       nterms : the number of terms in the multipole expansion.
c       mpole : the multipole expansions
c       nnodes : the number of plane wave nodes.
c       zs : the shifting coefficients.
c       comp : the precomputed coefficients for the mp->pw translation
c              operator.
c       ichildbox, irowbox, icolbox : the tree structure.
c
c     on output :
c       mexnall,mexn12,mexsall,mexs34,mexeall,mexe13,mexe1,
c       mexe3,mexwall,mexw24,mexw2,mexw4 : the plane wave expansion
c         coefficients.
c
c     subroutine called : expcoeff().
c
c***********************************************************************
c
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
c-----global variables.
c
      integer  ibox
      integer  nnodes
      integer  ichildbox(4,ibox)
      integer  icolbox(ibox)
      integer  irowbox(ibox), nterms
      complex *16 mpole(0:nterms,*)
      complex *16 zs(-3:3,-3:3,1:nnodesmax)
      complex *16 expn(1:nnodesmax),exps(1:nnodesmax)
      complex *16 expe(1:nnodesmax),expw(1:nnodesmax)
      complex *16 mexnall(1:50),mexn12(1:50),mexn123(50),mexn124(50)
      complex *16 mexsall(1:50),mexs34(1:50),mexs134(50),mexs234(50)
      complex *16 mexeall(1:nnodesmax),mexe13(1:nnodesmax)
      complex *16 mexe1(1:nnodesmax),mexe3(1:nnodesmax)
      complex *16 mexwall(1:nnodesmax),mexw24(1:nnodesmax)
      complex *16 mexw2(1:nnodesmax),mexw4(1:nnodesmax)
      real *8 comp(nnodesmax,-nterms:nterms)
c
c-----local variables.
c
      integer   ic, jj
      integer   ichild(4)
      complex *16 czero
c
      data czero/(0.0d0,0.0d0)/

      ichild(1) = ichildbox(1,ibox)
      ichild(2) = ichildbox(2,ibox)
      ichild(3) = ichildbox(3,ibox)
      ichild(4) = ichildbox(4,ibox)

c
c-----include contributions from child 1
c

      ic = ichild(4)
      if (ic .ne. -1) then
         call yexpcoeff(expe,expw,expn,exps,nnodes,mpole(0,ic),
     1        nterms, comp)
         do jj = 1,nnodes
            mexn12(jj) = expn(jj)
            mexsall(jj) = exps(jj)
            mexs134(jj) = exps(jj)
         enddo
         do jj = 1,nnodes
            mexe1(jj) = expe(jj)
            mexwall(jj) = expw(jj)
         enddo
      else
         do jj=1,nnodes
            mexn12(jj)=czero
            mexsall(jj)=czero
            mexs134(jj)=czero
         enddo
         do jj=1,nnodes
            mexe1(jj)=czero
            mexwall(jj)=czero
         enddo
      endif


c
c-----include contributions from child 2
c
      ic = ichild(3)
      if (ic .ne. -1) then
         call yexpcoeff(expe,expw,expn,exps,nnodes, mpole(0,ic),
     1        nterms, comp)
         do jj = 1,nnodes
            expn(jj) = expn(jj)*zs(0,1,jj)
            mexn12(jj) = mexn12(jj) + expn(jj)
            exps(jj) = exps(jj)*zs(0,1,jj)
            mexsall(jj) = mexsall(jj) + exps(jj)
            mexs234(jj) = exps(jj)
         enddo
         do jj = 1,nnodes
            mexeall(jj) = expe(jj)*zs(-1,0,jj)
            mexw2(jj) = expw(jj)*zs(1,0,jj)
         enddo
      else
         do jj = 1,nnodes
            expn(jj) = czero
c            mexn12(jj) = mexn12(jj) + czero
            exps(jj) = czero
c            mexsall(jj) = mexsall(jj) + czero
            mexs234(jj) = czero
         enddo
         do jj = 1,nnodes
            mexeall(jj) = czero
            mexw2(jj) = czero
         enddo
      endif
c
c-----include contributions from child 3
c
      ic = ichild(1)
      if (ic .ne. -1) then
         call yexpcoeff(expe,expw,expn,exps,nnodes, mpole(0,ic),
     1        nterms, comp)
c
         do jj = 1,nnodes
            expn(jj) = expn(jj)*zs(-1,0,jj)
            mexnall(jj) = expn(jj) + mexn12(jj)
            mexn123(jj) = mexnall(jj)
            exps(jj) = exps(jj)*zs(1,0,jj)
            mexs34(jj) = exps(jj)
         enddo
c
         do jj = 1,nnodes
            mexe3(jj) = expe(jj)*zs(0,1,jj)
            mexe13(jj) = mexe1(jj) + mexe3(jj)
            mexeall(jj) = mexeall(jj) + mexe13(jj)
            expw(jj) = expw(jj)*zs(0,1,jj)
            mexwall(jj) = mexwall(jj) + expw(jj)
         enddo
      else
         do jj = 1,nnodes
            expn(jj) = czero
            mexnall(jj) = czero + mexn12(jj)
            mexn123(jj) = mexnall(jj)
            exps(jj) = czero
            mexs34(jj) = czero
         enddo
c
         do jj = 1,nnodes
            mexe3(jj) = czero
            mexe13(jj) = mexe1(jj) + mexe3(jj)
            mexeall(jj) = mexeall(jj) + mexe13(jj)
            expw(jj) = czero
c            mexwall(jj) = mexwall(jj) + czero
         enddo
      endif

c
c-----include contributions from child 4
c
      ic = ichild(2)
      if (ic .ne. -1) then
         call yexpcoeff(expe,expw,expn,exps,nnodes, mpole(0,ic),
     1        nterms, comp)
c
         do jj = 1,nnodes
            expn(jj) = expn(jj)*zs(-1,1,jj)
            mexnall(jj) = mexnall(jj) + expn(jj)
            mexn124(jj) = mexn12(jj) + expn(jj)
            exps(jj) = exps(jj)*zs(1,1,jj)
            mexs34(jj) = mexs34(jj) + exps(jj)
            mexsall(jj) = mexsall(jj) + mexs34(jj)
            mexs134(jj) = mexs134(jj) + mexs34(jj)
            mexs234(jj) = mexs234(jj) + mexs34(jj)
         enddo
c
         do jj = 1,nnodes
            expe(jj) = expe(jj)*zs(-1,1,jj)
            mexeall(jj) = mexeall(jj) + expe(jj)
            mexw4(jj) = expw(jj)*zs(1,1,jj)
            mexw24(jj) = mexw2(jj) + mexw4(jj)
            mexwall(jj) = mexwall(jj) + mexw24(jj)
         enddo
      else
         do jj = 1,nnodes
            expn(jj) = czero
            mexnall(jj) = mexnall(jj) + czero
            mexn124(jj) = mexn12(jj) + czero
            exps(jj) = czero
            mexs34(jj) = mexs34(jj) + czero
            mexsall(jj) = mexsall(jj) + mexs34(jj)
            mexs134(jj) = mexs134(jj) + mexs34(jj)
            mexs234(jj) = mexs234(jj) + mexs34(jj)
         enddo
c
         do jj = 1,nnodes
            expe(jj) = czero
            mexeall(jj) = mexeall(jj) + czero
            mexw4(jj) = czero
            mexw24(jj) = mexw2(jj) + mexw4(jj)
            mexwall(jj) = mexwall(jj) + mexw24(jj)
         enddo
      endif
c
      return
      end



c**********************************************************************
      subroutine yexpcoeff(betae, betaw, betan, betas,
     1                    nnodes, mpole, nterms, comp)
c**********************************************************************
c     this subroutine generates four exponential expansions given
c     one multipole expansion as input.
c
c     on input,
c       nnodes : the number of terms in the exponential expansions.
c       mpole : the array of multipole coefficients.
c       nterms : the number of terms in the multipole expansion.
c       comp() : the precomputed coefficients for the multipole to
c                plane wave coefficients.
c
c     on output :
c       betae, betaw, betan, betas : the plane wave coefficients.
c
c**********************************************************************
c
      implicit none
c
      integer nnodesmax
      parameter (nnodesmax=50)
c
c-----global parameters.
c
      complex *16 betae(1:nnodesmax), betaw(1:nnodesmax)
      complex *16 betan(1:nnodesmax), betas(1:nnodesmax)
      integer  nnodes, nterms
      complex *16 mpole(0:nterms)
      real *8 comp(nnodesmax,-nterms:nterms)
c
c-----local parameters.
c
      integer  i, j
      complex *16 imag, conjmp(100)
      complex *16 sum10,sum11,sum2,sum30,sum31,sum4
      data imag/(0.0d0, 1.0d0)/
c
c-----note the difference between the four directions, so we need
c     4 different sum? for different (l mod 4).
c
      do j=1, nterms
        conjmp(j)=dconjg(mpole(j))
      enddo
c
      do i = 1, nnodes
        sum10=0.0d0
        sum11=0.0d0
        do j=1, nterms, 4
          sum10=sum10+mpole(j)*comp(i,j)
          sum11=sum11+conjmp(j)*comp(i,-j)
        enddo
c
        sum2=0.0d0
        do j=2, nterms, 4
          sum2=sum2+mpole(j)*comp(i,j)+conjmp(j)*comp(i,-j)
        enddo
c
        sum30=0.0d0
        sum31=0.0d0
        do j=3, nterms, 4
          sum30=sum30+mpole(j)*comp(i,j)
          sum31=sum31+conjmp(j)*comp(i,-j)
        enddo
c
        sum4=0.0d0
        do j=4, nterms, 4
          sum4=sum4+mpole(j)*comp(i,j)+conjmp(j)*comp(i,-j)
        enddo
        sum4=sum4+mpole(0)*comp(i,0)
c
        betae(i)=sum10+sum11+sum2+sum30+sum31+sum4
        betaw(i)=dconjg(sum2+sum4-sum10-sum11-sum30-sum31)
        betas(i)=sum4-sum2+imag*(-sum10+sum11+sum30-sum31)
        betan(i)=dconjg(sum4-sum2+imag*(sum10-sum11-sum30+sum31))
      end do
c
      return
      end



c**********************************************************************
      subroutine yexp4local(betahat,np, nnodes, betaw,
     1                     betas, betae, betan, comp)
c**********************************************************************
c     the following subroutine takes four multipole expansions (from four
c     different directions) as input after they have been shifted to a
c     common target point.  the four expansions are joined into one
c     common local expansion.
c
c     on input :
c       betae, betaw, betas, and betan : represent the exponential multipole
c         expansions.  the expansions are defined by the direction in which
c         they decay, so the east expansion came from the west, the west
c         expansion came from the east, the north expansion came from the
c         south, and the south expansoin came from the north.
c       p : is the number of terms in the local expansion.
c       nnodes : is the number of terms in the exponential expansion.
c       comp : is the precomputed matrix from the plane wave expansion
c              to multipole expansion.
c
c     on output :
c       betahat : represents the one local expansion.
c
c     none of the input parameters are altered.
c     the expansions are merged by counting through the loop mod 4.
c
c**********************************************************************
c
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
c-----global variables
c
      integer  nnodes, np
      real *8 comp(nnodesmax,-np:np)
      complex *16 betahat(0:np)
      complex *16 betan(1:nnodesmax), betas(1:nnodesmax)
      complex *16 betae(1:nnodesmax), betaw(1:nnodesmax)
c
c-----local variables
c
      integer  i, k
      complex *16 bsum1, bsum2, bsum3, bsum0
      complex *16 asum1, asum2, asum3, asum0
      complex *16 dconjw, dconjn
      complex *16 imag
      data imag/(0.0d0, 1.0d0)/
c
c-----the local coefficients are computed below.
c
      do i = 1, nnodes
         dconjw= dconjg(betaw(i))
         dconjn= dconjg(betan(i))
         bsum0 = betae(i) +dconjw
         bsum1 = -betae(i)+dconjw
         bsum2 = dconjn+betas(i)
         bsum3 = imag*(dconjn-betas(i))
         asum0=bsum0+bsum2
         asum1=bsum1+bsum3
         asum2=bsum0-bsum2
         asum3=bsum1-bsum3
         bsum0=dconjg(asum0)
         bsum1=dconjg(asum3)
         bsum2=dconjg(asum2)
         bsum3=dconjg(asum1)
c
         do k=0, np, 4
           betahat(k)=betahat(k)+comp(i,-k)*asum0+comp(i,k)*bsum0
         enddo
c
         do k=1,np,4
           betahat(k)=betahat(k)+comp(i,-k)*asum1+comp(i,k)*bsum1
         enddo
c
         do k=2,np,4
           betahat(k)=betahat(k)+comp(i,-k)*asum2+comp(i,k)*bsum2
         enddo
c
         do k=3,np,4
           betahat(k)=betahat(k)+comp(i,-k)*asum3+comp(i,k)*bsum3
         enddo
c
      end do
c
      return
      end



c**********************************************************************
      subroutine ymklists(ibox,inall,nnall,iynall,
     1    in12,nn12,iy12,in123,nn123,iy123,in124,nn124,iy124,
     2    isall,nsall,iysall,
     3    is34,ns34,iy34,is134,ns134,iy134,is234,ns234,iy234,
     4    ieall,neall,iyeall,
     5    ie13,ne13,iy13,iwall,nwall,iywall,iw24,nw24,iy24,
     5    iw2,iy2,nw2,iw4,iy4,nw4,ie1,iy1,ne1,ie3,iy3,ne3,
     7    icolleagbox,ichildbox,
     8    icolbox, irowbox,level)
c**********************************************************************
c     this subroutine is set up to compute all of the north, south,
c     east, and west interaction lists.  because this is for the
c     adaptive case, all of these lists have to be computed in one
c     loop.  north and south take precedence over the east and
c     west and this is reflected  in the else statements within
c     the routine.  processing is done at the parent level.
c     in the periodic case, the procedure is basically
c     the same, but if a column or row number lies outside the center
c     box, it is readjusted to account for the periodicty.
c     when the lists are actually processed, there is no difference
c     between free space and periodic case.
c
c     input:
c
c     ibox denotes the box being considered
c
c     level is the level of ibox
c
c     icolleagbox, irowbox, and icolbox define the tree
c
c     output:
c
c     the naming convention for the lists is is as follows:
c
c     inall is an array that denotes the boxes in the
c     north all list.
c
c     nnall is the number of boxes in the north all
c     list.
c
c     iynall represents the corresponding offsets of the boxes
c     in the north all list.
c
c     the same convention is used for the south all, east all, west all,
c     north12, south34, east13, west24, west2, west4, east1,
c     and east3 lists.
c
c**********************************************************************
c
      implicit none
c
c-----global variables
c
      integer  inall(4),nnall,iynall(4)
      integer  isall(4),nsall,iysall(4)
      integer  ieall(4),neall,iyeall(4)
      integer  iwall(4),nwall,iywall(4)
      integer  in12(4),nn12,iy12(4)
      integer  is34(4),ns34,iy34(4)
      integer  iw24(2),nw24,iy24(2)
      integer  ie13(2),ne13,iy13(2)
      integer  iw4(2),nw4,iy4(2)
      integer  iw2(2),nw2,iy2(2)
      integer  ie1(2),ne1,iy1(2)
      integer  ie3(2),ne3,iy3(2)
      integer  in123(1),nn123,iy123(1)
      integer  in124(1),nn124,iy124(1)
      integer  is134(1),ns134,iy134(1)
      integer  is234(1),ns234,iy234(1)
      integer  ibox
      integer  icolleagbox(9,1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  level
c
c-----local variables
c
      integer  i, j, iout
      integer  ichild, ncntr, scntr, ecntr, wcntr
      integer  n12cntr, s34cntr
      integer  w24cntr, e13cntr
      integer  w4cntr, w2cntr
      integer  e1cntr, e3cntr
      integer  n123cntr,n124cntr,s134cntr,s234cntr
      integer  icoltest, irowtest
      integer  icol, irow
c
c-----initially, set all list entries
c     and offsets to zero.
c
      do j = 1, 4
        inall(j)  = 0
        iynall(j) = 0
        isall(j)  = 0
        iysall(j) = 0
        ieall(j)  = 0
        iyeall(j) = 0
        iwall(j)  = 0
        iywall(j) = 0
        in12(j) = 0
        iy12(j) = 0
        is34(j) = 0
        iy34(j) = 0
      end do
c
      do j = 1, 2
        ie13(j) = 0
        iy13(j) = 0
        iw24(j) = 0
        iy24(j) = 0
        ie1(j) = 0
        iy1(j) = 0
        iw2(j) = 0
        iy2(j) = 0
        ie3(j) = 0
        iy3(j) = 0
        iw4(j) = 0
        iy4(j) = 0
      end do
c
      in123(1) = 0
      iy123(1) = 0
      in124(1) = 0
      iy124(1) = 0
      is134(1) = 0
      iy134(1) = 0
      is234(1) = 0
      iy234(1) = 0
c
c-----all of the offsets are set based from the
c     box in the lower left corner (child 4 in
c     out ordering convention).
c     icol and irow are the rows and columns of
c     the box whose list we are trying to generate.
c
c      icol = icolbox(ichildbox(4,ibox))
c      irow = irowbox(ichildbox(4,ibox))

      icol=2*icolbox(ibox)
      irow=2*irowbox(ibox)
      icol=icol-1
      irow=irow-1
c     the fourth child may not exist, so need to make sure that
c     icol and irow are set using a different method.  Done by 
c     computing the values icolbox and irowbox would be if it had a 
c     fourth child.


c
c-----first do the free space case:
c
c
c-------set all of the counters to 1
c
      ncntr   =  1
      scntr   =  1
      ecntr   =  1
      wcntr   =  1
      n12cntr =  1
      s34cntr =  1
      w24cntr =  1
      e13cntr =  1
      w4cntr  =  1
      w2cntr  =  1
      e1cntr  =  1
      e3cntr  =  1
      n123cntr=  1
      n124cntr=  1
      s134cntr=  1
      s234cntr=  1
c
c-------first scan through all nine of the boxes colleagues
c    


      do 100 i = 1, 9
         iout = icolleagbox(i,ibox)
c
c---------test to see if this colleague doesn't exist or is
c         childless, if so, skip it.
c
         if(iout .lt. 0)goto 100
c         if(ichildbox(1,iout) .lt. 0)goto 100
c
c---------scan all four of the colleagues children.
c
         do j = 1, 4
c
c-----------icoltest and irowtest represent the row and column
c           of the box being checked.
c
            ichild = ichildbox(j,iout)
            if (ichild .ne. -1) then
               icoltest = icolbox(ichild)
               irowtest = irowbox(ichild)
c     
               if(irowtest .eq. irow+3)then
                  if (icoltest .eq. icol-2) then
                     in123(n123cntr)=ichild
                     iy123(n123cntr)=icoltest-icol
                     n123cntr=n123cntr+1
                     iw4(w4cntr) = ichild
                     iy4(w4cntr) = irowtest - irow
                     w4cntr = w4cntr + 1
                  elseif (icoltest .eq. icol+3) then
                     in124(n124cntr)=ichild
                     iy124(n124cntr)=icoltest-icol
                     n124cntr=n124cntr+1
                     ie3(e3cntr) = ichild
                     iy3(e3cntr) = irowtest - irow
                     e3cntr = e3cntr + 1
                  else
                     inall(ncntr) = ichild
                     iynall(ncntr) = icoltest - icol
                     ncntr = ncntr + 1
                  endif
               elseif(irowtest .eq. irow-2)then
                  if (icoltest .eq. icol-2) then
                     is134(s134cntr)=ichild
                     iy134(s134cntr)=icoltest-icol
                     s134cntr=s134cntr+1
                     iw2(w2cntr) = ichild
                     iy2(w2cntr) = irowtest - irow
                     w2cntr = w2cntr + 1
                  elseif (icoltest .eq. icol+3) then
                     is234(s234cntr)=ichild
                     iy234(s234cntr)=icoltest-icol
                     s234cntr=s234cntr+1
                     ie1(e1cntr) = ichild
                     iy1(e1cntr) = irowtest - irow
                     e1cntr = e1cntr + 1
                  else
                     isall(scntr) = ichild
                     iysall(scntr) = icoltest - icol
                     scntr = scntr + 1
                  endif
               elseif(icoltest .eq. icol-2)then
                  iwall(wcntr) = ichild
                  iywall(wcntr) = irowtest - irow
                  wcntr = wcntr + 1
               elseif(icoltest .eq. icol+3)then
                  ieall(ecntr) = ichild
                  iyeall(ecntr) = irowtest - irow
                  ecntr = ecntr + 1
               elseif(irowtest .eq. irow+2)then
                  in12(n12cntr) = ichild
                  iy12(n12cntr) = icoltest - icol
                  n12cntr = n12cntr + 1
                  if(icoltest .eq. icol-1)then
                     iw4(w4cntr) = ichild
                     iy4(w4cntr) = irowtest - irow
                     w4cntr = w4cntr + 1
                  endif
                  if(icoltest .eq. icol+2)then
                     ie3(e3cntr) = ichild
                     iy3(e3cntr) = irowtest - irow
                     e3cntr = e3cntr + 1
                  endif
               elseif(irowtest .eq. irow-1)then
                  is34(s34cntr) = ichild
                  iy34(s34cntr) = icoltest - icol
                  s34cntr = s34cntr + 1
                  if(icoltest .eq. icol-1)then
                     iw2(w2cntr) = ichild
                     iy2(w2cntr) = irowtest - irow
                     w2cntr = w2cntr + 1
                  endif
                  if(icoltest .eq. icol+2)then
                     ie1(e1cntr) = ichild
                     iy1(e1cntr) = irowtest - irow
                     e1cntr = e1cntr + 1
                  endif
               elseif(icoltest .eq. icol-1)then
                  iw24(w24cntr) = ichild
                  iy24(w24cntr) = irowtest - irow
                  w24cntr = w24cntr + 1
               elseif(icoltest .eq. icol+2)then
                  ie13(e13cntr) = ichild
                  iy13(e13cntr) = irowtest - irow
                  e13cntr = e13cntr + 1
               endif
            endif
         enddo
 100  continue
      nnall = ncntr   -  1
      nsall = scntr   -  1
      neall = ecntr   -  1
      nwall = wcntr   -  1
      nn12  = n12cntr -  1
      ns34  = s34cntr -  1
      nw24  = w24cntr -  1
      ne13  = e13cntr -  1
      nw4   = w4cntr  -  1
      nw2   = w2cntr  -  1
      ne1   = e1cntr  -  1
      ne3   = e3cntr  -  1
      nn123 = n123cntr-  1
      nn124 = n124cntr-  1
      ns134 = s134cntr-  1
      ns234 = s234cntr-  1


      return
      end




c************************************************************************
      subroutine localexpansion(layers,locexp,icol,irow,level,z,
     $     pot,potx,poty,potxx,potxy,potyy,betascal,beta,nterms)
c     Evaulate the local expansion at the point z

c     Input:
c     layers - number of layer potentials to compute
c     locexp - local expansion coefficients
c     icol - column of finest box containing z
c     irow - row of finest box containing z
c     level - level of box containing z
c     z - location of target
c     pot - potential at z
c     potx - x-component of derivative of potential
c     poty - y-component of derivative of potential
c     betascal - scaled values of beta
c     beta - Yukawa parameter
c     nterms - Size of local expansion

c     Output
c     pot - updated potential at z
c     potx - updated x-component of derivative of potential
c     poty - updated y-component of derivative of potential
c     potxx - updated 2nd derivative of potential
c     potxy - updated 2nd derivative of potential
c     potyy - updated 2nd derivative of potential

      implicit none

      integer layers
      complex(8) locexp(0:nterms)
      complex(8) z,pot,potx,poty
      complex(8) potxx,potxy,potyy
      integer *4 icol,irow,level
      real *8 betascal,beta
      integer *4 nterms

      integer i
      real *8 r,theta
      complex(8) eye,cen
      complex(8) val0,val1,val2
      real *8 bi(0:100)
      real *8 sscale

      data eye/(0.0d0,1.0d0)/

      cen=dble(2*icol-1)+eye*dble(2*irow-1)
      cen=cen*(0.5d0)**(level+1)-(0.5d0,0.5d0)

      call polar_convert(z-cen,r,theta)
      call besseli(nterms+1,r*beta,bi,betascal)

      if (layers .eq. 1) then
        val0=locexp(0)*bi(0)
        pot=pot+val0
        do i=1,nterms
           val1=locexp(i)*bi(i)*exp(eye*i*theta)
           val2=dconjg(locexp(i))*bi(i)*exp(-1.0d0*eye*i*theta)
           pot=pot+val1+val2
        enddo
c       Evaluate the local expansion for the single layer potential
      elseif (layers .eq. 2) then
        val0=locexp(0)*bi(0)
        pot=pot+val0
        do i=1,nterms
           val1=locexp(i)*bi(i)*exp(eye*i*theta)
           val2=dconjg(locexp(i))*bi(i)*exp(-1.0d0*eye*i*theta)
           pot=pot+val1+val2
        enddo
c       Evaluate the local expansion for the single layer potential

        sscale = betascal/2.d0
        val0 = sscale*locexp(0)*bi(1)*beta*real(z-cen)/r
        potx = potx - val0
        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = locexp(i)*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + locexp(i)*bi(i)*exp(eye*i*theta)*
     $      (-eye*i*imag(z-cen)/r/r)
          val2 = dconjg(locexp(i))*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(-eye*i*theta)
          val2 = val2 + dconjg(locexp(i))*bi(i)*exp(-eye*i*theta)*
     $      (eye*i*imag(z-cen)/r/r)
          potx = potx - val1 - val2
        enddo
c     Evaluate the x componenet of the gradient of the local
c     expansion

        sscale = betascal/2.d0
        val0 = sscale*locexp(0)*bi(1)*beta*imag(z-cen)/r
        poty = poty - val0
        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = locexp(i)*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + locexp(i)*bi(i)*exp(eye*i*theta)*
     $      (eye*i*real(z-cen)/r/r)
          val2 = dconjg(locexp(i))*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(-eye*i*theta)
          val2 = val2 + dconjg(locexp(i))*bi(i)*exp(-eye*i*theta)*
     $     (-eye*i*real(z-cen)/r/r)
          poty = poty - val1 - val2
        enddo
c     Evaluate the y componenet of the gradient of the local
c     expansion


      elseif (layers .eq. 3) then
        val0=locexp(0)*bi(0)
        pot=pot+val0
        do i=1,nterms
           val1=locexp(i)*bi(i)*exp(eye*i*theta)
           val2=dconjg(locexp(i))*bi(i)*exp(-1.0d0*eye*i*theta)
           pot=pot+val1+val2
        enddo
c       Evaluate the local expansion for the single layer potential

        sscale = betascal/2.d0
        val0 = sscale*locexp(0)*bi(1)*beta*real(z-cen)/r
        potx = potx - val0
        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = locexp(i)*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + locexp(i)*bi(i)*exp(eye*i*theta)*
     $      (-eye*i*imag(z-cen)/r/r)
          val2 = dconjg(locexp(i))*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(-eye*i*theta)
          val2 = val2 + dconjg(locexp(i))*bi(i)*exp(-eye*i*theta)*
     $      (eye*i*imag(z-cen)/r/r)
          potx = potx - val1 - val2
        enddo
c     Evaluate the x componenet of the gradient of the local
c     expansion

        sscale = betascal/2.d0
        val0 = sscale*locexp(0)*bi(1)*beta*imag(z-cen)/r
        poty = poty - val0
        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = locexp(i)*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + locexp(i)*bi(i)*exp(eye*i*theta)*
     $      (eye*i*real(z-cen)/r/r)
          val2 = dconjg(locexp(i))*(sscale*bi(i+1)+i*bi(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(-eye*i*theta)
         val2 = val2 + dconjg(locexp(i))*bi(i)*exp(-eye*i*theta)*
     $     (-eye*i*real(z-cen)/r/r)
          poty = poty - val1 - val2
        enddo
c     Evaluate the y component of the gradient of the local
c     expansion

        sscale = betascal/2.d0
        val0 = -real(z-cen)**4.d0*bi(0) - 
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0*bi(0)+
     $      sscale*r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(1)
        val0 = val0*locexp(0)/r**4.d0
        potxx = potxx - val0

        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = (-real(z-cen)**4.d0-
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0-
     $      i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)+
     $      i*(real(z-cen)**2.d0-imag(z-cen)**2.d0))*bi(i)+
     $      r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*sscale*bi(i+1)+
     $      2.d0*eye*real(z-cen)*imag(z-cen)*i*
     $      (i*bi(i)-bi(i)+r*sscale*bi(i+1))
          val1 = val1*locexp(i)*exp(eye*i*theta)/r**4.d0
          sscale = 2.d0*dble(i)/betascal
          val2 = (-real(z-cen)**4.d0-
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0-
     $      i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)-
     $      i*(real(z-cen)**2.d0-imag(z-cen)**2.d0))*bi(i)+
     $      r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*sscale*bi(i-1)-
     $      2.d0*eye*real(z-cen)*imag(z-cen)*i*
     $      (-i*bi(i)-bi(i)+r*sscale*bi(i-1))
          val2 = val2*dconjg(locexp(i))*exp(-eye*i*theta)/r**4.d0
          potxx = potxx - val1 - val2
        enddo
c       Evaulate the second-derivative with respect to x


        sscale = betascal/2.d0
        val0 = real(z-cen)*imag(z-cen)*r**2.d0*bi(0) -
     $      2.d0*real(z-cen)*imag(z-cen)*r*sscale*bi(1)
        val0 = val0*locexp(0)/r**4.d0
        potxy = potxy + val0

        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = real(z-cen)*imag(z-cen)*r**2.d0*bi(i) + 
     $      eye*i*r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*
     $      sscale*bi(i+1)-
     $      eye*i*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(i) +
     $      eye*i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(i)-
     $      2.d0*real(z-cen)*imag(z-cen)*
     $      (sscale*r*bi(i+1)+i*bi(i)-i**2.d0*bi(i))
          val1 = val1*locexp(i)*exp(eye*i*theta)/r**4.d0
          sscale = 2.d0*dble(i)/betascal
          val2 = real(z-cen)*imag(z-cen)*r**2.d0*bi(i) - 
     $      eye*i*r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*
     $      sscale*bi(i-1)+
     $      eye*i*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(i) +
     $      eye*i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(i)-
     $      2.d0*real(z-cen)*imag(z-cen)*
     $      (sscale*r*bi(i-1)-i*bi(i)-i**2.d0*bi(i))
          val2 = val2*dconjg(locexp(i))*exp(-eye*i*theta)/r**4.d0
c          val2 = dconjg(val1)
          potxy = potxy + val1 + val2
        enddo
c       Evaulate the second-mixed-derivative 




        sscale = betascal/2.d0
        val0 = -imag(z-cen)**4.d0*bi(0) - 
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0*bi(0)-
     $      sscale*r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*bi(1)
        val0 = val0*locexp(0)/r**4.d0
        potyy = potyy - val0

        do i=1,nterms
          sscale = betascal/2.d0/dble(i+1)
          val1 = (-imag(z-cen)**4.d0-
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0+
     $      i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)-
     $      i*(real(z-cen)**2.d0-imag(z-cen)**2.d0))*bi(i)-
     $      r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*sscale*bi(i+1)-
     $      2.d0*eye*real(z-cen)*imag(z-cen)*i*
     $      (i*bi(i)-bi(i)+r*sscale*bi(i+1))
          val1 = val1*locexp(i)*exp(eye*i*theta)/r**4.d0
          sscale = 2.d0*dble(i)/betascal
          val2 = (-imag(z-cen)**4.d0-
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0+
     $      i**2.d0*(real(z-cen)**2.d0-imag(z-cen)**2.d0)+
     $      i*(real(z-cen)**2.d0-imag(z-cen)**2.d0))*bi(i)-
     $      r*(real(z-cen)**2.d0-imag(z-cen)**2.d0)*sscale*bi(i-1)+
     $      2.d0*eye*real(z-cen)*imag(z-cen)*i*
     $      (-i*bi(i)-bi(i)+r*sscale*bi(i-1))
          val2 = val2*dconjg(locexp(i))*exp(-eye*i*theta)/r**4.d0
          potyy = potyy - val1 - val2
        enddo
c       Evaulate the second-derivative with respect to y

      endif


      return
      end



c***************************************************************
      subroutine mpole_eval(layers,z,pot,potx,poty,potxx,potxy,potyy,
     $      mpole,icol,irow,level,beta,betascal,nterms)
c     Evaluate the multipole expansion.  This is how list 3 affects
c     a childless box

c     Input:
c     layers - number of layer potentials to compute
c     z - location of target strored as a complex variable
c     pot - potential at target point
c     potx - x-component of gradient of potential at target point
c     poty - y-component of gradient of potential at target point
c     potxx - x-component of second derivative of potential at target point
c     potxy - mixed component of second derivative of potential at target point
c     potyy - y-component of second derivative of potential at target point
c     mpole - multipole coefficients for box that is well-seperated
c       from z
c     icol - column of box whose multipole expansion we are evaulating
c     irow - row of box whose multipole expansion we are evaulating
c     level - level of box whose multipole expansion we are evaulating
c
c     Output
c     pot - updated potential at target point
c     potx - updated x-component of gradient of potential at 
c       target point
c     poty - updated y-component of gradient of potential at 
c       target point

      implicit none

      complex(8) z,pot,potx,poty,potxx,potxy,potyy
      complex(8) mpole(0:nterms)
      real *8 bk(0:100)
      integer icol,irow,level,nterms,layers
      real *8 beta,betascal

      integer i
      real *8 sscale
      real *8 r,theta
      complex *16 cen
      complex *16 val0,val1,val2
      complex *16 eye

      eye=(0.0d0,1.0d0)

      cen=dble(2*icol-1)+eye*dble(2*irow-1)
      cen=cen*(0.5d0)**(level+1)-(0.5d0,0.5d0)

      call polar_convert(z-cen,r,theta)
      call besselk(nterms+1,r*beta,bk,betascal)

      if (layers .eq. 1) then
        val0 = mpole(0)*bk(0)
        pot = pot + val0
        do i = 1,nterms
          val1 = mpole(i)*bk(i)*exp(eye*i*theta)
          val2 = dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)
          pot = pot+val1+val2
        enddo
c       Evaluate the multipole expansion for the single-layer potential

      elseif (layers .eq. 2) then
        val0 = mpole(0)*bk(0)
        pot = pot + val0
        do i = 1,nterms
          val1 = mpole(i)*bk(i)*exp(eye*i*theta)
          val2 = dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)
          pot = pot+val1+val2
        enddo
c       Evaluate the multipole expansion for the single-layer potential
        
        sscale = 2.d0/betascal
        val0 = -sscale*mpole(0)*bk(1)*beta*real(z-cen)/r
        potx = potx - val0
        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = mpole(i)*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + mpole(i)*bk(i)*exp(eye*i*theta)*
     $      (-eye*i*imag(z-cen)/r/r)
c          val2 = dconjg(mpole(i))*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
c     $      beta*real(z-cen)/r*exp(-eye*i*theta)
c          val2 = val2 + dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)*
c     $      (eye*i*imag(z-cen)/r/r)
          val2 = dconjg(val1)
          potx = potx - val1 - val2
        enddo
c     Evaluate the multipole expansion for the derivative w.r.t. the
c     first variable of the single-layer potential

        sscale = 2.d0/betascal
        val0 = -sscale*mpole(0)*bk(1)*beta*imag(z-cen)/r
        poty = poty - val0
        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = mpole(i)*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + mpole(i)*bk(i)*exp(eye*i*theta)*
     $      (eye*i*real(z-cen)/r/r)
c          val2 = dconjg(mpole(i))*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
c     $      beta*imag(z-cen)/r*exp(-eye*i*theta)
c          val2 = val2 + dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)*
c     $      (-eye*i*real(z-cen)/r/r)
          val2 = dconjg(val1)
          poty = poty - val1 - val2
        enddo
c     Evaluate the multipole expansion for the derivative w.r.t. the
c     second variable of the single-layer potential

      elseif (layers .eq. 3) then
        val0 = mpole(0)*bk(0)
        pot = pot + val0
        do i = 1,nterms
          val1 = mpole(i)*bk(i)*exp(eye*i*theta)
          val2 = dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)
          pot = pot+val1+val2
        enddo
c       Evaluate the multipole expansion for the single-layer potential
        
        sscale = 2.d0/betascal
        val0 = -sscale*mpole(0)*bk(1)*beta*real(z-cen)/r
        potx = potx - val0
        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = mpole(i)*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
     $      beta*real(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + mpole(i)*bk(i)*exp(eye*i*theta)*
     $      (-eye*i*imag(z-cen)/r/r)
c          val2 = dconjg(mpole(i))*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
c     $      beta*real(z-cen)/r*exp(-eye*i*theta)
c          val2 = val2 + dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)*
c     $      (eye*i*imag(z-cen)/r/r)
          val2 = dconjg(val1)
          potx = potx - val1 - val2
        enddo
c     Evaluate the multipole expansion for the derivative w.r.t. the
c     first variable of the single-layer potential

        sscale = 2.d0/betascal
        val0 = -sscale*mpole(0)*bk(1)*beta*imag(z-cen)/r
        poty = poty - val0
        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = mpole(i)*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
     $      beta*imag(z-cen)/r*exp(eye*i*theta)
          val1 = val1 + mpole(i)*bk(i)*exp(eye*i*theta)*
     $      (eye*i*real(z-cen)/r/r)
c          val2 = dconjg(mpole(i))*(-sscale*bk(i+1)+i*bk(i)/beta/r)*
c     $      beta*imag(z-cen)/r*exp(-eye*i*theta)
c          val2 = val2 + dconjg(mpole(i))*bk(i)*exp(-eye*i*theta)*
c     $      (-eye*i*real(z-cen)/r/r)
          val2 = dconjg(val1)
          poty = poty - val1 - val2
        enddo
c     Evaluate the multipole expansion for the derivative w.r.t. the
c     second variable of the single-layer potential

        sscale = 2.d0/betascal
        val0 = (real(z-cen)**4.d0 + 
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0)*bk(0) + 
     $      sscale*(real(z-cen)**2.d0 - imag(z-cen)**2.0d0)*r*bk(1)
        val0 = mpole(0)*val0/r**4.d0
        potxx = potxx + val0

        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = (real(z-cen)**4.d0 +
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0 + 
     $      dble(i)*(-real(z-cen)**2.d0 + imag(z-cen)**2.d0) + 
     $      dble(i*i)*(real(z-cen)**2.0d0 - imag(z-cen)**2.0d0))*bk(i)+
     $      (2.d0*real(z-cen)*imag(z-cen)*dble(i - i*i))*eye*bk(i) +
     $      (real(z-cen)**2.d0 - imag(z-cen)**2.d0 + 2.d0*eye*i*
     $      real(z-cen)*imag(z-cen))*r*sscale*bk(i+1)
          val1 = mpole(i) * val1*exp(eye*i*theta)/r**4.d0
          val2 = dconjg(val1)
c          val2 = (real(z-cen)**4.d0 +
c     $      real(z-cen)**2.d0*imag(z-cen)**2.d0 + 
c     $      dble(i)*(-real(z-cen)**2.d0 + imag(z-cen)**2.d0) + 
c     $      dble(i*i)*(real(z-cen)**2.0d0 - imag(z-cen)**2.0d0))*bk(i)+
c     $      (2.d0*real(z-cen)*imag(z-cen)*dble(-i + i*i))*eye*bk(i) +
c     $      (real(z-cen)**2.d0 - imag(z-cen)**2.d0 - 2.d0*eye*i*
c     $      real(z-cen)*imag(z-cen))*r*sscale*bk(i+1)
c          val2 = dconjg(mpole(i)) * val2*exp(-eye*i*theta)/r**4.d0
          potxx = potxx + val1 + val2
        enddo

       sscale = 2.d0/betascal
       val0 = (real(z-cen)**3.d0*imag(z-cen) +
     $    real(z-cen)*imag(z-cen)**3.d0)*bk(0) + 
     $    r*sscale*bk(1)*(2.d0*real(z-cen)*imag(z-cen))
        val0 = mpole(0)*val0/r**4.d0
        potxy = potxy + val0

        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal 
          val1 = (real(z-cen)**3.d0*imag(z-cen) + 
     $      real(z-cen)*imag(z-cen)**3.d0 + 
     $      2.d0*real(z-cen)*imag(z-cen)*(i*i-i))*bk(i) +
     $      (real(z-cen)**2.d0*(-i+i*i) + 
     $      imag(z-cen)**2.d0*(-i*i+i))*eye*bk(i) +
     $      (-eye*i*real(z-cen)**2.d0+2.d0*real(z-cen)*imag(z-cen)+
     $      eye*i*imag(z-cen)**2.d0)*r*sscale*bk(i+1)
          val1 = mpole(i) * val1*exp(eye*i*theta)/r**4.d0
          potxy = potxy + val1 + dconjg(val1)
        enddo




        sscale = 2.d0/betascal
        val0 = (imag(z-cen)**4.d0 + 
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0)*bk(0) + 
     $      sscale*(-real(z-cen)**2.d0 + imag(z-cen)**2.0d0)*r*bk(1)
        val0 = mpole(0)*val0/r**4.d0
        potyy = potyy + val0

        do i=1,nterms
          sscale = 2.d0*dble(i+1)/betascal
          val1 = (imag(z-cen)**4.d0 +
     $      real(z-cen)**2.d0*imag(z-cen)**2.d0 + 
     $      dble(i)*(real(z-cen)**2.d0 - imag(z-cen)**2.d0) + 
     $      dble(i*i)*(-real(z-cen)**2.0d0 + imag(z-cen)**2.0d0))*bk(i)+
     $      (2.d0*real(z-cen)*imag(z-cen)*dble(-i + i*i))*eye*bk(i) +
     $      (-real(z-cen)**2.d0 + imag(z-cen)**2.d0 - 2.d0*eye*i*
     $      real(z-cen)*imag(z-cen))*r*sscale*bk(i+1)
          val1 = mpole(i) * val1*exp(eye*i*theta)/r**4.d0
          val2 = dconjg(val1)
          potyy = potyy + val1 + val2
        enddo


      endif
  

      return
      end



c***************************************************************
      subroutine update_local_coe(z,q,icol,irow,level,
     $     betascal,nterms,local)
c     Update local coefficients of a box due to a particle.  This is
c     how a childless box affects a box in its list3

c     Input:
c     z - location of source point
c     q - charge at source point
c     icol - column of target box
c     irow - row of target box
c     level - level of target box
c     betascal - scaled values of beta
c     nterms - number of terms in local expansion
c     local - local coefficients
c
c     Output:
c     local - updated local coefficients
c

      implicit none

      complex *16 z,cen
      real *8 q
      integer icol,irow,level,nterms
      real *8 betascal(0:1)
      real *8 bk(0:100)
      complex *16 local(0:nterms)

      integer i
      real *8 r,theta,beta
      complex *16 eye

      data eye/(0.0d0,1.0d0)/

      beta=betascal(0)
      cen=dble(2*icol-1)+eye*dble(2*irow-1)
      cen=cen*(0.5d0)**(level+1)-(0.5d0,0.5d0)

      call polar_convert(z-cen,r,theta)
      call besselk(nterms+1,r*beta,bk,betascal(level))

      do i=0,nterms
       local(i) = local(i) + q*bk(i)*exp(-eye*i*theta)
      enddo
c     Update the local coefficients of the box corresponding
c     to icol,irow due to the source point z


      return
      end



c************************************************************************
      subroutine direct_eval(npts,layers,beta,dom,charge,permute,
     $    fpt1,fpt2,numptsbox1,numptsbox2,
     $    pot,potx,poty,potxx,potxy,potyy)
c     Direct evaluation for List 1

c     Input:
c     npts - total number of sources/targets
c     layers - number of layer potentials to compute
c     dom - set of sources/targets
c     charge - charge of sources
c     permute - ordering of points so that we can find them in the
c       tree structure
c     fpt1 - first point inside box of target points
c     fpt2 - first point inside box of source points
c     numptsbox1 - number of points inside box of target points
c     numptsbox2 - number of points inside box of source points
c     pot - potential at target points
c     potx - x-component of gradient of potential at target points
c     poty - y-component of gradient of potential at target points
c
c     Output:
c     pot - updated potential at target points
c     potx - updated x-component of gradient of potential at 
c       target points
c     poty - updated y-component of gradient of potential at 
c       target points
c     potxx - updated x-component of second derivative of potential at 
c       target points
c     potxy - updated mixed component of second derivative of 
c       potential at target points
c     potyy - updated y-component of second derivative of potential at 
c       target points

      implicit none
      
      integer *4 npts,layers
      real *8 charge(npts)
      complex *16 dom(npts)
      integer *4 permute(npts)
      integer *4 fpt1,fpt2
      integer *4 numptsbox1,numptsbox2
      complex *16 pot(npts)
      complex *16 potx(npts),poty(npts)
      complex *16 potxx(npts),potxy(npts),potyy(npts)
      real *8 beta

      integer i,j,index1,index2
      real *8 bk(0:1),dist
      complex *16 z1,z2


      if (layers .eq. 1) then
        do i=fpt1,fpt1+numptsbox1-1
c       Loop over target points
          index1 = permute(i)
          z1 = dom(index1)
          do j=fpt2,fpt2+numptsbox2-1
c       Loop over source points
            index2 = permute(j)
            z2 = dom(index2)
            if (index1 .ne. index2) then
              dist = abs(z1-z2)
              call besselk(1,beta*dist,bk,2.0d0)
              pot(index1)=pot(index1)+charge(index2)*bk(0)
            endif
          enddo
        enddo
c      update potential directly


      elseif (layers .eq. 2) then
        do i=fpt1,fpt1+numptsbox1-1
c       Loop over target points
          index1 = permute(i)
          z1 = dom(index1)
          do j=fpt2,fpt2+numptsbox2-1
c       Loop over source points
            index2 = permute(j)
            z2 = dom(index2)
            if (index1 .ne. index2) then
              dist = abs(z1-z2)
              call besselk(2,beta*dist,bk,2.0d0)
              pot(index1)=pot(index1)+charge(index2)*bk(0)
              potx(index1) = potx(index1) + charge(index2)*bk(1)*
     $            beta/dist*real(z1-z2)
              poty(index1) = poty(index1) + charge(index2)*bk(1)*
     $            beta/dist*imag(z1-z2)
            endif
          enddo
        enddo
c     Update potential and its gradient directly


      elseif (layers .eq. 3) then
        do i=fpt1,fpt1+numptsbox1-1
c       Loop over target points
          index1 = permute(i)
          z1 = dom(index1)
          do j=fpt2,fpt2+numptsbox2-1
c       Loop over source points
            index2 = permute(j)
            z2 = dom(index2)
            if (index1 .ne. index2) then
              dist = abs(z1-z2)
              call besselk(2,beta*dist,bk,2.0d0)
              pot(index1)=pot(index1)+charge(index2)*bk(0)
              potx(index1) = potx(index1) + charge(index2)*bk(1)*
     $            beta/dist*real(z1-z2)
              poty(index1) = poty(index1) + charge(index2)*bk(1)*
     $            beta/dist*imag(z1-z2)
              potxx(index1) = potxx(index1) + charge(index2)*
     $            (bk(0)/dist**2.d0*real(z1-z2)**2.d0 + 
     $            bk(1)/dist - 
     $            2.d0*imag(z1-z2)**2.d0*bk(1)/dist**3.d0)
              potxy(index1) = potxy(index1) + charge(index2)*
     $            (bk(0)/dist**2.d0 + bk(1)/dist**3.d0 *
     $            2.d0)*real(z1-z2)*imag(z1-z2)
              potyy(index1) = potyy(index1) + charge(index2)*
     $            (bk(0)/dist**2.d0*imag(z1-z2)**2.d0 + 
     $            bk(1)/dist - 
     $            2.d0*real(z1-z2)**2.d0*bk(1)/dist**3.d0)
            endif
          enddo
        enddo
c     Update potential, gradient, and second derivative directly
      endif



      return
      end



c********1*********2*********3*********4*********5*********6*********7**
      subroutine yprocessno(lexp1,inall,nnall,iynall,
     1           in12,nn12,iy12,in123,nn123,iy123,in124,nn124,iy124,
     2           mexnall,mexn12,mexn123,mexn124,zs,nnodes)
c***********************************************************************
c
c     this subroutine processes the north interaction lists.
c
c     input:
c
c     inall(nnall), iynall(nnall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other north lists are similarly defined (see ymklists).
c
c     mexnall is the exponential expansion for all eight children, etc.
c     zs are the shift coefficients computed by subroutine
c          mkexps.
c
c     output:
c
c     lexp1, which contains the local north expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
c
c***********************************************************************
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
      integer  nnall,inall(nnall),iynall(nnall)
      integer  nn12,in12(nn12),iy12(nn12)
      integer  nn123,in123(1),iy123(1)
      integer  nn124,in124(1),iy124(1)
      integer  i, jj
      integer  nnodes
c
      complex *16 lexp1(1:nnodesmax,*)
      complex *16 mexnall(1:nnodesmax)
      complex *16 mexn12(1:nnodesmax)
      complex *16 mexn123(1:nnodesmax)
      complex *16 mexn124(1:nnodesmax)
      complex *16 zs(-3:3,-3:3,nnodesmax)
      complex *16 zmul
c
      do i = 1,nnall
         do jj = 1,nnodes
            zmul = zs(3,-iynall(i),jj)
            lexp1(jj,inall(i)) = lexp1(jj,inall(i)) +
     1            mexnall(jj)*zmul
         enddo
      enddo
c
      do i = 1,nn12
         do jj = 1,nnodes
            zmul = zs(2,-iy12(i),jj)
            lexp1(jj,in12(i)) = lexp1(jj,in12(i)) +
     1            mexn12(jj)*zmul
         enddo
      enddo
c
      do i = 1,nn123
         do jj = 1,nnodes
            zmul = zs(3,-iy123(i),jj)
            lexp1(jj,in123(i)) = lexp1(jj,in123(i)) +
     1            mexn123(jj)*zmul
         enddo
      enddo
c
      do i = 1,nn124
         do jj = 1,nnodes
            zmul = zs(3,-iy124(i),jj)
            lexp1(jj,in124(i)) = lexp1(jj,in124(i)) +
     1            mexn124(jj)*zmul
         enddo
      enddo
c
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      subroutine yprocessso(lexp2,isall,nsall,iysall,
     1           is34,ns34,iy34,is134,ns134,iy134,is234,ns234,iy234,
     2           mexsall,mexs34,mexs134,mexs234,zs,nnodes)
c***********************************************************************
c
c     this subroutine processes the south interaction lists.
c
c     input:
c
c     isall(nsall), iysall(nsall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other south lists are similarly defined (see mksolist).
c
c     mexsall is the exponential expansion for all eight children, etc.
c     zs are the shift coefficients computed by subroutine
c          mkexps.
c
c     output:
c
c     lexp2, which contains the local north expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
c
c***********************************************************************
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
      integer  nsall, isall(nsall),iysall(nsall)
      integer  ns34,is34(ns34),iy34(ns34)
      integer  ns134,is134(1),iy134(1)
      integer  is234(1),ns234,iy234(1)
      integer  i, jj
      integer  nnodes
c
      complex *16 lexp2(1:nnodesmax,*)
      complex *16 mexsall(1:nnodesmax)
      complex *16 mexs34(1:nnodesmax)
      complex *16 mexs134(1:nnodesmax)
      complex *16 mexs234(1:nnodesmax)
      complex *16 zs(-3:3,-3:3,1:nnodesmax)
      complex *16 zmul
c
      do i = 1,nsall
         do jj = 1,nnodes
            zmul = zs(2,-iysall(i),jj)
            lexp2(jj,isall(i)) = lexp2(jj,isall(i)) +
     1            mexsall(jj)*zmul
         enddo
      enddo
c
      do i = 1,ns34
         do jj = 1,nnodes
            zmul = zs(1,-iy34(i),jj)
            lexp2(jj,is34(i)) = lexp2(jj,is34(i)) +
     1            mexs34(jj)*zmul
         enddo
      enddo
c
      do i = 1,ns134
         do jj = 1,nnodes
            zmul = zs(2,-iy134(i),jj)
            lexp2(jj,is134(i)) = lexp2(jj,is134(i)) +
     1            mexs134(jj)*zmul
         enddo
      enddo
c
      do i = 1,ns234
         do jj = 1,nnodes
            zmul = zs(2,-iy234(i),jj)
            lexp2(jj,is234(i)) = lexp2(jj,is234(i)) +
     1            mexs234(jj)*zmul
         enddo
      enddo
c
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      subroutine yprocessea(lexp1,ieall,neall,iyeall,
     2           ie13,ne13,iy13,ie1,ne1,iy1,ie3,ne3,iy3,
     5           mexeall,mexe13,mexe1,mexe3,zs,nnodes)
c***********************************************************************
c
c     this subroutine processes the east interaction lists.
c
c     input:
c
c     ieall(neall), iyeall(neall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other east lists are similarly defined (see mkealist).
c
c     mexeall is the exponential expansion for all eight children, etc.
c     zs are the shift coefficients computed by subroutine
c          mkexps.
c
c     output:
c
c     lexp1, which contains the local east expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
c
c***********************************************************************
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
      integer  neall,ieall(neall),iyeall(neall)
      integer  ne13,ie13(ne13),iy13(ne13)
      integer  ne1,ie1(ne1),iy1(ne1)
      integer  ne3,ie3(ne3),iy3(ne3)
      integer  i, jj
      integer  nnodes
c
      complex *16 lexp1(1:nnodesmax,*)
      complex *16 mexeall(1:nnodesmax)
      complex *16 mexe13(1:nnodesmax)
      complex *16 mexe1(1:nnodesmax)
      complex *16 mexe3(1:nnodesmax)
      complex *16 zs(-3:3,-3:3,1:nnodesmax)
      complex *16 zmul
c
      do i = 1,neall
         do jj = 1,nnodes
            zmul = zs(3,-iyeall(i),jj)
            lexp1(jj,ieall(i)) = lexp1(jj,ieall(i)) +
     1            mexeall(jj)*zmul
         enddo
      enddo
c
      do i = 1,ne13
         do jj = 1,nnodes
            zmul = zs(2,-iy13(i),jj)
            lexp1(jj,ie13(i)) = lexp1(jj,ie13(i)) +
     1            mexe13(jj)*zmul
         enddo
      enddo
c
      do i = 1,ne1
         do jj = 1,nnodes
            zmul = zs(-iy1(i)+1,-iy1(i),jj)
            lexp1(jj,ie1(i)) = lexp1(jj,ie1(i)) +
     1            mexe1(jj)*zmul
         enddo
      enddo
c
      do i = 1,ne3
         do jj = 1,nnodes
            zmul = zs(iy3(i),-iy3(i),jj)
            lexp1(jj,ie3(i)) = lexp1(jj,ie3(i)) +
     1            mexe3(jj)*zmul
         enddo
      enddo

c
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      subroutine yprocesswe(lexp2,iwall,nwall,iywall,
     2           iw24,nw24,iy24,iw2,nw2,iy2,iw4,nw4,iy4,
     5           mexwall,mexw24,mexw2,mexw4,zs,nnodes)
c***********************************************************************
c
c     this subroutine processes the west interaction lists.
c
c     input:
c
c     iwall(nwall), iywall(nwall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other west lists are similarly defined (see mkwelist).
c
c     mexeall is the exponential expansion for all eight children, etc.
c     zs are the shift coefficients computed by subroutine
c          mkexps.
c
c     output:
c
c     lexp1, which contains the local west expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
c
c***********************************************************************
      implicit none
      integer  nnodesmax
      parameter (nnodesmax=50)
c
      integer  nwall,iwall(nwall),iywall(nwall)
      integer  nw24,iw24(nw24),iy24(nw24)
      integer  nw2,iw2(nw2),iy2(nw2)
      integer  nw4,iw4(nw4),iy4(nw4)
      integer  i, jj
      integer  nnodes
c
      complex *16 lexp2(1:nnodesmax,*)
      complex *16 mexwall(1:nnodesmax)
      complex *16 mexw24(1:nnodesmax)
      complex *16 mexw2(1:nnodesmax)
      complex *16 mexw4(1:nnodesmax)
      complex *16 zs(-3:3,-3:3,1:nnodesmax)
      complex *16 zmul
c
      do i = 1,nwall
         do jj = 1,nnodes
            zmul = zs(2,-iywall(i),jj)
            lexp2(jj,iwall(i)) = lexp2(jj,iwall(i)) +
     1            mexwall(jj)*zmul
         enddo
      enddo
c
      do i = 1,nw24
         do jj = 1,nnodes
            zmul = zs(1,-iy24(i),jj)
            lexp2(jj,iw24(i)) = lexp2(jj,iw24(i)) +
     1            mexw24(jj)*zmul
         enddo
      enddo
c
      do i = 1,nw2
         do jj = 1,nnodes
            zmul = zs(-iy2(i),-iy2(i),jj)
            lexp2(jj,iw2(i)) = lexp2(jj,iw2(i)) +
     1            mexw2(jj)*zmul
         enddo
      enddo
c
      do i = 1,nw4
         do jj = 1,nnodes
            zmul = zs(iy4(i)-1,-iy4(i),jj)
            lexp2(jj,iw4(i)) = lexp2(jj,iw4(i)) +
     1            mexw4(jj)*zmul
         enddo
      enddo

      return
      end




c*********************************************************************
      subroutine ypwta_coef(nterms,nnodes,comp,comp2)
c*********************************************************************
c
c     this subroutine precomputes the necessary coefficients for the
c     plane wave to local expansions.
c
c     on input :
c       nterms : the number of multipole expansion terms.
c       nnodes : the number of plane wave expansion terms.
c       comp() : the precomputed coefficients from the
c                mp-> pw expansion.
c
c     on output :
c       comp2 : the pw->ta coefficients.
c
c*************************************************************************
      implicit none
c
      integer  nterms, nnodes, nnodesmax
      parameter(nnodesmax=50)
c
c-----global variables.
c
      real *8 comp(nnodesmax, -nterms:nterms)
      real *8 comp2(nnodesmax, -nterms:nterms)
c
c-----local variables.
c
      integer i, k
c
c-------now compute the coefficients comp().
c
      if ( nnodes.eq.1 ) then
        do i=1, nnodes
          do k=1, nterms
            comp2(i,k)=0.0d0
            comp2(i,-k)=0.0d0
          enddo
          comp2(i,0)=0.0d0
        enddo
      else
        do i=1, nnodes
          do k=1, nterms
            comp2(i,k)=comp(i,k)/comp(i,0)
            comp2(i,-k)=comp(i,-k)/comp(i,0)
          enddo
          comp2(i,0)=1.0d0
        enddo
      endif
c
      return
      end



c**********************************************************************
      subroutine polar_convert(z,r,theta)
c     Convert a complex number z to polar form
c
c     Input:
c     z - complex number
c
c     Output:
c     r - absolute value of z
c     theta - argument of z

      implicit none

      complex(8) z
      real *8 r,theta

      real *8 pi

      pi=4.0d0*datan(1.0d0)
      r=abs(z)
c     Most complex numbers fall in this case

      if(real(z) .ne. 0.0d0) then
         theta=atan(imag(z)/real(z))
         if (real(z) .lt. 0.0d0) then
            theta=theta+pi
         endif
      else
         if (imag(z) .gt. 0.0d0) then
            theta=pi/2.0d0
         else
            theta=3.0d0*pi/2.0d0
         endif
      endif
c     Need to worry about cases where z is on one of the axis


c     converts a complex variable z to (r,theta) so that
c     z=r*exp(eye*theta)

      return
      end


c*****************************************************************
      subroutine yformmp_coef(index,nterms,scale,betascal,levnow,wint)
c*****************************************************************
c     this subroutine precomputs the coefficients for forming
c     the multipole expansion for the necessary levels.
c     now everything is based on the table of table concepts.
c
c     on input :
c        index : the index which shows if this is the first time
c                this subroutine is excuted. if index=0, then first
c                time, otherwise, not the first time.
c        nterms : the number of terms in the multipole expansion.
c        scale : the length of the box for different levels.
c        betascale : the scale factor for beta.
c        levnow : the current level, a number from 0 to maxlevel.
c
c     on output :
c        wint : the coefficients for forming the multipole expansion,
c               it is the translation matrix from the polynomial
c               coefficients to multipole coefficients.
c
c     the current program is based on table of tables, it uses computed
c       chebyshev expansion of these translation operators as a
c       function of beta*h (see the notes), and extropolate at
c       other points using the 20 term chebyshev polynomials.
c
c     local parameters:
c       ndeg : the degree of the polynomial approximation.  (3 for 4th order)
c       maxlevel : the max level of the adaptive mesh structure. (17)
c       ntermmax : the max number of multipole expansions. (43/4+1=11)
c
c*******************************************************************
c
      implicit none
c
      integer  ndeg, nterms, maxlevel,maxbeta
      parameter (ndeg=3,maxlevel=17,maxbeta=6)
c
c-----global variables.
c
      integer  index, levnow
      real *8  rjudge(1:maxbeta,0:maxlevel)
      real *8 scale(0:maxlevel), betascal(0:maxlevel)
      real *8 wint(nterms/4+1,10)
      real *8 wintsave(1:maxbeta,0:maxlevel,11,10)
      real *8 www(1:4,1:21,0:10,10)
      save rjudge,wintsave,www
c
c-----local variables.
c
      integer  k,l,i,ind,nlegen,ndiff,indbeta
      real *8 rkc,xhalf,x(0:2*ndeg)
      real *8 dcsevl,rktemp,cheby(1:21)
c
      nlegen=5
      ndiff= nlegen*(nlegen+1)/2
c
      if (index .eq. 0) then
c
c-------the first time this subroutine is excuted.
c
        call formmp_cheby_tab(www)
        do i=0, maxlevel-1
          do k=1,maxbeta-1
            rjudge(k,i)=-1
          enddo
        enddo
        rjudge(1,maxlevel)=-1.0d0
        rjudge(2,maxlevel)=-1.0d0
        rjudge(3,maxlevel)=-1.0d0
        rjudge(4,maxlevel)=-1.0d0
        rjudge(5,maxlevel)=-1.0d0
        rjudge(6,maxlevel)=-1.0d0
        return
      endif
c
c-----now check if current current beta, and current level's coefficients
c     is already computed.
c
      if (dabs(rjudge(1,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=1
      elseif (dabs(rjudge(2,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=2
      elseif (dabs(rjudge(3,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=3
      elseif (dabs(rjudge(4,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=4
      elseif (dabs(rjudge(5,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=5
      elseif (dabs(rjudge(6,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=6
      elseif (rjudge(1,maxlevel) .le. 0) then
        indbeta=1
        rjudge(1,maxlevel)=betascal(0)
      elseif (rjudge(2,maxlevel) .le. 0) then
        indbeta=2
        rjudge(2,maxlevel)=betascal(0)
      elseif (rjudge(3,maxlevel) .le. 0) then
        indbeta=3
        rjudge(3,maxlevel)=betascal(0)
      elseif (rjudge(4,maxlevel) .le. 0) then
        indbeta=4
        rjudge(4,maxlevel)=betascal(0)
      elseif (rjudge(5,maxlevel) .le. 0) then
        indbeta=5
        rjudge(5,maxlevel)=betascal(0)
      elseif (rjudge(6,maxlevel) .le. 0) then
        indbeta=6
        rjudge(6,maxlevel)=betascal(0)
      else
        print *, 'too many beta, stop!', (rjudge(i,maxlevel),i=1,6)
        print *, 'difference', ((rjudge(i,maxlevel)-betascal(0)),i=1,6)
        stop
      endif
      
c
c-----now check if current level's coefficients is already computed.
c
      if (rjudge(indbeta,levnow).gt.0) then
        do l=1, nterms/4+1
          do i=1,10
            wint(l,i)=wintsave(indbeta,levnow, l,i)
          enddo
        enddo
c
      else
        rjudge(indbeta,levnow)=1
        xhalf=scale(levnow)/2.0d0
        x(0)=xhalf*xhalf
        do i=1,3
          x(i)=x(i-1)*xhalf
        enddo
c
c-------rkc is the constant before i_l( rkc*rk(k))
c
        rkc=betascal(0)*scale(levnow)
        do ind=1,nterms/4+1
          do k=1,10
            if (rkc.le. 8.0d0) then
              rktemp=(2.0d0*rkc-8.0d0)/8.0d0
              do i=1,21
                cheby(i)=www(1,i,ind-1,k)
              enddo
            elseif (rkc.le. 16.0d0) then
              rktemp=(2.0d0*rkc-24.0d0)/8.0d0
              do i=1,21
                cheby(i)=www(2,i,ind-1,k)
              enddo
            elseif (rkc.le. 24.0d0) then
              rktemp=(2.0d0*rkc-40.0d0)/8.0d0
              do i=1,21
                cheby(i)=www(3,i,ind-1,k)
              enddo
            elseif (rkc.le. 32.0d0) then
              rktemp=(2.0d0*rkc-56.0d0)/8.0d0
              do i=1,21
                cheby(i)=www(4,i,ind-1,k)
              enddo
            else
              rktemp=0.0d0
              do i=1,21
                cheby(i)=0.0d0
              enddo
            endif
            cheby(1)=2.0d0*cheby(1)
c
            wint(ind,k)=dcsevl(rktemp,cheby,21)
          enddo
c
          wint(ind,1 )=x(0)*wint(ind,1)
          wint(ind,2 )=x(2)*wint(ind,2)
          wint(ind,3 )=x(1)*wint(ind,3 )
          wint(ind,4 )=x(3)*wint(ind,4 )
          wint(ind,5 )=x(3)*wint(ind,5 )
          wint(ind,6 )=x(2)*wint(ind,6 )
          wint(ind,7 )=x(2)*wint(ind,7 )
          wint(ind,8 )=x(1)*wint(ind,8 )
          wint(ind,9 )=x(3)*wint(ind,9 )
          wint(ind,10)=x(3)*wint(ind,10)
        enddo
c
        do l=1, nterms/4+1
          do i=1, 10
            wintsave(indbeta,levnow,l,i)=wint(l,i)
          enddo
        enddo
      endif
c
      return
      end



c**********************************************************************
      subroutine ympmp_coef(index,nterms,scale,betascal,levnow,c)
c**********************************************************************
c
c     this subroutine precomputes the mp->mp shifting coefficients
c     for theta=pi/4, and all other thetas can be derived by taking
c     conjugate/minus operations and are not computed here. but
c     that is done in the subroutine
c
c     on input :
c        index : the index which shows if this is the first time
c                this subroutine is excuted. if index=0, then first
c                time, otherwise, not the first time.
c        nterms : the number of multipole expansion terms.
c        scale : the length of the box for different levels.
c        betascale : the scale factor for beta.
c        levnow : the current level, a number from 0 to maxlevel.
c
c     on output :
c       c : the computed coefficients.
c
c     the coefficients are scaled to avoid over/under flow.
c
c**********************************************************************
c
      implicit none
c
      integer  ntermsmax, maxlevel,maxbeta
      parameter(ntermsmax=42,maxlevel=17,maxbeta=6)
c
c-----global variables.
c
      integer  index, nterms, levnow
      real *8  rjudge(1:maxbeta,0:maxlevel)
      real *8 scale(0:maxlevel), betascal(0:maxlevel)
      complex(8) c(0:nterms, -nterms:nterms)
      complex(8) csave(1:maxbeta,0:maxlevel,0:ntermsmax,
     1                  -ntermsmax:ntermsmax)
      save rjudge,csave
c
c-----local variables.
c
      integer i,m,mabs,lmabs,l,indbeta,k
      real *8 bi(0:2*ntermsmax), power2(0:ntermsmax)
      real *8 rkc, betapower(-ntermsmax*2:ntermsmax*2)
      real *8 betahalf, betah
      real *8 fact(0:2*ntermsmax)
      real *8 pi
      complex *16 imag
      data imag/(0.0d0, 1.0d0)/
      save imag, pi, fact
c
      if (index.eq.0) then
c
c-------the first time this subroutine is called, initialize.
c
        pi=4.0d0*datan(1.0d0)
        fact(0)=1.0d0
        do i=1, 2*ntermsmax
          fact(i)=fact(i-1)*dble(i)
        enddo
c
        do i=0, maxlevel-1
          do k=1,maxbeta-1
            rjudge(k,i)=-1
          enddo
        enddo
        rjudge(1,maxlevel)=-1.0d0
        rjudge(2,maxlevel)=-1.0d0
        rjudge(3,maxlevel)=-1.0d0
        rjudge(4,maxlevel)=-1.0d0
        rjudge(5,maxlevel)=-1.0d0
        rjudge(6,maxlevel)=-1.0d0
        return
      endif
c
c
c-----now check if current current beta, and current level's coefficients
c     is already computed.
c
      if (dabs(rjudge(1,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=1
      elseif (dabs(rjudge(2,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=2
      elseif (dabs(rjudge(3,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=3
      elseif (dabs(rjudge(4,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=4
      elseif (dabs(rjudge(5,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=5
      elseif (dabs(rjudge(6,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=6
      elseif (rjudge(1,maxlevel) .le. 0) then
        indbeta=1
        rjudge(1,maxlevel)=betascal(0)
      elseif (rjudge(2,maxlevel) .le. 0) then
        indbeta=2
        rjudge(2,maxlevel)=betascal(0)
      elseif (rjudge(3,maxlevel) .le. 0) then
        indbeta=3
        rjudge(3,maxlevel)=betascal(0)
      elseif (rjudge(4,maxlevel) .le. 0) then
        indbeta=4
        rjudge(4,maxlevel)=betascal(0)
      elseif (rjudge(5,maxlevel) .le. 0) then
        indbeta=5
        rjudge(5,maxlevel)=betascal(0)
      elseif (rjudge(6,maxlevel) .le. 0) then
        indbeta=6
        rjudge(6,maxlevel)=betascal(0)
      else
        print *, 'too many beta, stop!', (rjudge(i,maxlevel),i=1,6)
        print *, 'difference', ((rjudge(i,maxlevel)-betascal(0)),i=1,6)
        stop
      endif
c
      if (rjudge(indbeta,levnow).le.0 ) then
        rjudge(indbeta,levnow)=1
c
c-------the scaling factor is betascal(0)*scale(levnow), and 2^{-l}
c
        betah=betascal(0)*scale(levnow)
        betahalf= betah/2.0d0
        betapower(0)=1.0d0
        do i=1,2*nterms
          betapower(i)=betapower(i-1)*betahalf
          betapower(-i)=betapower(-i+1)/betahalf
        enddo
c
        power2(0)=1.0d0
        do i=1, nterms
          power2(i)=power2(i-1)/2.0d0
        enddo
c
c-------rkc is the shifting distance, now compute bi(i).
c
        rkc=betah/dsqrt(2.0d0)
        call besseli(nterms*2, rkc, bi, betah )
c
c-------now compute the coefficients.
c
        do l=0, nterms
          do m=-nterms, nterms
            mabs=iabs(m)
            lmabs=iabs(l-m)
            c(l,m)=betapower(-l+mabs+lmabs)*power2(l)*fact(l)
     1             /fact(mabs)/fact(lmabs)*bi(lmabs)*
     2             exp(-imag*dble(l-m)*pi/4.0d0)
            csave(indbeta,levnow,l,m)=c(l,m)
          enddo
        enddo
      else
        do l=0, nterms
          do m=-nterms, nterms
            c(l,m)=csave(indbeta,levnow,l,m)
          enddo
        enddo
      endif
c
      return
      end


c**********************************************************************
      subroutine ytata_coef(index,nterms,scale,betascal,levnow,c)
c**********************************************************************
c
c     this subroutine precomputes the local->local shifting coefficients
c     for theta=pi/4, and all other thetas can be derived by taking
c     conjugate/minus operations and are not computed here. but
c     that is done in the subroutine
c
c     on input :
c        index : the index which shows if this is the first time
c                this subroutine is excuted. if index=0, then first
c                time, otherwise, not the first time.
c        nterms : the number of local expansion terms.
c        scale : the length of the box for different levels.
c        betascal : the scale factor for beta.
c        levnow : the current level, a number from 0 to maxlevel.
c
c     on output :
c       c : the computed coefficients.
c
c     the coefficients are scaled to avoid over/under flow.
c
c     function called :
c       besseli().
c
c**********************************************************************
c
      implicit none
c
      integer  ntermsmax, maxlevel,maxbeta
      parameter(ntermsmax=42,maxlevel=17,maxbeta=6)
c
c-----global variables.
c
      integer  index, nterms, levnow
      real *8  rjudge(1:maxbeta,0:maxlevel)
      real *8 scale(0:maxlevel), betascal(0:maxlevel)
      complex *16 c(0:nterms, -nterms:nterms)
      complex *16 csave(1:maxbeta,0:maxlevel,0:ntermsmax,
     1                  -ntermsmax:ntermsmax)
      save rjudge,csave
c
c-----local variables.
c
      integer i,m,mabs,lmabs,l,indbeta,k
      real *8 bi(0:2*ntermsmax), power2(0:ntermsmax)
      real *8 rkc, betapower(-ntermsmax*2:ntermsmax*2)
      real *8 betahalf,betah
      real *8 fact(0:2*ntermsmax)
      real *8 pi
      complex *16 imag
      data imag/(0.0d0, 1.0d0)/
      save imag, pi, fact
c
      if (index.eq.0) then
        pi=4.0d0*datan(1.0d0)
        fact(0)=1.0d0
        do i=1, 2*ntermsmax
          fact(i)=fact(i-1)*dble(i)
        enddo
c
        do i=0, maxlevel-1
          do k=1,maxbeta-1
            rjudge(k,i)=-1
          enddo
        enddo
        rjudge(1,maxlevel)=-1.0d0
        rjudge(2,maxlevel)=-1.0d0
        rjudge(3,maxlevel)=-1.0d0
        rjudge(4,maxlevel)=-1.0d0
        rjudge(5,maxlevel)=-1.0d0
        rjudge(6,maxlevel)=-1.0d0
        return
      endif
c
c-----now check if current current beta, and current level's coefficients
c     is already computed.
c
      if (dabs(rjudge(1,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=1
      elseif (dabs(rjudge(2,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=2
      elseif (dabs(rjudge(3,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=3
      elseif (dabs(rjudge(4,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=4
      elseif (dabs(rjudge(5,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=5
      elseif (dabs(rjudge(6,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=6
      elseif (rjudge(1,maxlevel) .le. 0) then
        indbeta=1
        rjudge(1,maxlevel)=betascal(0)
      elseif (rjudge(2,maxlevel) .le. 0) then
        indbeta=2
        rjudge(2,maxlevel)=betascal(0)
      elseif (rjudge(3,maxlevel) .le. 0) then
        indbeta=3
        rjudge(3,maxlevel)=betascal(0)
      elseif (rjudge(4,maxlevel) .le. 0) then
        indbeta=4
        rjudge(4,maxlevel)=betascal(0)
      elseif (rjudge(5,maxlevel) .le. 0) then
        indbeta=5
        rjudge(5,maxlevel)=betascal(0)
      elseif (rjudge(6,maxlevel) .le. 0) then
        indbeta=6
        rjudge(6,maxlevel)=betascal(0)
      else
        print *, 'too many beta, stop!', (rjudge(i,maxlevel),i=1,6)
        print *, 'difference', ((rjudge(i,maxlevel)-betascal(0)),i=1,6)
        stop
      endif
c
      if (rjudge(indbeta,levnow).le.0 ) then
        rjudge(indbeta,levnow)=1
c
c-------the scaling factor is betascal(0)*scale(levnow), and 2^{-l}
c
        betah=betascal(0)*scale(levnow)
c
c-------check if no local expansion needed.
c
        if (betah .le. 36.0d0) then
          betahalf= betah /2.0d0
          betapower(0)=1.0d0
          do i=1,2*nterms
            betapower(i)=betapower(i-1)*betahalf
            betapower(-i)=betapower(-i+1)/betahalf
          enddo
c
          power2(0)=1.0d0
          do i=1, nterms
            power2(i)=power2(i-1)/2.0d0
          enddo
c
c-------rkc is the shifting distance, now compute bi(i).
c
          rkc=betascal(0)*scale(levnow)*dsqrt(2.0d0)/4.0d0
          call besseli(nterms*2, rkc, bi, betah )
c
c-------now compute the coefficients.
c
          do l=0, nterms
            do m=-nterms, nterms
              mabs=iabs(m)
              lmabs=iabs(l-m)
              c(l,m)=betapower(+l-mabs+lmabs)*power2(l)/fact(l)
     1               *fact(mabs)/fact(lmabs)*bi(lmabs)*
     2                 exp(imag*dble(m-l)*pi/4.0d0)
              csave(indbeta,levnow,l,m)=c(l,m)
            enddo
          enddo
        else
          do l=0, nterms
            do m=-nterms, nterms
              c(l,m)=0.0d0
              csave(indbeta,levnow,l,m)=c(l,m)
            enddo
          enddo
        endif
      else
        do l=0, nterms
          do m=-nterms, nterms
            c(l,m)=csave(indbeta,levnow,l,m)
          enddo
        enddo
      endif
c
      return
      end
c
c**********************************************************************
      subroutine ymppw_coef(index,nterms,iprec,nnodes,scale,betascal,
     1                     levnow,xnodes, comp)
c**********************************************************************
c
c     this subroutine precomputes the necessary coefficients for the
c     multipole to plane wave and plane wave to local expansions.
c
c     on input :
c        index : the index which shows if this is the first time
c                this subroutine is excuted. if index=0, then first
c                time, otherwise, not the first time.
c        nterms : the number of multipole expansion terms.
c        scale : the length of the box for different levels.
c        betascal : the scale factor for beta.
c        levnow : the current level, a number from 0 to maxlevel.
c
c     on output :
c       wnodes : the weights for the gauss quadrature.
c       xnodes : the nodes for the quadrature.
c       comp : the coefficients so w(k)=\sum comp(l)*a_l, l=-nterms, nterms.
c
c     subroutine called : yukq2d() : generates the optimal quadratures
c       for different beta value. note the difference of the weights
c       defined in the subroutine and the weights defined in our
c       program.
c
c*************************************************************************
      implicit none
c
      integer  ntermsmax, maxlevel, nnodesmax, lwork,maxbeta
      parameter(ntermsmax=42,maxlevel=17,nnodesmax=50,lwork=100000)
      parameter(maxbeta=6)
c
c-----global variables.
c
      integer  index, nterms, nnodes, levnow
      integer  lnodes(1:maxbeta,0:maxlevel)
      real *8  rjudge(1:maxbeta,0:maxlevel)
      integer  ier, iprec, lused
      real *8 scale(0:maxlevel), betascal(0:maxlevel)
      real *8 xnodes(nnodesmax)
      real *8 comp(nnodesmax, -nterms:nterms)
      real *8 xnodesave(1:maxbeta,0:maxlevel,nnodesmax)
      real *8 compsave(1:maxbeta,0:maxlevel,nnodesmax,
     1                 -ntermsmax:ntermsmax)
      real *8 work(lwork)
      complex *16 zweights(nnodesmax)
      save rjudge,lnodes,xnodesave,compsave
c
c-----local variables.
c
      integer i,k,indbeta
      real *8 wnodes(nnodesmax)
      real *8 betahalf,betabox,beta0, betatemp
      real *8 temp0,pi
c

      pi=3.141592653589793d0
      if (index .eq. 0 ) then
        do i=0, maxlevel-1
          do k=1,maxbeta-1
            rjudge(k,i)=-1
          enddo
        enddo
        rjudge(1,maxlevel)=-1.0d0
        rjudge(2,maxlevel)=-1.0d0
        rjudge(3,maxlevel)=-1.0d0
        rjudge(4,maxlevel)=-1.0d0
        rjudge(5,maxlevel)=-1.0d0
        rjudge(6,maxlevel)=-1.0d0
        return
      endif
c
c-----now check if current current beta, and current level's coefficients
c     is already computed.
c
      if (dabs(rjudge(1,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=1
      elseif (dabs(rjudge(2,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=2
      elseif (dabs(rjudge(3,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=3
      elseif (dabs(rjudge(4,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=4
      elseif (dabs(rjudge(5,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=5
      elseif (dabs(rjudge(6,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=6
      elseif (rjudge(1,maxlevel) .le. 0) then
        indbeta=1
        rjudge(1,maxlevel)=betascal(0)
      elseif (rjudge(2,maxlevel) .le. 0) then
        indbeta=2
        rjudge(2,maxlevel)=betascal(0)
      elseif (rjudge(3,maxlevel) .le. 0) then
        indbeta=3
        rjudge(3,maxlevel)=betascal(0)
      elseif (rjudge(4,maxlevel) .le. 0) then
        indbeta=4
        rjudge(4,maxlevel)=betascal(0)
      elseif (rjudge(5,maxlevel) .le. 0) then
        indbeta=5
        rjudge(5,maxlevel)=betascal(0)
      elseif (rjudge(6,maxlevel) .le. 0) then
        indbeta=6
        rjudge(6,maxlevel)=betascal(0)
      else
        print *, 'too many beta, stop!', (rjudge(i,maxlevel),i=1,6)
        print *, 'difference', ((rjudge(i,maxlevel)-betascal(0)),i=1,6)
        stop
      endif
c
      if (rjudge(indbeta,levnow).le.0 ) then
        rjudge(indbeta,levnow)=1
c
c-------now compute the weights and nodes.
c       betabox : the scaled beta for the unit box.
c
        betabox=betascal(0)*scale(levnow)
        if (betabox .ge. 8.0d0*pi) then
          nnodes=1
          xnodes(1)=0.0d0
          zweights(1)=0.0d0
        else
c          call yukq2d(ier, betabox, iprec, nnodes, xnodes,zweights,
c     1                work, lwork, lused)
          call yukq2d(ier, betabox, nnodes, xnodes,zweights,
     1                work, lwork, lused)
        endif
c
c-------now rescale the weights and nodes.
c
        lnodes(indbeta,levnow)=nnodes
c
        if (nnodes .gt. 49) then
          print *, 'too many nodes needed, stoped.'
          stop
        endif
c
        do i=1, nnodes
          wnodes(i)=dble(zweights(i))*pi/dsqrt(xnodes(i)*xnodes(i)+
     1              betabox*betabox)/2.0d0
          xnodesave(indbeta,levnow,i)=xnodes(i)
        enddo
c
c-------now compute the coefficients comp().
c
        betahalf=betabox/2.0d0
        do i=1, nnodes
          betatemp=(dsqrt(betabox*betabox+xnodes(i)*xnodes(i))+
     1             xnodes(i))/betabox
          temp0=betahalf*betatemp
          beta0=betahalf/betatemp
c
          comp(i,0)=wnodes(i)
          compsave(indbeta,levnow,i,0)=comp(i,0)
c
          do k=1, nterms
            comp(i,k)=comp(i,k-1)*temp0/dble(k)
            compsave(indbeta,levnow,i,k)=comp(i,k)
            comp(i,-k)=comp(i,-k+1)*beta0/dble(k)
            compsave(indbeta,levnow,i,-k)=comp(i,-k)
c
          enddo
        enddo
c
      else
        nnodes=lnodes(indbeta,levnow)
        do i=1, nnodes
          xnodes(i)=xnodesave(indbeta,levnow,i)
        enddo
c
        do i=1, nnodes
          do k=-nterms,nterms
            comp(i,k)=compsave(indbeta,levnow, i,k)
          enddo
        enddo
      endif
c
      return
      end
c
c**********************************************************************
      subroutine ypwpw_coef(index,nnodes,scale,betascal,levnow,
     1        xnodes,zs)
c**********************************************************************
c
c     this subroutine computes the tables of exponentials needed
c     for translating exponential representations of harmonic
c     functions, discretized via norman's quadratures.
c
c     on input:
c
c     xlength : the length of the current box.
c     beta : the frequency in the equation \delta u - beta^2 u =f
c     betascal : the scaling factor for beta.
c     xnodes(nnodes)   discretization points in lambda integral
c     nnodes          number of discret. pts. in lambda integral
c
c     on output:
c
c     zs(n,m,k)   e^{- dsqrt(lambda_k^2+beta^2)*n}*e^{im lambda_k)}.
c
c**********************************************************************
      implicit none
      integer  nnodesmax, maxlevel,maxbeta
      parameter (nnodesmax=50,maxlevel=17,maxbeta=6)
c
      integer  nnodes
      integer  index, levnow
      real *8  rjudge(1:maxbeta,0:maxlevel)
      real *8 scale(0:maxlevel), betascal(0:maxlevel)
      real *8 xnodes(nnodesmax)
      complex *16 zs(-3:3,-3:3,1:nnodesmax)
      complex *16 zsave(1:maxbeta,0:maxlevel,-3:3,-3:3,1:nnodesmax)
      save rjudge,zsave
c
      integer   i,k,m,n,indbeta
      real *8 betabox, alpha
      complex *16 imag
      data imag/(0.0d0, 1.0d0)/
c
      if (index .eq. 0 ) then
        do i=0, maxlevel-1
          do k=1,maxbeta-1
            rjudge(k,i)=-1
          enddo
        enddo
        rjudge(1,maxlevel)=-1.0d0
        rjudge(2,maxlevel)=-1.0d0
        rjudge(3,maxlevel)=-1.0d0
        rjudge(4,maxlevel)=-1.0d0
        rjudge(5,maxlevel)=-1.0d0
        rjudge(6,maxlevel)=-1.0d0
        return
      endif
c
c-----now check if current current beta, and current level's coefficients
c     is already computed.
c
      if (dabs(rjudge(1,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=1
      elseif (dabs(rjudge(2,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=2
      elseif (dabs(rjudge(3,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=3
      elseif (dabs(rjudge(4,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=4
      elseif (dabs(rjudge(5,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=5
      elseif (dabs(rjudge(6,maxlevel)-betascal(0)).le.1.0e-14) then
        indbeta=6
      elseif (rjudge(1,maxlevel) .le. 0) then
        indbeta=1
        rjudge(1,maxlevel)=betascal(0)
      elseif (rjudge(2,maxlevel) .le. 0) then
        indbeta=2
        rjudge(2,maxlevel)=betascal(0)
      elseif (rjudge(3,maxlevel) .le. 0) then
        indbeta=3
        rjudge(3,maxlevel)=betascal(0)
      elseif (rjudge(4,maxlevel) .le. 0) then
        indbeta=4
        rjudge(4,maxlevel)=betascal(0)
      elseif (rjudge(5,maxlevel) .le. 0) then
        indbeta=5
        rjudge(5,maxlevel)=betascal(0)
      elseif (rjudge(6,maxlevel) .le. 0) then
        indbeta=6
        rjudge(6,maxlevel)=betascal(0)
      else
        print *, 'too many beta, stop!', (rjudge(i,maxlevel),i=1,6)
        print *, 'difference', ((rjudge(i,maxlevel)-betascal(0)),i=1,6)
        stop
      endif
c
      if (rjudge(indbeta,levnow).le.0 ) then
        rjudge(indbeta,levnow)=1
c
        betabox=scale(levnow)*betascal(0)
c
c-------loop over each lambda value
c
        if (betabox .le. 36.0d0) then
          do k = 1,nnodes
             alpha=dsqrt(xnodes(k)*xnodes(k)+betabox*betabox)
             do n = -3,3
                do m = -3,3
                   zs(n,m,k)=dexp(-alpha*dble(n))
     1                       *exp(-imag*xnodes(k)*dble(m))
                   zsave(indbeta,levnow,n,m,k)=zs(n,m,k)
                enddo
             enddo
          enddo
        else
          do k = 1,nnodes
             do n = -3,3
                do m = -3,3
                   zs(n,m,k)=0.0d0
                   zsave(indbeta,levnow,n,m,k)=zs(n,m,k)
                enddo
             enddo
          enddo
        endif
c
      else
        do k = 1,nnodes
          do n = -3,3
            do m = -3,3
              zs(n,m,k)=zsave(indbeta,levnow,n,m,k)
            enddo
          enddo
        enddo
      endif
c
      return
      end



c**********************************************************************
      subroutine besselk(nterms, x, bk, betascal)
c**********************************************************************
c
c     this subroutine computes the scaled bessel function
c     of first kind k_l(x)*(betascal/2)^l/gamma(l+1).
c     the reason for the scaling is to avoid under/over flow.
c
c     on input :
c       nterms, x : bk(0:nterms)(x) is given.
c       betascal : the scaling factor.
c
c     on output :
c       bk(0:nterms) : the bessel i function.
c
c     we will do the following in this subroutine :
c       1. for x < 1.0d-4, we use the two term approximation.
c       2. for x > 1.0d-4, we call the subroutine ribesl()
c            and then rescale the result.
c**********************************************************************
        implicit none
c
c-------global parameters.
c
        integer nterms
        real *8 x, bk(0:nterms), betascal
c
c-------local variables.
c
        integer ncalc, ize, ntermss, i
        parameter(ntermss=42)
        real *8 ensig, enmten
        real *8 alpha
        real *8 const
        real *8 fact(0:2*ntermss)
        data ensig, enmten/1.0d-4, 1.0d-300/
      data fact/1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     1          720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     2          3628800.0d0, 39916800.0d0, 479001600.0d0,
     3          6227020800.0d0, 87178291200.0d0, 1307674368000.0d0,
     4          20922789888000.0d0, 3.556874280960000d14,
     5          6.402373705728000d15, 1.216451004088320d17,
     6          2.432902008176640d18, 5.109094217170944d19,
     7          1.124000727777608d21, 2.585201673888498d22,
     8          6.204484017332394d23, 1.551121004333099d25,
     9          4.032914611266057d26, 1.088886945041835d28,
     a          3.048883446117138d29, 8.841761993739701d30,
     b          2.652528598121910d32, .822283865417792d+34,
     c          .263130836933694d+36, .868331761881189d+37,
     d          .295232799039604d+39, .103331479663861d+41,
     e          .371993326789901d+42, .137637530912263d+44,
     f          .523022617466601d+45, .203978820811974d+47,
     g          .815915283247898d+48, .334525266131638d+50,
     *          .140500611775288d+52, .604152630633738d+53,
     *          .265827157478845d+55, .119622220865480d+57,
     *          .550262215981209d+58, .258623241511168d+60,
     *          .124139155925361d+62, .608281864034268d+63,
     *          .304140932017134d+65, .155111875328738d+67,
     *          .806581751709439d+68, .427488328406003d+70,
     *          .230843697339241d+72, .126964033536583d+74,
     *          .710998587804863d+75, .405269195048772d+77,
     *          .235056133128288d+79, .138683118545690d+81,
     *          .832098711274139d+82, .507580213877225d+84,
     *          .314699732603879d+86, .198260831540444d+88,
     *          .126886932185884d+90, .824765059208247d+91,
     *          .544344939077443d+93, .364711109181887d+95,
     *          .248003554243683d+97, .171122452428141d+99,
     *          .119785716699699d101, .850478588567862d102,
     *          .612344583768861d104, .447011546151268d106,
     *          .330788544151939d108, .248091408113954d110,
     *          .188549470166605d112, .145183092028286d114,
     *          .113242811782063d116, .894618213078297d117,
     *          .715694570462638d119, .579712602074737d121,
     *          .475364333701284d123, .394552396972066d125,
     *          .331424013456535d127/
c
        if (x.lt. 0.0d0) then
          print *,' error in input for besseli()'
          stop
        elseif (x.gt. 700.0d0) then
          do i=0,nterms
            bk(i)=0.0d0
          enddo
          return
        endif
c
c-------we will call rkbesl() and then scale the bks'
c       by the scaling factor.
c
        const=1.0d0
c
c-------calculate k_(n).
c
        alpha=0.0d0
        ize=1
        call rkbesl(x, alpha, nterms+1, ize, bk, ncalc)
c
c---------now scale to i_n*scal**(-n)*fact(n).
c
        if (ncalc .lt. 0) then
          print *, 'error in besselk, stoped'
          stop
        endif
c
        do 2 i=0,ncalc-1
          bk(i)=bk(i)*const/fact(i)
          const=const*betascal/2.0d0
          if (dabs(bk(i)) .le. enmten ) then
            const=0.0d0
          endif
2       continue
c
        do i=ncalc, nterms
          bk(i)=bk(i-1)*bk(i)*betascal/2.0d0*fact(i-1)/fact(i)
        enddo
c
        return
        end

c**********************************************************************
      subroutine besseli(nterms, x, bi, betascal)
c**********************************************************************
c
c     this subroutine computes the scaled bessel function
c     of second kind i_l(x)/(betascal/2)^l*gamma(l+1).
c     the reason for the scaling is to avoid under/over flow.
c
c     on input :
c       nterms, x : bi(0:nterms)(x) is given.
c       betascal : the scaling factor.
c
c     on output :
c       bi(0:nterms) : the bessel i function.
c
c     we will do the following in this subroutine :
c       1. for x < 1.0d-4, we use the two term approximation.
c       2. for x > 1.0d-4, we call the subroutine ribesl()
c            and then rescale the result.
c**********************************************************************
        implicit none
c
c-------global parameters.
c
        integer nterms
        real *8 x, bi(0:nterms), betascal
c
c-------local variables.
c
        integer i,ncalc, ize, ntermss
        parameter(ntermss=42)
        real *8 ensig, enmten
        real *8 term1, term2, term3, term4, xscal, alpha
        real *8 const, dble
        real *8 fact(0:2*ntermss)
        data ensig, enmten/1.0d-4, 1.0d-300/
      data fact/1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     1          720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     2          3628800.0d0, 39916800.0d0, 479001600.0d0,
     3          6227020800.0d0, 87178291200.0d0, 1307674368000.0d0,
     4          20922789888000.0d0, 3.556874280960000d14,
     5          6.402373705728000d15, 1.216451004088320d17,
     6          2.432902008176640d18, 5.109094217170944d19,
     7          1.124000727777608d21, 2.585201673888498d22,
     8          6.204484017332394d23, 1.551121004333099d25,
     9          4.032914611266057d26, 1.088886945041835d28,
     a          3.048883446117138d29, 8.841761993739701d30,
     b          2.652528598121910d32, .822283865417792d+34,
     c          .263130836933694d+36, .868331761881189d+37,
     d          .295232799039604d+39, .103331479663861d+41,
     e          .371993326789901d+42, .137637530912263d+44,
     f          .523022617466601d+45, .203978820811974d+47,
     g          .815915283247898d+48, .334525266131638d+50,
     *          .140500611775288d+52, .604152630633738d+53,
     *          .265827157478845d+55, .119622220865480d+57,
     *          .550262215981209d+58, .258623241511168d+60,
     *          .124139155925361d+62, .608281864034268d+63,
     *          .304140932017134d+65, .155111875328738d+67,
     *          .806581751709439d+68, .427488328406003d+70,
     *          .230843697339241d+72, .126964033536583d+74,
     *          .710998587804863d+75, .405269195048772d+77,
     *          .235056133128288d+79, .138683118545690d+81,
     *          .832098711274139d+82, .507580213877225d+84,
     *          .314699732603879d+86, .198260831540444d+88,
     *          .126886932185884d+90, .824765059208247d+91,
     *          .544344939077443d+93, .364711109181887d+95,
     *          .248003554243683d+97, .171122452428141d+99,
     *          .119785716699699d101, .850478588567862d102,
     *          .612344583768861d104, .447011546151268d106,
     *          .330788544151939d108, .248091408113954d110,
     *          .188549470166605d112, .145183092028286d114,
     *          .113242811782063d116, .894618213078297d117,
     *          .715694570462638d119, .579712602074737d121,
     *          .475364333701284d123, .394552396972066d125,
     *          .331424013456535d127/
c
        if (x.lt. 0.0d0) then
          print *,' error in input for besseli()'
        endif
c
c-------if x.le. 10e-4, then use the 2 terms taylor expansion.
c
        if (x.le.ensig) then
          xscal=x/betascal
          term1=1.0d0
          term2=0.25d0*x*x
          bi(0)=term1*(1.0d0+term2)
          do 1 i=1, nterms
            term1=term1*xscal
            if (term1.le. enmten) then
              term1=0.0d0
            endif
            bi(i)=term1*(1.0d0+term2/dble(i+1) )
1         continue
c
c-------if   0.02 > x > 0.0001, then use four term approximation
c       to avoid over/underflow.
c
        elseif (x. le. 2.0d-2) then
          xscal=x/betascal
          term1=1.0d0
          term2=0.25d0*x*x
          term3=term2*term2/2.0d0
          term4=term2*term3/3.0d0
          bi(0)=term1*(1.0d0+term2+term3/2.0d0+term4/6.0d0)
          do i=1, nterms
            term1=term1*xscal
            if (term1.le. enmten) then
              term1=0.0d0
            endif
            bi(i)=term1*(1.0d0+(term2+(term3+term4/dble(i+3))
     1              /dble(i+2))/dble(i+1))
          enddo
c
c---------usually, scal should be less than one, however, in the
c         calculation of h_n, scal will be greater than one,
c         in this case, we should modify the program and replace every
c         small item with zero.
c
        else
c
c---------otherwise, we will call ribesl() and then scale the
c         in by the scaling factor.
c
          const=1.0d0
c
c---------calculate i_(n).
c
          alpha=0.0d0
          ize=1
          call ribesl(x, alpha, nterms+1, ize, bi, ncalc)
c
c---------now scale to i_n*scal**(-n)*fact(n).
c
          do 2 i=0,ncalc
            bi(i)=bi(i)/const*fact(i)
            const=const*betascal/2.0d0
2         continue
c
          do i=ncalc+1, nterms
            bi(i)=0.0d0
          enddo
        endif
        return
        end

c**********************************************************************
      subroutine inittable
c**********************************************************************
c     this subroutine initializes the required tables.
c
c**********************************************************************
c
      implicit none
      integer  nlev,nterms,iprec,nnodes
      real *8 betascal(0:10),scale(0:10)
      real *8 wint(5000)
      real *8 c(5000),xnodes(5000),comp(5000)
      real *8 zs(5000)
c 
      nlev=3
      nterms=21
c
      call yformmp_coef(0,nterms,scale,betascal,nlev,wint)
      call ympmp_coef(0,nterms,scale,betascal,nlev,c)
      call ymppw_coef(0,nterms,iprec,nnodes,scale,betascal,nlev,
     1               xnodes,comp)
      call ypwpw_coef(0,nnodes,scale,betascal,nlev,xnodes,zs)
      call ytata_coef(0,nterms,scale,betascal,nlev,c)

      return
      end






c**********************************************************************
      subroutine inf_norm(z,d)
      complex *16 z
      real*8 d

      d=max(dabs(dreal(z)),dabs(dimag(z)))

      return
      end



c***************************************************************
*deck dcsevl
      double precision function dcsevl (x, cs, n)
c***begin prologue  dcsevl
c***purpose  evaluate a chebyshev series.
c***library   slatec (fnlib)
c***category  c3a2
c***type      double precision (csevl-s, dcsevl-d)
c***keywords  chebyshev series, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  evaluate the n-term chebyshev series cs at x.  adapted from
c  a method presented in the paper by broucke referenced below.
c
c       input arguments --
c  x    value at which the series is to be evaluated.
c  cs   array of n terms of a chebyshev series.  in evaluating
c       cs, only half the first coefficient is summed.
c  n    number of terms in array cs.
c
c***references  r. broucke, ten subroutines for the manipulation of
c                 chebyshev series, algorithm 446, communications of
c                 the a.c.m. 16, (1973) pp. 254-256.
c               l. fox and i. b. parker, chebyshev polynomials in
c                 numerical analysis, oxford university press, 1968,
c                 page 56.
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900329  prologued revised extensively and code rewritten to allow
c           x to be slightly outside interval (-1,+1).  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dcsevl
      double precision b0, b1, b2, cs(*), onepl, twox, x, d1mach
      logical first
      save first, onepl
      data first /.true./
c***first executable statement  dcsevl
      if (first) onepl = 1.0d0 + d1mach(4)
      first = .false.
      if (n .lt. 1) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .le. 0', 2, 2)
      if (n .gt. 1000) call xermsg ('slatec', 'dcsevl',
     +   'number of terms .gt. 1000', 3, 2)
      if (dabs(x-1.0d0) .le. 1.0e-15) then
          x=1.0d0
      elseif (dabs(x+1.0d0) .le. 1.0e-15 ) then
          x=-1.0d0
      endif
c
      if (abs(x) .gt. onepl) call xermsg ('slatec', 'dcsevl',
     +   'x outside the interval (-1,+1)', 1, 1)
c
      b1 = 0.0d0
      b0 = 0.0d0
      twox = 2.0d0*x
      do 10 i = 1,n
         b2 = b1
         b1 = b0
         ni = n + 1 - i
         b0 = twox*b1 - b2 + cs(ni)
   10 continue
c
      dcsevl = 0.5d0*(b0-b2)
c
      return
      end
*deck d1mach
      double precision function d1mach (i)
c***begin prologue  d1mach
c***purpose  return floating point machine dependent constants.
c***library   slatec
c***category  r1
c***type      double precision (r1mach-s, d1mach-d)
c***keywords  machine constants
c***author  fox, p. a., (bell labs)
c           hall, a. d., (bell labs)
c           schryer, n. l., (bell labs)
c***description
c
c   d1mach can be used to obtain machine-dependent parameters for the
c   local machine environment.  it is a function subprogram with one
c   (input) argument, and can be referenced as follows:
c
c        d = d1mach(i)
c
c   where i=1,...,5.  the (output) value of d above is determined by
c   the (input) value of i.  the results for various values of i are
c   discussed below.
c
c   d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c   d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c   d1mach( 3) = b**(-t), the smallest relative spacing.
c   d1mach( 4) = b**(1-t), the largest relative spacing.
c   d1mach( 5) = log10(b)
c
c   assume double precision numbers are represented in the t-digit,
c   base-b form
c
c              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
c   emin .le. e .le. emax.
c
c   the values of b, t, emin and emax are provided in i1mach as
c   follows:
c   i1mach(10) = b, the base.
c   i1mach(14) = t, the number of base-b digits.
c   i1mach(15) = emin, the smallest exponent e.
c   i1mach(16) = emax, the largest exponent e.
c
c   to alter this function for a particular environment, the desired
c   set of data statements should be activated by removing the c from
c   column 1.  also, the values of d1mach(1) - d1mach(4) should be
c   checked for consistency with the local operating system.
c
c***references  p. a. fox, a. d. hall and n. l. schryer, framework for
c                 a portable library, acm transactions on mathematical
c                 software 4, 2 (june 1978), pp. 177-188.
c***routines called  xermsg
c***revision history  (yymmdd)
c   750101  date written
c   890213  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900618  added dec risc constants.  (wrb)
c   900723  added ibm rs 6000 constants.  (wrb)
c   900911  added sun 386i constants.  (wrb)
c   910710  added hp 730 constants.  (smr)
c   911114  added convex ieee constants.  (wrb)
c   920121  added sun -r8 compiler option constants.  (wrb)
c   920229  added touchstone delta i860 constants.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   920625  added convex -p8 and -pd8 compiler option constants.
c           (bks, wrb)
c   930201  added dec alpha and sgi constants.  (rwc and wrb)
c***end prologue  d1mach
c
      integer small(4)
      integer large(4)
      integer right(4)
      integer diver(4)
      integer log10(4)
c
      double precision dmach(5)
      save dmach
c
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
c
c     machine constants for the amiga
c     absoft fortran compiler using the 68020/68881 compiler option
c
c     data small(1), small(2) / z'00100000', z'00000000' /
c     data large(1), large(2) / z'7fefffff', z'ffffffff' /
c     data right(1), right(2) / z'3ca00000', z'00000000' /
c     data diver(1), diver(2) / z'3cb00000', z'00000000' /
c     data log10(1), log10(2) / z'3fd34413', z'509f79ff' /
c
c     machine constants for the amiga
c     absoft fortran compiler using software floating point
c
c     data small(1), small(2) / z'00100000', z'00000000' /
c     data large(1), large(2) / z'7fdfffff', z'ffffffff' /
c     data right(1), right(2) / z'3ca00000', z'00000000' /
c     data diver(1), diver(2) / z'3cb00000', z'00000000' /
c     data log10(1), log10(2) / z'3fd34413', z'509f79ff' /
c
c     machine constants for the apollo
c
c     data small(1), small(2) / 16#00100000, 16#00000000 /
c     data large(1), large(2) / 16#7fffffff, 16#ffffffff /
c     data right(1), right(2) / 16#3ca00000, 16#00000000 /
c     data diver(1), diver(2) / 16#3cb00000, 16#00000000 /
c     data log10(1), log10(2) / 16#3fd34413, 16#509f79ff /
c
c     machine constants for the burroughs 1700 system
c
c     data small(1) / zc00800000 /
c     data small(2) / z000000000 /
c     data large(1) / zdffffffff /
c     data large(2) / zfffffffff /
c     data right(1) / zcc5800000 /
c     data right(2) / z000000000 /
c     data diver(1) / zcc6800000 /
c     data diver(2) / z000000000 /
c     data log10(1) / zd00e730e7 /
c     data log10(2) / zc77800dc0 /
c
c     machine constants for the burroughs 5700 system
c
c     data small(1) / o1771000000000000 /
c     data small(2) / o0000000000000000 /
c     data large(1) / o0777777777777777 /
c     data large(2) / o0007777777777777 /
c     data right(1) / o1461000000000000 /
c     data right(2) / o0000000000000000 /
c     data diver(1) / o1451000000000000 /
c     data diver(2) / o0000000000000000 /
c     data log10(1) / o1157163034761674 /
c     data log10(2) / o0006677466732724 /
c
c     machine constants for the burroughs 6700/7700 systems
c
c     data small(1) / o1771000000000000 /
c     data small(2) / o7770000000000000 /
c     data large(1) / o0777777777777777 /
c     data large(2) / o7777777777777777 /
c     data right(1) / o1461000000000000 /
c     data right(2) / o0000000000000000 /
c     data diver(1) / o1451000000000000 /
c     data diver(2) / o0000000000000000 /
c     data log10(1) / o1157163034761674 /
c     data log10(2) / o0006677466732724 /
c
c     machine constants for the cdc 170/180 series using nos/ve
c
c     data small(1) / z"3001800000000000" /
c     data small(2) / z"3001000000000000" /
c     data large(1) / z"4ffefffffffffffe" /
c     data large(2) / z"4ffe000000000000" /
c     data right(1) / z"3fd2800000000000" /
c     data right(2) / z"3fd2000000000000" /
c     data diver(1) / z"3fd3800000000000" /
c     data diver(2) / z"3fd3000000000000" /
c     data log10(1) / z"3fff9a209a84fbcf" /
c     data log10(2) / z"3ffff7988f8959ac" /
c
c     machine constants for the cdc 6000/7000 series
c
c     data small(1) / 00564000000000000000b /
c     data small(2) / 00000000000000000000b /
c     data large(1) / 37757777777777777777b /
c     data large(2) / 37157777777777777777b /
c     data right(1) / 15624000000000000000b /
c     data right(2) / 00000000000000000000b /
c     data diver(1) / 15634000000000000000b /
c     data diver(2) / 00000000000000000000b /
c     data log10(1) / 17164642023241175717b /
c     data log10(2) / 16367571421742254654b /
c
c     machine constants for the celerity c1260
c
c     data small(1), small(2) / z'00100000', z'00000000' /
c     data large(1), large(2) / z'7fefffff', z'ffffffff' /
c     data right(1), right(2) / z'3ca00000', z'00000000' /
c     data diver(1), diver(2) / z'3cb00000', z'00000000' /
c     data log10(1), log10(2) / z'3fd34413', z'509f79ff' /
c
c     machine constants for the convex
c     using the -fn or -pd8 compiler option
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fffffffffffffff' /
c     data dmach(3) / z'3cc0000000000000' /
c     data dmach(4) / z'3cd0000000000000' /
c     data dmach(5) / z'3ff34413509f79ff' /
c
c     machine constants for the convex
c     using the -fi compiler option
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fefffffffffffff' /
c     data dmach(3) / z'3ca0000000000000' /
c     data dmach(4) / z'3cb0000000000000' /
c     data dmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the convex
c     using the -p8 compiler option
c
c     data dmach(1) / z'00010000000000000000000000000000' /
c     data dmach(2) / z'7fffffffffffffffffffffffffffffff' /
c     data dmach(3) / z'3f900000000000000000000000000000' /
c     data dmach(4) / z'3f910000000000000000000000000000' /
c     data dmach(5) / z'3fff34413509f79fef311f12b35816f9' /
c
c     machine constants for the cray
c
c     data small(1) / 201354000000000000000b /
c     data small(2) / 000000000000000000000b /
c     data large(1) / 577767777777777777777b /
c     data large(2) / 000007777777777777774b /
c     data right(1) / 376434000000000000000b /
c     data right(2) / 000000000000000000000b /
c     data diver(1) / 376444000000000000000b /
c     data diver(2) / 000000000000000000000b /
c     data log10(1) / 377774642023241175717b /
c     data log10(2) / 000007571421742254654b /
c
c     machine constants for the data general eclipse s/200
c     note - it may be appropriate to include the following card -
c     static dmach(5)
c
c     data small /    20k, 3*0 /
c     data large / 77777k, 3*177777k /
c     data right / 31420k, 3*0 /
c     data diver / 32020k, 3*0 /
c     data log10 / 40423k, 42023k, 50237k, 74776k /
c
c     machine constants for the dec alpha
c     using g_float
c
c     data dmach(1) / '0000000000000010'x /
c     data dmach(2) / 'ffffffffffff7fff'x /
c     data dmach(3) / '0000000000003cc0'x /
c     data dmach(4) / '0000000000003cd0'x /
c     data dmach(5) / '79ff509f44133ff3'x /
c
c     machine constants for the dec alpha
c     using ieee_format
c
c     data dmach(1) / '0010000000000000'x /
c     data dmach(2) / '7fefffffffffffff'x /
c     data dmach(3) / '3ca0000000000000'x /
c     data dmach(4) / '3cb0000000000000'x /
c     data dmach(5) / '3fd34413509f79ff'x /
c
c     machine constants for the dec risc
c
c     data small(1), small(2) / z'00000000', z'00100000'/
c     data large(1), large(2) / z'ffffffff', z'7fefffff'/
c     data right(1), right(2) / z'00000000', z'3ca00000'/
c     data diver(1), diver(2) / z'00000000', z'3cb00000'/
c     data log10(1), log10(2) / z'509f79ff', z'3fd34413'/
c
c     machine constants for the dec vax
c     using d_floating
c     (expressed in integer and hexadecimal)
c     the hex format below may not be suitable for unix systems
c     the integer format should be ok for unix systems
c
c     data small(1), small(2) /        128,           0 /
c     data large(1), large(2) /     -32769,          -1 /
c     data right(1), right(2) /       9344,           0 /
c     data diver(1), diver(2) /       9472,           0 /
c     data log10(1), log10(2) /  546979738,  -805796613 /
c
c     data small(1), small(2) / z00000080, z00000000 /
c     data large(1), large(2) / zffff7fff, zffffffff /
c     data right(1), right(2) / z00002480, z00000000 /
c     data diver(1), diver(2) / z00002500, z00000000 /
c     data log10(1), log10(2) / z209a3f9a, zcff884fb /
c
c     machine constants for the dec vax
c     using g_floating
c     (expressed in integer and hexadecimal)
c     the hex format below may not be suitable for unix systems
c     the integer format should be ok for unix systems
c
c     data small(1), small(2) /         16,           0 /
c     data large(1), large(2) /     -32769,          -1 /
c     data right(1), right(2) /      15552,           0 /
c     data diver(1), diver(2) /      15568,           0 /
c     data log10(1), log10(2) /  1142112243, 2046775455 /
c
c     data small(1), small(2) / z00000010, z00000000 /
c     data large(1), large(2) / zffff7fff, zffffffff /
c     data right(1), right(2) / z00003cc0, z00000000 /
c     data diver(1), diver(2) / z00003cd0, z00000000 /
c     data log10(1), log10(2) / z44133ff3, z79ff509f /
c
c     machine constants for the elxsi 6400
c     (assuming real*8 is the default double precision)
c
c     data small(1), small(2) / '00100000'x,'00000000'x /
c     data large(1), large(2) / '7fefffff'x,'ffffffff'x /
c     data right(1), right(2) / '3cb00000'x,'00000000'x /
c     data diver(1), diver(2) / '3cc00000'x,'00000000'x /
c     data log10(1), log10(2) / '3fd34413'x,'509f79ff'x /
c
c     machine constants for the harris 220
c
c     data small(1), small(2) / '20000000, '00000201 /
c     data large(1), large(2) / '37777777, '37777577 /
c     data right(1), right(2) / '20000000, '00000333 /
c     data diver(1), diver(2) / '20000000, '00000334 /
c     data log10(1), log10(2) / '23210115, '10237777 /
c
c     machine constants for the honeywell 600/6000 series
c
c     data small(1), small(2) / o402400000000, o000000000000 /
c     data large(1), large(2) / o376777777777, o777777777777 /
c     data right(1), right(2) / o604400000000, o000000000000 /
c     data diver(1), diver(2) / o606400000000, o000000000000 /
c     data log10(1), log10(2) / o776464202324, o117571775714 /
c
c     machine constants for the hp 730
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fefffffffffffff' /
c     data dmach(3) / z'3ca0000000000000' /
c     data dmach(4) / z'3cb0000000000000' /
c     data dmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the hp 2100
c     three word double precision option with ftn4
c
c     data small(1), small(2), small(3) / 40000b,       0,       1 /
c     data large(1), large(2), large(3) / 77777b, 177777b, 177776b /
c     data right(1), right(2), right(3) / 40000b,       0,    265b /
c     data diver(1), diver(2), diver(3) / 40000b,       0,    276b /
c     data log10(1), log10(2), log10(3) / 46420b,  46502b,  77777b /
c
c     machine constants for the hp 2100
c     four word double precision option with ftn4
c
c     data small(1), small(2) /  40000b,       0 /
c     data small(3), small(4) /       0,       1 /
c     data large(1), large(2) /  77777b, 177777b /
c     data large(3), large(4) / 177777b, 177776b /
c     data right(1), right(2) /  40000b,       0 /
c     data right(3), right(4) /       0,    225b /
c     data diver(1), diver(2) /  40000b,       0 /
c     data diver(3), diver(4) /       0,    227b /
c     data log10(1), log10(2) /  46420b,  46502b /
c     data log10(3), log10(4) /  76747b, 176377b /
c
c     machine constants for the hp 9000
c
c     data small(1), small(2) / 00040000000b, 00000000000b /
c     data large(1), large(2) / 17737777777b, 37777777777b /
c     data right(1), right(2) / 07454000000b, 00000000000b /
c     data diver(1), diver(2) / 07460000000b, 00000000000b /
c     data log10(1), log10(2) / 07764642023b, 12047674777b /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9, the sel systems 85/86, and
c     the perkin elmer (interdata) 7/32.
c
c     data small(1), small(2) / z00100000, z00000000 /
c     data large(1), large(2) / z7fffffff, zffffffff /
c     data right(1), right(2) / z33100000, z00000000 /
c     data diver(1), diver(2) / z34100000, z00000000 /
c     data log10(1), log10(2) / z41134413, z509f79ff /
c
c     machine constants for the ibm pc
c     assumes that all arithmetic is done in double precision
c     on 8088, i.e., not in 80 bit form for the 8087.
c
c     data small(1) / 2.23d-308  /
c     data large(1) / 1.79d+308  /
c     data right(1) / 1.11d-16   /
c     data diver(1) / 2.22d-16   /
c     data log10(1) / 0.301029995663981195d0 /
c
c     machine constants for the ibm rs 6000
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fefffffffffffff' /
c     data dmach(3) / z'3ca0000000000000' /
c     data dmach(4) / z'3cb0000000000000' /
c     data dmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the intel i860
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fefffffffffffff' /
c     data dmach(3) / z'3ca0000000000000' /
c     data dmach(4) / z'3cb0000000000000' /
c     data dmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the pdp-10 (ka processor)
c
c     data small(1), small(2) / "033400000000, "000000000000 /
c     data large(1), large(2) / "377777777777, "344777777777 /
c     data right(1), right(2) / "113400000000, "000000000000 /
c     data diver(1), diver(2) / "114400000000, "000000000000 /
c     data log10(1), log10(2) / "177464202324, "144117571776 /
c
c     machine constants for the pdp-10 (ki processor)
c
c     data small(1), small(2) / "000400000000, "000000000000 /
c     data large(1), large(2) / "377777777777, "377777777777 /
c     data right(1), right(2) / "103400000000, "000000000000 /
c     data diver(1), diver(2) / "104400000000, "000000000000 /
c     data log10(1), log10(2) / "177464202324, "476747767461 /
c
c     machine constants for pdp-11 fortran supporting
c     32-bit integers (expressed in integer and octal).
c
c     data small(1), small(2) /    8388608,           0 /
c     data large(1), large(2) / 2147483647,          -1 /
c     data right(1), right(2) /  612368384,           0 /
c     data diver(1), diver(2) /  620756992,           0 /
c     data log10(1), log10(2) / 1067065498, -2063872008 /
c
c     data small(1), small(2) / o00040000000, o00000000000 /
c     data large(1), large(2) / o17777777777, o37777777777 /
c     data right(1), right(2) / o04440000000, o00000000000 /
c     data diver(1), diver(2) / o04500000000, o00000000000 /
c     data log10(1), log10(2) / o07746420232, o20476747770 /
c
c     machine constants for pdp-11 fortran supporting
c     16-bit integers (expressed in integer and octal).
c
c     data small(1), small(2) /    128,      0 /
c     data small(3), small(4) /      0,      0 /
c     data large(1), large(2) /  32767,     -1 /
c     data large(3), large(4) /     -1,     -1 /
c     data right(1), right(2) /   9344,      0 /
c     data right(3), right(4) /      0,      0 /
c     data diver(1), diver(2) /   9472,      0 /
c     data diver(3), diver(4) /      0,      0 /
c     data log10(1), log10(2) /  16282,   8346 /
c     data log10(3), log10(4) / -31493, -12296 /
c
c     data small(1), small(2) / o000200, o000000 /
c     data small(3), small(4) / o000000, o000000 /
c     data large(1), large(2) / o077777, o177777 /
c     data large(3), large(4) / o177777, o177777 /
c     data right(1), right(2) / o022200, o000000 /
c     data right(3), right(4) / o000000, o000000 /
c     data diver(1), diver(2) / o022400, o000000 /
c     data diver(3), diver(4) / o000000, o000000 /
c     data log10(1), log10(2) / o037632, o020232 /
c     data log10(3), log10(4) / o102373, o147770 /
c
c     machine constants for the silicon graphics
c
c     data small(1), small(2) / z'00100000', z'00000000' /
c     data large(1), large(2) / z'7fefffff', z'ffffffff' /
c     data right(1), right(2) / z'3ca00000', z'00000000' /
c     data diver(1), diver(2) / z'3cb00000', z'00000000' /
c     data log10(1), log10(2) / z'3fd34413', z'509f79ff' /
c
c     machine constants for the sun
c
c     data dmach(1) / z'0010000000000000' /
c     data dmach(2) / z'7fefffffffffffff' /
c     data dmach(3) / z'3ca0000000000000' /
c     data dmach(4) / z'3cb0000000000000' /
c     data dmach(5) / z'3fd34413509f79ff' /
c
c     machine constants for the sun
c     using the -r8 compiler option
c
c     data dmach(1) / z'00010000000000000000000000000000' /
c     data dmach(2) / z'7ffeffffffffffffffffffffffffffff' /
c     data dmach(3) / z'3f8e0000000000000000000000000000' /
c     data dmach(4) / z'3f8f0000000000000000000000000000' /
c     data dmach(5) / z'3ffd34413509f79fef311f12b35816f9' /
c
c     machine constants for the sun 386i
c
c     data small(1), small(2) / z'fffffffd', z'000fffff' /
c     data large(1), large(2) / z'ffffffb0', z'7fefffff' /
c     data right(1), right(2) / z'000000b0', z'3ca00000' /
c     data diver(1), diver(2) / z'ffffffcb', z'3cafffff'
c     data log10(1), log10(2) / z'509f79e9', z'3fd34413' /
c
c     machine constants for the univac 1100 series ftn compiler
c
c     data small(1), small(2) / o000040000000, o000000000000 /
c     data large(1), large(2) / o377777777777, o777777777777 /
c     data right(1), right(2) / o170540000000, o000000000000 /
c     data diver(1), diver(2) / o170640000000, o000000000000 /
c     data log10(1), log10(2) / o177746420232, o411757177572 /
c
c***first executable statement  d1mach
      if (i .lt. 1 .or. i .gt. 5) call xermsg ('slatec', 'd1mach',
     +   'i out of bounds', 1, 2)
c
      d1mach = dmach(i)
      return
c
      end
*deck fdump
      subroutine fdump
c***begin prologue  fdump
c***purpose  symbolic dump (should be locally written).
c***library   slatec (xerror)
c***category  r3
c***type      all (fdump-a)
c***keywords  error, xermsg
c***author  jones, r. e., (snla)
c***description
c
c        ***note*** machine dependent routine
c        fdump is intended to be replaced by a locally written
c        version which produces a symbolic dump.  failing this,
c        it should be replaced by a version which prints the
c        subprogram nesting list.  note that this dump must be
c        printed on each of up to five files, as indicated by the
c        xgetua routine.  see xsetua and xgetua for details.
c
c     written by ron jones, with slatec common math library subcommittee
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  fdump
c***first executable statement  fdump
      return
      end
*deck i1mach
      integer function i1mach (i)
c***begin prologue  i1mach
c***purpose  return integer machine dependent constants.
c***library   slatec
c***category  r1
c***type      integer (i1mach-i)
c***keywords  machine constants
c***author  fox, p. a., (bell labs)
c           hall, a. d., (bell labs)
c           schryer, n. l., (bell labs)
c***description
c
c   i1mach can be used to obtain machine-dependent parameters for the
c   local machine environment.  it is a function subprogram with one
c   (input) argument and can be referenced as follows:
c
c        k = i1mach(i)
c
c   where i=1,...,16.  the (output) value of k above is determined by
c   the (input) value of i.  the results for various values of i are
c   discussed below.
c
c   i/o unit numbers:
c     i1mach( 1) = the standard input unit.
c     i1mach( 2) = the standard output unit.
c     i1mach( 3) = the standard punch unit.
c     i1mach( 4) = the standard error message unit.
c
c   words:
c     i1mach( 5) = the number of bits per integer storage unit.
c     i1mach( 6) = the number of characters per integer storage unit.
c
c   integers:
c     assume integers are represented in the s-digit, base-a form
c
c                sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c                where 0 .le. x(i) .lt. a for i=0,...,s-1.
c     i1mach( 7) = a, the base.
c     i1mach( 8) = s, the number of base-a digits.
c     i1mach( 9) = a**s - 1, the largest magnitude.
c
c   floating-point numbers:
c     assume floating-point numbers are represented in the t-digit,
c     base-b form
c                sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c                where 0 .le. x(i) .lt. b for i=1,...,t,
c                0 .lt. x(1), and emin .le. e .le. emax.
c     i1mach(10) = b, the base.
c
c   single-precision:
c     i1mach(11) = t, the number of base-b digits.
c     i1mach(12) = emin, the smallest exponent e.
c     i1mach(13) = emax, the largest exponent e.
c
c   double-precision:
c     i1mach(14) = t, the number of base-b digits.
c     i1mach(15) = emin, the smallest exponent e.
c     i1mach(16) = emax, the largest exponent e.
c
c   to alter this function for a particular environment, the desired
c   set of data statements should be activated by removing the c from
c   column 1.  also, the values of i1mach(1) - i1mach(4) should be
c   checked for consistency with the local operating system.
c
c***references  p. a. fox, a. d. hall and n. l. schryer, framework for
c                 a portable library, acm transactions on mathematical
c                 software 4, 2 (june 1978), pp. 177-188.
c***routines called  (none)
c***revision history  (yymmdd)
c   750101  date written
c   891012  added vax g-floating constants.  (wrb)
c   891012  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900618  added dec risc constants.  (wrb)
c   900723  added ibm rs 6000 constants.  (wrb)
c   901009  correct i1mach(7) for ibm mainframes. should be 2 not 16.
c           (rwc)
c   910710  added hp 730 constants.  (smr)
c   911114  added convex ieee constants.  (wrb)
c   920121  added sun -r8 compiler option constants.  (wrb)
c   920229  added touchstone delta i860 constants.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   920625  added convex -p8 and -pd8 compiler option constants.
c           (bks, wrb)
c   930201  added dec alpha and sgi constants.  (rwc and wrb)
c   930618  corrected i1mach(5) for convex -p8 and -pd8 compiler
c           options.  (dwl, rwc and wrb).
c***end prologue  i1mach
c
      integer imach(16),output
      save imach
      equivalence (imach(4),output)
c
c     machine constants for the amiga
c     absoft compiler
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the apollo
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        129 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1025 /
c
c     machine constants for the burroughs 1700 system
c
c     data imach( 1) /          7 /
c     data imach( 2) /          2 /
c     data imach( 3) /          2 /
c     data imach( 4) /          2 /
c     data imach( 5) /         36 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         33 /
c     data imach( 9) / z1ffffffff /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -256 /
c     data imach(13) /        255 /
c     data imach(14) /         60 /
c     data imach(15) /       -256 /
c     data imach(16) /        255 /
c
c     machine constants for the burroughs 5700 system
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         48 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /          8 /
c     data imach(11) /         13 /
c     data imach(12) /        -50 /
c     data imach(13) /         76 /
c     data imach(14) /         26 /
c     data imach(15) /        -50 /
c     data imach(16) /         76 /
c
c     machine constants for the burroughs 6700/7700 systems
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         48 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /          8 /
c     data imach(11) /         13 /
c     data imach(12) /        -50 /
c     data imach(13) /         76 /
c     data imach(14) /         26 /
c     data imach(15) /     -32754 /
c     data imach(16) /      32780 /
c
c     machine constants for the cdc 170/180 series using nos/ve
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -4095 /
c     data imach(13) /       4094 /
c     data imach(14) /         94 /
c     data imach(15) /      -4095 /
c     data imach(16) /       4094 /
c
c     machine constants for the cdc 6000/7000 series
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /    6loutput/
c     data imach( 5) /         60 /
c     data imach( 6) /         10 /
c     data imach( 7) /          2 /
c     data imach( 8) /         48 /
c     data imach( 9) / 00007777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /       -929 /
c     data imach(13) /       1070 /
c     data imach(14) /         94 /
c     data imach(15) /       -929 /
c     data imach(16) /       1069 /
c
c     machine constants for the celerity c1260
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / z'7fffffff' /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the convex
c     using the -fn compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the convex
c     using the -fi compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the convex
c     using the -p8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1023 /
c     data imach(13) /       1023 /
c     data imach(14) /        113 /
c     data imach(15) /     -16383 /
c     data imach(16) /      16383 /
c
c     machine constants for the convex
c     using the -pd8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         64 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 9223372036854775807 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1023 /
c     data imach(13) /       1023 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the cray
c     using the 46 bit integer compiler option
c
c     data imach( 1) /        100 /
c     data imach( 2) /        101 /
c     data imach( 3) /        102 /
c     data imach( 4) /        101 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         46 /
c     data imach( 9) / 1777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -8189 /
c     data imach(13) /       8190 /
c     data imach(14) /         94 /
c     data imach(15) /      -8099 /
c     data imach(16) /       8190 /
c
c     machine constants for the cray
c     using the 64 bit integer compiler option
c
c     data imach( 1) /        100 /
c     data imach( 2) /        101 /
c     data imach( 3) /        102 /
c     data imach( 4) /        101 /
c     data imach( 5) /         64 /
c     data imach( 6) /          8 /
c     data imach( 7) /          2 /
c     data imach( 8) /         63 /
c     data imach( 9) / 777777777777777777777b /
c     data imach(10) /          2 /
c     data imach(11) /         47 /
c     data imach(12) /      -8189 /
c     data imach(13) /       8190 /
c     data imach(14) /         94 /
c     data imach(15) /      -8099 /
c     data imach(16) /       8190 /
c
c     machine constants for the data general eclipse s/200
c
c     data imach( 1) /         11 /
c     data imach( 2) /         12 /
c     data imach( 3) /          8 /
c     data imach( 4) /         10 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /         16 /
c     data imach(11) /          6 /
c     data imach(12) /        -64 /
c     data imach(13) /         63 /
c     data imach(14) /         14 /
c     data imach(15) /        -64 /
c     data imach(16) /         63 /
c
c     machine constants for the dec alpha
c     using g_float
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the dec alpha
c     using ieee_float
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the dec risc
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the dec vax
c     using d_floating
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the dec vax
c     using g_floating
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1023 /
c     data imach(16) /       1023 /
c
c     machine constants for the elxsi 6400
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         32 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1022 /
c     data imach(16) /       1023 /
c
c     machine constants for the harris 220
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          0 /
c     data imach( 4) /          6 /
c     data imach( 5) /         24 /
c     data imach( 6) /          3 /
c     data imach( 7) /          2 /
c     data imach( 8) /         23 /
c     data imach( 9) /    8388607 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         38 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the honeywell 600/6000 series
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /         43 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          6 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         63 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 730
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the hp 2100
c     3 word double precision option with ftn4
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          4 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         39 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 2100
c     4 word double precision option with ftn4
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          4 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         23 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         55 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for the hp 9000
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          7 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         32 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -126 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1015 /
c     data imach(16) /       1017 /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9, the sel systems 85/86, and
c     the perkin elmer (interdata) 7/32.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          7 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) /  z7fffffff /
c     data imach(10) /         16 /
c     data imach(11) /          6 /
c     data imach(12) /        -64 /
c     data imach(13) /         63 /
c     data imach(14) /         14 /
c     data imach(15) /        -64 /
c     data imach(16) /         63 /
c
c     machine constants for the ibm pc
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          0 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        127 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1023 /
c
c     machine constants for the ibm rs 6000
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          0 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the intel i860
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the pdp-10 (ka processor)
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          5 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         54 /
c     data imach(15) /       -101 /
c     data imach(16) /        127 /
c
c     machine constants for the pdp-10 (ki processor)
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          5 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         62 /
c     data imach(15) /       -128 /
c     data imach(16) /        127 /
c
c     machine constants for pdp-11 fortran supporting
c     32-bit integer arithmetic.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for pdp-11 fortran supporting
c     16-bit integer arithmetic.
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          5 /
c     data imach( 4) /          6 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c     machine constants for the silicon graphics
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the sun
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -125 /
c     data imach(13) /        128 /
c     data imach(14) /         53 /
c     data imach(15) /      -1021 /
c     data imach(16) /       1024 /
c
c     machine constants for the sun
c     using the -r8 compiler option
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          6 /
c     data imach( 4) /          6 /
c     data imach( 5) /         32 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         31 /
c     data imach( 9) / 2147483647 /
c     data imach(10) /          2 /
c     data imach(11) /         53 /
c     data imach(12) /      -1021 /
c     data imach(13) /       1024 /
c     data imach(14) /        113 /
c     data imach(15) /     -16381 /
c     data imach(16) /      16384 /
c
c     machine constants for the univac 1100 series ftn compiler
c
c     data imach( 1) /          5 /
c     data imach( 2) /          6 /
c     data imach( 3) /          1 /
c     data imach( 4) /          6 /
c     data imach( 5) /         36 /
c     data imach( 6) /          4 /
c     data imach( 7) /          2 /
c     data imach( 8) /         35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /          2 /
c     data imach(11) /         27 /
c     data imach(12) /       -128 /
c     data imach(13) /        127 /
c     data imach(14) /         60 /
c     data imach(15) /      -1024 /
c     data imach(16) /       1023 /
c
c     machine constants for the z80 microprocessor
c
c     data imach( 1) /          1 /
c     data imach( 2) /          1 /
c     data imach( 3) /          0 /
c     data imach( 4) /          1 /
c     data imach( 5) /         16 /
c     data imach( 6) /          2 /
c     data imach( 7) /          2 /
c     data imach( 8) /         15 /
c     data imach( 9) /      32767 /
c     data imach(10) /          2 /
c     data imach(11) /         24 /
c     data imach(12) /       -127 /
c     data imach(13) /        127 /
c     data imach(14) /         56 /
c     data imach(15) /       -127 /
c     data imach(16) /        127 /
c
c***first executable statement  i1mach
      if (i .lt. 1  .or.  i .gt. 16) go to 10
c
      i1mach = imach(i)
      return
c
   10 continue
      write (unit = output, fmt = 9000)
 9000 format ('1error    1 in i1mach - i out of bounds')
c
c     call fdump
c
      stop
      end
*deck j4save
      function j4save (iwhich, ivalue, iset)
c***begin prologue  j4save
c***subsidiary
c***purpose  save or recall global variables needed by error
c            handling routines.
c***library   slatec (xerror)
c***type      integer (j4save-i)
c***keywords  error messages, error number, recall, save, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                = 5 refers to the total number of units to which
c                     each error message is to be written.
c                = 6 refers to the 2nd unit for error messages
c                = 7 refers to the 3rd unit for error messages
c                = 8 refers to the 4th unit for error messages
c                = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c***see also  xermsg
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900205  minor modifications to prologue.  (wrb)
c   900402  added type section.  (wrb)
c   910411  added keywords section.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  j4save
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
*deck xercnt
      subroutine xercnt (librar, subrou, messg, nerr, level, kontrl)
c***begin prologue  xercnt
c***subsidiary
c***purpose  allow user control over handling of errors.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xercnt-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xercnt.
c        if the user has provided his own version of xercnt, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        librar - the library that the routine is in.
c        subrou - the subroutine that xermsg is being called from
c        messg  - the first 20 characters of the error message.
c        nerr   - same as in the call to xermsg.
c        level  - same as in the call to xermsg.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900206  routine changed from user-callable to subsidiary.  (wrb)
c   900510  changed calling sequence to include library and subroutine
c           names, changed routine name from xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xercnt
      character*(*) librar, subrou, messg
c***first executable statement  xercnt
      return
      end
*deck xerhlt
      subroutine xerhlt (messg)
c***begin prologue  xerhlt
c***subsidiary
c***purpose  abort program execution and print error message.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerhlt-a)
c***keywords  abort program execution, error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        ***note*** machine dependent routine
c        xerhlt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg is as in xermsg.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  (none)
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900206  routine changed from user-callable to subsidiary.  (wrb)
c   900510  changed calling sequence to delete length of character
c           and changed routine name from xerabt to xerhlt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerhlt
      character*(*) messg
c***first executable statement  xerhlt
      stop
      end
*deck xermsg
      subroutine xermsg (librar, subrou, messg, nerr, level)
c***begin prologue  xermsg
c***purpose  process error messages for slatec and other libraries.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xermsg-a)
c***keywords  error message, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c   xermsg processes a diagnostic message in a manner determined by the
c   value of level and the current value of the library error control
c   flag, kontrl.  see subroutine xsetf for details.
c
c    librar   a character constant (or character variable) with the name
c             of the library.  this will be 'slatec' for the slatec
c             common math library.  the error handling package is
c             general enough to be used by many libraries
c             simultaneously, so it is desirable for the routine that
c             detects and reports an error to identify the library name
c             as well as the routine name.
c
c    subrou   a character constant (or character variable) with the name
c             of the routine that detected the error.  usually it is the
c             name of the routine that is calling xermsg.  there are
c             some instances where a user callable library routine calls
c             lower level subsidiary routines where the error is
c             detected.  in such cases it may be more informative to
c             supply the name of the routine the user called rather than
c             the name of the subsidiary routine that detected the
c             error.
c
c    messg    a character constant (or character variable) with the text
c             of the error or warning message.  in the example below,
c             the message is a character constant that contains a
c             generic message.
c
c                   call xermsg ('slatec', 'mmpy',
c                  *'the order of the matrix exceeds the row dimension',
c                  *3, 1)
c
c             it is possible (and is sometimes desirable) to generate a
c             specific message--e.g., one that contains actual numeric
c             values.  specific numeric values can be converted into
c             character strings using formatted write statements into
c             character variables.  this is called standard fortran
c             internal file i/o and is exemplified in the first three
c             lines of the following example.  you can also catenate
c             substrings of characters to construct the error message.
c             here is an example showing the use of both writing to
c             an internal file and catenating character strings.
c
c                   character*5 charn, charl
c                   write (charn,10) n
c                   write (charl,10) lda
c                10 format(i5)
c                   call xermsg ('slatec', 'mmpy', 'the order'//charn//
c                  *   ' of the matrix exceeds its row dimension of'//
c                  *   charl, 3, 1)
c
c             there are two subtleties worth mentioning.  one is that
c             the // for character catenation is used to construct the
c             error message so that no single character constant is
c             continued to the next line.  this avoids confusion as to
c             whether there are trailing blanks at the end of the line.
c             the second is that by catenating the parts of the message
c             as an actual argument rather than encoding the entire
c             message into one large character variable, we avoid
c             having to know how long the message will be in order to
c             declare an adequate length for that large character
c             variable.  xermsg calls xerprn to print the message using
c             multiple lines if necessary.  if the message is very long,
c             xerprn will break it into pieces of 72 characters (as
c             requested by xermsg) for printing on multiple lines.
c             also, xermsg asks xerprn to prefix each line with ' *  '
c             so that the total line length could be 76 characters.
c             note also that xerprn scans the error message backwards
c             to ignore trailing blanks.  another feature is that
c             the substring '$$' is treated as a new line sentinel
c             by xerprn.  if you want to construct a multiline
c             message without having to count out multiples of 72
c             characters, just use '$$' as a separator.  '$$'
c             obviously must occur within 72 characters of the
c             start of each line to have its intended effect since
c             xerprn is asked to wrap around at 72 characters in
c             addition to looking for '$$'.
c
c    nerr     an integer value that is chosen by the library routine's
c             author.  it must be in the range -99 to 999 (three
c             printable digits).  each distinct error should have its
c             own error number.  these error numbers should be described
c             in the machine readable documentation for the routine.
c             the error numbers need be unique only within each routine,
c             so it is reasonable for each routine to start enumerating
c             errors from 1 and proceeding to the next integer.
c
c    level    an integer value in the range 0 to 2 that indicates the
c             level (severity) of the error.  their meanings are
c
c            -1  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.  an attempt is made to only print this
c                message once.
c
c             0  a warning message.  this is used if it is not clear
c                that there really is an error, but the user's attention
c                may be needed.
c
c             1  a recoverable error.  this is used even if the error is
c                so serious that the routine cannot return any useful
c                answer.  if the user has told the error package to
c                return after recoverable errors, then xermsg will
c                return to the library routine which can then return to
c                the user's routine.  the user may also permit the error
c                package to terminate the program upon encountering a
c                recoverable error.
c
c             2  a fatal error.  xermsg will not return to its caller
c                after it receives a fatal error.  this level should
c                hardly ever be used; it is much better to allow the
c                user a chance to recover.  an example of one of the few
c                cases in which it is permissible to declare a level 2
c                error is a reverse communication library routine that
c                is likely to be called repeatedly until it integrates
c                across some interval.  if there is a serious error in
c                the input such that another step cannot be taken and
c                the library routine is called again without the input
c                error having been corrected by the caller, the library
c                routine will probably be called forever with improper
c                input.  in this case, it is reasonable to declare the
c                error to be fatal.
c
c    each of the arguments to xermsg is input; none will be modified by
c    xermsg.  a routine may make multiple calls to xermsg with warning
c    level messages; however, after a call to xermsg with a recoverable
c    error, the routine should return to the user.  do not try to call
c    xermsg with a second recoverable error after the first recoverable
c    error because the error package saves the error number.  the user
c    can retrieve this error number by calling another entry point in
c    the error handling package and then clear the error number when
c    recovering from the error.  calling xermsg in succession causes the
c    old error number to be overwritten by the latest error number.
c    this is considered harmless for error numbers associated with
c    warning messages but must not be done for error numbers of serious
c    errors.  after a call to xermsg with a recoverable error, the user
c    must be given a chance to call numxer or xerclr to retrieve or
c    clear the error number.
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  fdump, j4save, xercnt, xerhlt, xerprn, xersve
c***revision history  (yymmdd)
c   880101  date written
c   880621  revised as directed at slatec cml meeting of february 1988.
c           there are two basic changes.
c           1.  a new routine, xerprn, is used instead of xerprt to
c               print messages.  this routine will break long messages
c               into pieces for printing on multiple lines.  '$$' is
c               accepted as a new line sentinel.  a prefix can be
c               added to each line to be printed.  xermsg uses either
c               ' ***' or ' *  ' and long messages are broken every
c               72 characters (at most) so that the maximum line
c               length output can now be as great as 76.
c           2.  the text of all messages is now in upper case since the
c               fortran standard document does not admit the existence
c               of lower case.
c   880708  revised after the slatec cml meeting of june 29 and 30.
c           the principal changes are
c           1.  clarify comments in the prologues
c           2.  rename xrprnt to xerprn
c           3.  rework handling of '$$' in xerprn to handle blank lines
c               similar to the way format statements handle the /
c               character for new records.
c   890706  revised with the help of fred fritsch and reg clemens to
c           clean up the coding.
c   890721  revised to use new feature in xerprn to count characters in
c           prefix.
c   891013  revised to correct comments.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  changed test on nerr to be -9999999 < nerr < 99999999, but
c           nerr .ne. 0, and on level to be -2 < level < 3.  added
c           level=-1 logic, changed calls to xersav to xersve, and
c           xerctl to xercnt.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xermsg
      character*(*) librar, subrou, messg
      character*8 xlibr, xsubr
      character*72  temp
      character*20  lfirst
c***first executable statement  xermsg
      lkntrl = j4save (2, 0, .false.)
      maxmes = j4save (4, 0, .false.)
c
c       lkntrl is a local copy of the control flag kontrl.
c       maxmes is the maximum number of times any particular message
c          should be printed.
c
c       we print a fatal error message and terminate for an error in
c          calling xermsg.  the error number should be positive,
c          and the level should be between 0 and 2.
c
      if (nerr.lt.-9999999 .or. nerr.gt.99999999 .or. nerr.eq.0 .or.
     *   level.lt.-1 .or. level.gt.2) then
         call xerprn (' ***', -1, 'fatal error in...$$ ' //
     *      'xermsg -- invalid error number or level$$ '//
     *      'job abort due to fatal error.', 72)
         call xersve (' ', ' ', ' ', 0, 0, 0, kdummy)
         call xerhlt (' ***xermsg -- invalid input')
         return
      endif
c
c       record the message.
c
      i = j4save (1, nerr, .true.)
      call xersve (librar, subrou, messg, 1, nerr, level, kount)
c
c       handle print-once warning messages.
c
      if (level.eq.-1 .and. kount.gt.1) return
c
c       allow temporary user override of the control flag.
c
      xlibr  = librar
      xsubr  = subrou
      lfirst = messg
      lerr   = nerr
      llevel = level
      call xercnt (xlibr, xsubr, lfirst, lerr, llevel, lkntrl)
c
      lkntrl = max(-2, min(2,lkntrl))
      mkntrl = abs(lkntrl)
c
c       skip printing if the control flag value as reset in xercnt is
c       zero and the error is not fatal.
c
      if (level.lt.2 .and. lkntrl.eq.0) go to 30
      if (level.eq.0 .and. kount.gt.maxmes) go to 30
      if (level.eq.1 .and. kount.gt.maxmes .and. mkntrl.eq.1) go to 30
      if (level.eq.2 .and. kount.gt.max(1,maxmes)) go to 30
c
c       announce the names of the library and subroutine by building a
c       message in character variable temp (not exceeding 66 characters)
c       and sending it out via xerprn.  print only if control flag
c       is not zero.
c
      if (lkntrl .ne. 0) then
         temp(1:21) = 'message from routine '
         i = min(len(subrou), 16)
         temp(22:21+i) = subrou(1:i)
         temp(22+i:33+i) = ' in library '
         ltemp = 33 + i
         i = min(len(librar), 16)
         temp(ltemp+1:ltemp+i) = librar (1:i)
         temp(ltemp+i+1:ltemp+i+1) = '.'
         ltemp = ltemp + i + 1
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       if lkntrl is positive, print an introductory line before
c       printing the message.  the introductory line tells the choice
c       from each of the following three options.
c       1.  level of the message
c              'informative message'
c              'potentially recoverable error'
c              'fatal error'
c       2.  whether control flag will allow program to continue
c              'prog continues'
c              'prog aborted'
c       3.  whether or not a traceback was requested.  (the traceback
c           may not be implemented at some sites, so this only tells
c           what was requested, not what was delivered.)
c              'traceback requested'
c              'traceback not requested'
c       notice that the line including four prefix characters will not
c       exceed 74 characters.
c       we skip the next block if the introductory line is not needed.
c
      if (lkntrl .gt. 0) then
c
c       the first part of the message tells about the level.
c
         if (level .le. 0) then
            temp(1:20) = 'informative message,'
            ltemp = 20
         elseif (level .eq. 1) then
            temp(1:30) = 'potentially recoverable error,'
            ltemp = 30
         else
            temp(1:12) = 'fatal error,'
            ltemp = 12
         endif
c
c       then whether the program will continue.
c
         if ((mkntrl.eq.2 .and. level.ge.1) .or.
     *       (mkntrl.eq.1 .and. level.eq.2)) then
            temp(ltemp+1:ltemp+14) = ' prog aborted,'
            ltemp = ltemp + 14
         else
            temp(ltemp+1:ltemp+16) = ' prog continues,'
            ltemp = ltemp + 16
         endif
c
c       finally tell whether there should be a traceback.
c
         if (lkntrl .gt. 0) then
            temp(ltemp+1:ltemp+20) = ' traceback requested'
            ltemp = ltemp + 20
         else
            temp(ltemp+1:ltemp+24) = ' traceback not requested'
            ltemp = ltemp + 24
         endif
         call xerprn (' ***', -1, temp(1:ltemp), 72)
      endif
c
c       now send out the message.
c
      call xerprn (' *  ', -1, messg, 72)
c
c       if lkntrl is positive, write the error number and request a
c          traceback.
c
      if (lkntrl .gt. 0) then
         write (temp, '(''error number = '', i8)') nerr
         do 10 i=16,22
            if (temp(i:i) .ne. ' ') go to 20
   10    continue
c
   20    call xerprn (' *  ', -1, temp(1:15) // temp(i:23), 72)
         call fdump
      endif
c
c       if lkntrl is not zero, print a blank line and an end of message.
c
      if (lkntrl .ne. 0) then
         call xerprn (' *  ', -1, ' ', 72)
         call xerprn (' ***', -1, 'end of message', 72)
         call xerprn ('    ',  0, ' ', 72)
      endif
c
c       if the error is not fatal or the error is recoverable and the
c       control flag is set for recovery, then return.
c
   30 if (level.le.0 .or. (level.eq.1 .and. mkntrl.le.1)) return
c
c       the program will be stopped due to an unrecovered error or a
c       fatal error.  print the reason for the abort and the error
c       summary if the control flag and the maximum error count permit.
c
      if (lkntrl.gt.0 .and. kount.lt.max(1,maxmes)) then
         if (level .eq. 1) then
            call xerprn
     *         (' ***', -1, 'job abort due to unrecovered error.', 72)
         else
            call xerprn(' ***', -1, 'job abort due to fatal error.', 72)
         endif
         call xersve (' ', ' ', ' ', -1, 0, 0, kdummy)
         call xerhlt (' ')
      else
         call xerhlt (messg)
      endif
      return
      end
*deck xerprn
      subroutine xerprn (prefix, npref, messg, nwrap)
c***begin prologue  xerprn
c***subsidiary
c***purpose  print error messages processed by xermsg.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xerprn-a)
c***keywords  error messages, printing, xerror
c***author  fong, kirby, (nmfecc at llnl)
c***description
c
c this routine sends one or more lines to each of the (up to five)
c logical units to which error messages are to be sent.  this routine
c is called several times by xermsg, sometimes with a single line to
c print and sometimes with a (potentially very long) message that may
c wrap around into multiple lines.
c
c prefix  input argument of type character.  this argument contains
c         characters to be put at the beginning of each line before
c         the body of the message.  no more than 16 characters of
c         prefix will be used.
c
c npref   input argument of type integer.  this argument is the number
c         of characters to use from prefix.  if it is negative, the
c         intrinsic function len is used to determine its length.  if
c         it is zero, prefix is not used.  if it exceeds 16 or if
c         len(prefix) exceeds 16, only the first 16 characters will be
c         used.  if npref is positive and the length of prefix is less
c         than npref, a copy of prefix extended with blanks to length
c         npref will be used.
c
c messg   input argument of type character.  this is the text of a
c         message to be printed.  if it is a long message, it will be
c         broken into pieces for printing on multiple lines.  each line
c         will start with the appropriate prefix and be followed by a
c         piece of the message.  nwrap is the number of characters per
c         piece; that is, after each nwrap characters, we break and
c         start a new line.  in addition the characters '$$' embedded
c         in messg are a sentinel for a new line.  the counting of
c         characters up to nwrap starts over for each new line.  the
c         value of nwrap typically used by xermsg is 72 since many
c         older error messages in the slatec library are laid out to
c         rely on wrap-around every 72 characters.
c
c nwrap   input argument of type integer.  this gives the maximum size
c         piece into which to break messg for printing on multiple
c         lines.  an embedded '$$' ends a line, and the count restarts
c         at the following character.  if a line break does not occur
c         on a blank (it would split a word) that word is moved to the
c         next line.  values of nwrap less than 16 will be treated as
c         16.  values of nwrap greater than 132 will be treated as 132.
c         the actual line length will be npref + nwrap after npref has
c         been adjusted to fall between 0 and 16 and nwrap has been
c         adjusted to fall between 16 and 132.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  i1mach, xgetua
c***revision history  (yymmdd)
c   880621  date written
c   880708  revised after the slatec cml subcommittee meeting of
c           june 29 and 30 to change the name to xerprn and to rework
c           the handling of the new line sentinel to behave like the
c           slash character in format statements.
c   890706  revised with the help of fred fritsch and reg clemens to
c           streamline the coding and fix a bug that caused extra blank
c           lines to be printed.
c   890721  revised to add a new feature.  a negative value of npref
c           causes len(prefix) to be used as the length.
c   891013  revised to correct error in calculating prefix length.
c   891214  prologue converted to version 4.0 format.  (wrb)
c   900510  added code to break messages between words.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xerprn
      character*(*) prefix, messg
      integer npref, nwrap
      character*148 cbuff
      integer iu(5), nunit
      character*2 newlin
      parameter (newlin = '$$')
c***first executable statement  xerprn
      call xgetua(iu,nunit)
c
c       a zero value for a logical unit number means to use the standard
c       error message unit instead.  i1mach(4) retrieves the standard
c       error message unit.
c
      n = i1mach(4)
      do 10 i=1,nunit
         if (iu(i) .eq. 0) iu(i) = n
   10 continue
c
c       lpref is the length of the prefix.  the prefix is placed at the
c       beginning of cbuff, the character buffer, and kept there during
c       the rest of this routine.
c
      if ( npref .lt. 0 ) then
         lpref = len(prefix)
      else
         lpref = npref
      endif
      lpref = min(16, lpref)
      if (lpref .ne. 0) cbuff(1:lpref) = prefix
c
c       lwrap is the maximum number of characters we want to take at one
c       time from messg to print on one line.
c
      lwrap = max(16, min(132, nwrap))
c
c       set lenmsg to the length of messg, ignore any trailing blanks.
c
      lenmsg = len(messg)
      n = lenmsg
      do 20 i=1,n
         if (messg(lenmsg:lenmsg) .ne. ' ') go to 30
         lenmsg = lenmsg - 1
   20 continue
   30 continue
c
c       if the message is all blanks, then print one blank line.
c
      if (lenmsg .eq. 0) then
         cbuff(lpref+1:lpref+1) = ' '
         do 40 i=1,nunit
            write(iu(i), '(a)') cbuff(1:lpref+1)
   40    continue
         return
      endif
c
c       set nextc to the position in messg where the next substring
c       starts.  from this position we scan for the new line sentinel.
c       when nextc exceeds lenmsg, there is no more to print.
c       we loop back to label 50 until all pieces have been printed.
c
c       we look for the next occurrence of the new line sentinel.  the
c       index intrinsic function returns zero if there is no occurrence
c       or if the length of the first argument is less than the length
c       of the second argument.
c
c       there are several cases which should be checked for in the
c       following order.  we are attempting to set lpiece to the number
c       of characters that should be taken from messg starting at
c       position nextc.
c
c       lpiece .eq. 0   the new line sentinel does not occur in the
c                       remainder of the character string.  lpiece
c                       should be set to lwrap or lenmsg+1-nextc,
c                       whichever is less.
c
c       lpiece .eq. 1   the new line sentinel starts at messg(nextc:
c                       nextc).  lpiece is effectively zero, and we
c                       print nothing to avoid producing unnecessary
c                       blank lines.  this takes care of the situation
c                       where the library routine has a message of
c                       exactly 72 characters followed by a new line
c                       sentinel followed by more characters.  nextc
c                       should be incremented by 2.
c
c       lpiece .gt. lwrap+1  reduce lpiece to lwrap.
c
c       else            this last case means 2 .le. lpiece .le. lwrap+1
c                       reset lpiece = lpiece-1.  note that this
c                       properly handles the end case where lpiece .eq.
c                       lwrap+1.  that is, the sentinel falls exactly
c                       at the end of a line.
c
      nextc = 1
   50 lpiece = index(messg(nextc:lenmsg), newlin)
      if (lpiece .eq. 0) then
c
c       there was no new line sentinel found.
c
         idelta = 0
         lpiece = min(lwrap, lenmsg+1-nextc)
         if (lpiece .lt. lenmsg+1-nextc) then
            do 52 i=lpiece+1,2,-1
               if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
                  lpiece = i-1
                  idelta = 1
                  goto 54
               endif
   52       continue
         endif
   54    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      elseif (lpiece .eq. 1) then
c
c       we have a new line sentinel at messg(nextc:nextc+1).
c       don't print a blank line.
c
         nextc = nextc + 2
         go to 50
      elseif (lpiece .gt. lwrap+1) then
c
c       lpiece should be set down to lwrap.
c
         idelta = 0
         lpiece = lwrap
         do 56 i=lpiece+1,2,-1
            if (messg(nextc+i-1:nextc+i-1) .eq. ' ') then
               lpiece = i-1
               idelta = 1
               goto 58
            endif
   56    continue
   58    cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc = nextc + lpiece + idelta
      else
c
c       if we arrive here, it means 2 .le. lpiece .le. lwrap+1.
c       we should decrement lpiece by one.
c
         lpiece = lpiece - 1
         cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
         nextc  = nextc + lpiece + 2
      endif
c
c       print
c
      do 60 i=1,nunit
         write(iu(i), '(a)') cbuff(1:lpref+lpiece)
   60 continue
c
      if (nextc .le. lenmsg) go to 50
      return
      end
*deck xersve
      subroutine xersve (librar, subrou, messg, kflag, nerr, level,
     +   icount)
c***begin prologue  xersve
c***subsidiary
c***purpose  record that an error has occurred.
c***library   slatec (xerror)
c***category  r3
c***type      all (xersve-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c *usage:
c
c        integer  kflag, nerr, level, icount
c        character * (len) librar, subrou, messg
c
c        call xersve (librar, subrou, messg, kflag, nerr, level, icount)
c
c *arguments:
c
c        librar :in    is the library that the message is from.
c        subrou :in    is the subroutine that the message is from.
c        messg  :in    is the message to be saved.
c        kflag  :in    indicates the action to be performed.
c                      when kflag > 0, the message in messg is saved.
c                      when kflag=0 the tables will be dumped and
c                      cleared.
c                      when kflag < 0, the tables will be dumped and
c                      not cleared.
c        nerr   :in    is the error number.
c        level  :in    is the error severity.
c        icount :out   the number of times this message has been seen,
c                      or zero if the table has overflowed and does not
c                      contain this message specifically.  when kflag=0,
c                      icount will not be altered.
c
c *description:
c
c   record that this error occurred and possibly dump and clear the
c   tables.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  i1mach, xgetua
c***revision history  (yymmdd)
c   800319  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900413  routine modified to remove reference to kflag.  (wrb)
c   900510  changed to add library name and subroutine to calling
c           sequence, use if-then-else, make number of saved entries
c           easily changeable, changed routine name from xersav to
c           xersve.  (rwc)
c   910626  added libtab and subtab to save statement.  (bks)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xersve
      parameter (lentab=10)
      integer lun(5)
      character*(*) librar, subrou, messg
      character*8  libtab(lentab), subtab(lentab), lib, sub
      character*20 mestab(lentab), mes
      dimension nertab(lentab), levtab(lentab), kount(lentab)
      save libtab, subtab, mestab, nertab, levtab, kount, kountx, nmsg
      data kountx/0/, nmsg/0/
c***first executable statement  xersve
c
      if (kflag.le.0) then
c
c        dump the table.
c
         if (nmsg.eq.0) return
c
c        print to each unit.
c
         call xgetua (lun, nunit)
         do 20 kunit = 1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c
c           print the table header.
c
            write (iunit,9000)
c
c           print body of table.
c
            do 10 i = 1,nmsg
               write (iunit,9010) libtab(i), subtab(i), mestab(i),
     *            nertab(i),levtab(i),kount(i)
   10       continue
c
c           print number of other errors.
c
            if (kountx.ne.0) write (iunit,9020) kountx
            write (iunit,9030)
   20    continue
c
c        clear the error tables.
c
         if (kflag.eq.0) then
            nmsg = 0
            kountx = 0
         endif
      else
c
c        process a message...
c        search for this messg, or else an empty slot for this messg,
c        or else determine that the error table is full.
c
         lib = librar
         sub = subrou
         mes = messg
         do 30 i = 1,nmsg
            if (lib.eq.libtab(i) .and. sub.eq.subtab(i) .and.
     *         mes.eq.mestab(i) .and. nerr.eq.nertab(i) .and.
     *         level.eq.levtab(i)) then
                  kount(i) = kount(i) + 1
                  icount = kount(i)
                  return
            endif
   30    continue
c
         if (nmsg.lt.lentab) then
c
c           empty slot found for new message.
c
            nmsg = nmsg + 1
            libtab(i) = lib
            subtab(i) = sub
            mestab(i) = mes
            nertab(i) = nerr
            levtab(i) = level
            kount (i) = 1
            icount    = 1
         else
c
c           table is full.
c
            kountx = kountx+1
            icount = 0
         endif
      endif
      return
c
c     formats.
c
 9000 format ('0          error message summary' /
     +   ' library    subroutine message start             nerr',
     +   '     level     count')
 9010 format (1x,a,3x,a,3x,a,3i10)
 9020 format ('0other errors not individually tabulated = ', i10)
 9030 format (1x)
      end
*deck xgetua
      subroutine xgetua (iunita, n)
c***begin prologue  xgetua
c***purpose  return unit number(s) to which error messages are being
c            sent.
c***library   slatec (xerror)
c***category  r3c
c***type      all (xgetua-a)
c***keywords  error, xerror
c***author  jones, r. e., (snla)
c***description
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n .lt. 5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c***references  r. e. jones and d. k. kahaner, xerror, the slatec
c                 error-handling package, sand82-0800, sandia
c                 laboratories, 1982.
c***routines called  j4save
c***revision history  (yymmdd)
c   790801  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  xgetua
      dimension iunita(5)
c***first executable statement  xgetua
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue
      return
      end







