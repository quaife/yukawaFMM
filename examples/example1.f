      program example1
c********************************************************************
c  ****** DESCRIPTION:
c  Compute \sum q_{k}K_0(\beta |x_{j}-x_{k}|).  Can also compute
c  derivatives

c  ****** PARAMETERS:
c  npts     = total number of sources/targets
c  beta     = parameter in the argument of K_0
c  dom      = location of targets/sources
c  density  = charge densities
c  layers   = 1 only compute the single-layer potential
c           = 2 compute the single layer-potential and its gradient
c           = 3 compute the single layer-potential, its gradient
c               and its hessian
c maxptsperbox  = maxiumum number of points fmm is allowed to put 
c                in a childless box

c  ****** OUTPUT PARAMETERS:
c  potE   = single-layer potential using direct evaluation
c  potEx  = x-derivate of single-layer potential using direct evaluation
c  potEy  = y-derivate of single-layer potential using direct evaluation
c  potExx = xx-derivative of single-layer potential using 
c           direct evaluation
c  potExy = xx-derivative of single-layer potential using 
c           direct evaluation
c  potEyy = xx-derivative of single-layer potential using 
c           direct evaluation
c  potF   = single-layer potential using fmm
c  potFx  = x-derivate of single-layer potential using fmm
c  potFy  = y-derivate of single-layer potential using fmm
c  potFxx = xx-derivative of single-layer potential using fmm
c  potExy = xx-derivative of single-layer potential using fmm
c  potFyy = xx-derivative of single-layer potential using fmm
c  error = l2 error of potential between direct and fmm
c  errorx = l2 error of x-derivative of potential between direct
c           and fmm
c  errory = l2 error of y-derivative of potential between direct
c           and fmm
c  errorxx = l2 error of xx-derivative of potential between direct
c           and fmm
c  errorxy = l2 error of xy-derivative of potential between direct
c           and fmm
c  erroryy = l2 error of yy-derivative of potential between direct
c           and fmm


      implicit none

      integer *8 npts
c     number of points
      real*8 beta
c     Single parameter in the problem
      complex *16, allocatable :: dom(:)
c     The charge locations stored as complex variables ie. z = x+iy
      real *8, allocatable :: density(:) 
c     The charge density.  For now they must be real
      real *8, allocatable :: potE(:),potEx(:),potEy(:)
      real *8, allocatable :: potExx(:),potExy(:),potEyy(:)
c     Potential ie. Yukawa single-layer potential using the direct
c     method, ie. N^{2} algorithm
c     Gradient of single-layer potential using the direct method
c     Also may want the gradient of the double-layer potential
      complex *16, allocatable :: potF(:),potFx(:),potFy(:)
c     Potential ie. Yukawa single-layer potential using the FMM
      complex *16, allocatable :: potFxx(:),potFxy(:),potFyy(:)
c     Gradient of single-layer potential using the FMM
c     Also may want the gradient of the double-layer potential
      integer *4 layers
c     The number of layers you need to compute (1,2,3,etc).
c     For now can only do 2
      real *8 error,errorx,errory,errorxx,errorxy,erroryy
      real *8 potNorm,potxNorm,potyNorm,potxxNorm,potxyNorm,potyyNorm
      integer i,maxptsperbox
      real *4 timeep(2),ETIME
      real *4 t0_direct,t1_direct
      real *4 t0_fmm,t1_fmm
      integer stat,ierror

      layers = 3
c     Number of layer potentials to compute

      npts = 2**12
      allocate(dom(npts),stat = ierror);
      call emessage(stat)
      allocate(density(npts),stat = ierror);
      call emessage(stat)
      allocate(potE(npts),stat = ierror);
      call emessage(stat)
      allocate(potEx(npts),stat = ierror);
      call emessage(stat)
      allocate(potEy(npts),stat = ierror);
      call emessage(stat)
      allocate(potExx(npts),stat = ierror);
      call emessage(stat)
      allocate(potExy(npts),stat = ierror);
      call emessage(stat)
      allocate(potEyy(npts),stat = ierror);
      call emessage(stat)
      allocate(potF(npts),stat = ierror);
      call emessage(stat)
      allocate(potFx(npts),stat = ierror);
      call emessage(stat)
      allocate(potFy(npts),stat = ierror);
      call emessage(stat)
      allocate(potFxx(npts),stat = ierror);
      call emessage(stat)
      allocate(potFxy(npts),stat = ierror);
      call emessage(stat)
      allocate(potFyy(npts),stat = ierror);
      call emessage(stat)
c     allocate memeory for domain, density, and outputs


c     Set up the domain, density, and parameter beta
      call initData(npts,dom,density,beta)

c     Compute the potential and gradient directly
      t0_direct = ETIME(timeep)
      call potExact(npts,layers,dom,density,beta,potE,potEx,potEy,
     $    potExx,potExy,potEyy)
      t1_direct = ETIME(timeep)

      
c     Maximum number of points per box
      maxptsperbox = 10 
c     Compute the potential and gradient using the FMM
      t0_fmm = ETIME(timeep)
      call potFMM(npts,maxptsperbox,layers,dom,density,beta,
     $   potF,potFx,potFy,potFxx,potFxy,potFyy)
      t1_fmm = ETIME(timeep)


      open(unit=1,file='../output/pot.dat')
      do i=1,npts
        write(1,1030) potE(i)
      enddo
      close(unit=1)
      open(unit=1,file='../output/potx.dat')
      do i=1,npts
        write(1,1030) potEx(i)
      enddo
      close(unit=1)
      open(unit=1,file='../output/poty.dat')
      do i=1,npts
        write(1,1030) potEy(i)
      enddo
      close(unit=1)
      open(unit=1,file='../output/potxx.dat')
      do i=1,npts
        write(1,1030) potExx(i)
      enddo
      close(unit=1)
      open(unit=1,file='../output/potxy.dat')
      do i=1,npts
        write(1,1030) potExy(i)
      enddo
      close(unit=1)
      open(unit=1,file='../output/potyy.dat')
      do i=1,npts
        write(1,1030) potEyy(i)
      enddo
      close(unit=1)
c     Write all the potentials to files

c     Compute errors to check accuracy of FMM
      error = 0.d0
      errorx = 0.d0
      errory = 0.d0
      errorxx = 0.d0
      errorxy = 0.d0
      erroryy = 0.d0
c     To find the l2 errors
      potNorm = 0.d0
      potxNorm = 0.d0
      potyNorm = 0.d0
      potxxNorm = 0.d0
      potxyNorm = 0.d0
      potyyNorm = 0.d0
c     To find the l2 norm of the true layer potentials

      do i=1,npts
        error = error + (potE(i)-real(potF(i)))**2.d0
        errorx = errorx + (potEx(i)-real(potFx(i)))**2.d0
        errory = errory + (potEy(i)-real(potFy(i)))**2.d0
        errorxx = errorxx + (potExx(i)-real(potFxx(i)))**2.d0
        errorxy = errorxy + (potExy(i)-real(potFxy(i)))**2.d0
        erroryy = erroryy + (potEyy(i)-real(potFyy(i)))**2.d0

        potNorm = potNorm + potE(i)**2.d0
        potxNorm = potxNorm + potEx(i)**2.d0
        potyNorm = potyNorm + potEy(i)**2.d0
        potxxNorm = potxxNorm + potExx(i)**2.d0
        potxyNorm = potxyNorm + potExy(i)**2.d0
        potyyNorm = potyyNorm + potEyy(i)**2.d0
      enddo
c     compute absolve l2 errors and l2 norms of exact potentials
      error = sqrt(error)/sqrt(potNorm)
      errorx = sqrt(errorx)/sqrt(potxNorm)
      errory = sqrt(errory)/sqrt(potyNorm)
      errorxx = sqrt(errorxx)/sqrt(potxxNorm)
      errorxy = sqrt(errorxy)/sqrt(potxyNorm)
      erroryy = sqrt(erroryy)/sqrt(potyyNorm)
c     comptue relative l2 errors

      write(6,*) ''
      write(6,*) '******************ERRORS*******************'
      write(6,1010) 'l^{2} error in potential is ',error
      if (layers .gt. 1) then
        write(6,1010) 'l^{2} error in gradient is ',errorx
        write(6,1010) 'l^{2} error in gradient is ',errory
      endif
      if (layers .gt. 2) then
        write(6,1010) 'l^{2} error in second derivative is ',errorxx
        write(6,1010) 'l^{2} error in second derivative is ',errorxy
        write(6,1010) 'l^{2} error in second derivative is ',erroryy
      endif
      write(6,*) ''

      write(6,*) '******************TIMINGS******************'
      write(6,1020) 'Time to do Direct Evaluation ',
     $       t1_direct-t0_direct,' seconds' 
      write(6,1020) 'Time to do FMM Evaluation ',
     $       t1_fmm-t0_fmm,' seconds' 



 1010 format(A,ES8.2)
 1020 format(A,ES8.2,A)
 1030 format(F25.16)

      end



c********************************************************************
      subroutine initData(npts,dom,density,beta)
c  ****** DESCRIPTION:
c  build the source/target points and chose parameter beta 

c  ****** PARAMETERS:
c  npts     = total number of sources/targets

c  ****** OUTPUT PARAMETERS:
c  beta     = parameter in the argument of K_0
c  dom      = location of targets/sources
c  density  = charge densities

      implicit none

      integer npts
      real*8 beta
      complex *16 dom(npts)
      real*8 density(npts)

      integer i,ierror,stat
      real*8 theta,radius,twopi
      complex *16 eye

      twopi = 8.0d0*datan(1.0d0)
      eye = (0.d0,1.d0)

      beta = 1.d0

      do i=1,npts
        theta = dble(i-1)*twopi/dble(npts)
        density(i) = 1.d0/dble(npts)
        radius = 0.3d0+1.d-1*dcos(theta)+5.d-2*dsin(15.d0*theta)
        dom(i) = radius*exp(1.d0*eye*(theta+1.d-1))
      enddo
c     Assign a density function and target/source locations



      return
      end
      
c*******************************************************************
      subroutine potExact(npts,layers,dom,density,beta,potE,potEx,potEy,
     $      potExx,potExy,potEyy)
c  ****** DESCRIPTION:
c  build the potential using the N^2 direct algorithm.  Also computed
c  the requested number of derivatives

c  ****** PARAMETERS:
c  npts     = total number of sources/targets
c  layers   = number of derivatives to compute
c  dom      = complex discretization of the boundary
c  density  = charges/density function of sources
c  beta     = parameter in the kernel

c  ****** OUTPUT PARAMETERS:
c  potE     = potential
c  potEx    = x-derivative of the potential
c  potEy    = y-derivative of the potential
c  potExx   = xx-derivative of the potential
c  potExy   = xy-derivative of the potential
c  potEyy   = yy-derivative of the potential

      implicit none

      integer npts,layers
      real*8 beta
      complex *16 dom(npts)
      real*8 density(npts)
      real*8 potE(npts)
      real*8 potEx(npts),potEy(npts)
      real*8 potExx(npts),potExy(npts),potEyy(npts)

      integer i,j
      real*8 dist
      real*8 bk(0:1)
      integer ncalc

      
      do i=1,npts
        potE(i) = 0.d0
        potEx(i) = 0.d0
        potEy(i) = 0.d0
        potExx(i) = 0.d0
        potExy(i) = 0.d0
        potEyy(i) = 0.d0
      enddo
c     Initalize the potential to be zero 

      if (layers .eq. 1) then
        do i=1,npts
          do j=1,npts
            if (i .ne. j) then
              dist = abs(dom(i)-dom(j))
              call RKBESL(beta*dist,0.d0,1,1,bk,ncalc)
              potE(i) = potE(i) + density(j)*bk(0)
c             update the potential
            endif
          enddo
        enddo


      elseif (layers .eq. 2) then
        do i=1,npts
          do j=1,npts
            if (i .ne. j) then
              dist = abs(dom(i)-dom(j))
              call RKBESL(beta*dist,0.d0,2,1,bk,ncalc)
              potE(i) = potE(i) + density(j)*bk(0)
              potEx(i) = potEx(i) + density(j)*bk(1)/dist*beta*
     $           real(dom(i)-dom(j))
              potEy(i) = potEy(i) + density(j)*bk(1)/dist*beta*
     $           imag(dom(i)-dom(j))
            endif
          enddo
        enddo

      elseif (layers .eq. 3) then
        do i=1,npts
          do j=1,npts
            if (i .ne. j) then
              dist = abs(dom(i)-dom(j))
              call RKBESL(beta*dist,0.d0,2,1,bk,ncalc)
              potE(i) = potE(i) + density(j)*bk(0)
              potEx(i) = potEx(i) + density(j)*bk(1)/dist*beta*
     $           real(dom(i)-dom(j))
              potEy(i) = potEy(i) + density(j)*bk(1)/dist*beta*
     $           imag(dom(i)-dom(j))
              potExx(i) = potExx(i) + density(j)*
     $          (bk(0)/dist**2.d0 * real(dom(i)-dom(j))**2.d0 +
     $          bk(1)/dist-
     $          2.0d0*imag(dom(i)-dom(j))**2.d0*bk(1)/dist**3.d0)
              potExy(i) = potExy(i) + density(j)*
     $          (bk(0)/dist**2.d0 + bk(1)/dist**3.d0 * 
     $          2.d0)*real(dom(i)-dom(j))*imag(dom(i)-dom(j))
              potEyy(i) = potEyy(i) + density(j)*
     $          (bk(0)/dist**2.d0 * imag(dom(i)-dom(j))**2.d0 +
     $          bk(1)/dist-
     $          2.0d0*real(dom(i)-dom(j))**2.d0*bk(1)/dist**3.d0)
            endif
          enddo
        enddo


      endif


        

      return
      end
      


c********************************************************************
      subroutine emessage(stat)
c  ****** DESCRIPTION:
c  display an error message if alloc had trouble allocating memory
      implicit none

      integer stat

      if (stat .ne. 0) then
        print*,'PROBLEM ALLOCATING MEMORY'
      endif

      return
      end


