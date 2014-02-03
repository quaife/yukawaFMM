      program example2
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
c  iflagDirect = 1 if doing N-body directly.  0 otherwise
c  iflagFMM = 1 if doing N-body with FMM.  0 otherwise


      implicit none

      integer *8 npts,ntar
c     number of points on boundary and number of targets
      real*8 beta
c     Single parameter in the problem
      complex *16, allocatable :: dom(:)
c     domain locations stored as complex variables ie. z = x+iy
      complex *16, allocatable :: zz(:)
c     Target points where the volume potential is to be evaulated
      real *8, allocatable :: VLP(:)
c     Volume Potential 

      integer i
      real *4 timeep(2),ETIME
      real *4 t0,t1
      integer stat,ierror

      npts = 2**8
      ntar = 2**4
      allocate(dom(npts),stat = ierror);
      call emessage(stat)
      allocate(zz(ntar),stat = ierror);
      call emessage(stat)
      allocate(VLP(ntar),stat = ierror);
      call emessage(stat)
c     allocate memeory for domain, targets, and outputs


c     Set up the domain, density, and parameter beta
      call initData(npts,ntar,dom,beta,zz)

      t0 = ETIME(timeep)
      call forced(beta,npts,dom,ntar,zz,VLP)
      t1 = ETIME(timeep)

c      print*,VLP(1)


      write(6,*) '******************TIMINGS******************'
      write(6,1020) 'Evaluation time is ',
     $       t1 - t0,' seconds' 


 1020 format(A,ES8.2,A)

      end



c********************************************************************
      subroutine initData(npts,ntar,dom,beta,zz)
c  ****** DESCRIPTION:
c  build the source/target points and chose parameter beta 

c  ****** PARAMETERS:
c  npts     = total number of sources/targets

c  ****** OUTPUT PARAMETERS:
c  beta     = parameter in the argument of K_0
c  dom      = location of targets/sources

      implicit none

      integer npts,ntar
      real*8 beta
      complex *16 dom(npts)
      complex *16 zz(ntar)

      integer i,ierror,stat
      real*8 theta,radius,twopi
      complex *16 eye

      twopi = 8.0d0*datan(1.0d0)
      eye = (0.d0,1.d0)

      beta = 1.d0

      do i=1,npts
        theta = dble(i-1)*twopi/dble(npts)
        radius = 0.3d0+1.d-1*dcos(theta)+5.d-2*dsin(15.d0*theta)
        dom(i) = radius*exp(1.d0*eye*(theta+1.d-1))
      enddo
c     set up boundaries 

      do i=1,ntar
        zz(i) = (1.d0,0.d0)
      enddo


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


