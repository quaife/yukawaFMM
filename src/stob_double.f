c*********************************************************************
       SUBROUTINE STOB_DTOP_CHEBY_TAB(WINT,XWINT,YWINT)
c*********************************************************************
c
c      the precomputed table for small to big double layer top
c        interaction
c
c      output :
c        wint  : the precomputed chebyshev coefficients.
c        xwint : the precomputed chebyshev coefficients, x coefficient.
c        ywint : the precomputed chebyshev coefficients, y coefficient.
c
c********************************************************************
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l,m
       REAL *8 WINT(1:7,1:21,0:3,1:12,1:16)
       REAL *8 XWINT(1:7,1:21,0:3,1:12,1:16)
       REAL *8 YWINT(1:7,1:21,0:3,1:12,1:16)
c
       open(10,file='data/STOB_DTOP_CHEBY_TAB.dat',form='unformatted')
       do i=1,7
       do j=1,21
       do k=0,3
       do l=1,12
       do m=1,16
         read(10) wint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
cccccccccccccccthe following is to be added.
       return
ccccccccccccccc
c
       open(10,file='data/XSTOB_DTOP_CHEBY_TAB.dat',form='unformatted')
       do i=1,7
       do j=1,21
       do k=0,3
       do l=1,12
       do m=1,16
         read(10) xwint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10,file='data/YSTOB_DTOP_CHEBY_TAB.dat',form='unformatted')
       do i=1,7
       do j=1,21
       do k=0,3
       do l=1,12
       do m=1,16
         read(10) ywint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       RETURN
       END
c
       SUBROUTINE STOB_DT_CORRECT_TAB(WINT,XWINT,YWINT)
c*********************************************************************
c
c      the precomputed table for big to small double layer top interaction
c
c      output :
c        wint  : the precomputed chebyshev coefficients.
c        xwint : the precomputed chebyshev coefficients, x coefficient.
c        ywint : the precomputed chebyshev coefficients, y coefficient.
c
c********************************************************************
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l
       REAL *8 WINT(7,0:3,12,16)
       REAL *8 XWINT(7,0:3,12,16)
       REAL *8 YWINT(7,0:3,12,16)
c
       open(10,file='data/STOB_DT_CORRECT_TAB.dat',form='unformatted')
       do i=1,7
       do j=0,3
       do k=1,12
       do l=1,16
         read(10) wint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
ccccccccc
       return
cccccccccc
c
       open(10,file='data/XSTOB_DT_CORRECT_TAB.dat',form='unformatted')
       do i=1,7
       do j=0,3
       do k=1,12
       do l=1,16
         read(10) xwint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10,file='data/YSTOB_DT_CORRECT_TAB.dat',form='unformatted')
       do i=1,7
       do j=0,3
       do k=1,12
       do l=1,16
         read(10) ywint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       RETURN
       END
c

