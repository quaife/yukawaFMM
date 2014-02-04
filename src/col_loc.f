c*********************************************************************
       SUBROUTINE COL_CHEBY_TAB(WINT,XWINT,YWINT)
c*********************************************************************
c
c      the precomputed table for col to col box-box interaction
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
       REAL *8 WINT(8,21,10,9,16)
       REAL *8 XWINT(9,21,10,9,16)
       REAL *8 YWINT(9,21,10,9,16)
c
       open(10, file='data/COL_CHEBY_TAB.dat', form='unformatted')
       do i=1,8
       do j=1,21
       do k=1,10
       do l=1,9
       do m=1,16
         read(10) wint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10, file='data/XCOL_CHEBY_TAB.dat', form='unformatted')
       do i=1,9
       do j=1,21
       do k=1,10
       do l=1,9
       do m=1,16
         read(10) xwint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10, file='data/YCOL_CHEBY_TAB.dat', form='unformatted')
       do i=1,9
       do j=1,21
       do k=1,10
       do l=1,9
       do m=1,16
         read(10) ywint(i,j,k,l,m)
       enddo
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       return
       end
c
c*********************************************************************
       SUBROUTINE COL_CORRECT_TAB(WINT,XWINT,YWINT)
c*********************************************************************
c
c      the precomputed table for box-box col to col interactions.
c      Note: for the 1st interval, the following corrections
c        are necessary.
c
c      output :
c        wint : the precomputed chebyshev coefficients.
c        xwint : the precomputed chebyshev coefficients, x coefficient.
c        ywint : the precomputed chebyshev coefficients, y coefficient.
c
c********************************************************************
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l
       REAL *8 WINT(1:7,1:10,1:9,1:16)
       REAL *8 XWINT(1:7,1:10,1:9,1:16)
       REAL *8 YWINT(1:7,1:10,1:9,1:16)
c
       open(10, file='data/COL_CORRECT_TAB.dat', form='unformatted')
       do i=1,7
       do j=1,10
       do k=1,9
       do l=1,16
         read(10) wint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10, file='data/XCOL_CORRECT_TAB.dat', form='unformatted')
       do i=1,7
       do j=1,10
       do k=1,9
       do l=1,16
         read(10) xwint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       open(10, file='data/YCOL_CORRECT_TAB.dat', form='unformatted')
       do i=1,7
       do j=1,10
       do k=1,9
       do l=1,16
         read(10) ywint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       return
       end
