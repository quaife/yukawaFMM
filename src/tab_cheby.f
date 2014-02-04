c************************************************************************
       SUBROUTINE FORMMP_CHEBY_TAB(WINT)
c*********************************************************************
c
c      the precomputed table for forming multipole expansion from the
c      volume potentials.
c
c      output :
c        wint : the precomputed chebyshev coefficients.
c
c********************************************************************
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l
       REAL *8 WINT(1:4,1:21,0:10,10)
c
       open(10, file='../data/FORMMP_CHEBY_TAB.dat', form='unformatted')
       do i=1,4
       do j=1,21
       do k=0,10
       do l=1,10
         read(10) wint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       RETURN
       END
c*********************************************************************
       SUBROUTINE FORMPSLP_CHEBY_TAB(WINT)
c*********************************************************************
c
c      the precomputed table for forming multipole expansion from the
c      single layer potentials.
c
c      output :
c        wint : the precomputed chebyshev coefficients.
c
c********************************************************************
c
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l
       REAL *8 WINT(1:4,1:21,0:3,0:41)
c
       open(10, file='data/FORMPSLP_CHEBY_TAB.dat', form='unformatted')
       do i=1,4
       do j=1,21
       do k=0,3
       do l=0,41
         read(10) wint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       RETURN
       END
c
c*********************************************************************
       SUBROUTINE FORMPDLP_CHEBY_TAB(WINT)
c*********************************************************************
c
c      the precomputed table for forming multipole expansion from the
c      double layer potentials.
c
c      output :
c        wint : the precomputed chebyshev coefficients.
c
c********************************************************************
c
       IMPLICIT NONE
       INTEGER *4 i,j,k,l
       REAL *8 WINT(1:4,1:21,0:3,0:41)
c
       open(10, file='data/FORMPDLP_CHEBY_TAB.dat', form='unformatted')
       do i=1,4
       do j=1,21
       do k=0,3
       do l=0,41
         read(10) wint(i,j,k,l)
       enddo
       enddo
       enddo
       enddo
       close(10)
c
       RETURN
       END
c
