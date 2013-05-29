       SUBROUTINE RIBESL(X,ALPHA,NB,IZE,B,NCALC)
C-------------------------------------------------------------------
C
C  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA,
C  with or without exponential scaling.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         I's or exponentially scaled I's (I*EXP(-X))
C         are to be calculated.  If I's are to be calculated,
C         X must be less than EXPARG (see below).
C ALPHA - Working precision fractional part of order for which
C         I's or exponentially scaled I's (I*EXP(-X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
C         and 2 if exponentially scaled I's are to be calculated.
C B     - Working precision output vector of length NB.  If the routine
C         terminates normally (NCALC=NB), the vector B contains the 
C         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
C         corresponding exponentially scaled functions.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector B, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C 
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
C            in mind that if ABS(X)=N, then at least N iterations
C            of the backward recursion will be executed.  The value
C            of 10.0 ** 4 is used on every machine.
C   EXPARG = Largest working precision argument that the library
C            EXP routine can handle and upper limit on the
C            magnitude of X when IZE=1; approximately 
C            LOG(beta**maxexp)
C
C
C     Approximate values for some important machines are:
C
C                        beta       minexp      maxexp       it
C
C  CRAY-1        (S.P.)    2        -8193        8191        48
C  Cyber 180/855
C    under NOS   (S.P.)    2         -975        1070        48
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2         -126         128        24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2        -1022        1024        53
C  IBM 3033      (D.P.)   16          -65          63        14
C  VAX           (S.P.)    2         -128         127        24
C  VAX D-Format  (D.P.)    2         -128         127        56
C  VAX G-Format  (D.P.)    2        -1024        1023        53
C
C
C                        NSIG       ENTEN       ENSIG      RTNSIG
C
C CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
C Cyber 180/855
C   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
C IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
C VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
C VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
C
C
C                         ENMTEN      XLARGE   EXPARG 
C
C CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677 
C Cyber 180/855
C   under NOS   (S.P.)   1.25E-293    1.0E+4     741
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88  
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
C IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
C VAX           (S.P.)   1.17E-38     1.0E+4      88
C VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
C VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error,  NCALC .NE. NB, and not all I's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. 0:  An argument is out of range. For example,
C     NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG.
C     In this case, the B-vector is not calculated, and NCALC is
C     set to MIN0(NB,0)-1 so that NCALC .NE. NB.
C
C  NB .GT. NCALC .GT. 0: Not all requested function values could
C     be calculated accurately.  This usually occurs because NB is
C     much larger than ABS(X).  In this case, B(N) is calculated
C     to the desired accuracy for N .LE. NCALC, but precision
C     is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C     for N .GT. NCALC (because it is too small to be represented),
C     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C     significant figures of B(N) can be trusted.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by David J.
C  Sookne (2) that computes values of the Bessel functions J or
C  I of real argument and integer order.  Modifications include
C  the restriction of the computation to the I Bessel function
C  of non-negative real argument, the extension of the computation
C  to arbitrary positive order, the inclusion of optional
C  exponential scaling, and the elimination of most underflow.
C  An earlier version was published in (3).
C
C References: "A Note on Backward Recurrence Algorithms," Olver,
C              F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C              pp 941-947.
C
C             "Bessel Functions of Real Argument and Integer Order,"
C              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C              125-132.
C
C             "ALGORITHM 597, Sequence of Modified Bessel Functions
C              of the First Kind," Cody, W. J., Trans. Math. Soft.,
C              1983, pp. 242-245.
C
C  Latest modification: May 30, 1989
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER IZE,K,L,MAGX,N,NB,NBMX,NCALC,NEND,NSIG,NSTART
CS    REAL              GAMMA,
      DOUBLE PRECISION DGAMMA,             
     1 ALPHA,B,CONST,CONV,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG,
     2 ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL,
     3 RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
       DIMENSION B(NB)
      
C-------------------------------------------------------------------
C  Mathematical constants
C-------------------------------------------------------------------
CS    DATA ONE,TWO,ZERO,HALF,CONST/1.0E0,2.0E0,0.0E0,0.5E0,1.585E0/
       DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
C-------------------------------------------------------------------
C  Machine-dependent parameters
C-------------------------------------------------------------------
CS    DATA NSIG,XLARGE,EXPARG /8,1.0E4,88.0E0/
CS    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN/4.7E-38/
       DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
       DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
       DATA ENMTEN/8.9D-308/
C-------------------------------------------------------------------
C  Statement functions for conversion
C-------------------------------------------------------------------
CS    CONV(N) = REAL(N)
CS    FUNC(X) = GAMMA(X)
      CONV(N) = DBLE(N)
      FUNC(X) = DGAMMA(X)
C-------------------------------------------------------------------
C Check for X, NB, OR IZE out of range.
C-------------------------------------------------------------------
      IF ((NB.GT.0) .AND. (X .GE. ZERO) .AND.
     1    (ALPHA .GE. ZERO) .AND. (ALPHA .LT. ONE) .AND.
     2    (((IZE .EQ. 1) .AND. (X .LE. EXPARG)) .OR.
     3     ((IZE .EQ. 2) .AND. (X .LE. XLARGE)))) THEN
C-------------------------------------------------------------------
C Use 2-term ascending series for small X
C-------------------------------------------------------------------
            NCALC = NB
            MAGX = INT(X)
            IF (X .GE. RTNSIG) THEN
C-------------------------------------------------------------------
C Initialize the forward sweep, the P-sequence of Olver
C-------------------------------------------------------------------
                  NBMX = NB-MAGX
                  N = MAGX+1
                  EN = CONV(N+N) + (ALPHA+ALPHA)
                  PLAST = ONE
                  P = EN / X
C-------------------------------------------------------------------
C Calculate general significance test
C-------------------------------------------------------------------
                  TEST = ENSIG + ENSIG
                  IF (2*MAGX .GT. 5*NSIG) THEN
                        TEST = SQRT(TEST*P)
                     ELSE
                        TEST = TEST / CONST**MAGX
                  END IF
                  IF (NBMX .GE. 3) THEN
C-------------------------------------------------------------------
C Calculate P-sequence until N = NB-1.  Check for possible overflow.
C-------------------------------------------------------------------
                     TOVER = ENTEN / ENSIG
                     NSTART = MAGX+2
                     NEND = NB - 1
                     DO 100 K = NSTART, NEND
                        N = K
                        EN = EN + TWO
                        POLD = PLAST
                        PLAST = P
                        P = EN * PLAST/X + POLD
                        IF (P .GT. TOVER) THEN
C-------------------------------------------------------------------
C To avoid overflow, divide P-sequence by TOVER.  Calculate
C P-sequence until ABS(P) .GT. 1.
C-------------------------------------------------------------------
                           TOVER = ENTEN
                           P = P / TOVER
                           PLAST = PLAST / TOVER
                           PSAVE = P
                           PSAVEL = PLAST
                           NSTART = N + 1
   60                      N = N + 1
                              EN = EN + TWO
                              POLD = PLAST
                              PLAST = P
                              P = EN * PLAST/X + POLD
                           IF (P .LE. ONE) GO TO 60
                           TEMPB = EN / X
C-------------------------------------------------------------------
C Calculate backward test, and find NCALC, the highest N
C such that the test is passed.
C-------------------------------------------------------------------
                           TEST = POLD*PLAST / ENSIG
                           TEST = TEST*(HALF-HALF/(TEMPB*TEMPB))
                           P = PLAST * TOVER
                           N = N - 1
                           EN = EN - TWO
                           NEND = MIN0(NB,N)
                           DO 80 L = NSTART, NEND
                              NCALC = L
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN * PSAVEL/X + POLD
                              IF (PSAVE*PSAVEL .GT. TEST) GO TO 90
   80                      CONTINUE
                           NCALC = NEND + 1
   90                      NCALC = NCALC - 1
                           GO TO 120
                        END IF
  100                CONTINUE
                     N = NEND
                     EN = CONV(N+N) + (ALPHA+ALPHA)
C-------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C-------------------------------------------------------------------
                     TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
                  END IF
C-------------------------------------------------------------------
C Calculate P-sequence until significance test passed.
C-------------------------------------------------------------------
  110             N = N + 1
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN * PLAST/X + POLD
                  IF (P .LT. TEST) GO TO 110
C-------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C-------------------------------------------------------------------
  120             N = N + 1
                  EN = EN + TWO
                  TEMPB = ZERO
                  TEMPA = ONE / P
                  EM = CONV(N) - ONE
                  EMPAL = EM + ALPHA
                  EMP2AL = (EM - ONE) + (ALPHA + ALPHA)
                  SUM = TEMPA * EMPAL * EMP2AL / EM
                  NEND = N - NB
                  IF (NEND .LT. 0) THEN
C-------------------------------------------------------------------
C N .LT. NB, so store B(N) and set higher orders to zero.
C-------------------------------------------------------------------
                        B(N) = TEMPA
                        NEND = -NEND
                        DO 130 L = 1, NEND
  130                      B(N+L) = ZERO
                     ELSE
                        IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Recur backward via difference equation, calculating (but
C not storing) B(N), until N = NB.
C-------------------------------------------------------------------
                           DO 140 L = 1, NEND
                              N = N - 1
                              EN = EN - TWO
                              TEMPC = TEMPB
                              TEMPB = TEMPA
                              TEMPA = (EN*TEMPB) / X + TEMPC
                              EM = EM - ONE
                              EMP2AL = EMP2AL - ONE
                              IF (N .EQ. 1) GO TO 150
                              IF (N .EQ. 2) EMP2AL = ONE
                              EMPAL = EMPAL - ONE
                              SUM = (SUM + TEMPA*EMPAL) * EMP2AL / EM
  140                      CONTINUE
                        END IF
C-------------------------------------------------------------------
C Store B(NB)
C-------------------------------------------------------------------
  150                   B(N) = TEMPA
                        IF (NB .LE. 1) THEN
                           SUM = (SUM + SUM) + TEMPA
                           GO TO 230
                        END IF
C-------------------------------------------------------------------
C Calculate and Store B(NB-1)
C-------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N)  = (EN*TEMPA) / X + TEMPB
                        IF (N .EQ. 1) GO TO 220
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
                  END IF
                  NEND = N - 2
                  IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C-------------------------------------------------------------------
                     DO 200 L = 1, NEND
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*B(N+1)) / X +B(N+2)
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
  200                CONTINUE
                  END IF
C-------------------------------------------------------------------
C Calculate B(1)
C-------------------------------------------------------------------
                  B(1) = TWO*EMPAL*B(2) / X + B(3)
  220             SUM = (SUM + SUM) + B(1)
C-------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C-------------------------------------------------------------------
  230             IF (ALPHA .NE. ZERO)
     1               SUM = SUM * FUNC(ONE+ALPHA) * (X*HALF)**(-ALPHA)
                  IF (IZE .EQ. 1) SUM = SUM * EXP(-X)
                  TEMPA = ENMTEN
                  IF (SUM .GT. ONE) TEMPA = TEMPA * SUM
                  DO 260 N = 1, NB
                     IF (B(N) .LT. TEMPA) B(N) = ZERO
                     B(N) = B(N) / SUM
  260             CONTINUE
                  RETURN
C-------------------------------------------------------------------
C Two-term ascending series for small X.
C-------------------------------------------------------------------
               ELSE
                  TEMPA = ONE
                  EMPAL = ONE + ALPHA
                  HALFX = ZERO
                  IF (X .GT. ENMTEN) HALFX = HALF * X
                  IF (ALPHA .NE. ZERO) TEMPA = HALFX**ALPHA /FUNC(EMPAL)
                  IF (IZE .EQ. 2) TEMPA = TEMPA * EXP(-X)
                  TEMPB = ZERO
                  IF ((X+ONE) .GT. ONE) TEMPB = HALFX * HALFX
                  B(1) = TEMPA + TEMPA*TEMPB / EMPAL
                  IF ((X .NE. ZERO) .AND. (B(1) .EQ. ZERO)) NCALC = 0
                  IF (NB .GT. 1) THEN
                     IF (X .EQ. ZERO) THEN
                           DO 310 N = 2, NB
                              B(N) = ZERO
  310                      CONTINUE
                        ELSE
C-------------------------------------------------------------------
C Calculate higher-order functions.
C-------------------------------------------------------------------
                           TEMPC = HALFX
                           TOVER = (ENMTEN + ENMTEN) / X
                           IF (TEMPB .NE. ZERO) TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) .AND. (NCALC .GT. N))
     1                             NCALC = N-1
  340                      CONTINUE
                     END IF
                  END IF
            END IF
         ELSE
            NCALC = MIN0(NB,0)-1
      END IF
      RETURN
C---------- Last line of RIBESL ----------
      END

CS    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
CS    REAL 
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  Mathematical constants
C----------------------------------------------------------------------
CS    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
CS   1     SQRTPI/0.9189385332046727417803297E0/,
CS   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  Machine dependent parameters
C----------------------------------------------------------------------
CS    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
CS   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
     1     XINF/1.79D308/
C----------------------------------------------------------------------
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C----------------------------------------------------------------------
CS    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
CS   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
CS   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
CS   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
CS    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
CS   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
CS   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
CS   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  Coefficients for minimax approximation over (12, INF).
C----------------------------------------------------------------------
CS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
CS   1     -5.952379913043012E-04,7.93650793500350248E-04,
CS   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
CS   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C----------------------------------------------------------------------
C  Statement functions for conversion between integer and float
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C----------------------------------------------------------------------
C  Argument is negative
C----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Argument is positive
C----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
C----------------------------------------------------------------------
C  Argument .LT. EPS
C----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C----------------------------------------------------------------------
C  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
C----------------------------------------------------------------------
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C----------------------------------------------------------------------
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C----------------------------------------------------------------------
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
C----------------------------------------------------------------------
C  Evaluate for argument .GE. 12.0,
C----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C----------------------------------------------------------------------
C  Final adjustments and return
C----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
CS900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
C ---------- Last line of GAMMA ----------
      END
