! Made into .f90 standards and optimised - Elspeth Lee - May 2016
! Added updated BLAS + LAPACK routines from LAPACK 3.6.0 when possible
! Need to figure out how to link a machine optimised version to the main limex program
! Could save alot of time.

! c-----------------------------------------------------------------------
! c
! c                          LIMEX4_2_Auxiliaries
! c
! c-----------------------------------------------------------------------
! c
! c     This is a collection of subroutines for LIMEX, the extrapolation
! c     integrator for  the solution  of linearly-implicit differential-
! c     algebraic systems of the form
! c
! c          B (t,y) * y' (t) = f (t,y)
! c
! c     with B a (n,n)-matrix of rank less or equal n.
! c
! c     This file contains  auxiliary  routines from the BLAS and LAPACK
! c     libraries  for the LIMEX  version 4.2. The whole BLAS and LAPACK
! c     libraries are  available from the netlib repository  at
! c
! c        http://cm.bell-labs.com/netlib/master/readme.html
! c
! c     or one of its mirrors.
! c
! c-----------------------------------------------------------------------
! c
! c     Written by:
! c
! c     R. Ehrig, U. Nowak,
! c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
! c     Takustrasse 7
! c     D-14195 Berlin-Dahlem
! c
! c     phone : +49-30-84185-0
! c     fax   : +49-30-84185-125
! c     e-mail: ehrig@zib.de, nowak@zib.de
! c     URL   : http://www.zib.de
! c
! c-----------------------------------------------------------------------
! c
! c            **************************************************
! c            **                                              **
! c            **    This is version 4.2 of March, 17, 2000    **
! c            **                                              **
! c            **************************************************
! c
! c-----------------------------------------------------------------------
! c
! c     Contents:
! c
! c     daxpy   from BLAS1
! c     dcopy   from BLAS1
! c     ddot    from BLAS1
! c     dnrm2   from BLAS1
! c     dscal   from BLAS1
! c     dswap   from BLAS1
! c     idamax  from BLAS1
! c
! c     dgemv   from BLAS2
! c     dger    from BLAS2
! c     dtbsv   from BLAS2
! c
! c     dgemm   from BLAS3
! c     dtrsm   from BLAS3
! c
! c     dgbtf2  from LAPACK 3.0
! c     dgbtrf  from LAPACK 3.0
! c     dgbtrs  from LAPACK 3.0
! c     dgetf2  from LAPACK 3.0
! c     dgetrf  from LAPACK 3.0
! c     dgetrs  from LAPACK 3.0
! c     dlaswp  from LAPACK 3.0
! c     dtrtrs  from LAPACK 3.0
! c     ieeeck  from LAPACK 3.0
! c     ilaenv  from LAPACK 3.0
! c     lsame   from LAPACK 3.0
! c     xerbla  from LAPACK 3.0
! c
! c-----------------------------------------------------------------------
! c
! c     This is the original file daxpy.f from the BLAS1 routines.
! c
! c-----------------------------------------------------------------------
! c
SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION :: DA
  INTEGER :: INCX,INCY,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION :: DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  INTEGER :: I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     ..
  IF (N.LE.0) RETURN
  IF (DA.EQ.0.0d0) RETURN
  IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
    M = MOD(N,4)
    IF (M.NE.0) THEN
      DO I = 1,M
        DY(I) = DY(I) + DA*DX(I)
      END DO
    END IF
    IF (N.LT.4) RETURN
    MP1 = M + 1
    DO I = MP1,N,4
      DY(I) = DY(I) + DA*DX(I)
      DY(I+1) = DY(I+1) + DA*DX(I+1)
      DY(I+2) = DY(I+2) + DA*DX(I+2)
      DY(I+3) = DY(I+3) + DA*DX(I+3)
    END DO
  ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DY(IY) = DY(IY) + DA*DX(IX)
      IX = IX + INCX
      IY = IY + INCY
    END DO
  END IF

RETURN
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dcopy.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*),DY(*)
!    ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!    ..
  IF (N.LE.0) RETURN
  IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
    M = MOD(N,7)
    IF (M.NE.0) THEN
      DO I = 1,M
        DY(I) = DX(I)
      END DO
      IF (N.LT.7) RETURN
    END IF
    MP1 = M + 1
    DO I = MP1,N,7
      DY(I) = DX(I)
      DY(I+1) = DX(I+1)
      DY(I+2) = DX(I+2)
      DY(I+3) = DX(I+3)
      DY(I+4) = DX(I+4)
      DY(I+5) = DX(I+5)
      DY(I+6) = DX(I+6)
    END DO
  ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DY(IY) = DX(IX)
      IX = IX + INCX
      IY = IY + INCY
    END DO
  END IF

RETURN
END

!-----------------------------------------------------------------------
!
!     This is the original file ddot.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
  DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*),DY(*)
!    ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  DOUBLE PRECISION DTEMP
  INTEGER I,IX,IY,M,MP1
!    ..
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     ..
  DDOT = 0.0d0
  DTEMP = 0.0d0
  IF (N.LE.0) RETURN
  IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
    M = MOD(N,5)
    IF (M.NE.0) THEN
      DO I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
      END DO
      IF (N.LT.5) THEN
        DDOT=DTEMP
        RETURN
      END IF
    END IF
    MP1 = M + 1
    DO I = MP1,N,5
      DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
        & DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
    END DO
  ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DTEMP = DTEMP + DX(IX)*DY(IY)
      IX = IX + INCX
      IY = IY + INCY
    END DO
  END IF

  DDOT = DTEMP

RETURN
END

!-----------------------------------------------------------------------
!
!     This is the original file dnrm2.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER INCX,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION ONE,ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
  INTEGER IX
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC ABS,SQRT
!     ..
  IF (N.LT.1 .OR. INCX.LT.1) THEN
    NORM = ZERO
  ELSE IF (N.EQ.1) THEN
    NORM = ABS(X(1))
  ELSE
    SCALE = ZERO
    SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
    DO 10 IX = 1,1 + (N-1)*INCX,INCX
      IF (X(IX).NE.ZERO) THEN
        ABSXI = ABS(X(IX))
        IF (SCALE.LT.ABSXI) THEN
          SSQ = ONE + SSQ* (SCALE/ABSXI)**2
          SCALE = ABSXI
        ELSE
          SSQ = SSQ + (ABSXI/SCALE)**2
        END IF
      END IF
10  CONTINUE
    NORM = SCALE*SQRT(SSQ)
END IF

DNRM2 = NORM

RETURN
END

!-----------------------------------------------------------------------
!
!     This is the original file dscal.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION DA
  INTEGER INCX,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  INTEGER I,M,MP1,NINCX
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     ..
  IF (N.LE.0 .OR. INCX.LE.0) RETURN
  IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
    M = MOD(N,5)
    IF (M.NE.0) THEN
      DO I = 1,M
        DX(I) = DA*DX(I)
      END DO
      IF (N.LT.5) RETURN
    END IF
    MP1 = M + 1
    DO I = MP1,N,5
      DX(I) = DA*DX(I)
      DX(I+1) = DA*DX(I+1)
      DX(I+2) = DA*DX(I+2)
      DX(I+3) = DA*DX(I+3)
      DX(I+4) = DA*DX(I+4)
    END DO
  ELSE
!
!        code for increment not equal to 1
!
    NINCX = N*INCX
    DO I = 1,NINCX,INCX
      DX(I) = DA*DX(I)
    END DO
  END IF

RETURN
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dswap.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  DOUBLE PRECISION DTEMP
  INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     ..
  IF (N.LE.0) RETURN
  IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
    M = MOD(N,3)
    IF (M.NE.0) THEN
      DO I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
      END DO
      IF (N.LT.3) RETURN
    END IF
    MP1 = M + 1
    DO I = MP1,N,3
      DTEMP = DX(I)
      DX(I) = DY(I)
      DY(I) = DTEMP
      DTEMP = DX(I+1)
      DX(I+1) = DY(I+1)
      DY(I+1) = DTEMP
      DTEMP = DX(I+2)
      DX(I+2) = DY(I+2)
      DY(I+2) = DTEMP
    END DO
  ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
      DTEMP = DX(IX)
      DX(IX) = DY(IY)
      DY(IY) = DTEMP
      IX = IX + INCX
      IY = IY + INCY
    END DO
  END IF

RETURN
END

!
!-----------------------------------------------------------------------
!
!     This is the original file idamax.f from the BLAS1 routines.
!
!-----------------------------------------------------------------------
!
INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!  -- Reference BLAS level1 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  INTEGER INCX,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  DOUBLE PRECISION DMAX
  INTEGER I,IX
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC DABS
!     ..
  IDAMAX = 0
  IF (N.LT.1 .OR. INCX.LE.0) RETURN
  IDAMAX = 1
  IF (N.EQ.1) RETURN
  IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
    DMAX = DABS(DX(1))
    DO I = 2,N
      IF (DABS(DX(I)).GT.DMAX) THEN
        IDAMAX = I
        DMAX = DABS(DX(I))
      END IF
    END DO
  ELSE
!
!        code for increment not equal to 1
!
    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO I = 2,N
      IF (DABS(DX(IX)).GT.DMAX) THEN
        IDAMAX = I
        DMAX = DABS(DX(IX))
      END IF
      IX = IX + INCX
    END DO
  END IF

RETURN
END

!-----------------------------------------------------------------------
!
!     This is the original file dgemv.f from the BLAS2 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION ALPHA,BETA
  INTEGER INCX,INCY,LDA,M,N
  CHARACTER TRANS
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION ONE,ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
  INFO = 0
  IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
    &   .NOT.LSAME(TRANS,'C')) THEN
    INFO = 1
  ELSE IF (M.LT.0) THEN
    INFO = 2
  ELSE IF (N.LT.0) THEN
    INFO = 3
  ELSE IF (LDA.LT.MAX(1,M)) THEN
    INFO = 6
  ELSE IF (INCX.EQ.0) THEN
    INFO = 8
  ELSE IF (INCY.EQ.0) THEN
    INFO = 11
  END IF
  IF (INFO.NE.0) THEN
    CALL XERBLA('DGEMV ',INFO)
    RETURN
  END IF
!
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
    & ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
  IF (LSAME(TRANS,'N')) THEN
    LENX = N
    LENY = M
  ELSE
    LENX = M
    LENY = N
  END IF
  IF (INCX.GT.0) THEN
    KX = 1
  ELSE
    KX = 1 - (LENX-1)*INCX
  END IF
  IF (INCY.GT.0) THEN
    KY = 1
  ELSE
    KY = 1 - (LENY-1)*INCY
  END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
  IF (BETA.NE.ONE) THEN
    IF (INCY.EQ.1) THEN
      IF (BETA.EQ.ZERO) THEN
        DO 10 I = 1,LENY
          Y(I) = ZERO
10      CONTINUE
      ELSE
        DO 20 I = 1,LENY
          Y(I) = BETA*Y(I)
20      CONTINUE
      END IF
    ELSE
      IY = KY
      IF (BETA.EQ.ZERO) THEN
        DO 30 I = 1,LENY
          Y(IY) = ZERO
          IY = IY + INCY
30      CONTINUE
      ELSE
        DO 40 I = 1,LENY
          Y(IY) = BETA*Y(IY)
          IY = IY + INCY
40      CONTINUE
      END IF
    END IF
  END IF
  IF (ALPHA.EQ.ZERO) RETURN
  IF (LSAME(TRANS,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
    JX = KX
    IF (INCY.EQ.1) THEN
      DO 60 J = 1,N
        TEMP = ALPHA*X(JX)
        DO 50 I = 1,M
          Y(I) = Y(I) + TEMP*A(I,J)
50      CONTINUE
        JX = JX + INCX
60    CONTINUE
    ELSE
      DO 80 J = 1,N
        TEMP = ALPHA*X(JX)
        IY = KY
        DO 70 I = 1,M
          Y(IY) = Y(IY) + TEMP*A(I,J)
          IY = IY + INCY
70      CONTINUE
        JX = JX + INCX
80    CONTINUE
    END IF
  ELSE
!
!        Form  y := alpha*A**T*x + y.
!
    JY = KY
    IF (INCX.EQ.1) THEN
      DO 100 J = 1,N
        TEMP = ZERO
        DO 90 I = 1,M
          TEMP = TEMP + A(I,J)*X(I)
90      CONTINUE
        Y(JY) = Y(JY) + ALPHA*TEMP
        JY = JY + INCY
100   CONTINUE
    ELSE
      DO 120 J = 1,N
        TEMP = ZERO
        IX = KX
        DO 110 I = 1,M
          TEMP = TEMP + A(I,J)*X(IX)
          IX = IX + INCX
110     CONTINUE
        Y(JY) = Y(JY) + ALPHA*TEMP
        JY = JY + INCY
120   CONTINUE
    END IF
  END IF

!
!     End of DGEMV .
!
RETURN
END

!-----------------------------------------------------------------------
!
!     This is the original file dger.f from the BLAS2 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION ALPHA
  INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION ZERO
  PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
  INFO = 0
  IF (M.LT.0) THEN
    INFO = 1
  ELSE IF (N.LT.0) THEN
    INFO = 2
  ELSE IF (INCX.EQ.0) THEN
    INFO = 5
  ELSE IF (INCY.EQ.0) THEN
    INFO = 7
  ELSE IF (LDA.LT.MAX(1,M)) THEN
    INFO = 9
  END IF
  IF (INFO.NE.0) THEN
    CALL XERBLA('DGER  ',INFO)
    RETURN
  END IF
!
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  IF (INCY.GT.0) THEN
    JY = 1
  ELSE
    JY = 1 - (N-1)*INCY
  END IF
  IF (INCX.EQ.1) THEN
    DO 20 J = 1,N
      IF (Y(JY).NE.ZERO) THEN
        TEMP = ALPHA*Y(JY)
        DO 10 I = 1,M
          A(I,J) = A(I,J) + X(I)*TEMP
10      CONTINUE
      END IF
      JY = JY + INCY
20  CONTINUE
  ELSE
    IF (INCX.GT.0) THEN
      KX = 1
    ELSE
      KX = 1 - (M-1)*INCX
    END IF
    DO 40 J = 1,N
      IF (Y(JY).NE.ZERO) THEN
        TEMP = ALPHA*Y(JY)
        IX = KX
        DO 30 I = 1,M
          A(I,J) = A(I,J) + X(IX)*TEMP
          IX = IX + INCX
30      CONTINUE
      END IF
      JY = JY + INCY
40  CONTINUE
  END IF
!
!
!     End of DGER  .
!
RETURN
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dtbsv.f from the BLAS2 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER INCX,K,LDA,N
  CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION ZERO
  PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
  LOGICAL NOUNIT
!     ..
!     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MAX,MIN
!     ..
!
!     Test the input parameters.
!
  INFO = 0
  IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
    INFO = 1
  ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
    & .NOT.LSAME(TRANS,'C')) THEN
    INFO = 2
  ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
    INFO = 3
  ELSE IF (N.LT.0) THEN
    INFO = 4
  ELSE IF (K.LT.0) THEN
    INFO = 5
  ELSE IF (LDA.LT. (K+1)) THEN
    INFO = 7
  ELSE IF (INCX.EQ.0) THEN
    INFO = 9
  END IF
  IF (INFO.NE.0) THEN
    CALL XERBLA('DTBSV ',INFO)
    RETURN
  END IF
!
!     Quick return if possible.
!
  IF (N.EQ.0) RETURN
!
  NOUNIT = LSAME(DIAG,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  IF (INCX.LE.0) THEN
    KX = 1 - (N-1)*INCX
  ELSE IF (INCX.NE.1) THEN
    KX = 1
  END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
  IF (LSAME(TRANS,'N')) THEN
!
!        Form  x := inv( A )*x.
!
    IF (LSAME(UPLO,'U')) THEN
      KPLUS1 = K + 1
      IF (INCX.EQ.1) THEN
        DO 20 J = N,1,-1
          IF (X(J).NE.ZERO) THEN
            L = KPLUS1 - J
            IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
              TEMP = X(J)
              DO 10 I = J - 1,MAX(1,J-K),-1
                X(I) = X(I) - TEMP*A(L+I,J)
10            CONTINUE
            END IF
20      CONTINUE
          ELSE
            KX = KX + (N-1)*INCX
            JX = KX
            DO 40 J = N,1,-1
              KX = KX - INCX
              IF (X(JX).NE.ZERO) THEN
                IX = KX
                L = KPLUS1 - J
                IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                  TEMP = X(JX)
                  DO 30 I = J - 1,MAX(1,J-K),-1
                    X(IX) = X(IX) - TEMP*A(L+I,J)
                    IX = IX - INCX
30                CONTINUE
                END IF
              JX = JX - INCX
40          CONTINUE
          END IF
        ELSE
          IF (INCX.EQ.1) THEN
            DO 60 J = 1,N
              IF (X(J).NE.ZERO) THEN
                L = 1 - J
                IF (NOUNIT) X(J) = X(J)/A(1,J)
                  TEMP = X(J)
                  DO 50 I = J + 1,MIN(N,J+K)
                    X(I) = X(I) - TEMP*A(L+I,J)
50                CONTINUE
                END IF
60          CONTINUE
              ELSE
                JX = KX
                DO 80 J = 1,N
                  KX = KX + INCX
                  IF (X(JX).NE.ZERO) THEN
                    IX = KX
                    L = 1 - J
                    IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                    TEMP = X(JX)
                    DO 70 I = J + 1,MIN(N,J+K)
                      X(IX) = X(IX) - TEMP*A(L+I,J)
                      IX = IX + INCX
70                  CONTINUE
                  END IF
                  JX = JX + INCX
80              CONTINUE
              END IF
            END IF
          ELSE
!
!        Form  x := inv( A**T)*x.
!
          IF (LSAME(UPLO,'U')) THEN
            KPLUS1 = K + 1
            IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                TEMP = X(J)
                L = KPLUS1 - J
                DO 90 I = MAX(1,J-K),J - 1
                  TEMP = TEMP - A(L+I,J)*X(I)
90              CONTINUE
                IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                  X(J) = TEMP
100           CONTINUE
          ELSE
            JX = KX
            DO 120 J = 1,N
              TEMP = X(JX)
              IX = KX
              L = KPLUS1 - J
              DO 110 I = MAX(1,J-K),J - 1
                TEMP = TEMP - A(L+I,J)*X(IX)
                IX = IX + INCX
110           CONTINUE
                IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                  X(JX) = TEMP
                  JX = JX + INCX
                  IF (J.GT.K) KX = KX + INCX
120        CONTINUE
          END IF
        ELSE
          IF (INCX.EQ.1) THEN
            DO 140 J = N,1,-1
              TEMP = X(J)
              L = 1 - J
              DO 130 I = MIN(N,J+K),J + 1,-1
                TEMP = TEMP - A(L+I,J)*X(I)
130           CONTINUE
              IF (NOUNIT) TEMP = TEMP/A(1,J)
                X(J) = TEMP
140         CONTINUE
          ELSE
            KX = KX + (N-1)*INCX
            JX = KX
            DO 160 J = N,1,-1
              TEMP = X(JX)
              IX = KX
              L = 1 - J
              DO 150 I = MIN(N,J+K),J + 1,-1
                TEMP = TEMP - A(L+I,J)*X(IX)
                IX = IX - INCX
150           CONTINUE
              IF (NOUNIT) TEMP = TEMP/A(1,J)
                X(JX) = TEMP
                JX = JX - INCX
                IF ((N-J).GE.K) KX = KX - INCX
160         CONTINUE
          END IF
        END IF
      END IF
!
RETURN
!
!     End of DTBSV .
!
END
!
!-----------------------------------------------------------------------
!
!     This is the original file dgemm.f from the BLAS3 routines.
!
!-----------------------------------------------------------------------
!

SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  -- Reference BLAS level3 routine (version 3.6.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION ALPHA,BETA
  INTEGER K,LDA,LDB,LDC,M,N
  CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC MAX
!      ..
!     .. Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
  LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
  DOUBLE PRECISION ONE,ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
  NOTA = LSAME(TRANSA,'N')
  NOTB = LSAME(TRANSB,'N')
  IF (NOTA) THEN
    NROWA = M
    NCOLA = K
  ELSE
    NROWA = K
    NCOLA = M
  END IF
  IF (NOTB) THEN
    NROWB = K
  ELSE
    NROWB = N
  END IF
!
!     Test the input parameters.
!
  INFO = 0
  IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
    & (.NOT.LSAME(TRANSA,'T'))) THEN
    INFO = 1
  ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
    & (.NOT.LSAME(TRANSB,'T'))) THEN
    INFO = 2
  ELSE IF (M.LT.0) THEN
    INFO = 3
  ELSE IF (N.LT.0) THEN
    INFO = 4
  ELSE IF (K.LT.0) THEN
    INFO = 5
  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
    INFO = 8
  ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
    INFO = 10
  ELSE IF (LDC.LT.MAX(1,M)) THEN
    INFO = 13
  END IF
  IF (INFO.NE.0) THEN
    CALL XERBLA('DGEMM ',INFO)
    RETURN
  END IF
!
!     Quick return if possible.
!
  IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
    &  (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
  IF (ALPHA.EQ.ZERO) THEN
    IF (BETA.EQ.ZERO) THEN
      DO 20 J = 1,N
        DO 10 I = 1,M
          C(I,J) = ZERO
10      CONTINUE
20    CONTINUE
    ELSE
      DO 40 J = 1,N
        DO 30 I = 1,M
          C(I,J) = BETA*C(I,J)
30      CONTINUE
40    CONTINUE
    END IF
    RETURN
  END IF
!
!     Start the operations.
!
  IF (NOTB) THEN
    IF (NOTA) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
      DO 90 J = 1,N
        IF (BETA.EQ.ZERO) THEN
          DO 50 I = 1,M
            C(I,J) = ZERO
50        CONTINUE
        ELSE IF (BETA.NE.ONE) THEN
          DO 60 I = 1,M
            C(I,J) = BETA*C(I,J)
60        CONTINUE
        END IF
        DO 80 L = 1,K
          TEMP = ALPHA*B(L,J)
          DO 70 I = 1,M
            C(I,J) = C(I,J) + TEMP*A(I,L)
70        CONTINUE
80      CONTINUE
90    CONTINUE
    ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
      DO 120 J = 1,N
        DO 110 I = 1,M
          TEMP = ZERO
          DO 100 L = 1,K
            TEMP = TEMP + A(L,I)*B(L,J)
100       CONTINUE
          IF (BETA.EQ.ZERO) THEN
            C(I,J) = ALPHA*TEMP
          ELSE
            C(I,J) = ALPHA*TEMP + BETA*C(I,J)
          END IF
110     CONTINUE
120   CONTINUE
    END IF
  ELSE
    IF (NOTA) THEN
!
!           Form  C := alpha*A*B**T + beta*C
!
    DO 170 J = 1,N
      IF (BETA.EQ.ZERO) THEN
        DO 130 I = 1,M
          C(I,J) = ZERO
130     CONTINUE
      ELSE IF (BETA.NE.ONE) THEN
        DO 140 I = 1,M
          C(I,J) = BETA*C(I,J)
140     CONTINUE
      END IF
      DO 160 L = 1,K
        TEMP = ALPHA*B(J,L)
        DO 150 I = 1,M
          C(I,J) = C(I,J) + TEMP*A(I,L)
150     CONTINUE
160   CONTINUE
170 CONTINUE
    ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
    DO 200 J = 1,N
      DO 190 I = 1,M
        TEMP = ZERO
        DO 180 L = 1,K
          TEMP = TEMP + A(L,I)*B(J,L)
180     CONTINUE
        IF (BETA.EQ.ZERO) THEN
          C(I,J) = ALPHA*TEMP
        ELSE
          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
        END IF
190   CONTINUE
200 CONTINUE
    END IF
  END IF
!
RETURN
!
!     End of DGEMM .
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dtrsm.f from the BLAS3 routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine (version 3.4.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  DOUBLE PRECISION ALPHA
  INTEGER LDA,LDB,M,N
  CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL XERBLA
!    ..
!     .. Intrinsic Functions ..
  INTRINSIC MAX
!     ..
!     .. Local Scalars ..
  DOUBLE PRECISION TEMP
  INTEGER I,INFO,J,K,NROWA
  LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
  DOUBLE PRECISION ONE,ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!
!     Test the input parameters.
!
  LSIDE = LSAME(SIDE,'L')
  IF (LSIDE) THEN
    NROWA = M
  ELSE
    NROWA = N
  END IF
  NOUNIT = LSAME(DIAG,'N')
  UPPER = LSAME(UPLO,'U')
!
  INFO = 0
  IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
    INFO = 1
  ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
    INFO = 2
  ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
    & (.NOT.LSAME(TRANSA,'T')) .AND. &
    & (.NOT.LSAME(TRANSA,'C'))) THEN
    INFO = 3
  ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
    INFO = 4
  ELSE IF (M.LT.0) THEN
    INFO = 5
  ELSE IF (N.LT.0) THEN
    INFO = 6
  ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
    INFO = 9
  ELSE IF (LDB.LT.MAX(1,M)) THEN
    INFO = 11
  END IF
  IF (INFO.NE.0) THEN
    CALL XERBLA('DTRSM ',INFO)
    RETURN
  END IF
!
!     Quick return if possible.
!
  IF (M.EQ.0 .OR. N.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
  IF (ALPHA.EQ.ZERO) THEN
    DO 20 J = 1,N
      DO 10 I = 1,M
        B(I,J) = ZERO
10    CONTINUE
20  CONTINUE
    RETURN
  END IF
!
!     Start the operations.
!
  IF (LSIDE) THEN
    IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
    IF (UPPER) THEN
      DO 60 J = 1,N
        IF (ALPHA.NE.ONE) THEN
          DO 30 I = 1,M
            B(I,J) = ALPHA*B(I,J)
30        CONTINUE
        END IF
        DO 50 K = M,1,-1
          IF (B(K,J).NE.ZERO) THEN
            IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
              DO 40 I = 1,K - 1
                B(I,J) = B(I,J) - B(K,J)*A(I,K)
40            CONTINUE
            END IF
50      CONTINUE
60    CONTINUE
    ELSE
      DO 100 J = 1,N
        IF (ALPHA.NE.ONE) THEN
          DO 70 I = 1,M
            B(I,J) = ALPHA*B(I,J)
70        CONTINUE
        END IF
        DO 90 K = 1,M
          IF (B(K,J).NE.ZERO) THEN
            IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
              DO 80 I = K + 1,M
                B(I,J) = B(I,J) - B(K,J)*A(I,K)
80            CONTINUE
            END IF
90      CONTINUE
100  CONTINUE
    END IF
  ELSE
!
!           Form  B := alpha*inv( A**T )*B.
!
    IF (UPPER) THEN
      DO 130 J = 1,N
        DO 120 I = 1,M
          TEMP = ALPHA*B(I,J)
          DO 110 K = 1,I - 1
            TEMP = TEMP - A(K,I)*B(K,J)
110       CONTINUE
          IF (NOUNIT) TEMP = TEMP/A(I,I)
            B(I,J) = TEMP
120     CONTINUE
130   CONTINUE
          ELSE
            DO 160 J = 1,N
              DO 150 I = M,1,-1
                TEMP = ALPHA*B(I,J)
                DO 140 K = I + 1,M
                  TEMP = TEMP - A(K,I)*B(K,J)
140             CONTINUE
                IF (NOUNIT) TEMP = TEMP/A(I,I)
                  B(I,J) = TEMP
150           CONTINUE
160         CONTINUE
               END IF
          END IF
    ELSE
      IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
        IF (UPPER) THEN
          DO 210 J = 1,N
            IF (ALPHA.NE.ONE) THEN
              DO 170 I = 1,M
                B(I,J) = ALPHA*B(I,J)
170           CONTINUE
            END IF
            DO 190 K = 1,J - 1
              IF (A(K,J).NE.ZERO) THEN
              DO 180 I = 1,M
                B(I,J) = B(I,J) - A(K,J)*B(I,K)
180           CONTINUE
            END IF
190         CONTINUE
            IF (NOUNIT) THEN
              TEMP = ONE/A(J,J)
              DO 200 I = 1,M
                B(I,J) = TEMP*B(I,J)
200           CONTINUE
            END IF
210       CONTINUE
        ELSE
          DO 260 J = N,1,-1
            IF (ALPHA.NE.ONE) THEN
              DO 220 I = 1,M
                B(I,J) = ALPHA*B(I,J)
220           CONTINUE
            END IF
            DO 240 K = J + 1,N
              IF (A(K,J).NE.ZERO) THEN
                DO 230 I = 1,M
                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
230             CONTINUE
              END IF
240         CONTINUE
            IF (NOUNIT) THEN
              TEMP = ONE/A(J,J)
              DO 250 I = 1,M
                B(I,J) = TEMP*B(I,J)
250           CONTINUE
            END IF
260       CONTINUE
        END IF
      ELSE
!
!           Form  B := alpha*B*inv( A**T ).
!
        IF (UPPER) THEN
          DO 310 K = N,1,-1
            IF (NOUNIT) THEN
              TEMP = ONE/A(K,K)
              DO 270 I = 1,M
                B(I,K) = TEMP*B(I,K)
270           CONTINUE
            END IF
            DO 290 J = 1,K - 1
              IF (A(J,K).NE.ZERO) THEN
              TEMP = A(J,K)
              DO 280 I = 1,M
                B(I,J) = B(I,J) - TEMP*B(I,K)
280           CONTINUE
              END IF
290         CONTINUE
            IF (ALPHA.NE.ONE) THEN
              DO 300 I = 1,M
                B(I,K) = ALPHA*B(I,K)
300           CONTINUE
            END IF
310       CONTINUE
        ELSE
          DO 360 K = 1,N
            IF (NOUNIT) THEN
              TEMP = ONE/A(K,K)
              DO 320 I = 1,M
                B(I,K) = TEMP*B(I,K)
320           CONTINUE
            END IF
            DO 340 J = K + 1,N
              IF (A(J,K).NE.ZERO) THEN
                TEMP = A(J,K)
                DO 330 I = 1,M
                  B(I,J) = B(I,J) - TEMP*B(I,K)
330             CONTINUE
              END IF
340         CONTINUE
            IF (ALPHA.NE.ONE) THEN
              DO 350 I = 1,M
                B(I,K) = ALPHA*B(I,K)
350           CONTINUE
            END IF
360       CONTINUE
        END IF
      END IF
    END IF
!
RETURN
!
!     End of DTRSM .
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dgbtf2.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
  INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
  INTEGER            IDAMAX
  EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
  KV = KU + KL
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
   INFO = -1
  ELSE IF( N.LT.0 ) THEN
   INFO = -2
  ELSE IF( KL.LT.0 ) THEN
   INFO = -3
  ELSE IF( KU.LT.0 ) THEN
   INFO = -4
  ELSE IF( LDAB.LT.KL+KV+1 ) THEN
   INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
   CALL XERBLA( 'DGBTF2', -INFO )
   RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
  DO 20 J = KU + 2, MIN( KV, N )
    DO 10 I = KV - J + 2, KL
      AB( I, J ) = ZERO
10  CONTINUE
20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
JU = 1
!
  DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
    IF( J+KV.LE.N ) THEN
      DO 30 I = 1, KL
        AB( I, J+KV ) = ZERO
30    CONTINUE
    END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
    KM = MIN( KL, M-J )
    JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
    IPIV( J ) = JP + J - 1
    IF( AB( KV+JP, J ).NE.ZERO ) THEN
      JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
      IF( JP.NE.1 ) &
        & CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
        & AB( KV+1, J ), LDAB-1 )
!
      IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
        CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
        IF( JU.GT.J ) &
          &   CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
          &   AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
          &   LDAB-1 )
        END IF
      ELSE

!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
        IF( INFO.EQ.0 ) INFO = J
      END IF
40 CONTINUE

RETURN
!
!     End of DGBTF2
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dgbtrf.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
  INTEGER            NBMAX, LDWORK
  PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
  INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
    &  JU, K2, KM, KV, NB, NW
  DOUBLE PRECISION   TEMP
!     ..
!     .. Local Arrays ..
  DOUBLE PRECISION   WORK13( LDWORK, NBMAX ), &
    & WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
  INTEGER            IDAMAX, ILAENV
  EXTERNAL           IDAMAX, ILAENV
!     ..
!     .. External Subroutines ..
  EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL, &
    &                   DSWAP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
  KV = KU + KL
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
    INFO = -1
  ELSE IF( N.LT.0 ) THEN
    INFO = -2
  ELSE IF( KL.LT.0 ) THEN
    INFO = -3
  ELSE IF( KU.LT.0 ) THEN
    INFO = -4
  ELSE IF( LDAB.LT.KL+KV+1 ) THEN
    INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DGBTRF', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment
!
  NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
  NB = MIN( NB, NBMAX )
!
  IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!       Use unblocked code
!
    CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
  ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
    DO 20 J = 1, NB
      DO 10 I = 1, J - 1
      WORK13( I, J ) = ZERO
10    CONTINUE
20  CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
    DO 40 J = 1, NB
      DO 30 I = J + 1, NB
        WORK31( I, J ) = ZERO
30    CONTINUE
40  CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
    DO 60 J = KU + 2, MIN( KV, N )
      DO 50 I = KV - J + 2, KL
        AB( I, J ) = ZERO
50    CONTINUE
60  CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
    JU = 1
!
    DO 180 J = 1, MIN( M, N ), NB
      JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
      I2 = MIN( KL-JB, M-J-JB+1 )
      I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
      DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
        IF( JJ+KV.LE.N ) THEN
          DO 70 I = 1, KL
            AB( I, JJ+KV ) = ZERO
70        CONTINUE
        END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
        KM = MIN( KL, M-JJ )
        JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
        IPIV( JJ ) = JP + JJ - J
          IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
            IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
              IF( JP+JJ-1.LT.J+KL ) THEN

                CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                  & AB( KV+JP+JJ-J, J ), LDAB-1 )
              ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                  & WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                  & AB( KV+JP, JJ ), LDAB-1 )
              END IF
            END IF
!
!                 Compute multipliers
!
            CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
              & 1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
            JM = MIN( JU, J+JB-1 )
            IF( JM.GT.JJ ) &
              & CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
              & AB( KV, JJ+1 ), LDAB-1, &
              & AB( KV+1, JJ+1 ), LDAB-1 )
            ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
              IF( INFO.EQ.0 ) INFO = JJ
            END IF
!
!              Copy current column of A31 into the work array WORK31
!
            NW = MIN( JJ-J+1, I3 )
            IF( NW.GT.0 ) &
            & CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
            &          WORK31( 1, JJ-J+1 ), 1 )
80    CONTINUE
      IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
        J2 = MIN( JU-J+1, KV ) - JB
        J3 = MAX( 0, JU-J-KV+1 )
!
!              Use DLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
        CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
          & IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
        DO 90 I = J, J + JB - 1
          IPIV( I ) = IPIV( I ) + J - 1
90      CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
        K2 = J - 1 + JB + J2
        DO 110 I = 1, J3
          JJ = K2 + I
          DO 100 II = J + I - 1, J + JB - 1
            IP = IPIV( II )
            IF( IP.NE.II ) THEN
              TEMP = AB( KV+1+II-JJ, JJ )
              AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
              AB( KV+1+IP-JJ, JJ ) = TEMP
            END IF
100       CONTINUE
110     CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
         IF( J2.GT.0 ) THEN
!
!                 Update A12
!
           CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
              & JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
              &           AB( KV+1-JB, J+JB ), LDAB-1 )
!
            IF( I2.GT.0 ) THEN
!
!                   Update A22
!
              CALL DGEMM( 'No transpose', 'No transpose', I2, J2, &
              &            JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
              &           AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
             &          AB( KV+1, J+JB ), LDAB-1 )
            END IF
!
            IF( I3.GT.0 ) THEN
!
!                    Update A32
!
             CALL DGEMM( 'No transpose', 'No transpose', I3, J2, &
              &             JB, -ONE, WORK31, LDWORK, &
              &            AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
              &             AB( KV+KL+1-JB, J+JB ), LDAB-1 )
            END IF
         END IF
!
         IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
          DO 130 JJ = 1, J3
            DO 120 II = JJ, JB
              WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
120        CONTINUE
130       CONTINUE
!
!                 Update A13 in the work array
!
          CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
            & JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
            & WORK13, LDWORK )
!
            IF( I2.GT.0 ) THEN
!
!                    Update A23
!
          CALL DGEMM( 'No transpose', 'No transpose', I2, J3, &
           & JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
           & WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
           & LDAB-1 )
        END IF
!
        IF( I3.GT.0 ) THEN
!
!                    Update A33
!
          CALL DGEMM( 'No transpose', 'No transpose', I3, J3, &
           & JB, -ONE, WORK31, LDWORK, WORK13, &
           & LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
        END IF
!
!                 Copy the lower triangle of A13 back into place
!
        DO 150 JJ = 1, J3
          DO 140 II = JJ, JB
            AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
140       CONTINUE
150    CONTINUE
     END IF
    ELSE
!
!              Adjust the pivot indices.
!
      DO 160 I = J, J + JB - 1
        IPIV( I ) = IPIV( I ) + J - 1
160   CONTINUE
    END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
    DO 170 JJ = J + JB - 1, J, -1
      JP = IPIV( JJ ) - JJ + 1
      IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
        IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
          CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
            & AB( KV+JP+JJ-J, J ), LDAB-1 )
          ELSE
!
!                    The interchange does affect A31
!
          CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
            & WORK31( JP+JJ-J-KL, 1 ), LDWORK )
        END IF
      END IF
!
!              Copy the current column of A31 back into place
!
     NW = MIN( I3, JJ-J+1 )
        IF( NW.GT.0 ) &
         &  CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
         &  AB( KV+KL+1-JJ+J, JJ ), 1 )
170   CONTINUE
180  CONTINUE
  END IF
!
RETURN
!
!     End of DGBTRF
!
END


!
!-----------------------------------------------------------------------
!
!     This is the original file dgbtrs.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
&                   INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  CHARACTER          TRANS
  INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
!    ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LNOTI, NOTRAN
  INTEGER            I, J, KD, L, LM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
    &  LSAME( TRANS, 'C' ) ) THEN
    INFO = -1
  ELSE IF( N.LT.0 ) THEN
    INFO = -2
  ELSE IF( KL.LT.0 ) THEN
    INFO = -3
  ELSE IF( KU.LT.0 ) THEN
    INFO = -4
  ELSE IF( NRHS.LT.0 ) THEN
    INFO = -5
  ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
    INFO = -7
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
    INFO = -10
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DGBTRS', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
  KD = KU + KL + 1
  LNOTI = KL.GT.0
!
  IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
    IF( LNOTI ) THEN
      DO 10 J = 1, N - 1
        LM = MIN( KL, N-J )
        L = IPIV( J )
        IF( L.NE.J ) &
&            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
         CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
&                    LDB, B( J+1, 1 ), LDB )
10    CONTINUE
    END IF
!
    DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
      CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
&                 AB, LDAB, B( 1, I ), 1 )
20  CONTINUE
!
  ELSE
!
!        Solve A**T*X = B.
!
    DO 30 I = 1, NRHS
!
!           Solve U**T*X = B, overwriting B with X.
!
      CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
&                  LDAB, B( 1, I ), 1 )
30  CONTINUE
!
!        Solve L**T*X = B, overwriting B with X.
!
    IF( LNOTI ) THEN
      DO 40 J = N - 1, 1, -1
        LM = MIN( KL, N-J )
        CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
&                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
         L = IPIV( J )
         IF( L.NE.J ) &
&            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
40   CONTINUE
   END IF
END IF
RETURN
!
!     End of DGBTRS
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dgetrf.f from the LAPACK routines.
!
!-----------------------------------------------------------------------

SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DGETRF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
  INTEGER            ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
    INFO = -1
  ELSE IF( N.LT.0 ) THEN
    INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
    INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DGETRF', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
  IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
    CALL DGETRF2( M, N, A, LDA, IPIV, INFO )
  ELSE
!
!        Use blocked code.
!
    DO 20 J = 1, MIN( M, N ), NB
      JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
      CALL DGETRF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
      IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
&         INFO = IINFO + J - 1
      DO 10 I = J, MIN( M, J+JB-1 )
         IPIV( I ) = J - 1 + IPIV( I )
10    CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
      CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
      IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
        CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,&
&                      IPIV, 1 )
!
!              Compute block row of U.
!
        CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
&                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
&                     LDA )
        IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
          CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
&                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
&                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
&                        LDA )
        END IF
      END IF
20    CONTINUE
END IF
RETURN
!
!     End of DGETRF
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dgetrf2.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!

RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK computational routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * )
!    ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!    ..
!     .. Local Scalars ..
  DOUBLE PRECISION   SFMIN, TEMP
  INTEGER            I, IINFO, N1, N2
!     ..
!     .. External Functions ..
  DOUBLE PRECISION   DLAMCH
  INTEGER            IDAMAX
  EXTERNAL           DLAMCH, IDAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DSCAL, DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  IF( M.LT.0 ) THEN
    INFO = -1
  ELSE IF( N.LT.0 ) THEN
    INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
    INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DGETRF2', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

  IF ( M.EQ.1 ) THEN
!
!        Use unblocked code for one row case
!        Just need to handle IPIV and INFO
!
    IPIV( 1 ) = 1
    IF ( A(1,1).EQ.ZERO ) INFO = 1
!
    ELSE IF( N.EQ.1 ) THEN
!
!        Use unblocked code for one column case
!
!
!        Compute machine safe minimum
!
    SFMIN = DLAMCH('S')
!
!        Find pivot and test for singularity
!
    I = IDAMAX( M, A( 1, 1 ), 1 )
    IPIV( 1 ) = I
    IF( A( I, 1 ).NE.ZERO ) THEN
!
!           Apply the interchange
!
      IF( I.NE.1 ) THEN
        TEMP = A( 1, 1 )
        A( 1, 1 ) = A( I, 1 )
        A( I, 1 ) = TEMP
      END IF
!
!           Compute elements 2:M of the column
!
      IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
        CALL DSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
      ELSE
         DO 10 I = 1, M-1
           A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
10        CONTINUE
      END IF
!
    ELSE
      INFO = 1
    END IF
!
  ELSE
!
!        Use recursive code
!
    N1 = MIN( M, N ) / 2
    N2 = N-N1
!
!               [ A11 ]
!        Factor [ --- ]
!               [ A21 ]
!
    CALL DGETRF2( M, N1, A, LDA, IPIV, IINFO )

    IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO
!
!                              [ A12 ]
!        Apply interchanges to [ --- ]
!                              [ A22 ]
!
    CALL DLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
!
!        Solve A12
!
    CALL DTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, &
&               A( 1, N1+1 ), LDA )
!
!        Update A22
!
   CALL DGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, &
&               A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
!
!        Factor A22
!
   CALL DGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ),&
&                 IINFO )
!
!        Adjust INFO and the pivot indices
!
   IF ( INFO.EQ.0 .AND. IINFO.GT.0 )  INFO = IINFO + N1

   DO 20 I = N1+1, MIN( M, N )
      IPIV( I ) = IPIV( I ) + N1
20  CONTINUE
!
!        Apply interchanges to A21
!
  CALL DLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
!
END IF
RETURN
!
!     End of DGETRF2
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dlamch.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!

DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  CHARACTER          CMACH
!     ..
!
! =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!    ..
!     .. Local Scalars ..
  DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
&                   MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
  RND = ONE
!
  IF( ONE.EQ.RND ) THEN
    EPS = EPSILON(ZERO) * 0.5
  ELSE
    EPS = EPSILON(ZERO)
  END IF
!
  IF( LSAME( CMACH, 'E' ) ) THEN
    RMACH = EPS
  ELSE IF( LSAME( CMACH, 'S' ) ) THEN
    SFMIN = TINY(ZERO)
    SMALL = ONE / HUGE(ZERO)
    IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
      SFMIN = SMALL*( ONE+EPS )
    END IF
    RMACH = SFMIN
  ELSE IF( LSAME( CMACH, 'B' ) ) THEN
    RMACH = RADIX(ZERO)
  ELSE IF( LSAME( CMACH, 'P' ) ) THEN
    RMACH = EPS * RADIX(ZERO)
  ELSE IF( LSAME( CMACH, 'N' ) ) THEN
    RMACH = DIGITS(ZERO)
  ELSE IF( LSAME( CMACH, 'R' ) ) THEN
    RMACH = RND
  ELSE IF( LSAME( CMACH, 'M' ) ) THEN
    RMACH = MINEXPONENT(ZERO)
  ELSE IF( LSAME( CMACH, 'U' ) ) THEN
    RMACH = tiny(zero)
  ELSE IF( LSAME( CMACH, 'L' ) ) THEN
    RMACH = MAXEXPONENT(ZERO)
  ELSE IF( LSAME( CMACH, 'O' ) ) THEN
    RMACH = HUGE(ZERO)
  ELSE
    RMACH = ZERO
  END IF
!
DLAMCH = RMACH
RETURN
!
!     End of DLAMCH
!
END


!
!-----------------------------------------------------------------------
!
!     This is the original file dgetrs.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  CHARACTER          TRANS
  INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
    &    LSAME( TRANS, 'C' ) ) THEN
    INFO = -1
  ELSE IF( N.LT.0 ) THEN
    INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
    INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    INFO = -5
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
    INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DGETRS', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
  IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
    CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
    CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
      & ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
    CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
      & NRHS, ONE, A, LDA, B, LDB )
  ELSE
!
!        Solve A**T * X = B.
!
!        Solve U**T *X = B, overwriting B with X.
!
    CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
      & ONE, A, LDA, B, LDB )
!
!        Solve L**T *X = B, overwriting B with X.
!
    CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
      & A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
    CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
  END IF
!
RETURN
!
!     End of DGETRS
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dlaswp.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
  SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
  INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  DOUBLE PRECISION   A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
  INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
  DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
  IF( INCX.GT.0 ) THEN
    IX0 = K1
    I1 = K1
    I2 = K2
    INC = 1
  ELSE IF( INCX.LT.0 ) THEN
    IX0 = 1 + ( 1-K2 )*INCX
    I1 = K2
    I2 = K1
    INC = -1
  ELSE
    RETURN
  END IF
!
  N32 = ( N / 32 )*32
  IF( N32.NE.0 ) THEN
    DO 30 J = 1, N32, 32
      IX = IX0
      DO 20 I = I1, I2, INC
        IP = IPIV( IX )
        IF( IP.NE.I ) THEN
          DO 10 K = J, J + 31
            TEMP = A( I, K )
            A( I, K ) = A( IP, K )
            A( IP, K ) = TEMP
10        CONTINUE
        END IF
        IX = IX + INCX
20    CONTINUE
30 CONTINUE
  END IF

  IF( N32.NE.N ) THEN
    N32 = N32 + 1
    IX = IX0
    DO 50 I = I1, I2, INC
      IP = IPIV( IX )
      IF( IP.NE.I ) THEN
        DO 40 K = N32, N
          TEMP = A( I, K )
          A( I, K ) = A( IP, K )
          A( IP, K ) = TEMP
40      CONTINUE
      END IF
      IX = IX + INCX
50  CONTINUE
  END IF
!
RETURN
!
!     End of DLASWP
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file dtrtrs.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, &
  &  INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  CHARACTER          DIAG, TRANS, UPLO
  INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
  DOUBLE PRECISION   ZERO, ONE
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            NOUNIT
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NOUNIT = LSAME( DIAG, 'N' )
  IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
    INFO = -1
  ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. &
    & LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
    INFO = -2
  ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
    INFO = -3
  ELSE IF( N.LT.0 ) THEN
    INFO = -4
  ELSE IF( NRHS.LT.0 ) THEN
    INFO = -5
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    INFO = -7
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
    INFO = -9
  END IF
  IF( INFO.NE.0 ) THEN
    CALL XERBLA( 'DTRTRS', -INFO )
    RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) RETURN
!
!     Check for singularity.
!
  IF( NOUNIT ) THEN
    DO 10 INFO = 1, N
      IF( A( INFO, INFO ).EQ.ZERO ) RETURN
10  CONTINUE
  END IF
  INFO = 0
!
!     Solve A * x = b  or  A**T * x = b.
!
  CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, &
    & LDB )
!
RETURN
!
!     End of DTRTRS
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file ieeeck.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
  INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  INTEGER            ISPEC
  REAL               ONE, ZERO
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
    &                NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
  IEEECK = 1
!
  POSINF = ONE / ZERO
  IF( POSINF.LE.ONE ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  NEGINF = -ONE / ZERO
  IF( NEGINF.GE.ZERO ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  NEGZRO = ONE / ( NEGINF+ONE )
  IF( NEGZRO.NE.ZERO ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  NEGINF = ONE / NEGZRO
  IF( NEGINF.GE.ZERO ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  NEWZRO = NEGZRO + ZERO
  IF( NEWZRO.NE.ZERO ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  POSINF = ONE / NEWZRO
  IF( POSINF.LE.ONE ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  NEGINF = NEGINF*POSINF
  IF( NEGINF.GE.ZERO ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  POSINF = POSINF*POSINF
  IF( POSINF.LE.ONE ) THEN
    IEEECK = 0
    RETURN
  END IF

!
!     Return if we were only asked to check infinity arithmetic
!
  IF( ISPEC.EQ.0 ) RETURN
  !
  NAN1 = POSINF + NEGINF
  !
  NAN2 = POSINF / NEGINF
  !
  NAN3 = POSINF / POSINF
  !
  NAN4 = POSINF*ZERO
  !
  NAN5 = NEGINF*NEGZRO
  !
  NAN6 = NAN5*ZERO
!
  IF( NAN1.EQ.NAN1 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  IF( NAN2.EQ.NAN2 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  IF( NAN3.EQ.NAN3 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  IF( NAN4.EQ.NAN4 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  IF( NAN5.EQ.NAN5 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
  IF( NAN6.EQ.NAN6 ) THEN
    IEEECK = 0
    RETURN
  END IF
!
RETURN
END

!
!-----------------------------------------------------------------------
!
!     This is the original file ilaenv.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
  INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
  CHARACTER*( * )    NAME, OPTS
  INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
  INTEGER            I, IC, IZ, NB, NBMIN, NX
  LOGICAL            CNAME, SNAME
  CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
  INTEGER            IEEECK, IPARMQ
  EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
  GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
    & 130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
  ILAENV = -1
  RETURN
!
10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
  ILAENV = 1
  SUBNAM = NAME
  IC = ICHAR( SUBNAM( 1: 1 ) )
  IZ = ICHAR( 'Z' )
  IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
    IF( IC.GE.97 .AND. IC.LE.122 ) THEN
      SUBNAM( 1: 1 ) = CHAR( IC-32 )
      DO 20 I = 2, 6
        IC = ICHAR( SUBNAM( I: I ) )
        IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
20    CONTINUE
    END IF
!
  ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
    IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
      &   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
      &   ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
      SUBNAM( 1: 1 ) = CHAR( IC+64 )
      DO 30 I = 2, 6
        IC = ICHAR( SUBNAM( I: I ) )
        IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
           &  ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
           &  ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
30    CONTINUE
    END IF
!
  ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
    IF( IC.GE.225 .AND. IC.LE.250 ) THEN
      SUBNAM( 1: 1 ) = CHAR( IC-32 )
      DO 40 I = 2, 6
        IC = ICHAR( SUBNAM( I: I ) )
        IF( IC.GE.225 .AND. IC.LE.250 ) &
        &    SUBNAM( I: I ) = CHAR( IC-32 )
40    CONTINUE
    END IF
  END IF
!
  C1 = SUBNAM( 1: 1 )
  SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
  CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
  IF( .NOT.( CNAME .OR. SNAME ) ) RETURN
  C2 = SUBNAM( 2: 3 )
  C3 = SUBNAM( 4: 6 )
  C4 = C3( 2: 3 )
!
  GO TO ( 50, 60, 70 )ISPEC
!
50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
  NB = 1
!
  IF( C2.EQ.'GE' ) THEN
    IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
        NB = 64
      ELSE
        NB = 64
      END IF
    ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
    &  C3.EQ.'QLF' ) THEN
      IF( SNAME ) THEN
        NB = 32
      ELSE
        NB = 32
      END IF
    ELSE IF( C3.EQ.'HRD' ) THEN
      IF( SNAME ) THEN
         NB = 32
      ELSE
         NB = 32
      END IF
    ELSE IF( C3.EQ.'BRD' ) THEN
      IF( SNAME ) THEN
         NB = 32
      ELSE
         NB = 32
      END IF
    ELSE IF( C3.EQ.'TRI' ) THEN
      IF( SNAME ) THEN
         NB = 64
      ELSE
         NB = 64
      END IF
    END IF
  ELSE IF( C2.EQ.'PO' ) THEN
    IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
         NB = 64
      ELSE
         NB = 64
      END IF
   END IF
  ELSE IF( C2.EQ.'SY' ) THEN
    IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
         NB = 64
      ELSE
         NB = 64
      END IF
    ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
      NB = 32
    ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
      NB = 64
    END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
   IF( C3.EQ.'TRF' ) THEN
      NB = 64
   ELSE IF( C3.EQ.'TRD' ) THEN
      NB = 32
   ELSE IF( C3.EQ.'GST' ) THEN
      NB = 64
   END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&         'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&         THEN
         NB = 32
      END IF
   ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NB = 32
      END IF
   END IF
ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NB = 32
      END IF
   ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NB = 32
      END IF
   END IF
ELSE IF( C2.EQ.'GB' ) THEN
   IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
         IF( N4.LE.64 ) THEN
            NB = 1
         ELSE
            NB = 32
         END IF
      ELSE
         IF( N4.LE.64 ) THEN
            NB = 1
         ELSE
            NB = 32
         END IF
      END IF
   END IF
ELSE IF( C2.EQ.'PB' ) THEN
   IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
         IF( N2.LE.64 ) THEN
            NB = 1
         ELSE
            NB = 32
         END IF
      ELSE
         IF( N2.LE.64 ) THEN
            NB = 1
         ELSE
            NB = 32
         END IF
      END IF
   END IF
ELSE IF( C2.EQ.'TR' ) THEN
   IF( C3.EQ.'TRI' ) THEN
      IF( SNAME ) THEN
         NB = 64
      ELSE
         NB = 64
      END IF
   END IF
ELSE IF( C2.EQ.'LA' ) THEN
   IF( C3.EQ.'UUM' ) THEN
      IF( SNAME ) THEN
         NB = 64
      ELSE
         NB = 64
      END IF
   END IF
ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
   IF( C3.EQ.'EBZ' ) THEN
      NB = 1
   END IF
ELSE IF( C2.EQ.'GG' ) THEN
   NB = 32
   IF( C3.EQ.'HD3' ) THEN
      IF( SNAME ) THEN
         NB = 32
      ELSE
         NB = 32
      END IF
   END IF
END IF
ILAENV = NB
RETURN
!
60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
NBMIN = 2
IF( C2.EQ.'GE' ) THEN
   IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
&      'QLF' ) THEN
      IF( SNAME ) THEN
         NBMIN = 2
      ELSE
         NBMIN = 2
      END IF
   ELSE IF( C3.EQ.'HRD' ) THEN
      IF( SNAME ) THEN
         NBMIN = 2
      ELSE
         NBMIN = 2
      END IF
   ELSE IF( C3.EQ.'BRD' ) THEN
      IF( SNAME ) THEN
         NBMIN = 2
      ELSE
         NBMIN = 2
      END IF
   ELSE IF( C3.EQ.'TRI' ) THEN
      IF( SNAME ) THEN
         NBMIN = 2
      ELSE
         NBMIN = 2
      END IF
   END IF
ELSE IF( C2.EQ.'SY' ) THEN
   IF( C3.EQ.'TRF' ) THEN
      IF( SNAME ) THEN
         NBMIN = 8
      ELSE
         NBMIN = 8
      END IF
   ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
      NBMIN = 2
   END IF
ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
   IF( C3.EQ.'TRD' ) THEN
      NBMIN = 2
   END IF
ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NBMIN = 2
      END IF
   ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NBMIN = 2
      END IF
   END IF
ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NBMIN = 2
      END IF
   ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NBMIN = 2
      END IF
   END IF
ELSE IF( C2.EQ.'GG' ) THEN
   NBMIN = 2
   IF( C3.EQ.'HD3' ) THEN
      NBMIN = 2
   END IF
END IF
ILAENV = NBMIN
RETURN
!
70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
NX = 0
IF( C2.EQ.'GE' ) THEN
   IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
&      'QLF' ) THEN
      IF( SNAME ) THEN
         NX = 128
      ELSE
         NX = 128
      END IF
   ELSE IF( C3.EQ.'HRD' ) THEN
      IF( SNAME ) THEN
         NX = 128
      ELSE
         NX = 128
      END IF
   ELSE IF( C3.EQ.'BRD' ) THEN
      IF( SNAME ) THEN
         NX = 128
      ELSE
         NX = 128
      END IF
   END IF
ELSE IF( C2.EQ.'SY' ) THEN
   IF( SNAME .AND. C3.EQ.'TRD' ) THEN
      NX = 32
   END IF
ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
   IF( C3.EQ.'TRD' ) THEN
      NX = 32
   END IF
ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NX = 128
      END IF
   END IF
ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
   IF( C3( 1: 1 ).EQ.'G' ) THEN
      IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
&          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
&           THEN
         NX = 128
      END IF
   END IF
ELSE IF( C2.EQ.'GG' ) THEN
   NX = 128
   IF( C3.EQ.'HD3' ) THEN
      NX = 128
   END IF
END IF
ILAENV = NX
RETURN
!
80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
ILAENV = 6
RETURN
!
90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
ILAENV = 2
RETURN
!
100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
RETURN
!
110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
ILAENV = 1
RETURN
!
120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
ILAENV = 50
RETURN
!
130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
ILAENV = 25
RETURN
!
140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
ILAENV = 1
IF( ILAENV.EQ.1 ) THEN
   ILAENV = IEEECK( 1, 0.0, 1.0 )
END IF
RETURN
!
150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
ILAENV = 1
IF( ILAENV.EQ.1 ) THEN
   ILAENV = IEEECK( 0, 0.0, 1.0 )
END IF
RETURN
!
160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or related subroutines.
!
ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
RETURN
!
!     End of ILAENV
!
END

!
!-----------------------------------------------------------------------
!
!     This is the original file lsame.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32

      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
           &  INTA.GE.145 .AND. INTA.LE.153 .OR. &
           & INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
           &  INTB.GE.145 .AND. INTB.LE.153 .OR. &
           & INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
     RETURN
!
!     End of LSAME
!
      END
!
!-----------------------------------------------------------------------
!
!     This is the original file xerbla.f from the LAPACK routines.
!
!-----------------------------------------------------------------------
!
SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
  CHARACTER*(*)      SRNAME
  INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
  INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!
  WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO

  STOP
!
9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', 'an illegal value' )
!
!     End of XERBLA
!
END


INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.6.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2015
!
!     .. Scalar Arguments ..
INTEGER            IHI, ILO, ISPEC, LWORK, N
CHARACTER          NAME*( * ), OPTS*( * )
!
!  ================================================================
!     .. Parameters ..
INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
&                   ISHFTS = 15, IACC22 = 16 )
INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
&                  NIBBLE = 14, KNWSWP = 500 )
REAL               TWO
PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
INTEGER            NH, NS
INTEGER            I, IC, IZ
CHARACTER          SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. &
&    ( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
   NH = IHI - ILO + 1
   NS = 2
   IF( NH.GE.30 ) &
&      NS = 4
   IF( NH.GE.60 ) &
&      NS = 10
   IF( NH.GE.150 ) &
&      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
   IF( NH.GE.590 ) &
&      NS = 64
   IF( NH.GE.3000 ) &
&      NS = 128
   IF( NH.GE.6000 ) &
&      NS = 256
   NS = MAX( 2, NS-MOD( NS, 2 ) )
END IF
!
IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
   IPARMQ = NMIN
!
ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
   IPARMQ = NIBBLE
!
ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
   IPARMQ = NS
!
ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
   IF( NH.LE.KNWSWP ) THEN
      IPARMQ = NS
   ELSE
      IPARMQ = 3*NS / 2
   END IF
!
ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
!
!        Convert NAME to upper case if the first character is lower case.
!
   IPARMQ = 0
   SUBNAM = NAME
   IC = ICHAR( SUBNAM( 1: 1 ) )
   IZ = ICHAR( 'Z' )
   IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!           ASCII character set
!
      IF( IC.GE.97 .AND. IC.LE.122 ) THEN
         SUBNAM( 1: 1 ) = CHAR( IC-32 )
         DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.97 .AND. IC.LE.122 ) &
&               SUBNAM( I: I ) = CHAR( IC-32 )
         END DO
      END IF
!
   ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!           EBCDIC character set
!
      IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
&          ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
&          ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
         SUBNAM( 1: 1 ) = CHAR( IC+64 )
         DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
&                ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
&                ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: &
&                I ) = CHAR( IC+64 )
         END DO
      END IF
!
   ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!           Prime machines:  ASCII+128
!
      IF( IC.GE.225 .AND. IC.LE.250 ) THEN
         SUBNAM( 1: 1 ) = CHAR( IC-32 )
         DO I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.225 .AND. IC.LE.250 ) &
&               SUBNAM( I: I ) = CHAR( IC-32 )
         END DO
      END IF
   END IF
!
   IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR. &
&       SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
      IPARMQ = 1
      IF( NH.GE.K22MIN ) &
&         IPARMQ = 2
   ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
      IF( NH.GE.KACMIN ) &
&        IPARMQ = 1
      IF( NH.GE.K22MIN ) &
&         IPARMQ = 2
   ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR. &
&             SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
      IF( NS.GE.KACMIN ) &
&         IPARMQ = 1
      IF( NS.GE.K22MIN ) &
&         IPARMQ = 2
   END IF
!
ELSE
!        ===== invalid value of ispec =====
   IPARMQ = -1
!
END IF
!
!     ==== End of IPARMQ ====
!
END
