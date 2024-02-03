      SUBROUTINE ROS4(N,FCN,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,DFX,IDFX,
     &                  MAS ,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
C     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4  
C     (WITH STEP SIZE CONTROL).
C     C.F. SECTION IV.7
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
C     
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
C         SPRINGER-VERLAG (1990)               
C      
C     VERSION OF NOVEMBER 17, 1992
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F)
C                    REAL*8 X,Y(N),F(N)
C                    F(1)=...   ETC.
C
C     IFCN        GIVES INFORMATION ON FCN:
C                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS)
C                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS)
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, 
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY
C                 ADAPTS ITS STEP SIZE. STUDY THE CHOSEN VALUES FOR A FEW
C                 STEPS IN SUBROUTINE "SOLOUT", WHEN YOU ARE NOT SURE.
C                 (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,LDFY)
C                    REAL*8 X,Y(N),DFY(LDFY,N)
C                    DFY(1,1)= ...
C                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN 
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=N.
C
C     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X
C                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1;
C                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0).
C                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE DFX(N,X,Y,FX)
C                    REAL*8 X,Y(N),FX(N)
C                    FX(1)= ...
C                
C     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX:
C                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED.
C                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX.
C
C     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
C                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
C                 MATRIX AND NEEDS NOT TO BE DEFINED;
C                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
C                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
C                    SUBROUTINE MAS(N,AM,LMAS)
C                    REAL*8 AM(LMAS,N)
C                    AM(1,1)= ....
C                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
C                    AS FULL MATRIX LIKE
C                         AM(I,J) = M(I,J)
C                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
C                    DIAGONAL-WISE AS
C                         AM(I-J+MUMAS+1,J) = M(I,J).
C
C     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
C                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
C                       MATRIX, MAS IS NEVER CALLED.
C                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
C
C     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
C                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
C
C     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLMAS=N.
C                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,IRTRN)
C                    REAL*8 X,Y(N)
C                    ....  
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, ROS4 RETURNS TO THE CALLING PROGRAM.
C           
C     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                             N*(LJAC+LMAS+LE1+8)+5
C                 WHERE
C                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
C                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
C                 AND                  
C                    LMAS=0              IF IMAS=0
C                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
C                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
C                 AND
C                    LE1=N               IF MLJAC=N (FULL JACOBIAN)
C                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.). 
c
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
C                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
C                 STORAGE REQUIREMENT IS 
C                             LWORK = 2*N*N+8*N+5.
C
C     LWORK       DECLARED LENGHT OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 "LIWORK" MUST BE AT LEAST N+2.
C
C     LIWORK      DECLARED LENGHT OF ARRAY "IWORK".
C
C ----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(5)
C              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
C
C    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
C              IF IWORK(2).EQ.1  METHOD OF SHAMPINE
C              IF IWORK(2).EQ.2  METHOD GRK4T OF KAPS-RENTROP
C              IF IWORK(2).EQ.3  METHOD GRK4A OF KAPS-RENTROP
C              IF IWORK(2).EQ.4  METHOD OF VAN VELDHUIZEN (GAMMA=1/2)
C              IF IWORK(2).EQ.5  METHOD OF VAN VELDHUIZEN ("D-STABLE")
C              IF IWORK(2).EQ.6  AN L-STABLE METHOD
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=2.
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(3) <= HNEW/HOLD <= WORK(4)
C              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0
C
C    WORK(5)   AVOID THE HUMP: AFTER TWO CONSECUTIVE STEP REJECTIONS
C              THE STEP SIZE IS MULTIPLIED BY WORK(5)
C              DEFAULT VALUES: WORK(5)=0.1D0
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE WHERE THE SOLUTION IS COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND)
C
C     Y(N)        SOLUTION AT X
C  
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID=1  COMPUTATION SUCCESSFUL,
C                   IDID=-1 COMPUTATION UNSUCCESSFUL.
C
C --------------------------------------------------------- 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N),ATOL(1),RTOL(1),WORK(LWORK),IWORK(LIWORK)
      LOGICAL AUTNMS,IMPLCT,JBAND,ARRET
      EXTERNAL FCN,JAC,DFX,MAS,SOLOUT
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C -----------------------------------------------------
C --- COMMON STAT CAN BE USED FOR STATISTICS
C ---    NFCN      NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                  EVALUATION OF THE JACOBIAN ARE NOT COUNTED)  
C ---    NJAC      NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                  OR NUMERICALLY)
C ---    NSTEP     NUMBER OF COMPUTED STEPS
C ---    NACCPT    NUMBER OF ACCEPTED STEPS
C ---    NREJCT    NUMBER OF REJECTED STEPS (AFTER AT LEAST ONE STEP
C                  HAS BEEN ACCEPTED)
C ---    NDEC      NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX)
C ---    NSOL      NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS
C -----------------------------------------------------
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
      NFCN=0
      NJAC=0
      NSTEP=0
      NACCPT=0
      NREJCT=0
      NDEC=0
      NSOL=0
      ARRET=.FALSE.
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWORK(1).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(1)
         IF(NMAX.LE.0)THEN
            WRITE(6,*)' WRONG INPUT IWORK(1)=',IWORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- METH   COEFFICIENTS OF THE METHOD
      IF(IWORK(2).EQ.0)THEN
         METH=2
      ELSE
         METH=IWORK(2)
         IF(METH.LE.0.OR.METH.GE.7)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  
      IF(WORK(1).EQ.0.D0)THEN
         UROUND=1.D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-14.OR.UROUND.GE.1.D0)THEN
            WRITE(6,*)' COEFFICIENTS HAVE 16 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WORK(2).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(2)
      END IF
C -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(3).EQ.0.D0)THEN
         FAC1=5.D0
      ELSE
         FAC1=1.D0/WORK(3)
      END IF
      IF(WORK(4).EQ.0.D0)THEN
         FAC2=1.D0/6.0D0
      ELSE
         FAC2=1.D0/WORK(4)
      END IF
C -------  FACREJ    FOR THE HUMP
      IF(WORK(5).EQ.0.D0)THEN
         FACREJ=0.1D0
      ELSE
         FACREJ=WORK(5)
      END IF
C --------- CHECK IF TOLERANCES ARE O.K.
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
      ELSE
          DO 15 I=1,N
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          END IF
  15      CONTINUE
      END IF
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C         COMPUTATION OF ARRAY ENTRIES
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ?
      AUTNMS=IFCN.EQ.0
      IMPLCT=IMAS.NE.0
      JBAND=MLJAC.NE.N
      ARRET=.FALSE.
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
C -- JACOBIAN 
      IF(JBAND)THEN
         LDJAC=MLJAC+MUJAC+1
      ELSE
         LDJAC=N
      END IF
C -- MATRIX E FOR LINEAR ALGEBRA
      IF(JBAND)THEN
         LDE=2*MLJAC+MUJAC+1
      ELSE
         LDE=N
      END IF
C -- MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.N) THEN
              LDMAS=MLMAS+MUMAS+1
          ELSE
              LDMAS=N
          END IF
C ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
              WRITE (6,*) 'BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF
     & "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
      END IF
      LDMAS2=MAX(1,LDMAS)
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEYNEW=6
      IEDY1=IEYNEW+N
      IEDY=IEDY1+N
      IEAK1=IEDY+N
      IEAK2=IEAK1+N
      IEAK3=IEAK2+N
      IEAK4=IEAK3+N
      IEFX =IEAK4+N
      IEJAC=IEFX +N
      IEMAS=IEJAC+N*LDJAC
      IEE  =IEMAS+N*LDMAS
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE+N*LDE-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP=3
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIP+N-1
      IF(ISTORE.GT.LIWORK)THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL RO4COR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,IJAC,
     &   MLJAC,MUJAC,DFX,IDFX,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,METH,FAC1,FAC2,FACREJ,AUTNMS,IMPLCT,JBAND,
     &   LDJAC,LDE,LDMAS2,WORK(IEYNEW),WORK(IEDY1),WORK(IEDY),
     &   WORK(IEAK1),WORK(IEAK2),WORK(IEAK3),WORK(IEAK4),
     &   WORK(IEFX),WORK(IEJAC),WORK(IEE),WORK(IEMAS),IWORK(IEIP))
C ----------- RETURN -----------
      RETURN
      END
C
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE RO4COR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,
     &  IJAC,MLJAC,MUJAC,DFX,IDFX,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &  NMAX,UROUND,METH,FAC1,FAC2,FACREJ,AUTNMS,IMPLCT,BANDED,
     &  LFJAC,LE,LDMAS,YNEW,DY1,DY,AK1,AK2,AK3,AK4,FX,FJAC,E,FMAS,IP)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR ROS4
C     PARAMETERS SAME AS IN ROS4 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(N),YNEW(N),DY1(N),DY(N),AK1(N),
     *  AK2(N),AK3(N),AK4(N),FX(N),
     *  FJAC(LFJAC,N),E(LE,N),FMAS(LDMAS,N),ATOL(1),RTOL(1)
      INTEGER IP(N)
      LOGICAL REJECT,RJECT2,AUTNMS,IMPLCT,BANDED
      COMMON/STAT/NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) CALL MAS(N,FMAS,LDMAS)
C ---- PREPARE BANDWIDTHS -----
      IF (BANDED) THEN
          MLE=MLJAC
          MUE=MUJAC
          MBJAC=MLJAC+MUJAC+1
          MBB=MLMAS+MUMAS+1
          MDIAG=MLE+MUE+1
          MBDIAG=MUMAS+1
          MDIFF=MLE+MUE-MUMAS
      END IF
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
      IF (METH.EQ.1) CALL SHAMP (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.2) CALL GRK4T (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.3) CALL GRK4A (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.4) CALL VELDS (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.5) CALL VELDD (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IF (METH.EQ.6) CALL LSTAB (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
C --- INITIAL PREPARATIONS
      POSNEG=SIGN(1.D0,XEND-X)
      HMAXN=MIN(ABS(HMAX),ABS(XEND-X))
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAXN) 
      H=SIGN(H,POSNEG) 
      REJECT=.FALSE.
      NSING=0
      IRTRN=1 
      XXOLD=X
      IF (IOUT.NE.0) CALL SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN)
      IF (IRTRN.LT.0) GOTO 79
C --- BASIC INTEGRATION STEP  
   1  IF (NSTEP.GT.NMAX.OR.X+.1D0*H.EQ.X.OR.ABS(H).LE.UROUND) GOTO 79
      IF ((X-XEND)*POSNEG+UROUND.GT.0.D0) THEN
          H=HOPT
          IDID=1
          RETURN
      END IF
      HOPT=H
      IF ((X+H-XEND)*POSNEG.GT.0.D0) H=XEND-X
      CALL FCN(N,X,Y,DY1)
      NFCN=NFCN+1
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
      IF (IJAC.EQ.0) THEN
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY
          IF (BANDED) THEN
C --- JACOBIAN IS BANDED
              MUJACP=MUJAC+1
              MD=MIN(MBJAC,N)
              DO 16 K=1,MD
              J=K
 12           AK2(J)=Y(J)
              AK3(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
              Y(J)=Y(J)+AK3(J)
              J=J+MD
              IF (J.LE.N) GOTO 12 
              CALL FCN(N,X,Y,AK1)
              J=K
              LBEG=MAX(1,J-MUJAC)
 14           LEND=MIN(N,J+MLJAC)
              Y(J)=AK2(J)
              MUJACJ=MUJACP-J
              DO 15 L=LBEG,LEND
 15           FJAC(L+MUJACJ,J)=(AK1(L)-DY1(L))/AK3(J) 
              J=J+MD
              LBEG=LEND+1
              IF (J.LE.N) GOTO 14
 16           CONTINUE
          ELSE
C --- JACOBIAN IS FULL
              DO 18 I=1,N
              YSAFE=Y(I)
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
              Y(I)=YSAFE+DELT
              CALL FCN(N,X,Y,AK1)
              DO 17 J=1,N
  17          FJAC(J,I)=(AK1(J)-DY1(J))/DELT
  18          Y(I)=YSAFE
              MLJAC=N
          END IF
      ELSE
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
          CALL JAC(N,X,Y,FJAC,LFJAC)
      END IF
      IF (.NOT.AUTNMS) THEN
          IF (IDFX.EQ.0) THEN
C --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(X)))
              XDELT=X+DELT
              CALL FCN(N,XDELT,Y,AK1)
              DO 19 J=1,N
  19          FX(J)=(AK1(J)-DY1(J))/DELT
          ELSE
C --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X
              CALL DFX(N,X,Y,FX)
          END IF
      END IF
   2  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTE THE STAGES
C *** *** *** *** *** *** ***
      NDEC=NDEC+1
      HC21=C21/H
      HC31=C31/H
      HC32=C32/H
      HC41=C41/H
      HC42=C42/H
      HC43=C43/H
      FAC=1.D0/(H*GAMMA)
      IF (IMPLCT) THEN
          IF (BANDED) THEN
C --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX)
              DO 101 J=1,N
              I1=MAX0(1,MUJAC+2-J)
              I2=MIN0(MBJAC,N+MUJAC+1-J)
              DO 101 I=I1,I2
  101         E(I+MLE,J)=-FJAC(I,J)
              DO 102 J=1,N
              I1=MAX0(1,MUMAS+2-J)
              I2=MIN0(MBB,N+MUMAS+1-J)
              DO 102 I=I1,I2
              IB=I+MDIFF
  102         E(IB,J)=E(IB,J)+FAC*FMAS(I,J)
              CALL DECB(N,LE,E,MLE,MUE,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
C ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 103 I=1,N
  103             AK1(I)=DY1(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 110 I=1,N
  110             YNEW(I)=Y(I)+A21*AK1(I) 
                      CALL FCN(N,X,YNEW,DY)
                  DO 111 I=1,N
  111             YNEW(I)=HC21*AK1(I)
                  DO 114 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 113 J=J1,J2
  113             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  114             AK2(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 120 I=1,N
  120             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                      CALL FCN(N,X,YNEW,DY)
                  DO 121 I=1,N
  121             YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                  DO 124 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 123 J=J1,J2
  123             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  124             AK3(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 131 I=1,N
  131             YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                  DO 134 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 133 J=J1,J2
  133             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  134             AK4(I)=SUM+DY(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              ELSE
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B AND THE JACOBIAN OF F ARE BANDED
C ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 203 I=1,N
  203             AK1(I)=DY1(I)+HD1*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 210 I=1,N
  210             YNEW(I)=Y(I)+A21*AK1(I) 
                      CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 211 I=1,N
  211             YNEW(I)=HC21*AK1(I)
                  DO 214 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 213 J=J1,J2
  213             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  214             AK2(I)=SUM+DY(I)+HD2*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 220 I=1,N
  220             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                      CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 221 I=1,N
  221             YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                  DO 224 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 223 J=J1,J2
  223             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  224             AK3(I)=SUM+DY(I)+HD3*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 231 I=1,N
  231             YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                  DO 234 I=1,N
                  SUM=0.D0
                  J1=MAX0(1,I-MLMAS)
                  J2=MIN0(N,I+MUMAS)
                  DO 233 J=J1,J2
  233             SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
  234             AK4(I)=SUM+DY(I)+HD4*FX(I)
                      CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              END IF
          ELSE
              IF (MLMAS.NE.N) THEN
C --- THE MATRIX E (B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX)
                  MADD=MUMAS+1
                  DO 301 J=1,N
                  DO 301 I=1,N
  301             E(I,J)=-FJAC(I,J)
                  DO 302 J=1,N
                  I1=MAX0(1,J-MUMAS)
                  I2=MIN0(N,J+MLMAS)
                  DO 302 I=I1,I2
  302             E(I,J)=E(I,J)+FAC*FMAS(I-J+MADD,J)
                  CALL DEC(N,LE,E,IP,INFO)
                  IF (INFO.NE.0) GOTO 80
                  IF (AUTNMS) THEN
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
C ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                      DO 303 I=1,N
  303                 AK1(I)=DY1(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 310 I=1,N
  310                 YNEW(I)=Y(I)+A21*AK1(I) 
                          CALL FCN(N,X,YNEW,DY)
                      DO 311 I=1,N
  311                 YNEW(I)=HC21*AK1(I)
                      DO 314 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 313 J=J1,J2
  313                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  314                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 320 I=1,N
  320                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                           CALL FCN(N,X,YNEW,DY)
                      DO 321 I=1,N
  321                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 324 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 323 J=J1,J2
  323                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  324                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 331 I=1,N
  331                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                      DO 334 I=1,N
                      SUM=DY(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 333 J=J1,J2
  333                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  334                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  ELSE
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B IS BANDED BUT THE JACOBIAN OF F IS NOT
C ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                      HD1=H*D1
                      HD2=H*D2
                      HD3=H*D3
                      HD4=H*D4
                      DO 353 I=1,N
  353                 AK1(I)=DY1(I)+HD1*FX(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 360 I=1,N
  360                 YNEW(I)=Y(I)+A21*AK1(I) 
                          CALL FCN(N,X+C2*H,YNEW,DY)
                      DO 361 I=1,N
  361                 YNEW(I)=HC21*AK1(I)
                      DO 364 I=1,N
                      SUM=DY(I)+HD2*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 363 J=J1,J2
  363                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  364                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 370 I=1,N
  370                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                           CALL FCN(N,X+C3*H,YNEW,DY)
                      DO 371 I=1,N
  371                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 374 I=1,N
                      SUM=DY(I)+HD3*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 373 J=J1,J2
  373                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  374                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 381 I=1,N
  381                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                      DO 384 I=1,N
                      SUM=DY(I)+HD4*FX(I)
                      J1=MAX0(1,I-MLMAS)
                      J2=MIN0(N,I+MUMAS)
                      DO 383 J=J1,J2
  383                 SUM=SUM+FMAS(I-J+MADD,J)*YNEW(J)
  384                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  END IF
              ELSE
C --- THE MATRIX E (B IS A FULL MATRIX, JACOBIAN A FULL OR BANDED MATRIX)
                  IF (MLJAC.EQ.N) THEN
                      DO 401 J=1,N
                      DO 401 I=1,N
  401                 E(I,J)=FMAS(I,J)*FAC-FJAC(I,J)
                  ELSE
                      MADD=MUJAC+1
                      DO 405 J=1,N
                      DO 405 I=1,N
  405                 E(I,J)=FMAS(I,J)*FAC
                      DO 406 J=1,N
                      I1=MAX0(1,J-MUJAC)
                      I2=MIN0(N,J+MLJAC)
                      DO 406 I=I1,I2
  406                 E(I,J)=E(I,J)-FJAC(I-J+MADD,J)
                  END IF
                  CALL DEC(N,LE,E,IP,INFO)
                  IF (INFO.NE.0) GOTO 80
                  IF (AUTNMS) THEN
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B IS NOT BANDED
C ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                      DO 403 I=1,N
  403                 AK1(I)=DY1(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 410 I=1,N
  410                 YNEW(I)=Y(I)+A21*AK1(I) 
                          CALL FCN(N,X,YNEW,DY)
                      DO 411 I=1,N
  411                 YNEW(I)=HC21*AK1(I)
                      DO 414 I=1,N
                      SUM=DY(I)
                      DO 413 J=1,N
  413                 SUM=SUM+FMAS(I,J)*YNEW(J)
  414                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 420 I=1,N
  420                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                           CALL FCN(N,X,YNEW,DY)
                      DO 421 I=1,N
  421                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 424 I=1,N
                      SUM=DY(I)
                      DO 423 J=1,N
  423                 SUM=SUM+FMAS(I,J)*YNEW(J)
  424                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 431 I=1,N
  431                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                      DO 434 I=1,N
                      SUM=DY(I)
                      DO 433 J=1,N
  433                 SUM=SUM+FMAS(I,J)*YNEW(J)
  434                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  ELSE
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
C ---   2) THE MATRIX B IS NOT BANDED
C ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                      HD1=H*D1
                      HD2=H*D2
                      HD3=H*D3
                      HD4=H*D4
                      DO 503 I=1,N
  503                 AK1(I)=DY1(I)+HD1*FX(I)
                          CALL SOL(N,LE,E,AK1,IP)
                      DO 510 I=1,N
  510                 YNEW(I)=Y(I)+A21*AK1(I) 
                          CALL FCN(N,X+C2*H,YNEW,DY)
                      DO 511 I=1,N
  511                 YNEW(I)=HC21*AK1(I)
                      DO 514 I=1,N
                      SUM=DY(I)+HD2*FX(I)
                      DO 513 J=1,N
  513                 SUM=SUM+FMAS(I,J)*YNEW(J)
  514                 AK2(I)=SUM
                          CALL SOL(N,LE,E,AK2,IP)
                      DO 520 I=1,N
  520                 YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                           CALL FCN(N,X+C3*H,YNEW,DY)
                      DO 521 I=1,N
  521                 YNEW(I)=HC31*AK1(I)+HC32*AK2(I)
                      DO 524 I=1,N
                      SUM=DY(I)+HD3*FX(I)
                      DO 523 J=1,N
  523                 SUM=SUM+FMAS(I,J)*YNEW(J)
  524                 AK3(I)=SUM
                          CALL SOL(N,LE,E,AK3,IP)
                      DO 531 I=1,N
  531                 YNEW(I)=HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                      DO 534 I=1,N
                      SUM=DY(I)+HD4*FX(I)
                      DO 533 J=1,N
  533                 SUM=SUM+FMAS(I,J)*YNEW(J)
  534                 AK4(I)=SUM
                          CALL SOL(N,LE,E,AK4,IP)
                  END IF
              END IF
          END IF
      ELSE
          IF (BANDED) THEN
C --- THE MATRIX E (B=IDENTITY, JACOBIAN A BANDED MATRIX)
              DO 601 J=1,N
              I1=MAX0(1,MUJAC+2-J)
              I2=MIN0(MBJAC,N+MUJAC+1-J)
              DO 600 I=I1,I2
  600         E(I+MLE,J)=-FJAC(I,J)
  601         E(MDIAG,J)=E(MDIAG,J)+FAC
              CALL DECB(N,LE,E,MLE,MUE,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
C ---   2) THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
C ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 603 I=1,N
  603             AK1(I)=DY1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 610 I=1,N
  610             YNEW(I)=Y(I)+A21*AK1(I) 
                  CALL FCN(N,X,YNEW,DY)
                  DO 611 I=1,N
  611             AK2(I)=DY(I)+HC21*AK1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 620 I=1,N
  620             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                  CALL FCN(N,X,YNEW,DY)
                  DO 621 I=1,N
  621             AK3(I)=DY(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 631 I=1,N
  631             AK4(I)=DY(I)+HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                  CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              ELSE
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
C ---   2) THE JACOBIAN OF THE PROBLEM IS A BANDED MATRIX
C ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 703 I=1,N
  703             AK1(I)=DY1(I)+HD1*FX(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK1,IP)
                  DO 710 I=1,N
  710             YNEW(I)=Y(I)+A21*AK1(I) 
                  CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 711 I=1,N
  711             AK2(I)=DY(I)+HD2*FX(I)+HC21*AK1(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK2,IP)
                  DO 720 I=1,N
  720             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                  CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 721 I=1,N
  721             AK3(I)=DY(I)+HD3*FX(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOLB(N,LE,E,MLE,MUE,AK3,IP)
                  DO 731 I=1,N
  731             AK4(I)=DY(I)+HD4*FX(I)+HC41*AK1(I)+HC42*AK2(I)
     &                  +HC43*AK3(I) 
                  CALL SOLB(N,LE,E,MLE,MUE,AK4,IP)
              END IF
          ELSE
C --- THE MATRIX E (B=IDENTITY, JACOBIAN A FULL MATRIX)
              DO 801 J=1,N
              DO 800 I=1,N
  800         E(I,J)=-FJAC(I,J)
  801         E(J,J)=E(J,J)+FAC
              CALL DEC(N,LE,E,IP,INFO)
              IF (INFO.NE.0) GOTO 80
              IF (AUTNMS) THEN
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
C ---   2) THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
C ---   3) THE DIFFERENTIAL EQUATION IS AUTONOMOUS
                  DO 803 I=1,N
  803             AK1(I)=DY1(I)
                  CALL SOL(N,LE,E,AK1,IP)
                  DO 810 I=1,N
  810             YNEW(I)=Y(I)+A21*AK1(I) 
                  CALL FCN(N,X,YNEW,DY)
                  DO 811 I=1,N
  811             AK2(I)=DY(I)+HC21*AK1(I)
                  CALL SOL(N,LE,E,AK2,IP)
                  DO 820 I=1,N
  820             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                  CALL FCN(N,X,YNEW,DY)
                  DO 821 I=1,N
  821             AK3(I)=DY(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOL(N,LE,E,AK3,IP)
                  DO 831 I=1,N
  831             AK4(I)=DY(I)+HC41*AK1(I)+HC42*AK2(I)+HC43*AK3(I) 
                  CALL SOL(N,LE,E,AK4,IP)
              ELSE
C --- THIS PART COMPUTES THE STAGES IN THE CASE WHERE
C ---   1) THE DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
C ---   2) THE JACOBIAN OF THE PROBLEM IS A FULL MATRIX
C ---   3) THE DIFFERENTIAL EQUATION IS NON-AUTONOMOUS
                  HD1=H*D1
                  HD2=H*D2
                  HD3=H*D3
                  HD4=H*D4
                  DO 903 I=1,N
  903             AK1(I)=DY1(I)+HD1*FX(I)
                  CALL SOL(N,LE,E,AK1,IP)
                  DO 910 I=1,N
  910             YNEW(I)=Y(I)+A21*AK1(I) 
                  CALL FCN(N,X+C2*H,YNEW,DY)
                  DO 911 I=1,N
  911             AK2(I)=DY(I)+HD2*FX(I)+HC21*AK1(I)
                  CALL SOL(N,LE,E,AK2,IP)
                  DO 920 I=1,N
  920             YNEW(I)=Y(I)+A31*AK1(I)+A32*AK2(I)  
                  CALL FCN(N,X+C3*H,YNEW,DY)
                  DO 921 I=1,N
  921             AK3(I)=DY(I)+HD3*FX(I)+HC31*AK1(I)+HC32*AK2(I)
                  CALL SOL(N,LE,E,AK3,IP)
                  DO 931 I=1,N
  931             AK4(I)=DY(I)+HD4*FX(I)+HC41*AK1(I)+HC42*AK2(I)
     &                  +HC43*AK3(I) 
                  CALL SOL(N,LE,E,AK4,IP)
              END IF
          END IF
      END IF
      NSOL=NSOL+4
      NFCN=NFCN+2
C *** *** *** *** *** *** ***
C  ERROR ESTIMATION  
C *** *** *** *** *** *** ***
      NSTEP=NSTEP+1
C ------------ NEW SOLUTION ---------------
      DO 240 I=1,N
  240 YNEW(I)=Y(I)+B1*AK1(I)+B2*AK2(I)+B3*AK3(I)+B4*AK4(I)  
C ------------ COMPUTE ERROR ESTIMATION ----------------
      ERR=0.D0
      DO 300 I=1,N
      S=E1*AK1(I)+E2*AK2(I)+E3*AK3(I)+E4*AK4(I) 
      IF (ITOL.EQ.0) THEN
         SK=ATOL(1)+RTOL(1)*DMAX1(DABS(Y(I)),DABS(YNEW(I)))
      ELSE
         SK=ATOL(I)+RTOL(I)*DMAX1(DABS(Y(I)),DABS(YNEW(I)))
      END IF
  300 ERR=ERR+(S/SK)**2
      ERR=DSQRT(ERR/N)
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .2<=HNEW/H<=6.
      FAC=DMAX1(FAC2,DMIN1(FAC1,(ERR)**.25D0/.9D0))
      HNEW=H/FAC  
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LE.1.D0) THEN
C --- STEP IS ACCEPTED  
         NACCPT=NACCPT+1
         DO 44 I=1,N 
  44     Y(I)=YNEW(I) 
         XXOLD=X 
         X=X+H
         IF (IOUT.NE.0) CALL SOLOUT(NACCPT+1,XXOLD,X,Y,N,IRTRN)
         IF (IRTRN.LT.0) GOTO 79
         IF (DABS(HNEW).GT.HMAXN) HNEW=POSNEG*HMAXN
         IF (REJECT) HNEW=POSNEG*DMIN1(DABS(HNEW),DABS(H)) 
         REJECT=.FALSE.
         RJECT2=.FALSE.
         H=HNEW
         GOTO 1
      ELSE
C --- STEP IS REJECTED  
         IF (RJECT2) HNEW=H*FACREJ
         IF (REJECT) RJECT2=.TRUE.
         REJECT=.TRUE.
         H=HNEW
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         GOTO 2
      END IF
C --- EXIT
  80  WRITE (6,*) ' MATRIX E IS SINGULAR, INFO = ',INFO
      NSING=NSING+1
      IF (NSING.GE.5) GOTO 79
      H=H*0.5D0
      GOTO 2
  79  WRITE(6,979)X,H
 979  FORMAT(' EXIT OF ROS4 AT X=',D16.7,'   H=',D16.7)
      IDID=-1
      RETURN
      END
C
      SUBROUTINE SHAMP (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
         A21=2.D0
         A31=48.D0/25.D0
         A32=6.D0/25.D0
         C21=-8.D0
         C31=372.D0/25.D0
         C32=12.D0/5.D0
         C41=-112.D0/125.D0
         C42=-54.D0/125.D0
         C43=-2.D0/5.D0
         B1=19.D0/9.D0
         B2=1.D0/2.D0
         B3=25.D0/108.D0
         B4=125.D0/108.D0
         E1=17.D0/54.D0
         E2=7.D0/36.D0
         E3=0.D0
         E4=125.D0/108.D0
         GAMMA=.5D0
         C2= 0.1000000000000000D+01
         C3= 0.6000000000000000D+00
         D1= 0.5000000000000000D+00
         D2=-0.1500000000000000D+01
         D3= 0.2420000000000000D+01
         D4= 0.1160000000000000D+00
      RETURN
      END
C
      SUBROUTINE GRK4A (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.1108860759493671D+01
       A31= 0.2377085261983360D+01
       A32= 0.1850114988899692D+00
       C21=-0.4920188402397641D+01
       C31= 0.1055588686048583D+01
       C32= 0.3351817267668938D+01
       C41= 0.3846869007049313D+01
       C42= 0.3427109241268180D+01
       C43=-0.2162408848753263D+01
       B1= 0.1845683240405840D+01
       B2= 0.1369796894360503D+00
       B3= 0.7129097783291559D+00
       B4= 0.6329113924050632D+00
       E1= 0.4831870177201765D-01
       E2=-0.6471108651049505D+00
       E3= 0.2186876660500240D+00
       E4=-0.6329113924050632D+00
       GAMMA= 0.3950000000000000D+00
       C2= 0.4380000000000000D+00
       C3= 0.8700000000000000D+00
       D1= 0.3950000000000000D+00
       D2=-0.3726723954840920D+00
       D3= 0.6629196544571492D-01
       D4= 0.4340946962568634D+00
      RETURN
      END
C
      SUBROUTINE GRK4T (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.4524708207373116D+01
       A32= 0.4163528788597648D+01
       C21=-0.5071675338776316D+01
       C31= 0.6020152728650786D+01
       C32= 0.1597506846727117D+00
       C41=-0.1856343618686113D+01
       C42=-0.8505380858179826D+01
       C43=-0.2084075136023187D+01
       B1= 0.3957503746640777D+01
       B2= 0.4624892388363313D+01
       B3= 0.6174772638750108D+00
       B4= 0.1282612945269037D+01
       E1= 0.2302155402932996D+01
       E2= 0.3073634485392623D+01
       E3=-0.8732808018045032D+00
       E4=-0.1282612945269037D+01
       GAMMA= 0.2310000000000000D+00
       C2= 0.4620000000000000D+00
       C3= 0.8802083333333334D+00
       D1= 0.2310000000000000D+00
       D2=-0.3962966775244303D-01
       D3= 0.5507789395789127D+00
       D4=-0.5535098457052764D-01
      RETURN
      END
C
      SUBROUTINE VELDS (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
C --- METHOD GIVEN BY VAN VELDHUIZEN
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.1750000000000000D+01
       A32= 0.2500000000000000D+00
       C21=-0.8000000000000000D+01
       C31=-0.8000000000000000D+01
       C32=-0.1000000000000000D+01
       C41= 0.5000000000000000D+00
       C42=-0.5000000000000000D+00
       C43= 0.2000000000000000D+01
       B1= 0.1333333333333333D+01
       B2= 0.6666666666666667D+00
       B3=-0.1333333333333333D+01
       B4= 0.1333333333333333D+01
       E1=-0.3333333333333333D+00
       E2=-0.3333333333333333D+00
       E3=-0.0000000000000000D+00
       E4=-0.1333333333333333D+01
       GAMMA= 0.5000000000000000D+00
       C2= 0.1000000000000000D+01
       C3= 0.5000000000000000D+00
       D1= 0.5000000000000000D+00
       D2=-0.1500000000000000D+01
       D3=-0.7500000000000000D+00
       D4= 0.2500000000000000D+00
      RETURN
      END
C
      SUBROUTINE VELDD (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
C --- METHOD GIVEN BY VAN VELDHUIZEN
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.4812234362695436D+01
       A32= 0.4578146956747842D+01
       C21=-0.5333333333333331D+01
       C31= 0.6100529678848254D+01
       C32= 0.1804736797378427D+01
       C41=-0.2540515456634749D+01
       C42=-0.9443746328915205D+01
       C43=-0.1988471753215993D+01
       B1= 0.4289339254654537D+01
       B2= 0.5036098482851414D+01
       B3= 0.6085736420673917D+00
       B4= 0.1355958941201148D+01
       E1= 0.2175672787531755D+01
       E2= 0.2950911222575741D+01
       E3=-0.7859744544887430D+00
       E4=-0.1355958941201148D+01
       GAMMA= 0.2257081148225682D+00
       C2= 0.4514162296451364D+00
       C3= 0.8755928946018455D+00
       D1= 0.2257081148225682D+00
       D2=-0.4599403502680582D-01
       D3= 0.5177590504944076D+00
       D4=-0.3805623938054428D-01
      RETURN
      END 
C
      SUBROUTINE LSTAB (A21,A31,A32,C21,C31,C32,C41,C42,C43,
     &          B1,B2,B3,B4,E1,E2,E3,E4,GAMMA,C2,C3,D1,D2,D3,D4)
C --- AN L-STABLE METHOD
      IMPLICIT REAL*8 (A-H,O-Z)
       A21= 0.2000000000000000D+01
       A31= 0.1867943637803922D+01
       A32= 0.2344449711399156D+00
       C21=-0.7137615036412310D+01
       C31= 0.2580708087951457D+01
       C32= 0.6515950076447975D+00
       C41=-0.2137148994382534D+01
       C42=-0.3214669691237626D+00
       C43=-0.6949742501781779D+00
       B1= 0.2255570073418735D+01
       B2= 0.2870493262186792D+00
       B3= 0.4353179431840180D+00
       B4= 0.1093502252409163D+01
       E1=-0.2815431932141155D+00
       E2=-0.7276199124938920D-01
       E3=-0.1082196201495311D+00
       E4=-0.1093502252409163D+01
       GAMMA= 0.5728200000000000D+00
       C2= 0.1145640000000000D+01
       C3= 0.6552168638155900D+00
       D1= 0.5728200000000000D+00
       D2=-0.1769193891319233D+01
       D3= 0.7592633437920482D+00
       D4=-0.1049021087100450D+00
      RETURN
      END 

