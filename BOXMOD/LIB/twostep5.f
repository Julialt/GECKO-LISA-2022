C
      SUBROUTINE TWOSTEP5
     1 (maxsp, maxre, numre, mxleft, mxright, mself,
     2  numsp, numstoi, idrestoi, idpdstoi, nself,idselfreac,
     2  t, te, dtmin, dtmax, anrat, restoicf, pdstoicf,
     3  numit, atol, rtol, rem, rdep, rex, rdil,
     4  mtr, ntr, idtr, trprod, trloss, conid,ncons,
     1  y, noxfix, sumnox, idno, idno2, idno3,
     1  lpmap)

!JMLT!
! loops condensed where possible
! yields additional terms (rates / s):
!        trprod(T)=DT*YP(T,idtr) &
!        trloss(T)=DT*YL(T,idtr)*Y(T,idtr)
! !! OR !! (tracer prod / loss, organic rxns only)
!        trprod(T)=DT*TP(T,idtr) &
!        trloss(T)=DT*TL(T,idtr)*Y(T,idtr)
! Allows concentrations to be fixed at iteration level
C
C **********************************************************************
C PROGRAMMES APPELANTS : CHIMV4, CHIMPUFF
C PROGRAMMES APPELES   : ITER, FIT, NEWDT
C **********************************************************************
C
C     AUTHOR: JAN VERWER, CENTRE FOR MATHEMATICS AND COMPUTER SCIENCE
C     (CWI), KRUISLAAN 413, 1098 SJ AMSTERDAM, THE NETHERLANDS
C
C     EMAIL ADDRESS: JANV@CWI.NL
C     VERSION: NOVEMBER 1994 / 2
C
C     PURPOSE: TWOSTEP IS DESIGNED FOR THE NUMERICAL SOLUTION OF
C     STIFF ODE SYSTEMS
C
C        DY(T)/DT = YP(T,Y) - YL(Y,T)*Y(T)
C
C     ORIGINATING FROM ATMOSPHERIC CHEMISTRY. ITS UNDERLYING INTEGRATION
C     METHOD IS THE IMPLICIT, 2-ND ORDER, 2-STEP BDF FORMULA. A SIMPLE
C     EXPLICIT GAUSS-SEIDEL TECHNIQUE IS USED FOR APPROXIMATING THE
C     IMPLICITLY DEFINED BDF SOLUTIONS, RATHER THAN THE USUAL MODIFIED
C     NEWTON METHOD [1,2]. BY THIS APPROACH TWOSTEP IS EXPLICIT, WHEREAS
C     THE EXCELLENT STABILITY OF THE BDF METHOD IS MAINTAINED. ALSO THE
C     JACOBI TECHNIQUE MAY BE APPROPRIATE. IN GENERAL, HOWEVER, ITS USE
C     GENERALLY DECREASES THE STABILITY OF TWOSTEP AND HENCE ALSO THE
C     EFFICIENCY.
C
C     WHILE IN THE PROTOTYPE SOLVER DISCUSSED IN [1] THE NUMBER OF
C     ITERATIONS IS DETERMINED BY A CONVERGENCE CRITERION, TWOSTEP
C     WORKS WITH A FIXED, A PRIORI GIVEN NUMBER OF ITERATIONS. THIS
C     LEADS TO A SOMEWHAT SIMPLER CODE AND IN VARIOUS EXPERIMENTS
C     THIS MORE SIMPLE STRATEGY HAS TURNED OUT TO BE AS EFFICIENT AS
C     THE ITERATION STRATEGY [2]. IN FACT, TWOSTEP MAY WORK WELL WITH
C     ONLY A FEW ITERATIONS PER INTEGRATION STEP. IN THIS SITUATION
C     TWOSTEP HAS A WORKLOAD COMPARABLE TO EXPLICIT SOLVERS BASED
C     ON THE QUASI-STEADY-STATE-APPROACH (QSSA) [3]. NOTE THAT WHEN
C     A FEW GAUSS-SEIDEL ITERATIONS ARE USED, THE ORDER OF THE COMPONENTS
C     GENERALLY WILL INFLUENCE THE ACCURACY. TO OUR EXPERIENCE THIS
C     INFLUENCE IS MINOR.
C
C     [1] J.G. VERWER, GAUSS-SEIDEL ITERATION FOR STIFF ODES FROM
C     CHEMICAL KINETICS, SIAM JOURNAL ON SCIENTIFIC COMPUTING 15,
C     1243 - 1250, 1994.
C
C     [2] J.G. VERWER & D. SIMPSON, A COMPARISON BETWEEN TWO EXPLICIT
C     METHODS FOR STIFF ODES FROM ATMOSPHERIC CHEMISTRY. REPORT
C     NM-R9414, CWI, AMSTERDAM, 1994 (TO APPEAR IN APPL. NUMER. MATH.).
C
C     [3] J.G. VERWER & M. VAN LOON,  AN EVALUATION OF EXPLICIT PSEUDO-
C     STEADY-STATE APPROXIMATION SCHEMES FOR STIFF ODE SYSTEMS FROM
C     CHEMICAL KINETICS, J. COMPUT. PHYS. 113, 347 - 352, 1994.
C
C----------------------------------------------------------------------
C
C     MEANING OF PARAMETERS:
C
C     NUMSP   - INTEGER. NUMBER OF COMPONENTS.
C     T       - REAL. THE INDEPENDENT VARIABLE TIME.
C     TE      - REAL. THE ENDPOINT OF TIME.
C     DT      - REAL. THE STEPSIZE.
C     DTMIN   - REAL. A MINIMUM FOR DT.
C     DTMAX   - REAL. A MAXIMUM FOR DT.
C     YOLD    - REAL ARRAY (N). SOLUTION AT PREVIOUS TIME POINT.
C     Y       - REAL ARRAY (N). SOLUTION AT CURRENT TIME POINT.
C     YNEW    - REAL ARRAY (N). SOLUTION AT FORWARD TIME POINT.
C     YP      - REAL ARRAY (N). STORAGE FOR THE PRODUCTION TERM.
C     YL      - REAL ARRAY (N). STORAGE FOR THE LOSS TERM.
C     YL_ORG  - REAL ARRAY (N). STORAGE FOR THE LOSS TERM due to ORGANIC
C     REACTIONS ONLY.
C     tot     - REAL ARRAY (N). WORKARRAY FOR THE BDF2 METHOD.
C     ATOL    - REAL ARRAY (N). ABSOLUTE TOLERANCES.
C     RTOL    - REAL ARRAY (N). RELATIVE TOLERANCES.
C     NUMIT   - INTEGER. NUMBER OF GAUSS-SEIDEL OR JACOBI ITERATIONS.
C     NFCN    - INTEGER. THE TOTAL NUMBER OF (FUNCTION) CALLS OF ITER.
C     NACCPT  - INTEGER. THE NUMBER OF ACCEPTED INTEGRATION STEPS.
C     NREJEC  - INTEGER. THE NUMBER OF REJECTED INTEGRATION STEPS.
C     NSTART  - INTEGER. THE NUMBER OF RESTARTS + 1.
C
C----------------------------------------------------------------------
C
C     STORAGE: 9 REAL ARRAYS OF DIMENSION N.
C
C----------------------------------------------------------------------
C
C     INPUT:
C
C     NUMSP   - NUMBER OF COMPONENTS.
C     T       - INITIAL TIME; T IS CHANGED.
C     TE      - THE ENDPOINT OF TIME.
C     DTMIN   - THE MINIMAL STEPSIZE THAT TWOSTEP IS ALLOWED TO USE.
C     DTMAX   - THE MAXIMAL STEPSIZE THAT TWOSTEP IS ALLOWED TO USE.
C
C               IF ON INPUT DTMIN = DTMAX, THEN DT:=DTMIN AND DT IS KEPT
C               FIXED THROUGHOUT THE INTEGRATION, POSSIBLY EXCEPT FOR
C               THE FINAL STEP WHERE DT MAY BE ADJUSTED TO PRECISELY HIT
C               THE END POINT TE. SO, THE LENGTH OF THE INTERVAL NEED NOT BE
C               AN INTEGER MULTIPLE OF THE SELECTED STEPSIZE.
C
C               IF DTMIN = DTMAX, THE STEPSIZE CONTROL IS SWITCHED OFF.
C               THE USER THUS SHOULD ASCERTAIN THAT THE INTEGRATION
C               PROCESS REMAINS STABLE FOR THE SELECTED STEPSIZE.
C
C               IF DTMAX > DTMIN, THEN STEPSIZE CONTROL IS CARRIED OUT
C               AND TWOSTEP DETERMINES AN INITIAL STEPSIZE ITSELF. HOWEVER,
C               THE CONTROL AND THE INITIAL STEPSIZE SELECTION CAN BE
C               OVERRULED TO SATISFY THE CONSTRAINT DTMIN <= DT <= DTMAX.
C
C     Y       - INITIAL SOLUTION VECTOR; Y IS CHANGED.
C     ATOL    - ABSOLUTE TOLERANCES.
C     RTOL    - RELATIVE TOLERANCES.
C
C     NUMIT   - THE NUMBER OF GAUSS-SEIDEL (OR JACOBI ITERATIONS) USED
C               PER TIME STEP. THIS NUMBER THUS IS A PRIORI DESCRIBED
C               FOR THE WHOLE INTEGRATION. A LOW NUMBER IS RECOMMENDED
C               FOR GAUSS-SEIDEL ITERATION [1,2].
C
C               NOTE THAT FOR THE GAUSS-SEIDEL TECHNIQUE THE ORDER OF
C               THE COMPONENTS WITHIN ITER WILL PLAY A ROLE. IT IS
C               RECOMMENDED TO ORDER THE COMPONENTS IN DECREASING
C               OF THE LOSS RATES FOR THE VERY SHORT LIVING SPECIES
C               LIKE THE RADICALS.
C
C----------------------------------------------------------------------
C
C     OUTPUT:
C
C     T       - T = TE.
C     DTOLD   - THE LAST STEPSIZE VALUE USED WHEN WORKING WITH
C               VARIABLE STEPSIZES.
C     DT      - THE LAST STEPSIZE USED TO HIT THE ENDPOINT TE
C               WHEN WORKING WITH CONSTANT STEPSIZES.
C     Y       - THE COMPUTED SOLUTION AT T = TE.
C     NFCN, NACCPT, NREJEC, NSTART, STARTDT - SEE MEANING OF PARAMETERS.
C
C----------------------------------------------------------------------
C
C     TRACER INFO:
C
C     ntr : number of tracers for which information is required
C     idtr(mtr) : i.d. numbers of these tracers
C     trprod(maxsp) : total production rate of each tracer (molec/cc/s)
C     trloss(maxsp) : total loss ratw of each tracer (molec/cc/s)
C
C----------------------------------------------------------------------
C
C     SUBROUTINES:
C
C     TWOSTEP CALLS THREE SUBROUTINES, VIZ. NEWDT, FIT AND ITER. ONLY
C     SUBROUTINE ITER IS TO BE DEFINED BY THE USER.
C
C     NEWDT  - COMPUTES THE NEW STEPSIZE. NEWDT ITSELF CALLS FIT.
C     FIT    - MAY ADJUST DT TO GUARANTEE THAT THE REMAINDER
C              OF THE INTEGRATION INTERVAL IS AN INTEGER MULTIPLE OF
C              THE CURRENT STEPSIZE. THE ADJUSTMENT IS CARRIED AS SOON
C              AS (TE-T)/DT <= 10.0. HENCE THE ADJUSTMENT MAY LEAD TO
C              A STEPSIZE SMALLER THAN DTMIN FOR APPROXIMATELY
C              TEN INTEGRATION STEPS.
C     ITER   - A USER DEFINED ROUTINE FOR THE ODE SYSTEM. WITHIN ITER
C              ALSO THE GAUSS-SEIDEL OR JACOBI TECHNIQUE IS TO BE
C              IMPLEMENTED BY THE USER, AS EXEMPLIFIED BELOW:
C
C----------------------------------------------------------------------
C

       USE prodloss_module
       IMPLICIT NONE
C
C ----------------------------------------------------------------------
C 1.1 VARIABLES D INTERFACE
C ----------------------------------------------------------------------
C
C n             I    INPUT     NB D ESPECES
C maxsp         I    INPUT     NB MAX D ESPECES
C maxre         I    INPUT     NB MAX DE REACTIONS
C numre         I    INPUT     NB DE REACTIONS
C MAXSTOI       I    INPUT     NB MAX DE COEFF. STOECHIO
C numsp         I    INPUT     NB D ESPECES
C t             R    INPUT     TEMPS DE DEBUT
C te            R    INPUT
C dtmin, dtmax  R    INPUT     PAS DE TEMPS MIM, MAX
C anrat         R    INPUT
C stoicf        R    INPUT     STOICHIO. COEFFICIENTS FOR CVAR REACTIONS
C numit         I    INPUT     NB D ITERATIONS POUR SEIDEL OU JACOBI (= 2)
C startdt       R    INPUT     PAS DE TEMPS MIN
C atol, rtol    R    INPUT     TOLERANCE ABSOLUTE, RELATIVE
C mtr           I    INPUT     NB MAX DE TRACERS
C ntr           I    INPUT     NB DE TRACERS
C idtr          I    INPUT     I.D. DE TRACERS
C y             R    IN/OUTPUT CONCENTRATION
C dt            R    OUTPUT
C trloss,trprod R    OUTPUT    TRACER loss / production rates

C tot           R    LOCAL
C yold          R    LOCAL
C ynew          R    LOCAL
C yl            R    LOCAL
C yp            R    LOCAL
C
C
       INTEGER maxsp, maxre, numsp, numre
       INTEGER numit
       INTEGER mxleft, mxright
       INTEGER mself
       INTEGER numstoi(maxre,2)
       INTEGER idrestoi(maxre,mxleft)
       INTEGER idpdstoi(maxre,mxright)
       INTEGER nself,idselfreac(mself,2)
* arrays for tracer production/loss rate tracking
      INTEGER  mtr,ntr
      INTEGER  idtr(mtr)
      REAL     trprod(maxsp),trloss(maxsp)
* constrained concentrations
      INTEGER  ncons
      INTEGER  conid(ncons)
      INTEGER  noxfix, idno, idno2, idno3
      REAL     sumnox
! index map
      TYPE(spec_reac_map),TARGET, intent(in)  :: lpmap(maxsp)

c       REAL atol(maxsp), rtol(maxsp)
       REAL atol, rtol
       REAL t, te, dtmin, dtmax, dt
       REAL anrat(maxre), y(maxsp)
       REAL restoicf(maxre,mxleft)
       REAL pdstoicf(maxre,mxright)
       REAL tot(maxsp), yold(maxsp)
       REAL yl(maxsp), yp(maxsp)
       REAL tl(maxsp), tp(maxsp)
       REAL rem(maxsp),rdep(maxsp)
       REAL rex(maxsp),rdil
       REAL errt,memerrt
       INTEGER memid
       REAL ynew(maxsp)

C
C ----------------------------------------------------------------------
C 1.2 VARIABLES LOCALES
C ----------------------------------------------------------------------
C
      INTEGER i, j, k
      INTEGER nfcn, naccpt, nrejec
      INTEGER nstart
! JMLT: max orders of magnitude available for DT calculation
      INTEGER,PARAMETER :: MACH = 38
C
      REAL dtold, errlte, dtg
      REAL ratio, ytol, dy
      REAL a1, a2, c, cp1
      REAL zero
C
      LOGICAL accept, restart, failer
C
C ***************************************************************************
C ETAPE 2 :
C ***************************************************************************
!      WRITE(6,*) 'starting twostep'
C --------------------------------------------------------------------------
C     INITIALIZATION OF COUNTERS, ETC.
C -------------------------------------------------------------------------
C
      NACCPT=0
      NREJEC=0
      NFCN=0
      NSTART=0
      FAILER=.FALSE.
      RESTART=.FALSE.
      ACCEPT=.TRUE.

      trprod = 0.
      trloss = 0.

C
C --------------------------------------------------------------------------
C     INITIAL STEPSIZE COMPUTATION.
C --------------------------------------------------------------------------
C
      IF (DTMIN.EQ.DTMAX) THEN
       NSTART=1
       DT=MIN(DTMIN,(TE-T)/2)
       GOTO 28
      ENDIF
C
      tot =Y
      CALL iter4
     1   (maxsp, maxre, mxleft, mxright, mself,
     2    numsp, numre, numstoi, idrestoi, idpdstoi, nself,idselfreac,
     3    restoicf,pdstoicf,
     3    anrat, tot, 0.,
     4    rem, rdep, rex, rdil,
     5    y, yp, yl, tp, tl, conid,ncons,
     6    noxfix, sumnox, idno, idno2, idno3,
     7    lpmap)

      NFCN=NFCN+1
      DT=TE-T
C
      DO 20 I=1,numsp
c       YTOL=ATOL(I)+RTOL(I)*ABS(Y(I))
       YTOL=ATOL+RTOL*ABS(Y(I))
       DY=YP(I)-Y(I)*YL(I)
       IF (DY.NE.0.0) THEN
!JMLT added condition to prevent FPEs for very large DT values
!      PRINT*,YTOL,ABS(LOG10(YTOL))+ABS(LOG10(ABS(DY))),DT
         IF(ABS(LOG10(YTOL))+ABS(LOG10(ABS(DY)))
     &      .GT.MACH)THEN
           DT = DT
         ELSE
           DT=MIN(DT,YTOL/ABS(DY))
         ENDIF
       ENDIF
   20 CONTINUE
C
   25 CONTINUE
      NSTART=NSTART+1
      IF (RESTART) DT=DT/10.0
      RESTART=.TRUE.
      DT=MAX(DTMIN,MIN(DT,DTMAX))
C
      CALL FIT (T,TE,DT)
C
      DT=MIN(DT,(TE-T)/2)
C
C --------------------------------------------------------------------------
C     THE STARTING STEP IS CARRIED OUT, USING THE IMPLICIT EULER METHOD.
C --------------------------------------------------------------------------
C
   28 CONTINUE
!      DO 30 I=1,numsp
!       YNEW(I)=Y(I)
!       YOLD(I)=Y(I)
!       tot(I)=Y(I)
!   30 CONTINUE
!JMLT!
      YNEW(1:numsp)=Y(1:numsp)
      YOLD(1:numsp)=Y(1:numsp)
      tot(1:numsp)=Y(1:numsp)
C
      DO 40 I=1,NUMIT
C
      CALL iter4
     1   (maxsp, maxre, mxleft, mxright, mself,
     2    numsp, numre, numstoi, idrestoi, idpdstoi, nself,idselfreac,
     3    restoicf,pdstoicf,
     3    anrat, tot, dt,
     4    rem, rdep, rex, rdil,
     5    ynew, yp, yl, tp, tl, conid,ncons,
     6    noxfix, sumnox, idno, idno2, idno3,
     7    lpmap)
C
       NFCN=NFCN+1
   40 CONTINUE
C
      NACCPT=NACCPT+1
      T=T+DT
C
!      DO 50 J=1,numsp
!       Y(J)=YNEW(J)
!   50 CONTINUE
!JMLT!
      Y(1:numsp)=YNEW(1:numsp)
C
C --------------------------------------------------------------------------
C     SUBSEQUENT STEPS ARE CARRIED OUT WITH THE TWO-STEP BDF METHOD.
C --------------------------------------------------------------------------
C
      DTOLD=DT
      RATIO=1.0
   60 CONTINUE
      C=1.0/RATIO
      CP1=C+1.0
      A1=((C+1.0)**2)/(C*C+2.0*C)
      A2=-1.0/(C*C+2.0*C)
      DTG=DT*(1.0+C)/(2.0+C)
      ZERO = 0.
C
!      DO 70 J=1,numsp
!       tot(J)=A1*Y(J)+A2*YOLD(J)
!       YNEW(J)=MAX(ZERO,Y(J)+RATIO*(Y(J)-YOLD(J)))
!   70 CONTINUE
!JMLT!
      tot(1:numsp)=A1*Y(1:numsp)+A2*YOLD(1:numsp)
       YNEW(1:numsp)
     & = MAX(ZERO,Y(1:numsp)+RATIO*(Y(1:numsp)-YOLD(1:numsp)))
C
      DO 80 I=1,NUMIT
C
      CALL iter4
     1   (maxsp, maxre, mxleft, mxright, mself,
     2    numsp, numre, numstoi, idrestoi, idpdstoi, nself,idselfreac,
     3    restoicf,pdstoicf,
     3    anrat, tot, dtg,
     4    rem, rdep, rex, rdil,
     5    ynew, yp, yl, tp, tl, conid,ncons,
     6    noxfix, sumnox, idno, idno2, idno3,
     7    lpmap)
C
       NFCN=NFCN+1
C
   80 CONTINUE
C
C
C --------------------------------------------------------------------------
C     IF STEPSIZES SHOULD REMAIN EQUAL, STEPSIZE CONTROL IS OMITTED.
C --------------------------------------------------------------------------
C
      IF (DTMIN.EQ.DTMAX) THEN
C
       T=T+DTOLD
       NACCPT=NACCPT+1
C
!       DO 85 J=1,numsp
!        YOLD(J)=Y(J)
!        Y(J)=YNEW(J)
!   85  CONTINUE
!JMLT!
        YOLD(1:numsp)=Y(1:numsp)
        Y(1:numsp)=YNEW(1:numsp)
C
       IF (DT.NE.DTOLD) THEN
        T=T-DTOLD+DT
        GOTO 120
       ENDIF
C
       DT=MIN(DTOLD,TE-T)
       RATIO=DT/DTOLD
       IF (T.GE.TE) GOTO 120
       GOTO 60
      ENDIF
C
C --------------------------------------------------------------------------
C     OTHERWISE STEPSIZE CONTROL IS CARRIED OUT.
C --------------------------------------------------------------------------
C
      ERRLTE=0.0
C
cba      DO 90 I=1,numsp
cba       YTOL=ATOL(I)+RTOL(I)*ABS(Y(I))
cba       ERRLTE=MAX(ERRLTE,ABS(C*YNEW(I)-CP1*Y(I)+YOLD(I))/YTOL)
cba   90 CONTINUE
      DO 90 I=1,numsp
c       YTOL=ATOL(I)+RTOL(I)*ABS(Y(I))
       YTOL=ATOL+RTOL*ABS(Y(I))
       errt=ABS(C*YNEW(I)-CP1*Y(I)+YOLD(I))/YTOL
       IF (errt.gt.errlte) THEN
         errlte=errt
         memerrt=errt
         memid=i
       ENDIF
   90 CONTINUE
! TROUBLESHOOTING !
!      write(55,'(3(E13.5),i7)') t,dt,memerrt,memid

C
      ERRLTE=2.0*ERRLTE/(C+C*C)
C
      CALL NEWDT (T, TE, DTOLD, ERRLTE, DTMIN, DTMAX,
     1           DT, RATIO, ACCEPT)
C
C --------------------------------------------------------------------------
C     HERE THE STEP HAS BEEN ACCEPTED.
C --------------------------------------------------------------------------
C
      IF (ACCEPT) THEN
C
       FAILER=.FALSE.
       RESTART=.FALSE.
       T=T+DTOLD
       NACCPT=NACCPT+1
C
!       DO 100 J=1,numsp
!        YOLD(J)=Y(J)
!        Y(J)=YNEW(J)
!  100  CONTINUE
!JMLT!
        YOLD(1:numsp)=Y(1:numsp)
        Y(1:numsp)=YNEW(1:numsp)
C
!JMLT: tracer rates!
        DO k=1,ntr
! all rates
!          trprod(idtr(k))=YP(idtr(k))
!          trloss(idtr(k))=YL(idtr(k))*Y(idtr(k))
! rates for organic reactions only
          trprod(idtr(k))=TP(idtr(k))
          trloss(idtr(k))=TL(idtr(k))*Y(idtr(k))
        ENDDO
C
       IF (T.GE.TE) GOTO 120
       GOTO 60
      ENDIF
C
C --------------------------------------------------------------------------
C     A RESTART CHECK IS CARRIED OUT.
C --------------------------------------------------------------------------
C
      IF (FAILER) THEN
C
       NREJEC=NREJEC+1
       FAILER=.FALSE.
       NACCPT=NACCPT-1
       T=T-DTOLD
       WRITE(6,*) 'going back in time'
C
!       DO 110 J=1,numsp
!        Y(J)=YOLD(J)
!  110  CONTINUE
!JMLT!
        Y(1:numsp)=YOLD(1:numsp)
C
       GOTO 25
      ENDIF
C
C --------------------------------------------------------------------------
C     HERE THE STEP HAS BEEN REJECTED.
C --------------------------------------------------------------------------
C
      NREJEC=NREJEC+1
      FAILER=.TRUE.
      GOTO 60
C
  120 CONTINUE
      !PRINT*,Y
C
      END
