C
      SUBROUTINE NEWDT (T, TE, DTOLD, ERRLTE, DTMIN, DTMAX,
     1                 DT, RATIO, ACCEPT)
C
C **********************************************************************
C PROGRAMMES APPELANTS : TWOSTEP
C PROGRAMMES APPELES   : FIT 
C **********************************************************************
C
C
C **********************************************************************
C ETAPE 1 : DECLARATION
C **********************************************************************
C
      IMPLICIT NONE
C
C ---------------------------------------------------------------------
C 1.1 VARIABLES D INTERFACE
C --------------------------------------------------------------------
C
C ERRLTE      R     INPUT
C DTMIN       R     INPUT
C DTMAX       R     INPUT
C TE          R     INPUT
C T           R     INPUT
C DTOLD       R     INPUT
C DT          R     INPUT/OUTPUT
C RATIO       R     OUTPUT
C ACCEPT      L     OUTPUT
C
       REAL ERRLTE, DTMIN, T, TE, DTMAX, DTOLD
       REAL DT, RATIO
       LOGICAL ACCEPT
C
C --------------------------------------------------------------------
C 1.2 VARIABLES LOCALES
C --------------------------------------------------------------------
C
      REAL TS
      REAL ZERO5, DEUX, ZERO8
C
C *************************************************************************
C ETAPE 2 :
C ************************************************************************
C
      IF (ERRLTE.GT.1.0.AND.DT.GT.DTMIN) THEN
       ACCEPT=.FALSE.
       TS=T
      ELSE
       ACCEPT=.TRUE.
       DTOLD=DT
       TS=T+DTOLD
      ENDIF
C
      ZERO5 = 0.5
      DEUX = 2.
      ZERO8 = 0.8
      DT=MAX(ZERO5,MIN(DEUX,ZERO8/SQRT(ERRLTE)))*DT
      DT=MAX(DTMIN,MIN(DT,DTMAX))
C
      CALL FIT (TS,TE,DT)
C
      RATIO=DT/DTOLD
C
      END
