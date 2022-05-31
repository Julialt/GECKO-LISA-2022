* INCA : subroutine chkreac
*
* Purpose : 
* Check the reaction in the chemical scheme, especially look
* for trivial and duplicate reaction). TAKE CARE ABOUT NCONTROL BELOW
*
* INPUT :
*  lout      : unit number for the output file
*  maxre     : maximum number of reactions
*  maxstoi   : maximum number of species in a given reaction
*  num_m     : number of +M reactions
*  numextra  : number of EXTRA reactions
*  numhv     : number of HV reactions
*  numre     : total number of reactions
*  numcvar   : number of cvar reactions
*  numo2     : number of O2 reaction
*  numiso    : number of isomerisation reaction
*  itype(i)  : type gor the ith reaction (3=FO,6=HV,8=CVAR,7=EXTRA)
*  idrestoi(i,k) : id number of the kth species involved in reaction i
*                  at the reactant side.
*  idpdstoi(i,k) : id number of the kth species involved in reaction i
*                  at the product side.
*
* OUTPUT
*  lostop    : logical to stop program if error
*  locheck(i): logical check fot ith reaction
*  num_n     : total number of simple thermal reaction
*  id_n(i)   : reaction number for the ith thermal reaction
**********************************************************************

      SUBROUTINE chkreac (
     1    lout, 
     2    maxre, mxleft, mxright, mxrpero, maxro2,numiso,
     3    num_n, num_m, numfo, numhv, numre, numo2, nummeO2, nrpero,
     3    mxrdimer,nrdimer,maxdimer,
     4    numcvar, numextra, itype, id_n,
     4    numain, numaou, numwin, numwou,
     5    numstoi, idrestoi, idpdstoi,
     6    lostop, locheck)

      IMPLICIT NONE

* input/output
      INTEGER lout
      INTEGER maxre, mxleft, mxright, mxrpero, maxro2,mxrdimer,maxdimer
      INTEGER num_n, num_m, numfo, numhv, numo2, nummeo2
      INTEGER numcvar, numextra, numre,numiso
      INTEGER numain, numaou, numwin, numwou
      INTEGER itype(maxre), id_n(maxre)
      INTEGER numstoi(maxre,2)
c      INTEGER idstoi(maxre,maxstoi,2)
      INTEGER idrestoi(maxre,mxleft)
      INTEGER idpdstoi(maxre,mxright)
      !INTEGER nrpero(mxrpero),nrdimer(mxrdimer)
      INTEGER nrpero(maxro2),nrdimer(maxdimer)
      LOGICAL lostop,locheck(maxre)

* local
      INTEGER i, j, k, l, ncontrol
      LOGICAL loequal


* initialize
      loequal = .false.

* to save time in the interpretation of hugge chemical
* scheme, the following control has been ignored.
* turn the ncontrol variable to 0 to performe the tests

      ncontrol=1
      IF (ncontrol.eq.1) THEN
        WRITE(lout,*) 'THE SET OF REACTION WAS NOT CHECKED'
        GOTO 7000
      ENDIF

* -----------------------------------------------
* SEARCH TRIVIAL REACTION (REACTANTS = PRODUCTS)
* -----------------------------------------------

      DO 5200 i=1,numre
        IF (locheck(i)) GOTO 5200
        IF (numstoi(i,1).EQ.numstoi(i,2)) THEN
          DO j=1,numstoi(i,1)
            IF (idrestoi(i,j).ne.idpdstoi(i,j)) GOTO 5200
          ENDDO 
        ELSE
          GOTO 5200
        ENDIF

        WRITE(lout,*)
        WRITE(lout,*)'   --note--   equation number ',i,' is trivial'
        WRITE(lout,*)
        locheck(i)=.true.
5200  CONTINUE

* ----------------------------------------
* SEARCH FOR IDENTICAL REACTION EQUATIONS
* ----------------------------------------

      DO 5300 i=1,numre-1
        IF (locheck(i)) GOTO 5300
        DO 5290 j=i+1,numre
          IF (locheck(j)) GOTO 5290

* check if reaction i is identical to reaction j. 
          IF (numstoi(i,1).eq.numstoi(j,1).and.
     &        numstoi(i,2).eq.numstoi(j,2))     THEN
            DO l=1,numstoi(i,1)
              IF (idrestoi(i,l).ne.idrestoi(j,l)) GOTO 5250
            ENDDO
            DO l=1,numstoi(i,2)
              IF (idpdstoi(i,l).ne.idpdstoi(j,l)) GOTO 5250
            ENDDO
            loequal=.true.
            GOTO 5280
          ENDIF

* check if reaction i is the revers of reaction j. Remove comment
* if this kind of check is wanted - do the change for idstoi, since
* idstoi is not used anymore
5250      CONTINUE
C          IF (numstoi(i,1).eq.numstoi(j,2).and.
C     &        numstoi(i,2).eq.numstoi(j,1)) THEN
C            DO k=1,2
C              DO l=1,numstoi(i,k)
C                IF (idstoi(i,l,k).ne.idstoi(j,l,3-k)) GOTO 5290
C              ENDDO
C            ENDDO
C          ELSE
            GOTO 5290
C          ENDIF

* check if identical reactions belong to the same reaction type. 
5280      CONTINUE
          IF (loequal) THEN
            loequal=.false.
            IF (itype(i).eq.itype(j)) THEN
               IF (itype(i).eq.7) THEN
                 WRITE(lout,*)
                 WRITE(lout,*)'   --note--    equations number ',i,
     &                     ' and ',j
                 WRITE(lout,*)'               are identical except for',
     &                     '                  "extra" number'
                 GOTO 5290
               ENDIF
              WRITE(lout,*)
c2s              WRITE(lout,*)'   --error--   equations number ',i,
              WRITE(lout,*)'   --note--   equations number ',i,
     &                     ' and ',j
              WRITE(lout,*)'               are identical'
              WRITE(lout,*)
c              lostop=.true.

* in a previous version, if 2 reactions were identical but not of the
* same then a warning was given. Remove the comment if the warning is
* wanted 
C            ELSE
C              WRITE(lout,*)
C              WRITE(lout,*)'   --note--    equations number ',i,
C     &                     ' and ',j
C              WRITE(lout,*)'               are identical except for',
C     &                     ' "m", "hv", "cvar", "extra"'
C              WRITE(lout,*)'               or fall-off type'
C              WRITE(lout,*)
            ENDIF

* following statement are for the case i is the reverse of j. Remove
* the comment if this kind of check is wanted (don't forget to remove
* the comment above).
C          ELSE
C            IF(itype(i).eq.itype(j))THEN
C              WRITE(lout,*)
C              WRITE(lout,*)'   --note--    equation number ',i,
C     &                     ' is the revers of'
C              WRITE(lout,*)'               equation number ',j
C              WRITE(lout,*)
C            ELSE
C              WRITE(lout,*)
C              WRITE(lout,*)'   --note--    equation number ',i,
C     &                     ' is the revers of'
C              WRITE(lout,*)'               equation number ',j,
C     &                 ' except for "m", "hv", "cvar", "extra"'
C              WRITE(lout,*)'               or fall-off type'
C              WRITE(lout,*)
C            ENDIF
          ENDIF

5290    CONTINUE
5300  CONTINUE

* -----------------------
* FINAL CHECK 
* -----------------------

* get the number of "simple" thermal reaction (store data)
7000  CONTINUE
      num_n=0
      DO i=1,numre
        IF (itype(i).eq.0) THEN
          num_n=num_n+1
          id_n(num_n)=i
        ENDIF
      ENDDO


* check that nothing was forgotten ...
      WRITE(lout,*)
      WRITE(lout,*)'   num_n=',num_n
      WRITE(lout,*)'   num_m=',num_m
      WRITE(lout,*)'   numfo=',numfo
      WRITE(lout,*)'   numhv=',numhv
      WRITE(lout,*)' numcvar=',numcvar
      WRITE(lout,*)'numextra=',numextra
      WRITE(lout,*)'   numo2=',numo2
      WRITE(lout,*)'  numiso=',numiso
      WRITE(lout,*)' nummeo2=',nummeo2
      DO i=1,maxro2
        WRITE(lout,*)'numpero',i,'=',nrpero(i)
      ENDDO
      DO i=1,maxdimer
        WRITE(lout,*)'numdimer',i,'=',nrdimer(i)
      ENDDO
      WRITE(lout,*)'numain=',numain
      WRITE(lout,*)'numaou=',numaou
      WRITE(lout,*)'numwin=',numwin
      WRITE(lout,*)'numwou=',numwou
      
      WRITE(lout,*)'------------------------------------------------'
      WRITE(lout,*)'   numre=',numre
      WRITE(lout,*)

      ncontrol=num_n+num_m+numfo+numhv+numcvar+numextra+numo2+nummeo2
     &         +numiso
      DO i=1,maxro2
        ncontrol=ncontrol+nrpero(i)
      ENDDO
      DO i=1,maxdimer
        ncontrol=ncontrol+nrdimer(i)
      ENDDO
      ncontrol=ncontrol+numain+numaou+numwin+numwou
      IF (ncontrol .NE. numre) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   in the equation count'
        WRITE(lout,*)
        lostop=.true.
      ENDIF

      END
