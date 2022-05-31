************************************************************************
* MASTER MECHANISM - ROUTINE NAME : gettrack                           *
*                                                                      *
* PURPOSE :                                                            *
*  Browse the bond matrix to find all possible node at a given         *
*  position (npos), starting from a given node (top). Ring are allowed *
* (track is stopped when a full loop along the circle is made).        *
*                                                                      *
* The subroutine call "gettrack" which give all the possible track,    *
* starting from top.                                                   *
*  track(*,2) give all nodes in alpha position regarding node top      *
*  track(*,3) give all nodes in beta position regarding node top       *
*  track(*,4) give all nodes in gamma position regarding node top      *
*  etc...                                                              *
*                                                                      *
*                                                                      *
* INPUT:                                                               *
* - bond(i,j)    : carbon-carbon bond matrix                           *
* - nca          : number of group                                     *
* - top          : starting node number                                *
* - npos         : position to be searched for (n=2 is alpha position  *
*                  relative to top)                                    *
*                                                                      *
* OUTPUT:                                                              *
* - nnod         : number of distinct node in position "npos" relative *
*                  to top                                              * 
* - tnod(i)      : table of all nodes at position "npos"               *
************************************************************************
      SUBROUTINE nodmap(bond,top,nca,npos,nnod,tnod)
      IMPLICIT NONE
      INCLUDE 'general.h'

* input
      INTEGER  top
      INTEGER  nca
      INTEGER  bond(mca,mca)
      INTEGER  npos

* output:
      INTEGER  nnod, tnod(mca)

* internal:
      INTEGER  track(mco,mca)
      INTEGER  trlen(mco)
      INTEGER  ntr

      INTEGER  i,j,k

* initialize
* ----------
      nnod=0
      DO i=1,mca
        tnod(i)=0
      ENDDO

      CALL gettrack(bond,top,nca,ntr,track,trlen)

* avoid duplicate - set 0 to duplicate nodes
      IF (ntr.gt.1) THEN
        DO 10 i=1,ntr-1
          IF (track(i,npos).eq.0) GOTO 10
          DO j=i+1,ntr
            IF (track(j,npos).ne.0) THEN
              IF (track(i,npos).eq.track(j,npos)) THEN
                track(j,npos)=0
              ENDIF
            ENDIF  
          ENDDO
10        CONTINUE
      ENDIF
       
* get nodes
      DO i=1,ntr
        IF (track(i,npos).ne.0) THEN
          nnod=nnod+1
          tnod(nnod)=track(i,npos)
        ENDIF
      ENDDO

      END
