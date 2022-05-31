************************************************************************
* MASTER MECHANISM - ROUTINE NAME : lntree                             *
*                                                                      *
* PURPOSE :                                                            *
*  Browse the bond matrix to find all possible tracks starting from    *
*  a given node (top). Ring are allowed (track is stopped when a full  *
*  loop along the circle is made).                                     *
*                                                                      *
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
*                                                                      *
* OUTPUT:                                                              *
* - ntr          : number of distinct track found starting at "top"    *
* - track(i,j)   : possible track starting at top. Index i is for the  *
*                  for the ith track. Index j are the nodes for the    *
*                  track i.                                            *
* - trlen(mco)   : length of ith track found                           *
************************************************************************
      SUBROUTINE gettrack(bond,top,nca,ntr,track,trlen)
      IMPLICIT NONE
      INCLUDE 'general.h'

* input
      INTEGER  top
      INTEGER  nca
      INTEGER  bond(mca,mca)

* output:
      INTEGER  track(mco,mca)
      INTEGER  trlen(mco)
      INTEGER  ntr

* internal:
      INTEGER  ptr,niv,nod
      INTEGER  i,j,k
      INTEGER  memo(mca)
      INTEGER  slope

* -----------
* initialize
* -----------
      ptr = 0
      DO i=1,mco
        trlen(i)=0
        DO j=1,mca
          track(i,j)=0
        ENDDO
      ENDDO
      DO i=1,mca
        memo(i)=0
      ENDDO

* initialize parameters to find the tracks 
      track(1,1)=top
      memo(1)=top
      niv=1     ! number of node since starting node (i.e. top)
      ptr=0     ! current "search" pointer (i.e. to find next node in the track)
      ntr=1    ! number of track found - current track
      nod=top   ! current node along the track
      slope=1   ! equal 1 when going forward along the track, otherwise 0

* -----------
* start loop
* -----------

* get next - reentry point
10    CONTINUE
      ptr=ptr+1

* end of line reached - must go backward or exit
      IF (ptr.gt.nca) THEN      
        IF (niv.eq.1) THEN       ! all the possible track are found - exit
          GOTO 100
        ELSE                     ! go backward
          ptr=memo(niv)          ! set pointer to previous memo
          memo(niv)=0            ! reset memo   
          niv=niv-1              ! decrease niv (go backward)
          nod=memo(niv)          ! set new current node 

          IF (slope.ne.0) THEN   ! make a new track (if required)
            ntr=ntr+1
            IF (ntr.gt.mco) THEN
              STOP '--error-- in geattrack. ntr is greater than mco'
            ENDIF
            DO i=1,niv
                 track(ntr,i)=track(ntr-1,i)
            ENDDO
          ENDIF
             track(ntr,niv+1)=0     ! remove previous track

          slope=0                ! remember ... I am now going backward
          GOTO 10
        ENDIF
      ENDIF

* no bond between current node (nod) and next possible node (ptr)
      IF (bond(ptr,nod).eq.0) GOTO 10   ! get next

* new bond found - must go forward
      IF (bond(ptr,nod).ne.0) THEN
        DO k=1,niv
          IF (ptr.eq.memo(k)) GOTO 10   ! end circle or previous node
        ENDDO
        niv=niv+1                       ! increase niv (go forward)
        track(ntr,niv)=ptr              ! keep track
        nod=ptr                         ! set current node to the new node
        memo(niv)=nod                   ! remember which track was used
        ptr=0                           ! set pointer to 0 (to find next node)
        slope=1                         ! remember ... I am going forward
        GOTO 10                         ! go find next node
      ENDIF

100   CONTINUE

* remove last "case" (fill with top only) and clean track
      DO i=1,mca
        track(ntr,i)=0
      ENDDO
      ntr=ntr-1

* compute the length of each track 
      DO 200 i=1,ntr
        DO j=mca,1,-1
          IF (track(i,j).ne.0) THEN
            trlen(i)=j
            GOTO 200
          ENDIF
        ENDDO
200   CONTINUE

* end of gettrack
      RETURN
      END
