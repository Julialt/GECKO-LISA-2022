************************************************************************
* MASTER MECHANISM - ROUTINE NAME : grbond                             *
*                                                                      *
* PURPOSE: Contruction of the group matrix and the carbon-carbon       *
*          bond matrix of a molecule given as input.                   *
*                                                                      *
* INPUT:                                                               *
* - chem         : chemical formula                                    *
* - nc           : number of characters in chem                        *
*                                                                      *
* OUTPUT:                                                              *
* - group(i)     : groups at position (carbon) i                       *
* - bond(i,j)    : carbon-carbon bond matrix of CHEM                   *
* - dbflg        : flag for double-bond indication                     *
* - nring        : number of separate rings in CHEM                    *
*                                                                      *
* First groups are constructed. A group starts with a carbon "C", an   *
* aromatic carbon "c", or an ether bond "-O-",                         *
* and ends with either a next "C" or "c", or a                         *
* double-bond "=", or an ether bond "-O-", or  blank "NUL". In each    *
* group the trailing parentheses are deleted if it is an opening       *
* with no closing and vice-versa.                                      * 
* Next the carbon skeleton is built with carbon, oxygen (if ether)     *
* and remaining parentheses and double-bonds.                          *
* The bond matrix is then built based on the skeleton.                 *
* The program checks the valence for each carbon center.               *
*                                                                      *
************************************************************************
      SUBROUTINE grbond(chem,nc,group,bond,dbflg,nring)
      IMPLICIT NONE
      INCLUDE 'general.h'
      INCLUDE 'common.h'

* input:
      CHARACTER*(lfo) chem
      INTEGER         nc

* output:
      CHARACTER*(lgr) group(mca)
      INTEGER         bond(mca,mca), dbflg
      INTEGER         ring(mca,mca), nring

* internal:
      CHARACTER*(lfo) skelet
      INTEGER         ig,ng,ik,istop,nb,tp,p
      INTEGER         i,j,k,l,nok,ii,jj,kk,ir,n,nn
      INTEGER         dbbrk ! flag for chem where SMILES contains 'C=1'
      INTEGER         branch(mca)
      LOGICAL         loeter

* ----------
* INITIALIZE
* ----------
      !print*,'*grbond*',chem
      ig = 0
      ik = 0
      dbflg = 0
      nring = 0
      dbbrk = 0
      istop = 0
      l = 0
      skelet = ' '
      loeter=.false.

      group(:) = ' '

      bond(:,:) = 0
      ring(:,:) = 0

      IF (chem(1:1).NE.'C' .AND. chem(1:1).NE.'c' .AND. 
     &    chem(1:2).NE.'-O')THEN
        IF (chem(1:4).EQ.'=Cd1') THEN
          dbbrk=1
        ELSE
          WRITE(6,*)'WARNING : input to grbond must begin with C or c'
          WRITE(6,*)'begin with C/c:',chem
          STOP
        ENDIF
      ENDIF

* -------------------------------------------------------------
* CONSTRUCT CARBON-CENTERED GROUPS, SKELETON AND DIAGONAL TERM
* -------------------------------------------------------------
      DO 10 i=1,nc

* identify group
* --------------

* find where the group starts
        nok=0
        IF (chem(i:i).EQ.'C'.OR.chem(i:i).EQ.'c') nok=1
        IF (chem(i:i+1).EQ.'Cl') nok=0
        IF (chem(i:i+1).EQ.'-O') nok=1
        IF (nok.EQ.0) GO TO 10  ! cycle

        ig = ig + 1

* find where the group ends
        DO j=i+1,i+lgr
          nok=0
          IF (chem(j:j).EQ.'C'.OR.chem(j:j).EQ.'c')    nok=1
          IF (chem(j:j+1).EQ.'Cl')   nok=0
          IF (chem(j:j)  .EQ.'=')    nok=1
          IF (chem(j:j+1).EQ.'-O')  nok=1
          IF (chem(j:j)  .EQ.' ')    nok=1
          IF (nok.EQ.1) THEN
            l = j - i 
            istop = j - 1
            GO TO 20  ! exit
          ENDIF
        ENDDO

* write out group               
20      group(ig) = chem(i:istop)

* delete trailing parentheses from group
100     CONTINUE
        IF (group(ig)(l:l).EQ.'(') THEN
           group(ig)(l:l) = ' '
           l = l - 1
           GOTO 100
        ELSE IF (group(ig)(l:l) .EQ. ')') THEN
           p = 0
           DO j=1,l
             if (group(ig)(j:j) .eq. '(') p = p + 1
             if (group(ig)(j:j) .eq. ')') p = p - 1
           ENDDO
           IF (p .LT. 0) THEN
              group(ig)(l:l) = ' '
              l = l - 1
              GOTO 100
           ENDIF
        ENDIF

* build skeleton
* --------------
        ik = ik + 1
* skeletal oxygen
        IF (group(ig)(1:2).eq.'-O') THEN
          skelet(ik:ik) = 'O'
          j=3
* ring-joining oxygen (Oxygen is never a multiple center)
          DO n=1,mri
            IF(group(ig)(3:3).EQ.digit(n))THEN
              nring = nring+1
              ik=ik+1
              skelet(ik:ik)=digit(n)
              j=4
            ENDIF
          ENDDO
          group(ig)(j:j) = '-'
          loeter=.true.
        ELSE
* carbon (either flavor)
          IF(group(ig)(1:1).EQ.'c'.OR.group(ig)(1:1).EQ.'C')THEN
            skelet(ik:ik) = group(ig)(1:1)
            jj=2
            IF(group(ig)(2:2).EQ.'d') jj=3
*  ring-joining carbon
            DO 110 j=jj,l
              DO n=1,mri
                IF(group(ig)(j:j).EQ.digit(n))THEN
                  ik=ik+1
                  skelet(ik:ik)=group(ig)(j:j)
                  nring = nring+1
                  GOTO 110 ! try next character (in case multiple center)
                ENDIF
              ENDDO
! If reach this point, no more ring-joiners on group  
              GOTO 111  
110           CONTINUE
111         CONTINUE
          ENDIF

        ENDIF
        
* put parenthesis if any
        DO j=i+l,istop
          ik = ik + 1
          skelet(ik:ik) = chem(j:j)
        ENDDO

        IF (chem(istop+1:istop+1) .EQ. '=') THEN
          ik = ik + 1
          skelet(ik:ik) = '='
        ENDIF

* get the number of bonds at each node
* => construct diagonal of bond matrix
* (aromaticity is regarded as valence 3 for bookkeeping)
* --------------------------------------------
        nb = 0
        p = 0

        DO 120 j=2,l
* number of functional group inside parenthesis
          IF (group(ig)(j:j).EQ.'(' .AND. p .EQ. 0) nb=nb+1
          IF (group(ig)(j:j).EQ.'(') p = p + 1
          IF (group(ig)(j:j).EQ.')') p = p - 1

* H-,O-atom and radical attached to a carbon:
          IF (p .NE. 0) GOTO 120 ! cycle
          IF (group(ig)(j:j).EQ.'H')  nb = nb + 1
          IF (group(ig)(j-1:j-1).EQ.'H'.AND.
     &        group(ig)(j:j).EQ.'3')  nb = nb + 2
          IF (group(ig)(j-1:j-1).EQ.'H'.AND.
     &        group(ig)(j:j).EQ.'2')  nb = nb + 1
          IF (group(ig)(j:j).EQ.'O')  nb = nb + 2
          IF (group(ig)(j:j).EQ.'.')  nb = nb + 1 
120     CONTINUE

* assign to the diagonal term of the bond matrix
        bond(ig,ig) = nb

10    CONTINUE              

* keep number of groups:
      ng = ig

! troubleshooting
!      WRITE(6,'(15X,A11,A50)') 'skeleton = ',SKELET(1:50)
!      print*,' diagonal' 
!      DO j=1,ng
!        write(6,'(a10,30(i2))')group(j),(bond(i,j),i=1,ng)
!      ENDDO
!      print*,' '

      IF(MOD(nring,2).NE.0)THEN
        print*,'ERROR in grbond: RING NOT CLOSED'
        STOP
      ENDIF
      nring=nring/2

* ----------------------------------------------------------
* FIND TERMS OF BOND MATRIX (C-C, c-c, c-C, AND C-O-C BONDS)
* ----------------------------------------------------------

      ig = 0
      DO 50 i=1,ik-1
         nok=0
         IF (skelet(i:i).EQ.'C'.OR.
     &       skelet(i:i).EQ.'c'.OR.
     &       skelet(i:i).EQ.'O') nok=1
         IF (nok.EQ.0) GO TO 50  ! cycle

         ig = ig + 1
         k  = ig
         p  = 0
         tp = 0
         nb = 1
         DO 55 j=i+1,ik
           IF (skelet(j:j).EQ.'(') THEN
              tp = 1
              p  = p + 1
           !print*,'(    ',ig,k,nb,bond(ig,k)
           ELSE IF (skelet(j:j).EQ.')') THEN
              tp = 0
              p  = p - 1
           !print*,')    ',ig,k,nb,bond(ig,k)
              IF (p.LT.0) GOTO 50  ! exit
           ELSE IF (skelet(j:j).EQ.'C'.OR.skelet(j:j).EQ.'c') THEN
              k = k + 1
              IF (p.EQ.0 .AND. tp.EQ.0) THEN
                 bond(ig,k) = nb
                 bond(k,ig) = nb
           !print*,'C n()',ig,k,nb,bond(ig,k)
                 GOTO 50  ! exit
              ELSE IF (p.EQ.1 .AND. tp.EQ.1) THEN
                 bond(ig,k) = nb
                 bond(k,ig) = nb
           !print*,'C y()',ig,k,nb,bond(ig,k)
                 nb = 1
                 tp = 0 
                 GOTO 55  ! cycle 
! new section in case of double bonds up branches (2 lines of code new)
              ELSE IF (p.NE.tp) THEN
                nb = 1
           !print*,'p<>tp',ig,k,nb,bond(ig,k)
              ENDIF                       
           ELSE IF (skelet(j:j).EQ.'O') THEN
              k = k + 1
              IF (p.EQ.0 .AND. tp.EQ.0) THEN
                 bond(ig,k) = nb
                 bond(k,ig) = nb
                 GOTO 50  ! exit
              ELSE IF (p.EQ.1 .AND. tp.EQ.1) THEN
                 bond(ig,k) = nb
                 bond(k,ig) = nb
                 nb = 1
                 tp = 0 
                 GOTO 55  ! cycle 
              ENDIF                       
           ELSE IF (skelet(j-1:j).EQ.'(='.AND. p.EQ.1) THEN
             nb    = 2
             dbflg = 1
           ELSE IF (skelet(j:j).EQ.'=' .AND. p.EQ.0) THEN
             nb    = 2
             dbflg = 1
           ENDIF
55       CONTINUE
50    CONTINUE

      DO i=1,ng
        nb = 0
        DO j=1,mca
           nb = nb + bond(i,j)
        ENDDO

      ENDDO

!troubleshooting
      !print*,' basic skeleton' 
      !DO j=1,ng
      !  write(6,'(a10,30(i2))')group(j),(bond(i,j),i=1,ng)
      !ENDDO
      !print*,' '

! ----------------------------------------------------------
! FIND BOND MATRIX TERMS FOR RING-JOINING CARBONS
! ----------------------------------------------------------
! i  = group index of current node
! ii = group index of subsequent node 
! j  = skeleton index of current node
! jj = skeleton index of subsequent node
! k  = skeleton index of current ring-join character
! kk = skeleton index of subsequent ring-join character
! ----------------------------------------------------------
      IF(nring.GT.0)THEN
        i=0
        ii=0
        j=0
        jj=0
        k=0
        kk=0

        DO j=1,ik-1  ! loop skeleton
          IF (skelet(j:j).EQ.'C'.OR.
     &        skelet(j:j).EQ.'c'.OR.
     &        skelet(j:j).EQ.'O')THEN  ! node
            k=j
            i=i+1
            nb = 0
60          CONTINUE
            k=k+1
            DO n=1,nring  ! loop digits
              IF (skelet(k:k).EQ.digit(n)) THEN  ! ring start
                IF (ring(n,1).EQ.1) THEN ! ring already dealt with
                  GO TO 60  ! test for multi-ring center 
                ENDIF
                ii = i  
                nb = 1
                DO jj=k+1,ik  ! loop skeleton forwards
                   IF (skelet(jj:jj).EQ.'C'.OR.
     &                skelet(jj:jj).EQ.'c'.OR.
     &                skelet(jj:jj).EQ.'O')THEN  ! node
                    ii=ii+1 
                    kk = jj+1
                    IF(skelet(kk:kk).EQ.digit(n))THEN ! ring joined
                      !print*,"  ring ",digit(n)," joins at nodes", i,ii
                      ring(n,1)=1
                      bond(i,ii) = nb
                      bond(ii,i) = nb
                      IF(digit(n).EQ.'1'.AND.dbbrk.EQ.1)THEN
                        bond(i,ii) = bond(i,ii)+1
                        bond(ii,i) = bond(ii,i)+1
                      ENDIF
                      GO TO 60  ! test for multi-ring center on 1st node
                    ELSE ! test for multi-ring center on 2nd node
                      DO nn=1,nring
                        IF (skelet(kk:kk).EQ.digit(nn)) THEN  ! ring char
                          kk=jj+2
                          IF(skelet(kk:kk).EQ.digit(n))THEN ! ring joined
                      !print*,"  ring ",digit(n)," joins at nodes", i,ii
                            ring(n,1)=1
                            bond(i,ii) = nb
                            bond(ii,i) = nb
                            GO TO 60  ! test for multi-ring center on 1st node
                          ENDIF  ! ring joined
                        ENDIF  ! ring char
                      ENDDO  ! loop digits for ring char to ignore
                    ENDIF  ! ring joined
                  ENDIF  ! node
                ENDDO  ! loop skeleton forwards
              ENDIF  ! ring start
            ENDDO  ! loop digits
          ENDIF  ! node
        ENDDO  ! loop skeleton
      ENDIF  ! ring(s) present

! troubleshooting
      !print*,' after ring-joining' 
      !DO j=1,ng
      !  write(6,'(a10,30(i2))')group(j),(bond(i,j),i=1,ng)
      !ENDDO
      !print*,' '

* -----------
* FINAL CHECK
* -----------

* verify valence = 4 on each carbon, or '3' if aromatic carbon:
      DO i=1,ng
        nb = 0
        DO j=1,mca
           nb = nb + bond(i,j)
        ENDDO

        IF (group(i)(1:1).EQ.'C'.AND.nb.NE.4 .OR.
     &      group(i)(1:1).EQ.'c'.AND.nb.NE.3) THEN
          WRITE(6,'(a)') '--error-- '
          WRITE(6,'(a)') 'from SGMM routine : grbond'
          WRITE(6,'(a)') 'incorrect valence for species:'
          WRITE(6,'(a)') chem
          WRITE(6,'(15X,A11,A50)') 'skeleton = ',SKELET(1:50)
          WRITE(6,'(10X,A16)') 'bond - matrix = '
          DO k=1,ng
            WRITE(6,'(30(i2))') (bond(k,l),l=1,ng)
          ENDDO
          STOP
        ENDIF
        bond(i,i) = 0
      ENDDO
      !print*,'after zeroing main diagonal'
      !DO j=1,ng
      !  write(6,'(a10,#0(i2))')group(j),(bond(i,j),i=1,ng)
      !ENDDO
      !print*,' '

* check that no 'Cd' group exists without a Cd=Cd structure
      dbflg=0
      DO i=1,ng
        IF(INDEX(group(i),'Cd').NE.0)THEN
          DO j=1,mca
            IF(bond(i,j).EQ.2)dbflg=1
          ENDDO
          IF(dbflg.NE.1)THEN
            WRITE(6,'(a)') '--error-- '
            WRITE(6,'(a)') 'from SGMM routine : grbond'
            WRITE(6,'(a)') 'Cd but no >C=C< structure in species:'
            WRITE(6,'(a)') chem
            WRITE(6,'(15X,A11,A50)') 'skeleton = ',SKELET(1:50)
            WRITE(6,'(10X,A16)') 'bond - matrix = '
            DO k=1,ng
              WRITE(6,'(30(i2))') (bond(k,l),l=1,ng)
            ENDDO
            STOP
          ENDIF              
        ENDIF              
      ENDDO

* check that the -O- function is not a terminal position
      IF (loeter) THEN
        DO i=1,ng
          IF (group(i)(1:2).EQ.'-O') THEN
            nok=0
            DO j=1,ng
              nok=nok+bond(i,j)
            ENDDO
            IF (nok.NE.2) THEN
              WRITE(6,'(a)') '--error-- '
              WRITE(6,'(a)') 'from SGMM routine : grbond'
              WRITE(6,'(a)') 'wrong number of bond at an -O- site'
              WRITE(6,'(a)') 'valence should be 2 at that site.'
              WRITE(6,'(a)') 'chemical is :'
              WRITE(6,'(a)') chem
              WRITE(6,'(15X,A11,A50)') 'skeleton = ',SKELET(1:50)
              WRITE(6,'(10X,A16)') 'bond - matrix = '
              DO k=1,mca
                WRITE(6,'(30(1X,I2))') (bond(k,l),l=1,mca)
              ENDDO
              STOP
            ENDIF
          ENDIF
        ENDDO   

* C-O-C bonds turned to 3 (instead of 1) to be recognized later. 
        DO i=1,ng-1
          IF (group(i)(1:2).EQ.'-O') THEN
            DO j=1,ng
              IF (bond(i,j).EQ.1) THEN
                bond(i,j)=3
                bond(j,i)=3
              ENDIF
            ENDDO
          ENDIF
        ENDDO  
      ENDIF 

! troubleshooting
      !print*,'from grbond'
      !print*,'ng=',ng
      !DO j=1,ng
      !  write(6,'(a10,14(i2))')group(j),(bond(i,j),i=1,ng)
      !ENDDO
      !print*,' '

1000  RETURN
      END
