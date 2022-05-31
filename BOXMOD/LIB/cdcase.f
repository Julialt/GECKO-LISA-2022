************************************************************************
* MASTER MECHANISM - ROUTINE NAME : cdcase                             *
*                                                                      *
*                                                                      *
* PURPOSE : the parameterization used to assign a rate constant for    *
*           OH addition to C=C bonds depends on whether the C=C bond   *
*           is conjugated or not with a C=O bond (i.e. structure of    *
*           type -C=C-C=O). The subroutine returns the "case" to which *
*           the species belongs. Six cases are considered :            *
*                                                                      *
* CASE 1 : regular Cd molecule : only >C=C< and >C=C-C=C< bonds in     *
*          the molecule (i.e, without conjugated C=C-C=O)              *
* CASE 2 : is for structures containing the >C=C-C=O structure         *
*          but no C=C-C=C (or C=C=O)                                   *
* CASE 3 : is for the -CO-C=C-C=C-C=O structure only (i.e. containing  *
*          carbonyl at both sides of the conjugated C=C-C=C)           *
* CASE 4 : Two double bonds non conjugated (i.e. C=C-C-C=C) but one    *
*          containing at least one C=C-C=O                             *
* CASE 5 : is for the -CO-C=C-C=C< structure (i.e. containing          *
*          carbonyl at only one side of the conjugated C=C-C=C)        *
* CASE 6 : -C=C=O e.g. ketene (only case for this group)               *
* CASE 7 : -C=C-O- => vinyl ether chemistry                            *
* INPUT:                                                               *
* - chem         : chemical formula                                    *
* - group(i)     : groups at position (carbon) i                       *
* - bond(i,j)    : carbon-carbon bond matrix of chem                   *
*                                                                      *
* OUTPUT:                                                              *
* - cdcase       : case to which chem belongs                          *
* - ncd          : number of "Cd" carbons in chem                      *
* - conjug       : =1 if conjugated Cd (C=C-C=C), otherwise =0         *
* - cdtable(i)   : carbon number bearing a "Cd"                        *
* - tcdcase(i)   : carbon number at terminal position in C=C-C=C       *
* - cdsub(i)     : number of -C- substitutents (including -CO-) bonded *
*                  to the Cd corresponding to cdtable(i)               *
* - cdcarbo(i)   : number of -CO- substitutents bonded to the Cd       *
*                  corresponding to cdtable(i)                         *
*                                                                      *
************************************************************************
      SUBROUTINE cdcase(chem,bond,group,
     &                  ncd,ncdcase,conjug,
     &                  cdtable,tcdtable,cdsub,cdcarbo,cdeth)
      IMPLICIT NONE
      INCLUDE 'general.h'
      INCLUDE 'organic.h'

* input:
      CHARACTER(lfo) chem
      INTEGER         bond(mca,mca)
      CHARACTER(lgr) group(mca)

* output :
      INTEGER         ncd
      INTEGER         ncdcase, cdtable(4),tcdtable(4), cdsub(4)
      INTEGER         cdcarbo(4,2)
      INTEGER         cdeth(4,2)

* internal
      INTEGER         conjug, l, i , j, nb,o

      !WRITE(6,*)'*cdcase*'

* initialize
      ncdcase=0
      ncd=0
      l=0

      DO i=1,4
        cdtable(i)=0
        cdsub(i)=0
        tcdtable(i)=0
        DO j=1,2
          cdcarbo(i,j)=0
          cdeth(i,j)=0
        ENDDO
      ENDDO

* count number of Cd in the molecule and put Cd number in cdtable
      DO i=1,mca
        IF (INDEX(group(i),'Cd').NE.0)  THEN
          ncd=ncd+1
          cdtable(ncd)=i      
        ENDIF    
      ENDDO

* check that the molecule is allowed. 
* -----------------------------------
      IF (ncd.eq.0) THEN
        RETURN
      ENDIF

* more than 2 double bounds are not treated here 
      IF (ncd.gt.4) THEN
        WRITE(6,'(a)') '--error--, in cdcase subroutine'
        WRITE(6,'(a)') 'number of double bonded carbon'
        WRITE(6,'(a)') 'is greater 4'
        WRITE(6,'(a)') 'for the species :'
        WRITE(6,'(a)') chem
        WRITE(6,'(a)') 'ncd=',ncd
        WRITE(99,*)  'cdcase',chem 
        RETURN !STOP
      ENDIF
       
* check >C=C=C< structure
      IF (ncd.eq.3) THEN
        WRITE(6,'(a)') '--error--, in cdcase subroutine'
        WRITE(6,'(a)') '>C=C=C< structure not allowed'
        WRITE(6,'(a)') 'for the species :'
        WRITE(6,'(a)') chem
        WRITE(6,'(a)') 'ncd=',ncd
        WRITE(99,*)  'cdcase',chem 
        RETURN !STOP
      ENDIF

* check >C=CR-OH or >C=CR-ONO2 (not available yet)
      DO i=1,4
        IF(cdtable(i).GT.0)THEN
          IF (INDEX(group(cdtable(i)),'(OH)').NE.0)THEN
            WRITE(6,'(a)') '--error--, in cdcase subroutine'
            WRITE(6,'(a)') '>C=CR-OH structure not allowed'
            WRITE(6,'(a)') 'for the species :'
            WRITE(6,'(a)') chem
            WRITE(99,*)  'cdcase',chem 
            RETURN !STOP
          ELSE IF (INDEX(group(cdtable(i)),'(ONO2)').NE.0)THEN
            WRITE(6,'(a)') '--error--, in cdcase subroutine'
            WRITE(6,'(a)') '>C=CR-ONO2 structure not allowed'
            WRITE(6,'(a)') 'for the species :'
            WRITE(6,'(a)') chem
            WRITE(99,*)  'cdcase',chem 
            RETURN !STOP
          ENDIF
        ENDIF
      ENDDO            

* check if the molecule is a >C=C-C=C< structure. At the end, tcdtable
* should only contain the "terminal Cd's". This is necessary for
* further checks
      conjug=0
      IF (ncd.eq.4) THEN
        DO i=1,4
          tcdtable(i)=cdtable(i)
        ENDDO
        DO i=1,3
          DO j=i+1,4
            IF (bond(cdtable(i),cdtable(j)).eq.1) THEN
              conjug=1
              tcdtable(i)=0
              tcdtable(j)=0
              GOTO 457
            ENDIF
          ENDDO
        ENDDO               
      ENDIF
457   CONTINUE

* count the number of carbons or -O- bonded to each Cd, except the carbon 
* involved in the C=C bond (for which bond(i,j)=2)
      DO i = 1, 4
        IF (cdtable(i).ne.0) THEN

* if C=C=O present, default to case 6 (jlt)
          IF (group(cdtable(i))(1:3).EQ.'CdO') GOTO 600
* end default to case 6

          nb=0
          DO j=1,mca
            IF ((bond(cdtable(i),j).eq.1).OR.
     &          (bond(cdtable(i),j).eq.3)) nb=nb+1
          ENDDO
          cdsub(i)=nb

! Ludo: check the number of ether functions on C=C
          o=0
          DO j=1,mca
             IF (bond(cdtable(i),j).EQ.3) THEN
               IF (o.LT.1) THEN
                  o=o+1
                  cdeth(i,o)=j
               ELSE
                  WRITE(6,'(a)') '--error--, in cdcase subroutine'
                  WRITE(6,'(a)') 'more than one ether function on a Cd'
                  WRITE(6,'(a)') 'for the species :'
                  WRITE(6,'(a)') chem
                  WRITE(99,*)  'cdcase',chem 
                  STOP
               ENDIF                   
             ENDIF
          ENDDO
        ENDIF
      ENDDO

! Ludo: if C=C-O- structure, case 7
      DO i=1,3,2
         IF (((cdeth(i,1)).NE.0).AND.((cdeth(i+1,1)).NE.0)) THEN
            WRITE(6,'(a)') '--error--, in cdcase subroutine'
            WRITE(6,'(a)') 'the two C of a double bond each have'
            WRITE(6,'(a)') 'an ether function'
            WRITE(6,'(a)') 'for the species :'
            WRITE(6,'(a)') chem
            WRITE(99,*)  'cdcase',chem 
            STOP
         ELSE IF ((((cdeth(i,1)).NE.0).AND.((cdeth(i+1,1)).EQ.0)).OR.
     &           (((cdeth(i,1)).EQ.0).AND.((cdeth(i+1,1)).NE.0))) THEN
            ncdcase=7
            WRITE(32,*) 'case 7',chem(1:100)
            GOTO 700
         ENDIF
      ENDDO

* count the number of carbonyls conjugated to each Cd
      DO i=1,4
        l=0
        IF (cdtable(i).ne.0) THEN
          DO j=1,mca
            IF (bond(cdtable(i),j).eq.1) THEN
              IF (group(j).eq.'CHO') THEN
                l=l+1
                cdcarbo(i,l)=j
              ELSE IF (group(j)(1:2).eq.'CO') THEN
                l=l+1
                cdcarbo(i,l)=j
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        IF (l.gt.2) THEN
          WRITE(6,'(a)') '--error--, in cdcase subroutine'
          WRITE(6,'(a)') 'more than 2 carbonyl bounded to a Cd'
          WRITE(6,'(a)') 'for the species :'
          WRITE(6,'(a)') chem
          WRITE(99,*)  'cdcase',chem 
          RETURN !STOP
        ENDIF
      ENDDO

* check if the molecule has a -CO-C=C-C=C-CO- structure (case 3)
      IF (conjug.eq.1) THEN
        l=0
        DO 345 i=1,4
          IF (tcdtable(i).ne.0) THEN
            DO j=1,mca
              IF (bond(tcdtable(i),j).eq.1) THEN
                IF (group(j)(1:3).EQ.'CHO') THEN
                  l=l+1
c                  Ci=tcdtable(i)
                  GOTO 345
                ENDIF
                IF (group(j)(1:2).EQ.'CO') THEN
                  l=l+1
c                 Ci=tcdtable(i)
                  GOTO 345
                ENDIF
              ENDIF
            ENDDO
          ENDIF
345     CONTINUE
        IF (l.ge.2) THEN
          ncdcase=3
          GOTO 700
        ENDIF
      ENDIF

* -C=C-C=C-CO- is substituted with at least one carbonyl (case 5)      
      IF (conjug.eq.1) THEN
        DO i=1,4
          IF (cdtable(i).ne.0) THEN
            DO j=1,mca
              IF (bond(cdtable(i),j).eq.1) THEN
                IF (group(j)(1:3).EQ.'CHO') THEN
                  IF ((cdcarbo(2,1).NE.0).AND.
     &               (cdcarbo(3,1).NE.0)) THEN
                    ncdcase=3
                    GOTO 700
                  ELSE
                    ncdcase=5
                    GOTO 700
                  ENDIF
                ENDIF
                IF (group(j)(1:2).EQ.'CO') THEN
                  IF ((cdcarbo(2,1).NE.0).AND.
     &               (cdcarbo(3,1).NE.0)) THEN
                    ncdcase=3
                    GOTO 700
                  ELSE
                    ncdcase=5
                    GOTO 700
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

* check for C=C-C=O structure (case 2 or 4)
      DO i=1,4
        IF (cdtable(i).ne.0) THEN
          DO j=1,mca
            IF (bond(cdtable(i),j).eq.1) THEN
              IF (group(j)(1:3).EQ.'CHO') THEN
                IF (ncd.eq.2) THEN
                  ncdcase=2
                  GOTO 700
                ELSE IF (ncd.eq.4) THEN
                  ncdcase=4
                  GOTO 700
                ENDIF    
              ENDIF
              IF (group(j)(1:2).EQ.'CO') THEN
                IF (ncd.eq.2) THEN
                  ncdcase=2
                  GOTO 700
                ELSE IF (ncd.eq.4) THEN
                  ncdcase=4
                  GOTO 700
                ENDIF    
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

600   CONTINUE

* check for C=C=O structure (case 6)
* count the number of carbons bonded to each Cd, except the carbon 
* involved in the C=C bond (for which bond(i,j)=2)
      DO i = 1, 4
        IF (cdtable(i).ne.0) THEN
          IF (INDEX(group(cdtable(i)),'CdO').NE.0) THEN
            ncdcase=6
            GOTO 700
          ENDIF
        ENDIF
      ENDDO

* otherwise case 1 (if this point is reached)
      ncdcase=1

700   CONTINUE
      
      END
