!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input :
! --------
!  chem : chemical formula
!  node : number of nodes is chem
!  group : table of groups
!  bond : bond matrix
!
! outputs :
! ---------
!  ngrp : number of functional group in chem
!  nodetype : table of character for type node
!             'y' is for carbonyl
!             'r' is for aromatic
!             'o' is for -O- node
!             'd' is for Cd
!             'n' is for the others (i.e. normal)
! alifun(i) : number of group of type "i" on a aliphatic carbon
! cdfun(i) : number of group of type "i" on a cd carbon
! arofun(i) : number of group of type "i" on a aromatic carbon
!   the index in those tables are :
!     index  1 -  5 : 'ROH  ','RNO2 ','RONO2','ROOH ','RF   '  
!     index  6 - 10 : 'RCl  ','RBr  ','RI   ','RCHO ','RCOR '
!     index 11 - 15 : 'RCOOH','COOOH','PAN  ','ROR  ','RCOOR'
!     index 16 - 20 : 'HCOOR','RCOF ','RCOCl','RCOBr','RCOI '
!     index 21      : 'CO(ONO2)'
! mapfun(a,b,c) : provide the number of function of type 'c'
!                 at position (node) 'a'. index 'b' if for node
!                 type with 1=aliphatic, 2=cd and 3=aromatic
!                 for example, the molecule CH2(OH)CdH=CdHCHO 
!                 should have non zero values at the positions :
!                 mapfun(1,1,1)=1
!                 mapfun(4,2,9)=1
! funflg(a)     : get the number of functional group at 
!                 node a. For the example above, non-zero
!                 values are found at position 1 and 4, where
!                 it is set to 1.  
! tabester : provide the position of ester "couple" (i.e
!            the O and CO nodes. For example, the molecule
!            CH3CO-O-CH2-O-COCH3 has the following values
!            tabester(1,1)=3,tabester(1,2)=2 
!            tabester(2,1)=5,tabester(2,2)=6
! ierr     : if not set to 0, then an error occured in the
!             subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE chemmap(chem,node,group,bond,ngrp,nodetype,
     &            alifun,cdfun,arofun,mapfun,funflg,aroflg,
     &            arobig,tabester,ierr)
      IMPLICIT NONE
      INCLUDE 'general.h'


      CHARACTER(lfo) chem

      INTEGER         ngrp
      CHARACTER(lgr) group(mca)
      INTEGER         bond(mca,mca),node
      CHARACTER(1)     nodetype(mca)
!      REAL            alifun(20),cdfun(20),arofun(20)
!      REAL            mapfun(mca,3,20)
      REAL            alifun(21),cdfun(21),arofun(21)
      REAL            mapfun(mca,3,21)
      INTEGER         funflg(mca),aroflg(mca),arobig(mca)
      INTEGER         tabester(4,2)  ! 1= -O- side, 2= CO side
      INTEGER         ierr

      INTEGER         i,j,k
      INTEGER         nnod,tnod(mca)
      INTEGER         nf,ialpha,ialpha2,iy,rflg,dflg,yflg
      INTEGER         ichecko, ichecky
      INTEGER         ytab(2)
      INTEGER         nester
c      CHARACTER(5)     funnam(20)

c      DATA  funnam/'ROH  ','RNO2 ','RONO2','ROOH ','RF   ',  
c     &             'RCl  ','RBr  ','RI   ','RCHO ','RCOR ',
c     &             'RCOOH','COOOH','PAN  ','ROR  ','RCOOR',
c     &             'HCOOR','RCOF ','RCOCl','RCOBr','RCOI '/

! ---------
! Initialize
! ---------
      ierr=0
      DO i=1,21
        alifun(i)=0
        cdfun(i)=0
        arofun(i)=0
        DO j=1,3
          DO k=1,mca
            mapfun(k,j,i)=0
          ENDDO
        ENDDO
      ENDDO
      DO i=1,mca
        nodetype(i)=' '
        funflg(i)=0
        aroflg(i)=0
        arobig(i)=0
      ENDDO
      nester=0
      DO i=1,2
        DO j=1,4
          tabester(j,i)=0
        ENDDO
      ENDDO
    

! get the type of nodes in groups
      DO i=1,node
        IF (group(i)(1:2).eq.'CO') THEN
          nodetype(i)='y'
        ELSE IF (group(i)(1:3).eq.'CHO') THEN
          nodetype(i)='y'
        ElSE IF (group(i)(1:1).eq.'c') THEN
          nodetype(i)='r'
        ELSE IF (group(i)(1:3).eq.'-O-') THEN
          nodetype(i)='o'
        ELSE IF (group(i)(1:2).eq.'Cd') THEN
          nodetype(i)='d'
        ELSE 
          nodetype(i)='n'
        ENDIF 
      ENDDO

! feed table of functions 
! --------------------------
!  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
!  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
! 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= R-COO-H 
! 16 = HCO-O-R; 17= -CO(F) ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
! 21 = R-CO(ONO2)


* ------- Alkohols (index 1) and Acids (index 11) ----------
      IF (INDEX(chem,'(OH)').eq.0) GOTO 100
      DO i=1, node
        IF (INDEX(group(i),'(OH)').ne.0) THEN
          nf=0            
          DO j=1,lgr-3
            IF (group(i)(j:j+3).eq.'(OH)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! alcohol aliphatic
            alifun(1)=alifun(1)+nf
            mapfun(i,1,1)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! phenols 
            arofun(1)=arofun(1)+nf
            mapfun(i,3,1)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-OH (not expected, since enol) 
            cdfun(1)=cdfun(1)+nf
            mapfun(i,2,1)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! Carboxylic acid
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a carboxylic group'
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(OH) aromatic 
              arofun(11)=arofun(11)+1
              mapfun(i,3,11)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(OH) on Cd
              cdfun(11)=cdfun(11)+1
c              alifun(11)=alifun(11)+1
              mapfun(i,2,11)=1
c              mapfun(i,1,11)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(OH) aliphatic
              alifun(11)=alifun(11)+1
              mapfun(i,1,11)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
100   CONTINUE

* ----------- Nitro (index 2) -----------------
      IF (INDEX(chem,'(NO2)').eq.0) GOTO 110
      DO i = 1, node 
        IF (INDEX(group(i),'(NO2)').ne.0) THEN
          nf=0
          DO j=1,lgr-4
            IF (group(i)(j:j+4).eq.'(NO2)') nf=nf+1 
          ENDDO

          IF (nodetype(i).eq.'n') THEN     ! NO2 aliphatic
            alifun(2)=alifun(2)+nf
            mapfun(i,1,2)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN     ! NO2 on Cd 
            cdfun(2)=cdfun(2)+nf
            mapfun(i,2,2)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN     ! NO2 on aromatic 
            arofun(2)=arofun(2)+nf
            mapfun(i,3,2)=nf
            ngrp=ngrp+nf
          ELSE
            WRITE(6,*) '-- error --, in '
            WRITE(6,*) 'a nitro group (NO2) is borne by an '
            WRITE(6,*) 'unexpected group in chem :'
            WRITE(6,*) chem(1:70)
            STOP
          ENDIF
        ENDIF
      ENDDO
110   CONTINUE

* ----------- Nitrate (index 3) -----------------
      IF (INDEX(chem,'(ONO2)').eq.0) GOTO 120
      DO i = 1, node 
        IF (INDEX(group(i),'(ONO2)').ne.0) THEN
          IF (INDEX(group(i),'CO(ONO2)').NE.0) GOTO 120
          nf=0
          DO j=1,lgr-5
            IF (group(i)(j:j+5).eq.'(ONO2)') nf=nf+1 
          ENDDO

          IF (nodetype(i).eq.'n') THEN          ! ONO2 aliphatic
            alifun(3)=alifun(3)+nf
            mapfun(i,1,3)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN     ! ONO2 on Cd 
            cdfun(3)=cdfun(3)+nf
            mapfun(i,2,3)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN     ! ONO2 on aromatic 
            arofun(3)=arofun(3)+nf
            mapfun(i,3,3)=nf
            ngrp=ngrp+nf
          ELSE
            WRITE(6,*) '-- error --, in '
            WRITE(6,*) 'a nitro group (ONO2) is borne by an '
            WRITE(6,*) 'unexpected group in chem :'
            WRITE(6,*) chem(1:70)
            STOP
          ENDIF
        ENDIF
      ENDDO
120   CONTINUE

* ------- hydroperoxydes (index 4) and peracids (index 12) ----------
      IF (INDEX(chem,'(OOH)').eq.0) GOTO 130
      DO i=1, node
        IF (INDEX(group(i),'(OOH)').ne.0) THEN
          nf=0            
          DO j=1,lgr-4
            IF (group(i)(j:j+4).eq.'(OOH)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! hydroperoxyde aliphatic
            alifun(4)=alifun(4)+nf
            mapfun(i,1,4)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! aromaric -OOH 
            arofun(4)=arofun(4)+nf
            mapfun(i,3,4)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-OOH (not expected) 
            cdfun(4)=cdfun(4)+nf
            mapfun(i,2,4)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! peracid acid
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a peracid group'
              WRITE(6,*) chem(1:70)
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(OOH) aromatic 
              arofun(12)=arofun(12)+1
              mapfun(i,3,12)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(OOH) on Cd
              cdfun(12)=cdfun(12)+1
c              alifun(12)=alifun(12)+1
              mapfun(i,2,12)=1
c              mapfun(i,1,12)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(OOH) aliphatic
              alifun(12)=alifun(12)+1
              mapfun(i,1,12)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
130   CONTINUE

* ------- fluroro (index 5) and fluoro acyl (index 17) ----------
      IF (INDEX(chem,'(F)').eq.0) GOTO 140
      DO i=1, node
        IF (INDEX(group(i),'(F)').ne.0) THEN
          nf=0            
          DO j=1,lgr-2
            IF (group(i)(j:j+2).eq.'(F)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! F aliphatic
            alifun(5)=alifun(5)+nf
            mapfun(i,1,5)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! aromaric F 
            arofun(5)=arofun(5)+nf
            mapfun(i,3,5)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-F 
            cdfun(5)=cdfun(5)+nf
            mapfun(i,2,5)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! fluoro acyl
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a fluoro acyl group'
              WRITE(6,*) chem(1:70)
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(F) aromatic 
              arofun(17)=arofun(17)+1
              mapfun(i,3,17)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(F) on Cd
              cdfun(17)=cdfun(17)+1
              mapfun(i,2,17)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(F) aliphatic
              alifun(17)=alifun(17)+1
              mapfun(i,1,17)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
140   CONTINUE

* ------- chloro (index 6) and chloro acyl (index 18) ----------
      IF (INDEX(chem,'(Cl)').eq.0) GOTO 150
      DO i=1, node
        IF (INDEX(group(i),'(Cl)').ne.0) THEN
          nf=0            
          DO j=1,lgr-3
            IF (group(i)(j:j+3).eq.'(Cl)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! Cl aliphatic
            alifun(6)=alifun(6)+nf
            mapfun(i,1,6)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! aromaric Cl
            arofun(6)=arofun(6)+nf
            mapfun(i,3,6)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-Cl  
            cdfun(6)=cdfun(6)+nf
            mapfun(i,2,6)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! Chloro acyl
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a chloro acyl group'
              WRITE(6,*) chem(1:70)
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(Cl) aromatic 
              arofun(18)=arofun(18)+1
              mapfun(i,3,18)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(Cl) on Cd
              cdfun(18)=cdfun(18)+1
              mapfun(i,2,18)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(Cl) aliphatic
              alifun(18)=alifun(18)+1
              mapfun(i,1,18)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
150   CONTINUE

* ------- Bromo (index 7) and bromo acyl (index 19) ----------
      IF (INDEX(chem,'(Br)').eq.0) GOTO 160
      DO i=1, node
        IF (INDEX(group(i),'(Br)').ne.0) THEN
          nf=0            
          DO j=1,lgr-3
            IF (group(i)(j:j+3).eq.'(Br)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! Br aliphatic
            alifun(7)=alifun(7)+nf
            mapfun(i,1,7)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! aromatic Br
            arofun(7)=arofun(7)+nf
            mapfun(i,3,7)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-Br  
            cdfun(7)=cdfun(7)+nf
            mapfun(i,2,7)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! Bromo acyl
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a bromo acyl group'
              WRITE(6,*) chem(1:70)
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(Br) aromatic 
              arofun(19)=arofun(19)+1
              mapfun(i,3,19)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(Br) on Cd
              cdfun(19)=cdfun(19)+1
              mapfun(i,2,19)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(Br) aliphatic
              alifun(19)=alifun(19)+1
              mapfun(i,1,19)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
160   CONTINUE

* ------- iodo (index 8) and iodo acyl (index 20) ----------
      IF (INDEX(chem,'(I)').eq.0) GOTO 170
      DO i=1, node
        IF (INDEX(group(i),'(I)').ne.0) THEN
          nf=0            
          DO j=1,lgr-2
            IF (group(i)(j:j+2).eq.'(I)') nf=nf+1 
          ENDDO
          IF (nodetype(i).eq.'n') THEN       ! I aliphatic
            alifun(8)=alifun(8)+nf
            mapfun(i,1,8)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'r') THEN  ! aromaric I 
            arofun(8)=arofun(8)+nf
            mapfun(i,3,8)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'d') THEN  ! =Cd-I 
            cdfun(8)=cdfun(8)+nf
            mapfun(i,2,8)=nf
            ngrp=ngrp+nf
          ELSE IF (nodetype(i).eq.'y') THEN  ! fluoro acyl
            CALL nodmap(bond,i,node,2,nnod,tnod)
            IF (nnod.gt.1) THEN 
              DO j=1,node
                write(*,*) 'group',j,'-',group(j)
              ENDDO
              write(*,*) 'nnod=',nnod
              write(*,*) 'tnod=',(tnod(j),j=1,10)
              WRITE(6,*) '-- error --, a unique C is expected in'
              WRITE(6,*) 'alpha position of a Iodo acyl group'
              WRITE(6,*) chem(1:70)
              STOP
            ENDIF
            ialpha=tnod(1)
            IF (nodetype(ialpha).EQ.'r')THEN         ! CO(I) aromatic 
              arofun(20)=arofun(20)+1
              mapfun(i,3,20)=1
              ngrp=ngrp+1
            ELSE IF (nodetype(ialpha).EQ.'d')THEN    ! CO(I) on Cd
              cdfun(20)=cdfun(20)+1
              mapfun(i,2,20)=1
              ngrp=ngrp+1
            ELSE                                    ! CO(I) aliphatic
              alifun(20)=alifun(20)+1
              mapfun(i,1,20)=1
              ngrp=ngrp+1
            ENDIF 
          ENDIF
        ENDIF
      ENDDO 
170   CONTINUE

*---------- ester (index 15 and 16) -------------
      IF (INDEX(chem,'-O').eq.0) GOTO 221
      DO 220 i = 1, node
        IF (group(i)(1:3).EQ.'-O-') THEN
          CALL nodmap(bond,i,node,2,nnod,tnod)
          rflg=0
          dflg=0
          yflg=0
          DO j=1,nnod
            ialpha=tnod(j)
            IF (nodetype(ialpha).EQ.'y') THEN   ! RCO-O-R function
              ichecky=0
              DO k=1,4
                IF (tabester(k,2).eq.ialpha) ichecky=1   ! carbonyl already used
              ENDDO
              IF (ichecky.eq.0) THEN
                yflg=yflg+1  
                ytab(yflg)=ialpha
              ENDIF
            ENDIF
c            IF (nodetype(ialpha).EQ.'r') rflg=rflg+1       
c            IF (nodetype(ialpha).EQ.'d') dflg=dflg+1       
          ENDDO
          IF (yflg.eq.0) GOTO 220                      ! simple ether

! simple ester 
          IF (yflg.ne.0) THEN
            nester=nester+1
            IF (nester.gt.4) STOP 'stop in chemmap. > 4'
            iy=ytab(1)
            tabester(nester,1)=i
            tabester(nester,2)=iy
            DO j=1,2
              IF (tnod(j).ne.ytab(1)) THEN
                ialpha=tnod(j)   ! second side of the -O- (first is iy)
              ENDIF
            ENDDO
            IF (group(iy)(1:3).eq.'CHO') THEN    ! HCO-O-R
              IF (nodetype(ialpha).EQ.'r') THEN       ! aromatic CHO-O-R
                arofun(16)=arofun(16)+1
                mapfun(i,3,16)=1
                mapfun(iy,3,16)=1
                ngrp=ngrp+1
              ELSE IF (nodetype(ialpha).EQ.'d') THEN  ! =C-O-CHO
                cdfun(16)=cdfun(16)+1
                mapfun(i,2,16)=1
                mapfun(iy,2,16)=1
                ngrp=ngrp+1
              ELSE                                    ! R-O-CHO
                alifun(16)=alifun(16)+1
                mapfun(i,1,16)=1
                mapfun(iy,1,16)=1
                ngrp=ngrp+1
              ENDIF
            ELSE IF (group(iy)(1:2).eq.'CO') THEN    ! RCO-O-R
              CALL nodmap(bond,iy,node,2,nnod,tnod)
              DO j=1,2
                IF (tnod(j).ne.i) ialpha2=tnod(j)  ! ialpha2 is the node next to the
CO of the ester
              ENDDO
              rflg=0
              dflg=0

!              IF ( (nodetype(ialpha).EQ.'y').OR.
!     &             (nodetype(ialpha2).EQ.'y')     ) THEN   ! R-CO-O-CO-R function
!                  alifun(10)=alifun(10)+1
!                  mapfun(i,1,10)=1
!                  ngrp=ngrp+1                  
!cc                WRITE(29,*) chem(1:50)
!cc                ierr=1
!cc                RETURN
!              ENDIF   
! structure is ialpha2-CO-O-ialpha
              IF (nodetype(ialpha).EQ.'r') rflg=rflg+1       
              IF (nodetype(ialpha).EQ.'d') dflg=dflg+1       
              IF (nodetype(ialpha2).EQ.'r') rflg=rflg+1       
              IF (nodetype(ialpha2).EQ.'d') dflg=dflg+1       

              IF (rflg.ne.0) THEN                    ! aromatic ester
c                arofun(15)=arofun(15)+real(rflg)/2.
c                mapfun(i,3,15)=real(rflg)/2.
c                mapfun(iy,3,15)=real(rflg)/2.
                arofun(15)=arofun(15)+1.
                mapfun(i,3,15)=1.
                mapfun(iy,3,15)=1.
              ELSE IF (dflg.ne.0) THEN                    ! =C-CO-O-
c                cdfun(15)=cdfun(15)+real(dflg)/2.
c                mapfun(i,2,15)=real(dflg)/2.
c                mapfun(iy,2,15)=real(dflg)/2.
                cdfun(15)=cdfun(15)+1.
                mapfun(i,2,15)=1.
                mapfun(iy,2,15)=1.
c              nf=2-rflg-dflg
c              nf=1-rflg-dflgaroflg(i)
              ELSE                       ! R-COO-R
c                alifun(15)=alifun(15)+real(nf)/2.
c                mapfun(i,1,15)=real(dflg)/2.
c                mapfun(iy,1,15)=real(dflg)/2.
                alifun(15)=alifun(15)+1.
                mapfun(i,1,15)=1.
                mapfun(iy,1,15)=1.
              ENDIF
              ngrp=ngrp+1
            ENDIF
          ENDIF

* -CO-O-CO- 
cc          IF (yflg.eq.2) THEN
cc            WRITE(29,*) chem(1:50)
cc            ierr=1
cc            RETURN
cc          ENDIF

        ENDIF 
          
220   CONTINUE
221   CONTINUE 



*------------Aldehydes (index 9) --------- 
      IF (INDEX(chem,'CHO').eq.0) GOTO 181
      DO 180 i = 1, node
        IF (group(i)(1:3).EQ.'CHO') THEN
        CALL nodmap(bond,i,node,2,nnod,tnod)
        IF (nnod.ne.1) THEN
           WRITE(6,*) '-- error --, a unique C is expected in'
           WRITE(6,*) 'alpha position of a CHO  group'
           WRITE(6,*) chem(1:70)
           STOP
        ENDIF
        ialpha=tnod(1)
        IF (nodetype(ialpha).EQ.'o') THEN   ! HCO-O-R function
           ichecko=0 ! check if the ether is already involved in an ester
           ichecky=0 ! check if the carbonyl is already involved in an ester
           DO k=1,4
             IF (tabester(k,1).eq.ialpha) ichecko=1   ! ether already used
             IF (tabester(k,2).eq.i) ichecky=1   ! carbonyl already used
           ENDDO
           IF (ichecky.eq.1) GOTO 180  ! ether already involved
           IF (ichecko.eq.0) GOTO 180  ! carbonyl that must be an ester
        ENDIF  ! if that point is reached then must be counted as aldehyde

        IF (nodetype(ialpha).EQ.'r') THEN       ! aromatic aldehyde
           arofun(9)=arofun(9)+1
           mapfun(i,3,9)=1
           ngrp=ngrp+1
        ELSE IF (nodetype(ialpha).EQ.'d') THEN  ! =C-CHO
           cdfun(9)=cdfun(9)+1
c           alifun(9)=alifun(9)+1
           mapfun(i,2,9)=1
c           mapfun(i,1,9)=1
           ngrp=ngrp+1
        ELSE                                    ! R-CHO
           alifun(9)=alifun(9)+1
           mapfun(i,1,9)=1
           ngrp=ngrp+1
        ENDIF
       ENDIF
180   CONTINUE
181   CONTINUE 

*---------- ketones (index 10) -------------
      IF (INDEX(chem,'CO').eq.0) GOTO 191
      DO 190 i = 1, node
        IF (group(i)(1:3).EQ.'CO ') THEN
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod.ne.2) THEN
            WRITE(6,*) '-- error --, only 2 C is expected in'
            WRITE(6,*) 'alpha position of a -CO-  group'
            WRITE(6,*) chem(1:70)
            STOP
          ENDIF
          rflg=0
          dflg=0
          DO j=1,nnod
            ialpha=tnod(j)
            IF (nodetype(ialpha).EQ.'o') THEN   ! RCO-O-R function
              ichecko=0 ! check if the ether is already involved in an ester
              ichecky=0 ! check if the carbonyl is already involved in an ester
              DO k=1,4
                IF (tabester(k,1).eq.ialpha) ichecko=1   ! ether already used
                IF (tabester(k,2).eq.i) ichecky=1   ! carbonyl already used
              ENDDO
              IF (ichecky.eq.1) GOTO 190  ! carbonyl already involved
              IF (ichecko.eq.0) GOTO 190  ! ether that must be an ester
            ENDIF  ! if that point is reached then must be counted as ketone
            IF (nodetype(ialpha).EQ.'r') rflg=rflg+1       
            IF (nodetype(ialpha).EQ.'d') dflg=dflg+1       
          ENDDO
          IF (rflg.ne.0) THEN                    ! aromatic ketone
c            arofun(10)=arofun(10)+real(rflg)/2.
c            mapfun(i,3,10)=real(rflg)/2.
            arofun(10)=arofun(10)+1.
            mapfun(i,3,10)=1.
          ELSE IF (dflg.ne.0) THEN                    ! =C-CO-R
c            cdfun(10)=cdfun(10)+real(dflg)/2.
            cdfun(10)=cdfun(10)+1.
c            alifun(10)=alifun(10)+real(dflg)/2.
c            mapfun(i,2,10)=real(dflg)/2.
            mapfun(i,2,10)=1.
c            mapfun(i,1,10)=real(dflg)/2.
c          nf=2-rflg-dflg
c          nf=1-rflg-dflg
          ELSE                       ! R-CO-R
c            alifun(10)=alifun(10)+real(nf)/2.
c            mapfun(i,1,10)=real(dflg)/2.
            alifun(10)=alifun(10)+1.
            mapfun(i,1,10)=1.
          ENDIF
          ngrp=ngrp+1
        ENDIF
190   CONTINUE
191   CONTINUE 


* feed table of functions 
* --------------------------
*  1= -OH    ;  2= -NO2     ;  3= -ONO2   ; 4= -OOH    ;  5= -F     
*  6= -Cl    ;  7=  -Br     ;  8= -I      ; 9= -CHO    ; 10= -CO-   
* 11= -COOH  ; 12= -CO(OOH) ; 13= -PAN    ; 14= -O-    ; 15= -COO- 
* 16=HCO-O-  ; 17= -CO(F)   ; 18= -CO(Cl) ;  19= -CO(Br) ; 20= -CO(I)
* 21= -CH3
*------------- PAN (index 13) -----------------
      IF (INDEX(chem,'CO(OONO2').eq.0) GOTO 200
      DO i = 1, node
        IF (group(i)(1:9).eq.'CO(OONO2)') THEN
         CALL nodmap(bond,i,node,2,nnod,tnod)
         IF (nnod.ne.1) THEN
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a CO(OONO2)  group'
            WRITE(6,*) chem(1:70)
            STOP
         ENDIF
         ialpha=tnod(1)
         IF (nodetype(ialpha).EQ.'o') THEN   ! R-O-CO(OONO2) function
            WRITE(6,*) '-- error --,  in'
            WRITE(6,*) '-O-CO(OONO2) group is unexpected'
            WRITE(6,*) chem(1:70)
            STOP
         ENDIF
         IF (nodetype(ialpha).EQ.'r') THEN       ! aromatic PAN
            arofun(13)=arofun(13)+1
            mapfun(i,3,13)=1
            ngrp=ngrp+1
         ELSE IF (nodetype(ialpha).EQ.'d') THEN  ! =C-CO(OONO2)
            cdfun(13)=cdfun(13)+1
c            alifun(13)=alifun(13)+1
            mapfun(i,2,13)=1
c            mapfun(i,1,13)=1
            ngrp=ngrp+1
         ELSE                                    ! R-CO(OONO2)
            alifun(13)=alifun(13)+1
            mapfun(i,1,13)=1
            ngrp=ngrp+1
         ENDIF
        ENDIF
      ENDDO
200   CONTINUE 

*---------- ether (index 14) -------------
      IF (INDEX(chem,'-O').eq.0) GOTO 211
      DO 210 i = 1, node
        IF (group(i)(1:3).EQ.'-O-') THEN
          CALL nodmap(bond,i,node,2,nnod,tnod)
          IF (nnod.ne.2) THEN
            WRITE(6,*) '-- error --, only 2 C is expected in'
            WRITE(6,*) 'alpha position of a -O-  group'
            WRITE(6,*) chem(1:70)
            STOP
          ENDIF
          rflg=0
          dflg=0
          DO j=1,nnod
            ialpha=tnod(j)
            IF (nodetype(ialpha).EQ.'y') THEN   ! RCO-O-R function
              ichecko=0 ! check if the ether is already involved in an ester
              ichecky=0 ! check if the carbonyl is already involved in an ester
              DO k=1,4
                IF (tabester(k,1).eq.i) ichecko=1   ! ether already used
                IF (tabester(k,2).eq.ialpha) ichecky=1   ! carbonyl already used
              ENDDO
              IF (ichecko.eq.1) GOTO 210  ! ether already involved in an ester
              IF (ichecky.eq.0) GOTO 210  ! carbonyl that must be an ester
            ENDIF  ! if that point is reached then must be counted as ether

            IF (nodetype(ialpha).EQ.'r') rflg=rflg+1       
            IF (nodetype(ialpha).EQ.'d') dflg=dflg+1       
          ENDDO
          IF (rflg.ne.0) THEN                    ! aromatic ether
c            arofun(14)=arofun(14)+real(rflg)/2.
c            mapfun(i,3,14)=real(rflg)/2.
            arofun(14)=arofun(14)+1.
            mapfun(i,3,14)=1.
            IF (rflg.gt.1) THEN
              WRITE(29,*) chem(1:50)
              ierr=1
              RETURN
            ENDIF
          ELSE IF (dflg.ne.0) THEN                    ! =C-O-R
c            cdfun(14)=cdfun(14)+real(dflg)/2.
c            mapfun(i,2,14)=real(dflg)/2.
            cdfun(14)=cdfun(14)+1.
            mapfun(i,2,14)=1.
c            IF (dflg.gt.1) THEN
c              WRITE(29,*) chem(1:50)
c              ierr=1
c              RETURN
c            ENDIF
c          nf=2-rflg-dflg
c          nf=1-rflg-dflg
          ELSE                      ! R-O-R
c            alifun(14)=alifun(14)+real(nf)/2.
c            mapfun(i,1,14)=real(dflg)/2.
            alifun(14)=alifun(14)+1.
            mapfun(i,1,14)=1.
          ENDIF
          ngrp=ngrp+1
        ENDIF
210   CONTINUE
211   CONTINUE 

*------------- CO(ONO2) (index 21) --------------
      IF (INDEX(chem,'CO(ONO2').eq.0) GOTO 230
      DO i = 1, node
        IF (group(i)(1:8).eq.'CO(ONO2)') THEN
         CALL nodmap(bond,i,node,2,nnod,tnod)
         IF (nnod.ne.1) THEN
            WRITE(6,*) '-- error --, a unique C is expected in'
            WRITE(6,*) 'alpha position of a CO(ONO2)  group'
            WRITE(6,*) chem(1:70)
            STOP
         ENDIF
         ialpha=tnod(1)
         IF (nodetype(ialpha).EQ.'o') THEN   ! R-O-CO(ONO2) function
            WRITE(6,*) '-- error --,  in'
            WRITE(6,*) '-O-CO(ONO2) group is unexpected'
            WRITE(6,*) chem(1:70)
            STOP
         ENDIF
         IF (nodetype(ialpha).EQ.'r') THEN       ! aromatic
            arofun(21)=arofun(21)+1
            mapfun(i,3,21)=1
            ngrp=ngrp+1
         ELSE IF (nodetype(ialpha).EQ.'d') THEN  ! =C-CO(ONO2)
            cdfun(21)=cdfun(21)+1
c            alifun(13)=alifun(13)+1
            mapfun(i,2,21)=1
c            mapfun(i,1,13)=1
            ngrp=ngrp+1
         ELSE                                    ! R-CO(ONO2)
            alifun(21)=alifun(21)+1
            mapfun(i,1,21)=1
            ngrp=ngrp+1
         ENDIF
        ENDIF
      ENDDO
230   CONTINUE



* set the table telling if a "functional group" is available at a given node
      DO i=1,node
        DO j=1,3
c        DO j=1,2
c          DO k=1,20
          DO k=1,21
            IF (mapfun(i,j,k).ne.0.) THEN
              IF (mapfun(i,j,k).lt.1.) THEN
                funflg(i)=funflg(i)+1
              ELSE
                funflg(i)=funflg(i)+INT(mapfun(i,j,k)) 
              ENDIF
              IF (j.eq.3) THEN
                aroflg(i)=aroflg(i)+1
                arobig(i)=1
                IF (k.eq.1) arobig(i)=0
                IF (k.eq.5) arobig(i)=0
                IF (k.eq.6) arobig(i)=0
                IF (k.eq.7) arobig(i)=0
                IF (k.eq.8) arobig(i)=0
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

* END 
      END
