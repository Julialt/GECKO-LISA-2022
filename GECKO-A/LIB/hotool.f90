MODULE hotool
IMPLICIT NONE
CONTAINS

!=======================================================================
! PURPOSE :                                                           
!   Find rate constants for H-atom abstraction from VOC by HO based on
!   Jenkin et al., ACP, 2018
!=======================================================================
SUBROUTINE rabsoh(tbond,tgroup,ig,arrhc,nring)
  USE keyparameter, ONLY: mxcp,mxlfo,saru
  USE keyflag, ONLY: sar_info
  USE ho_aromtool, ONLY: arom_data
  USE ringtool, ONLY: ring_data
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: tbond(:,:) ! node matrix for the current species
  CHARACTER(LEN=*),INTENT(IN):: tgroup(:)  ! groups for the current species
  INTEGER,INTENT(IN)         :: ig         ! group bearing the leaving H 
  INTEGER,INTENT(IN)         :: nring      ! # of rings
  REAL,INTENT(OUT)           :: arrhc(:)   ! arrh. coef. for H abstraction
!  CHARACTER(LEN=*),INTENT(INOUT) :: com(:)     ! comment about the rate constant 

  INTEGER  :: i,j,k
  REAL     :: mult
  INTEGER  :: nca,nether
  INTEGER  :: nbF                          ! # of "correction factors" applied
  INTEGER  :: o_sub(2),m_sub(2),p_sub(2)   ! branches in o,m,p position - if any
  INTEGER  :: nCd                          ! # of Cd bounded to ig
  INTEGER  :: alf(3)                       ! index of Cd neighbors
  INTEGER  :: Cdtable(4)                   ! Cd nodes (bounded to ig)
  LOGICAL  :: rmFc

! initialize
! -----------
  arrhc(:)=0. 
  nca=COUNT(tgroup/=' ')
  IF (ig>nca) STOP "ig > nca in rabsoh" 
      
  IF (sar_info==1) THEN
    WRITE(saru,*) '   '                                  !! debug
    WRITE(saru,*) '============= RABSOH ======'          !! debug
    DO i=1,nca ;  WRITE(saru,*) i,tgroup(i) ; ENDDO      !! debug
    WRITE(saru,*) 'ig,tgroup(ig)',ig,tgroup(ig)          !! debug
  ENDIF
  
! add the references
!  DO i=1,SIZE(com)-1
!    IF (com(i)==' ') THEN
!      com(i)='HABSAR' ; com(i+1)='MJ18KMV000'
!    ENDIF
!  ENDDO     
            
! --------------------
! --- FIND K(0) VALUE
! ------------------ -
  IF (tgroup(ig)(1:3)=='CH3') THEN
    arrhc(1)=2.90E-12 ; arrhc(3)=925.

  ELSE IF(tgroup(ig)(1:3)=='CH2') THEN
    arrhc(1)=4.95E-12 ; arrhc(3)=555.

  ELSE IF(tgroup(ig)(1:3)=='CHO') THEN
    DO i=1,nca
      IF (tbond(ig,i)==3) THEN                     ! -O-CHO group
        arrhc(1)=1.70E-12 ; arrhc(3)=910.
        IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ELSE IF (tbond(ig,i)==1) THEN
        DO j=1,nca
          IF ((tbond(j,i)==1).AND.(j/=ig)) THEN    ! C(OH)-C(polar)-CHO
            IF ((INDEX(tgroup(i),'O')==0).AND.(INDEX(tgroup(j),'(OH)')/=0)) THEN   
              arrhc(1)=11.7E-12 ; arrhc(3)=78.
              IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
            ENDIF
          ENDIF
        ENDDO
        IF (tgroup(i)(1:3)=='CH3')       THEN ; arrhc(1)=4.60E-12 ; arrhc(3)=-350. ! CH3CHO
        ELSE IF (tgroup(i)(1:2)=='Cd')   THEN ; arrhc(1)=1.30E-11 ; arrhc(3)=0.    ! C=CCHO
        ELSE IF (tgroup(i)(1:4)=='CH2 ') THEN ; arrhc(1)=5.08E-12 ; arrhc(3)=-420. ! -CH2CHO
        ELSE IF (tgroup(i)(1:4)=='CHO ') THEN ; arrhc(1)=1.55E-12 ; arrhc(3)=-340. ! CHOCHO
        ELSE IF (tgroup(i)(1:5)=='CH2(O')THEN ; arrhc(1)=11.7E-12 ; arrhc(3)=140.  ! CH2(polar)CHO
        ELSE IF (tgroup(i)(1:2)=='CO')   THEN ; arrhc(1)=1.78E-12 ; arrhc(3)=-590. ! -COCHO           
        ELSE IF ((tgroup(i)(1:4)=='CH(O').OR.(tgroup(i)(1:3)=='C(O')) THEN         ! >C(polar)CHO
          arrhc(1)=11.7E-12 ; arrhc(3)=-25.            
        ELSE IF (tgroup(i)(1:1)=='c') THEN                                         ! phi-CHO
          arrhc(1)=1.21E-11 ; arrhc(3)=0.  
          CALL arom_data(i,tgroup,tbond,nca,o_sub,m_sub,p_sub)
          DO k=1,2
            IF (o_sub(k)/=0) THEN 
              IF (INDEX(tgroup(o_sub(k)),'O')==0) arrhc(3)=arrhc(3)-115.
            ENDIF
            IF (m_sub(k)/=0) THEN 
              IF (INDEX(tgroup(m_sub(k)),'O')==0) arrhc(3)=arrhc(3)-78.
            ENDIF
            IF (p_sub(k)/=0) THEN 
              IF (INDEX(tgroup(p_sub(k)),'O')==0) arrhc(3)=arrhc(3)-78.
            ENDIF
          ENDDO
        ELSE IF (tgroup(i)(1:1)=='C')   THEN ; arrhc(1)=5.22E-12 ;  arrhc(3)=-490. ! CCHO           
        ELSE
          WRITE(6,*) 'No constant found in rabsoh for group: ', tgroup(i)
          STOP "in rabsoh, expected reaction not found"
        ENDIF
        IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ENDIF
    ENDDO
    
  ELSE IF(tgroup(ig)(1:2)=='CH') THEN
    arrhc(1)=3.17E-12 ; arrhc(3)=225.

  ELSE
    WRITE(6,*) '--error--, in rabsoh. No reaction found for group: ', tgroup(ig)
    STOP "in rabsoh"
  ENDIF
      
! K0 for alpha ether
! -----------------
  eloop: DO i=1,nca
    IF (tbond(ig,i)==3) THEN
      DO j=1,nca                                  ! skip -CO-O-CH-
        IF ( (tbond(i,j)==3).AND.(j/=ig).AND.&   
             ((tgroup(j)(1:3)=='CO ').OR.(tgroup(j)(1:3)=='CHO')) ) THEN        
          CYCLE eloop
        ENDIF
      ENDDO

      IF (tgroup(ig)(1:3)=='CH3') THEN
        arrhc(1)=2.22E-12 ; arrhc(3)=160.         ! -O-CH3
        IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ELSE IF (tgroup(ig)(1:2)=='CH') THEN        ! -O-CH<
        DO j=1,nca
          IF (tbond(ig,j)==1) THEN
            DO k=1,nca
              IF ((tbond(j,k)==3).AND.(nring==0)) THEN  ! -O-CH-C-O- 
                arrhc(1)=1.17E-12 ;  arrhc(3)=-760.    
                IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        arrhc(1)=1.20E-12 ; arrhc(3)=-460.        ! -O-CH< "regular" 
        IF (nring/=0) THEN
          CALL ringfac(ig,nca,tbond,tgroup,mult)
          arrhc(3)=arrhc(3)-298.*log(mult)   
        ENDIF
        IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc) ; RETURN !--- Out
      ENDIF
    ENDIF
  ENDDO eloop
  
! K0 for Aromatic branches (Fph1(1)=8.6, Fph1(3)=345, Fph2(1)=7.0, Fph2(3)=580)
  aloop: DO i=1,nca
    IF (tbond(ig,i)/=1) CYCLE
    IF (tgroup(i)(1:1)=='c') THEN
      IF ((tgroup(ig)(1:4)=='CH3 ').OR.(tgroup(ig)(1:4)=='CH2(')) THEN
        arrhc(1)=arrhc(1)*8.6 ; arrhc(3)=arrhc(3)+345.
      ELSE IF (tgroup(ig)(1:2)=='CH') THEN
        arrhc(1)=arrhc(1)*7.0 ; arrhc(3)=arrhc(3)+580.
      ENDIF
! scaling factor due to alkyl group in orto and para positions
      CALL arom_data(i,tgroup,tbond,nca,o_sub,m_sub,p_sub)
      DO k=1,2
        IF (o_sub(k)/=0) THEN 
          IF (INDEX(tgroup(o_sub(k)),'O')==0) arrhc(3)=arrhc(3)-140.
        ENDIF
        IF (p_sub(k)/=0) THEN 
          IF (INDEX(tgroup(p_sub(k)),'O')==0) arrhc(3)=arrhc(3)-140.
        ENDIF
      ENDDO
      EXIT aloop
    ENDIF
  ENDDO aloop
      
  IF (sar_info==1) WRITE(saru,*) ' 1) rate constant =',arrhc(1:3)

      
! -----------------------
! --- Ea CORRECTION 
! -----------------------

! FIND SUBSTITUENTS ON SAME CARBON
  mult=1.
  IF (INDEX(tgroup(ig),'(OH)')/=0)   mult=mult*3.6 
  IF (INDEX(tgroup(ig),'(ONO2)')/=0) mult=mult*0.16
  IF (INDEX(tgroup(ig),'(OOH)') /=0) mult=mult*3.6 
  IF (INDEX(tgroup(ig),'(OOOH)')/=0) mult=mult*3.6 
  IF (mult/=1.) arrhc(3)=arrhc(3)-298.*LOG(mult)

  IF (sar_info==1) THEN
    WRITE(saru,*) ' 2) fact on carbon bearing leaving H:',mult !! debug
    WRITE(saru,*) ' 3) Substituents factors:'                  !! debug
  ENDIF
            
! C=C neighbors - Seek the 2 neighbors (if any) in case of superallyl 
! resonant structures
  alf(:)=0 ; Cdtable(:)=0 ; nCd=0
  DO i=1,nca
    IF ((tgroup(i)(1:2)=='Cd').AND.(tbond(ig,i)==1)) THEN
      nCd=nCd+1 ; alf(nCd)=i
    ENDIF
  ENDDO

  IF (nCd/=0) THEN           ! CH-Cd=Cd
    Cdtable(1)=alf(1)        ! look for the Cd on the chain
    DO i=1,nca
      IF (tbond(Cdtable(1),i)==2) THEN ; Cdtable(2)=i ; EXIT ; ENDIF
    ENDDO
    DO i=1,nca
      IF ((tbond(Cdtable(2),i)==1).AND.(tgroup(i)(1:2)=='Cd')) THEN
        Cdtable(3)=i ; EXIT
      ENDIF
    ENDDO
    IF (Cdtable(3)/=0) THEN  ! CH-Cd=Cd-Cd=Cd-
      DO i=1,nca
        IF (tbond(Cdtable(3),i)==2) THEN ; Cdtable(4)=i ; EXIT ; ENDIF
      ENDDO
    ENDIF
  ENDIF

  IF (nCd==2) THEN           ! C=C at both side of the C-H node (overwrite if used above)   
    Cdtable(3)=alf(2)        ! look for the Cd on the chain
    DO i=1,nca
      IF (tbond(Cdtable(3),i)==2) THEN ; Cdtable(4)=i ; EXIT; ENDIF
    ENDDO
  ENDIF

  IF (nCd>2) STOP "unexpected 3 Cd attached to the C-H node"      

! FIND SUBSTITUENTS ON ALPHA CARBONS:
  alphaloop: DO i=1,nca
    mult=1. ; nbF=0 
    IF (tbond(ig,i)/=0) THEN

      IF (sar_info==1) THEN
        WRITE(saru,*) ' 3.1) alpha group factors'           !! debug
        WRITE(saru,*) '  >>> group_neighbor=',i,tgroup(i)   !! debug
      ENDIF

! simple alkyl:
      IF (tgroup(i)(1:3)=='CH3') THEN ; mult=1.   ; nbF=nbF+1 ; ENDIF
      IF (tgroup(i)(1:4)=='CH2 ')THEN ; mult=1.35 ; nbF=nbF+1 ; ENDIF
      IF (tgroup(i)(1:3)=='CH ') THEN ; mult=1.35 ; nbF=nbF+1 ; ENDIF
      IF (tgroup(i)(1:2)=='C ')  THEN ; mult=1.35 ; nbF=nbF+1 ; ENDIF

! contribution of functional groups

! alcohol
      IF ((tgroup(i)(1:2)/='CO').AND.(INDEX(tgroup(i),'(OH)')/=0)) THEN
        mult=mult*2.7 ; nbF=nbF+1
      ENDIF

! carbonyl          
      IF (tgroup(i)(1:2)=='CO') THEN
        DO j=1,nca
          IF (tbond(i,j)==3) THEN                  ! ester -O-CO-CH-
            mult=mult*0.4 ; nbF=nbF+1
            EXIT
          ENDIF
        ENDDO
        DO j=1,nca
          IF ((tbond(i,j)/=0).AND.(j/=ig)) THEN    ! CO(OH)-CO-CH-
            IF (tgroup(j)=='CO(OH)') THEN          
              mult=0.00000001 ; nbF=1                      ! reaction only occurs on -COOH side
              GOTO 10
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF (tgroup(i)(1:6)=='CO(OH)')    THEN; mult=mult*0.4 ; nbF=nbF+1; ENDIF
      IF (INDEX(tgroup(i),'(ONO2)')/=0)THEN; mult=mult*0.34; nbF=nbF+1; ENDIF
      IF (INDEX(tgroup(i),'(NO2)')/=0) THEN; mult=mult*0.31; nbF=nbF+1; ENDIF
      IF (tgroup(i)(1:10)=='CO(OONO2)')THEN; mult=0.1      ; nbF=nbF+1; ENDIF  ! BABABA skip ?????

! Esters
      IF (INDEX(tgroup(i),'-O-')/=0) THEN          
        DO j=1,nca
          IF ((tbond(i,j)==3).AND.(j/=ig)) THEN
            IF (tgroup(j)=='CHO') THEN              ! -CH-O-CHO
              mult=0.85 ; nbF=nbF+1
              GOTO 10
            ELSE IF (tgroup(j)=='CO') THEN          ! -CH-O-CO-
              mult=2.20 ; nbF=nbF+1
              GOTO 10
            ENDIF
          ENDIF                 
        ENDDO
      ENDIF

! FIND SUBSTITUENTS ON BETA CARBONS: 
! for -CH2CO- (but not -CO-O-), use MULT=3.4
! for -CH2CO-O- , use MULT=2.7 
! for -CH2-O-, >CH-O-,>C(-O-R)- , use MULT=3.5 (but not for -CH2-O-CO-)
! must remove the 1.35 contribution if there is an alkyl in alpha
      nether=0 ; rmFc=.FALSE.
      betaloop: DO j=1,nca
        IF ((tbond(i,j)/=0 ).AND.(j/=ig)) THEN
          IF (tgroup(j)(1:6)=='CO(OH)') THEN             ! CH-C-COOH
            IF ( ((tgroup(i)=='CH2 ').OR.(tgroup(i)=='CH ').OR.(tgroup(i)=='C ')).AND.&
                 (.NOT.rmFc) ) THEN
              mult=mult/1.35 ; nbF=nbF-1 ; rmFc=.TRUE.   ! rm contribution from aliphatic alpha C
            ENDIF
            mult=mult*2.7 ; nbF=nbF+1
            CYCLE betaloop
          
          ELSE IF (tgroup(j)(1:2)=='CO') THEN            ! CH-C-CO
            DO k=1,nca
              IF ((tbond(j,k)==1).AND.(k/=i)) THEN       ! CH-C-CO-C
                IF ( ((tgroup(i)=='CH2 ').OR.(tgroup(i)=='CH ').OR.(tgroup(i)=='C ')).AND.&
                     (.NOT.rmFc) )  THEN
                  mult= mult/1.35 ; nbF=nbF-1 ; rmFc=.TRUE.
                ENDIF
                IF (tgroup(i)(1:3)/='CO ') THEN          ! CH-CO-CO-C
                  mult=mult*3.4 ; nbF=nbF+1              
                  CYCLE betaloop
                ENDIF
              ELSE IF (tbond(j,k)==3) THEN               ! CH-C-CO-O
                IF ( ((tgroup(i)=='CH2 ').OR.(tgroup(i)=='CH ').OR.(tgroup(i)=='C ')).AND.&
                     (.NOT.rmFc) )  THEN
                  mult=mult/1.35 ; nbF=nbF-1 ; rmFc=.TRUE.
                ENDIF
                mult=mult*2.7 ; nbF=nbF+1
                CYCLE betaloop
              ENDIF
            ENDDO    
          
          ELSE IF (tgroup(j)(1:3)=='CHO') THEN           ! CH-C-CHO
            IF ( ((tgroup(i)=='CH2 ').OR.(tgroup(i)=='CH ').OR.(tgroup(i)=='C ')).AND.&
                 (.NOT.rmFc) )  THEN
                mult= mult/1.35 ; nbF=nbF-1 ; rmFc=.TRUE.
            ENDIF
            mult= mult*3.4 ; nbF=nbF+1
            CYCLE betaloop               
          
          ELSE IF (tgroup(j)=='-O-') THEN                ! CH-C-O-
            DO k=1,nca
              IF ((tbond(j,k)==3).AND.(k/=i)) THEN
                IF ((tgroup(i)(1:3)/='CO ').AND.(tgroup(k)(1:3)/='CO ').AND. & ! BABABABA IF impossible ???
                    (tgroup(k)(1:3)/='CHO')) THEN
                  nether=nether+1                        ! CH-CO-O-CO
                  IF (((tgroup(i)=='CH2 ').OR.(tgroup(i)=='CH ').OR. &
                      (tgroup(i)=='C ')).AND.(nether==1).AND.(.NOT.rmFc)) THEN
                    mult=mult/1.35 ; nbF=nbF-1 ; rmFc=.TRUE.
                  ENDIF                
                  mult=mult*3.5 ; nbF=nbF+1
                  CYCLE betaloop
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO betaloop
        
10    CONTINUE

! select the right factor
! first IF : only count contribution if the neighbor is the 1st Cd
! to avoid to count twice the contribution for Cd=Cd-CH-Cd=Cd structures
      IF (i==Cdtable(1)) THEN 
        IF (Cdtable(3)==0) THEN     ! only one C=C double bond
          IF (tgroup(Cdtable(2))(1:4)=='CdH2') THEN
            mult=mult*2.5 ; nbF=nbF+1
          ELSE 
            mult=mult*6.2 ; nbF=nbF+1
          ENDIF
        ELSE                          ! two C=C double bonds
          IF ((tgroup(Cdtable(2))(1:4)=='CdH2').AND.(tgroup(Cdtable(4))(1:4)=='CdH2')) THEN
            mult=mult*5.0  ; nbF=nbF+1
          ELSE IF ((tgroup(Cdtable(2))(1:4)=='CdH2').OR.(tgroup(Cdtable(4))(1:4)=='CdH2')) THEN
            mult=mult*8.7  ; nbF=nbF+1
          ELSE 
            mult=mult*12.4 ; nbF=nbF+1
          ENDIF
        ENDIF
      ENDIF

! rescale factors (for multiple interaction)
!      IF (nbF/=0) mult=mult**(1./nbF)
      IF (nbF>1) mult=mult**(1./nbF)
      arrhc(3)=arrhc(3)-298.*LOG(mult)

      IF (sar_info==1) THEN
        IF (nCd>0) WRITE(saru,*) 'nCd=',nCd              !! debug
        IF (nCd>0) WRITE(saru,*) 'Cdtable=',Cdtable(:)   !! debug
        IF (nbF>1) WRITE(saru,*) ' nbF=',nbF             !! debug
        WRITE(saru,*) ' F_neighbor=',mult                !! debug
      ENDIF
    ENDIF
  ENDDO alphaloop

! ADD RING FACTOR 
  IF (nring/=0) THEN
    CALL ringfac(ig,nca,tbond,tgroup,mult)
    arrhc(3)=arrhc(3)-298.*log(mult)   
  ENDIF

  IF (sar_info==1) CALL wrtkabs(saru,tgroup(ig),arrhc)

END SUBROUTINE rabsoh

!=======================================================================
! Purpose: compute the ring factor (ring strength) for H-abstraction 
! on a node belonging to rings. Data are from Jenkin et al., ACPs, 2018 
!=======================================================================
SUBROUTINE ringfac(ig,nca,tbond,tgroup,mfac)
  USE keyparameter, ONLY: saru
  USE ringtool, ONLY: ring_data
  USE keyflag, ONLY: sar_info
  IMPLICIT NONE

  INTEGER,INTENT(IN)         :: ig         ! group bearing the leaving H 
  INTEGER,INTENT(IN)         :: nca        ! # of nodes 
  INTEGER,INTENT(IN)         :: tbond(:,:) ! node matrix for the current species
  CHARACTER(LEN=*),INTENT(IN):: tgroup(:)  ! groups for the current species
  REAL,INTENT(OUT)           :: mfac       ! ring factor (@298 K) to be applied on Ea

  INTEGER,PARAMETER :: mxirg=6             ! max # of distinct rings 
  INTEGER  :: nring_ind                    ! # of disctinct rings
  INTEGER  :: trackrg(mxirg,SIZE(tgroup))  ! (a,:)== track (node #) belonging ring a
  LOGICAL  :: ring_ind(mxirg,SIZE(tgroup)) ! (a,b)==true if node b belong to ring a
  INTEGER  :: i,k
  REAL     :: mult
  INTEGER  :: rgord,nCO,nether

  mfac=1.
!nring_ind,ring_ind,trackrg) 
  CALL ring_data(ig,nca,tbond,tgroup,nring_ind,ring_ind,trackrg) 
  IF (sar_info==1) WRITE(saru,*) '# of distinct cycles:',nring_ind
  ringloop: DO i=1,nring_ind
    mult=1.
    rgord=COUNT(ring_ind(i,:).EQV..TRUE.)   ! get ring size                    

    SELECT CASE (rgord)
      CASE (3)
        mult=0.018 
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))=='-O-') mult=0.0079
        ENDDO
        
      CASE (4)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO
        IF ((nCO==0).AND.(nether==0))      THEN ; mult=0.41            ! factor for cycloalkane
        ELSE IF ((nCO/=0).AND.(nether==1)) THEN ; mult=(0.08*0.5)**0.5 ! square of the product of factors
        ELSE IF (nCO/=0)                   THEN ; mult=0.08
        ELSE IF (nether==1)                THEN ; mult=0.5  
        ENDIF                 
        
      CASE (5)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO
  
        IF ((nCO==0).AND.(nether==0))      THEN ; mult=0.69             ! factor for cycloalkane
        ELSE IF ((nCO/=0).AND.(nether==2)) THEN ; mult=(0.32*0.59)**0.5 ! square of the product of factors
        ELSE IF (nCO/=0)                   THEN ; mult=0.32
        ELSE IF (nether==2)                THEN ; mult=0.59             ! BABABA what if nether =1 ???
        ENDIF                       
        
      CASE (6)
        nether=0 ; nCO=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
          IF (tgroup(trackrg(i,k))(1:3)=='CO ') nCO=nCO+1
        ENDDO   
        IF ((nCO==0).AND.(nether==0))      THEN ; mult=0.95             ! factor for cycloalkane
        ELSE IF ((nCO/=0).AND.(nether==2)) THEN ; mult=(0.61*0.42)**0.5 ! square of the product of factors
        ELSE IF ((nCO/=0).AND.(nether==1)) THEN ; mult=(0.61*0.60)**0.5 ! square of the product of factors
        ELSE IF (nCO/=0)                   THEN ; mult=0.61
        ELSE IF (nether==1)                THEN ; mult=0.60
        ELSE IF (nether>1)                 THEN ; mult=0.42
        ENDIF        
                       
      CASE (7)
        mult=1. ; nether=0
        DO k=1,rgord
          IF (tgroup(trackrg(i,k))(1:3)=='-O-') nether=nether+1
        ENDDO
        IF (nether==1) mult=0.85
        IF (nether==2) mult=0.58
  
      CASE (8)
        mult=1.
    END SELECT
    
    IF (sar_info==1) THEN
      WRITE(saru,*) 'CYCLE of size:',rgord
      WRITE(saru,*) 'Fring=',mult
    ENDIF
!    arrhc(3)=arrhc(3)-298.*log(mult)
    mfac=mfac*mult   
  ENDDO ringloop
END SUBROUTINE ringfac 

!=======================================================================
! PURPOSE: small routine to write basic info
!=======================================================================
SUBROUTINE wrtkabs(lout,grp,arrhc)
  IMPLICIT NONE

  INTEGER :: lout
  CHARACTER(LEN=*) :: grp
  REAL :: arrhc(:)

  WRITE(lout,*)' '
  WRITE(lout,*)' 4) TOTAL RATE for abstraction at:',TRIM(grp)
  WRITE(lout,*)'  H abs by OH rate 298:',arrhc(1)*EXP(-arrhc(3)/298.)
END SUBROUTINE wrtkabs

END MODULE hotool
