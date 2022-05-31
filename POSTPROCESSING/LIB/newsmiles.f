!***********************************************************************
* MASTER MECHANISM - ROUTINE NAME :      SMILES                        *
*                                                                      *
* PURPOSE:  -  convert standard chemical formula to SMILES format      *
*                                                                      *
* INPUT                                                                *
* - chem         : standard formula of species                         *
*                                                                      *
* OUTPUT                                                               *
* - smile        : SMILES-notation formula in same orientation as chem *
*
************************************************************************
      SUBROUTINE SMILES(chem,smile)
      USE PARAMETERS
      IMPLICIT NONE
* input:
      CHARACTER*(lfo) chem
* ouput:
      CHARACTER*(lfo) smile

* local
      CHARACTER*(lgr) grpnew,grpold

      CHARACTER*(lfo) tchem
      INTEGER         i,p
      INTEGER         nc 
! for standard code (not CU)
      INTEGER         nca,nce
      INTEGER         cnum, onum
      INTEGER         ncbn,nhyd,noxy,nnit

* ------------------------------------------------------
* INITIALIZE
* ------------------------------------------------------
!      PRINT*,'subroutine smiles'
!      PRINT*,'in smiles: ',TRIM(chem)

      smile(1:lfo)=" "

      p = INDEX(chem," ")
      chem(p:lfo)=" "

! for standard code (not CU)
      ncbn = 0
      nhyd = 0
      noxy = 0
      nnit = 0

* ------------------------------------------------------
* CHECK INPUT FORMULA
* ------------------------------------------------------
* CHECK FOR RADICALS (not represented)
      IF(INDEX(chem,'.').NE.0)THEN
        smile = 'N/A:radical'
        RETURN
      ENDIF

* CHECK FOR '#mm' notation
      IF(INDEX(chem,'#mm').NE.0.OR.INDEX(chem,'#bb').NE.0)THEN
        smile = 'aromatic/bioburning compound'
        !PRINT*,TRIM(smile)
        tchem(1:lfo-3)=chem(4:lfo)//'   '
        chem=tchem
        !PRINT*,TRIM(chem)
      ELSE
        IF(INDEX(chem,'#').NE.0)THEN
          smile = 'handwritten compound'
          !PRINT*,TRIM(smile)
          tchem(1:lfo-1)=chem(2:lfo)//'   '
          chem=tchem
        ENDIF
      ENDIF

!! for standard code (not CU)
* get the number of >C< , >c< (nca) and -O- (nce) in the molecule
* this also filters out named species
      nc  = INDEX(chem,' ') - 1
      ncbn = cnum(chem,nc)
      nce = onum(chem,nc)
      nca = ncbn + nce

      !PRINT*,"nca = ",nca
* only multi-node molecules checked:
      IF (nca.LE.1) THEN
        IF(INDEX(chem(1:9),'CH3(ONO2)').EQ.1) THEN
          smile = 'CO[N](=O)[O]'
        ELSEIF(INDEX(chem(1:8),'CH3(OOH)').EQ.1) THEN
          smile = 'COO'
        ELSEIF(INDEX(chem(1:7),'CHO(OH)').EQ.1) THEN
          smile = 'OC=O'
        ELSEIF(INDEX(chem(1:4),'CH2O').EQ.1) THEN
          smile = 'C=O'
        ELSE
          smile = 'C1/inorg: NO_SMILES_ASSESSED'
        ENDIF
        RETURN
      ENDIF
! End non-CU section

* ------------------------------------------------------
* convert individual groups, counting atoms as you go.
* ------------------------------------------------------

      smile = chem

! a) remove 'd'
      grpold = "d"
      grpnew = ""
      CALL chemtosmile(grpold,grpnew,smile)

! b) remove CH* hydrogens
      grpold = "H3"
      grpnew = ""
      CALL chemtosmile(grpold,grpnew,smile)

      grpold = "H2"
      grpnew = ""
      CALL chemtosmile(grpold,grpnew,smile)
!new
      grpold = "CHO"
      grpnew = "C(=O)"
      CALL chemtosmile(grpold,grpnew,smile)
!new
      grpold = "C(=O)  "
      grpnew = "C=O "
      CALL chemtosmile(grpold,grpnew,smile)

      grpold = "H"
      grpnew = ""
      CALL chemtosmile(grpold,grpnew,smile)

! c) convert NO2 to N(=O)=O {N,P,V}
      grpold = "NO2"
      grpnew = "N(=O)=O"
      CALL chemtosmile(grpold,grpnew,smile)

! d) convert CO to C(=O) {A,G,P} 
!new
      grpold = "CO("
      grpnew = "C(=O)("
      CALL chemtosmile(grpold,grpnew,smile)

      grpold = "CO"
      grpnew = "C(=O)"
      CALL chemtosmile(grpold,grpnew,smile)
      grpold = "C1O"
      grpnew = "C1(=O)"
      CALL chemtosmile(grpold,grpnew,smile)
      grpold = "C2O"
      grpnew = "C2(=O)"
      CALL chemtosmile(grpold,grpnew,smile)

! e) convert OH to O {A,G,H,O} (probably redundant, given 'b')
      grpold = "OH"
      grpnew = "O"
      CALL chemtosmile(grpold,grpnew,smile)

! f) convert -O- to O {E} 
      grpold = "-O"
      grpnew = "O"
      CALL chemtosmile(grpold,grpnew,smile)
      grpold = "-"
      grpnew = ""
      CALL chemtosmile(grpold,grpnew,smile)

! The SMILE is non-standard - we would have to create a routine similar
! to prioty.f to align with the UNIQUE SMILES standard
! and use stdchm-type routines to identify terminal groups.
!------------------------------------------

* ------------------------------------------------------
* RETURN THE STANDARDIZED SMILES FORMULA 
* ------------------------------------------------------

      IF(INDEX(smile," ").EQ.1)THEN
        smile(1:lfo-1)=smile(2:lfo)
        smile(lfo:lfo)=" "
      ENDIF
      
! JMLT test lines 
      !PRINT*,"C_H_O_N_",ncbn,nhyd,noxy,nnit
      !PRINT*,"end smile: ",TRIM(smile)

! skip the parentheses check
!      RETURN

* check parenthesis: if open parentheses,then stop the run.
      nc = INDEX(smile,' ') - 1
      p = 0
      DO i=1,nc
        IF(smile(i:i).EQ.'(') p = p + 1
        IF(smile(i:i).EQ.')') p = p - 1
      ENDDO

      IF(p.NE.0) THEN
        WRITE(6,'(a)') '--error-- 1'
        WRITE(6,'(a)') 'from routine : smile'
        WRITE(6,'(a)') 'parentheses mismatch in output SMILES:'
        WRITE(6,'(a)') chem
        WRITE(6,'(a)') smile
        STOP
      ENDIF
     
      END
!========================================================
      SUBROUTINE chemtosmile(grpold,grpnew,smile)

      USE PARAMETERS

* input:
      CHARACTER*(lgr) grpnew,grpold
* in/output:
      CHARACTER*(lfo) smile
* internal:
      INTEGER :: i,j,k,tend
      CHARACTER*(lfo) tchem

! -------------------------------------
      tchem(1:lfo)=" "

      j=LEN_TRIM(grpold)
      k=LEN_TRIM(grpnew)
      tend = INDEX(smile," ")
      !PRINT*,smile(1:tend)," ",grpold(1:j)," ",grpnew(1:k)

      DO WHILE(INDEX(smile(1:tend),grpold(1:j)).GT.0)
        i = INDEX(smile,grpold(1:j))
        IF(i.EQ.0) RETURN
        tend = INDEX(smile," ")
        tchem(1:i-1)=smile(1:i-1)
        tchem(i+k:tend+k)=smile(i+j:tend+j)
        IF(k.GT.0) tchem(i:i+k-1)=grpnew(1:k)
        smile = tchem
!        PRINT*,smile(1:LEN_TRIM(smile))
!new
        tend = INDEX(smile," ")
      ENDDO

      END SUBROUTINE chemtosmile
