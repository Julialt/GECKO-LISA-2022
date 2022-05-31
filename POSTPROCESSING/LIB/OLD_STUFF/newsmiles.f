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
      INTEGER         cnum
      INTEGER         onum
      INTEGER         dbflg
      INTEGER         bond(mca,mca)
      INTEGER         path(mca,mca,mca)
      INTEGER         clngth(mca,mca)
      INTEGER         rjg(mri,2)
      CHARACTER*(lgr) group(mca),tgrp(mca)
      CHARACTER*(lgr) grpnew,grpold
      LOGICAL         lobond(mca,mca)

      CHARACTER*(lfo) tsmile,tchem,copy(mco)
      INTEGER         nc, nca, nce, nring 
      INTEGER         i, j, k
      INTEGER         kn,ko,kg
      INTEGER         rank(mca)
      INTEGER         ncbn,nhyd,noxy,nnit
      INTEGER         p, ptr

* ------------------------------------------------------
* INITIALIZE
* ------------------------------------------------------
!      PRINT*,'subroutine smiles'
!      PRINT*,'in smiles: ',trim(chem)

      ncbn = 0
      nhyd = 0
      noxy = 0
      nnit = 0
      smile(1:lfo)=" "
      tsmile(1:lfo)=" "

      DO i=1,mco
        copy(i) = ' '
      ENDDO

      DO i=1,mca
        DO j=1,mca
          clngth(i,j) = 0
          lobond(i,j) = .false.
          DO k=1,mca
            path(i,j,k) = 0
          ENDDO
        ENDDO
      ENDDO


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
        PRINT*,trim(smile)
        tchem(1:lfo-3)=chem(4:lfo)//'   '
        chem=tchem
      ELSE
        IF(INDEX(chem,'#').NE.0)THEN
          smile = 'handwritten compound'
          PRINT*,trim(smile)
          tchem(1:lfo-1)=chem(2:lfo)//'   '
          chem=tchem
        ENDIF
      ENDIF

* get the number of >C< , >c< (nca) and -O- (nce) in the molecule
* this also filters out named species
      nc  = INDEX(chem,' ') - 1
      ncbn = cnum(chem,nc)
      nce = onum(chem,nc)
      nca = ncbn + nce

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

* ------------------------------------------------------
* build the bond and group matrix
* ------------------------------------------------------
      CALL stdchm(chem)
      !PRINT*,"stdchm op: ",trim(chem)

      CALL grbond(chem,nc,group,bond,dbflg,nring)

      CALL ratings(nca,group,bond,nring,rank)

      IF(nring.GT.0) THEN
      !  smile = 'ERROR:bug_in_SMILES_for_cyclics'
      !  RETURN
        CALL uniqring(nring,nca,group,bond,rank,rjg)
        DO k=1,nring
          i=rjg(k,1)
          j=rjg(k,2)
          bond(i,j) = 0
          bond(j,i) = 0
        ENDDO
        CALL rjgrm(nring,group,rjg)
      ENDIF

      !DO i=1,nca
      !  PRINT*,i,group(i),bond(i,1:nca)
      !ENDDO
      !PRINT*,'----------'

* ------------------------------------------------------
* convert individual groups, counting atoms as you go.
* ------------------------------------------------------
      DO i=1,nca
        tgrp(i)=group(i)

* a) remove CH* hydrogens
        IF(INDEX(tgrp(i),'CH').NE.0.OR.INDEX(tgrp(i),'cH').NE.0)THEN
          ptr=MAX(INDEX(tgrp(i),'CH')+1,INDEX(tgrp(i),'cH')+1)
          IF(INDEX(tgrp(i),'CH2').NE.0)THEN
            k=2
            nhyd = nhyd+2
          ELSE IF(INDEX(tgrp(i),'CH3').NE.0)THEN
            k=2
            nhyd = nhyd+3
          ELSE
            k=1
            nhyd = nhyd+1
          ENDIF
          tgrp(i)(ptr:mca-k)=tgrp(i)(ptr+k:mca)
          !PRINT*,"a)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF

* b) remove CdH hydrogens
        IF(INDEX(tgrp(i),'CdH').NE.0)THEN
          ptr=INDEX(tgrp(i),'CdH')+2
          IF(INDEX(tgrp(i),'CdH2').NE.0)THEN
            k=2
            nhyd = nhyd+2
          ELSE
            k=1
            nhyd = nhyd+1
          ENDIF
          tgrp(i)(ptr:mca-k)=tgrp(i)(ptr+k:mca)
          !PRINT*,"b)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF

* c) convert NO2 to N(=O)=O {N,P,V}
        IF(INDEX(tgrp(i),'ONO2').NE.0)THEN
          noxy = noxy + 1
        ENDIF
        IF(INDEX(tgrp(i),'NO2').NE.0)THEN
          noxy = noxy + 2
          nnit = nnit + 1
          grpold = 'NO2'
          grpnew = 'N(=O)=O'
          ko=LEN(TRIM(grpold))
          kn=LEN(TRIM(grpnew))

          ptr = INDEX(tgrp(i),TRIM(grpold))
          kg = LEN(tgrp(i))
          tgrp(i)(ptr+kn:kg+kn-ko)=tgrp(i)(ptr+ko:kg)
          tgrp(i)(ptr:ptr+kn-1)=grpnew
          !PRINT*,"c)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF

* d) convert CO to C(=O) {A,G,P} 
        IF(INDEX(tgrp(i),'CO').NE.0)THEN
          noxy = noxy + 1
          grpold = 'CO'
          grpnew = 'C(=O)'
          ko=LEN(TRIM(grpold))
          kn=LEN(TRIM(grpnew))

          ptr = INDEX(tgrp(i),TRIM(grpold))
          kg = LEN(tgrp(i))
          tgrp(i)(ptr+kn:kg+kn-ko)=tgrp(i)(ptr+ko:kg)
          tgrp(i)(ptr:ptr+kn-1)=grpnew
          !PRINT*,"d)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF

* e) convert OH to O {A,G,H,O} 
        IF(INDEX(tgrp(i),'(OOH)').NE.0)THEN
          noxy = noxy + 2
          nhyd = nhyd + 1
        ENDIF
        IF(INDEX(tgrp(i),'(OH)').NE.0)THEN
          noxy = noxy + 1
          nhyd = nhyd + 1
        ENDIF
        IF(INDEX(tgrp(i),'OH').NE.0)THEN
          grpold = 'OH'
          grpnew = 'O'
          ko=LEN(TRIM(grpold))
          kn=LEN(TRIM(grpnew))

          ptr = INDEX(tgrp(i),TRIM(grpold))
          kg = LEN(tgrp(i))
          tgrp(i)(ptr+kn:kg+kn-ko)=tgrp(i)(ptr+ko:kg)
          tgrp(i)(ptr:ptr+kn-1)=grpnew
          !PRINT*,"e)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF
        
* f) convert -O- to O {E} 
        IF(INDEX(tgrp(i),'-O-').NE.0)THEN
          noxy = noxy + 1
          grpold = '-O-'
          grpnew = 'O'
          ko=LEN(TRIM(grpold))
          kn=LEN(TRIM(grpnew))

          ptr = INDEX(tgrp(i),TRIM(grpold))
          kg = LEN(tgrp(i))
          tgrp(i)(ptr+kn:kg+kn-ko)=tgrp(i)(ptr+ko:kg)
          tgrp(i)(ptr:ptr+kn-1)=grpnew
          !PRINT*,"f)",tgrp(i)(1:INDEX(tgrp(i)," ")),ncbn,nhyd,noxy,nnit
        ENDIF
        
!2  new lines: string groups together in order presented w/o any parentheses
!        iblank=INDEX(tsmile," ") 
!        tsmile=TRIM(tsmile(1:iblank))//trim(ADJUSTL(tgrp(i)))
      ENDDO
      !PRINT*,'----------'
* ------------------------------------------------------
* add ring chars back in, rebuild chem as a SMILE.
      IF(nring.GT.0) CALL rjgadd(nring,tgrp,rjg)

      CALL rebond(bond,tgrp,tsmile,nring)
* The SMILE is non-standard - we would have to create a routine similar
* to prioty.f to align with the UNIQUE SMILES standard

!------------------------------------------
* remove CdH hydrogens
* remove 'd'

60    IF(INDEX(tsmile,'d').NE.0) THEN
        ptr=INDEX(tsmile,'d')
        tsmile(ptr:lfo-1)=tsmile(ptr+1:lfo)
      ENDIF
      IF(INDEX(tsmile,'d').NE.0) GOTO 60

* ------------------------------------------------------
* RETURN THE STANDARDIZED SMILES FORMULA 
* ------------------------------------------------------

      smile = tsmile
      IF(INDEX(smile," ").EQ.1)THEN
        smile(1:lfo-1)=smile(2:lfo)
        smile(lfo:lfo)=" "
      ENDIF
       
      !PRINT*,"C_H_O_N_",ncbn,nhyd,noxy,nnit
      !PRINT*,"end smile: ",trim(smile)

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
