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
      IMPLICIT NONE
      INCLUDE 'general.h'

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
      INTEGER         leaf, last, ptr, p, i, j, k, ncp
      INTEGER         ig, pg, ng
      INTEGER         kn,ko,kg
      INTEGER         rank(mca)
      INTEGER         ncbn,nhyd,noxy,nnit

* ------------------------------------------------------
* INITIALIZE
* ------------------------------------------------------
!      print*,'subroutine smiles'

      ncbn = 0
      nhyd = 0
      noxy = 0
      nnit = 0

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
        smile = 'N/A'
        RETURN
      ENDIF

* CHECK FOR '#mm' notation
      IF(INDEX(chem,'#mm').NE.0)THEN
        tchem=chem(4:lfo)//'   '
        chem=tchem
!        print*,chem
      ELSE
        IF(INDEX(chem,'#').NE.0)THEN
          PRINT*,"species ",chem," cannot be parsed by smiles"
          RETURN
        ENDIF
      ENDIF

* get the number of >C< , >c< (nca) and -O- (nce) in the molecule
      nc  = INDEX(chem,' ') - 1
      ncbn = cnum(chem,nc)
      nce = onum(chem,nc)
      nca = ncbn + nce
* only multi-node molecules checked:
      IF (nca.LE.1) THEN
        smile = 'N/A'
        RETURN
      ENDIF

* ------------------------------------------------------
* build the bond and group matrix
* ------------------------------------------------------
      CALL grbond(chem,nc,group,bond,dbflg,nring)
      CALL ratings(nca,group,bond,nring,rank)
      IF(nring.GT.0) THEN
        CALL uniqring(nring,nca,group,bond,rank,rjg)
        DO k=1,nring
          i=rjg(k,1)
          j=rjg(k,2)
          bond(i,j) = 0
          bond(j,i) = 0
        ENDDO
      ENDIF


      !DO i=1,nca
      !  PRINT*,group(i)      
      !ENDDO
      !PRINT*,'----------'

      IF(nring.GT.0) CALL rjgrm(nring,group,rjg)
 
* ------------------------------------------------------
* convert individual groups, counting atoms as you go.
* ------------------------------------------------------
      DO i=1,nca
        !PRINT*,group(i)      
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
        

      ENDDO
* ------------------------------------------------------
* add ring chars back in
      IF(nring.GT.0) CALL rjgadd(nring,tgrp,rjg)

* ------------------------------------------------------
* make a logical copy of the bond matrix
      DO i=1,nca
        DO j=1,nca
          IF (bond(i,j).NE.0) lobond(i,j)= .true.
        ENDDO
      ENDDO

* find longest tree, top-down, starting with the first group
      CALL lntree(bond,dbflg,1,2,nca,clngth,path)

* look down-top for the very longest tree ...
      ncp = 0
      DO 50 i=1,nca
        IF (clngth(1,i).NE.0) THEN
           leaf = path(1,i,clngth(1,i))
           last = path(1,i,clngth(1,i)-1)
           CALL lntree(bond,dbflg,leaf,last,nca,clngth,path)

* write the SMILES of species according to the longest tree 
* in path. Each group are written in the formula with all the
* branches by the mkcopy subroutine. OMIT revers subroutine. (Should be
* in standard format already.) 
           DO 55 j=1,nca
             IF (clngth(leaf,j).NE.0) THEN
                ncp=ncp+1
                ptr = 1
                DO k=1,nca
                  ig = path(leaf,j,k)
                  IF (ig.NE.0) THEN
                     IF (k.GT.1)   pg = path(leaf,j,k-1)
                     IF (k.LT.nca) ng = path(leaf,j,k+1)
        CALL mkcopy (lobond,tgrp,nca,rank,nring,ig,pg,ng,ptr,copy(ncp))
                  ENDIF
                ENDDO
             ENDIF
55         CONTINUE
        ENDIF
50    CONTINUE

* just use the first available version - SMILES doesn't care
      tsmile = copy(ncp)

* write double bonds
      IF (dbflg.NE.0) CALL dwrite(tsmile)

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
!      PRINT*,"C_H_O_N_",ncbn,nhyd,noxy,nnit

! skip the parentheses check
      RETURN

* check parenthesis:if open parentheses,then stop the run.
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
!        STOP
      ENDIF
     
      END
