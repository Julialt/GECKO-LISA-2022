! **********************************************************************
! PROGRAMMES APPELANTS :
! PROGRAMMES APPELES   : NONE
!
!  SUBROUTINE CALCULANT LES TRANSFORMATIONS
!  CHIMIQUES DES ESPECES PLUS
!  GAUSS-SEIDEL ITERATIONS POUR TWOSTEP.
!
! SIMPLIFIED NOTATION:
! nr = numre
! ns = numstoi
! nf = nself
! idr = idrestoi
! idp = idpdstoi
! idf = idselfreac
! cfr = restoicf
! cfp = pdstoicf
!
! 2016.05.06 VERSION iter4 combines Paris version iter3dyn
!            with Boulder version iter3_cfix
! 2016.05.28 version nodyn removes dyn variables for OFR & amb
! 2021.09.16 converted to f90 to implement Paris OMP fix (JMLT)
! **********************************************************************
      SUBROUTINE iter4                             
     $   (maxsp, maxre, mxleft, mxright, mself,     
     $    numsp, nr, ns, idr, idp, nf,idf,         
     $    cfr,cfp, rk, tot, gdt,                    
     $    rem, rdep, rex, rdil,                     
     $    c, xfr, xfl, trr, trl, conid,ncons,       
     $    noxfix, sumnox, idno, idno2, idno3, lpmap)

      !$ use OMP_LIB

      USE prodloss_module
      USE module_data_gecko_main, ONLY: small
      IMPLICIT NONE

      INTEGER maxsp, maxre, mxleft, mxright, mself

      INTEGER numsp, nr
      INTEGER ns(maxre,2)
      INTEGER idr(maxre,mxleft)
      INTEGER idp(maxre,mxright)
      INTEGER nf,idf(mself,2)
      REAL    cfr(maxre,mxleft)
      REAL    cfp(maxre,mxright)

! index map
      TYPE(spec_reac_map), TARGET, intent(in) :: lpmap(maxsp)
      TYPE(spec_reac_map), POINTER:: p_lpmap
      REAL    stpd

      REAL    GDT
      REAL    C(maxsp),RK(maxre),XFR(maxsp)
      REAL    XFL(maxsp),TOT(maxsp)
      REAL    TRL(maxsp),TRR(maxsp)
      INTEGER  ncons,conid(ncons)
      INTEGER  noxfix, idno, idno2, idno3
      REAL     sumnox


      REAL    rem(maxsp),rdep(maxsp)
      REAL    rex(maxsp),rdil

! local variables
      REAL    rate(maxre)
      INTEGER i, j, ire, ire2
      REAL    ctotaer
      !REAL,PARAMETER :: small = 1.0E-30
      REAL    cfix(ncons)

! new variables (for do-loop elimination)
      INTEGER    idarr(mxleft)

! initialisation
!$OMP PARALLEL DO private(i)
      DO i=1,numsp
        xfr(i)=0.
        xfl(i)=0.
        IF (c(i) .LT. small) c(i) = small
      ENDDO
!$OMP END PARALLEL DO

      cfix = 0.

!$OMP PARALLEL DO private(i)
      DO i=1,nr
        rate(i) = rk(i)
      ENDDO
!$OMP END PARALLEL DO

!************************************************************
! compute the reaction rate. The stoi. coef. of the reactants is 1,
! except for the species given in idf where it is 2.
! Do self reaction calculation second: mathematically identical to
! original code which required splitting the 1,nr loop.
! edited 210913 per issue #129
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!=====loop over 1,nr
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ire,i) SHARED(rate,ns,c,idr)
      DO ire=1,nr
        DO i=1,ns(ire,1)
         rate(ire)=rate(ire)*c(idr(ire,i))
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DO i=1,nf
        ire=idf(i,1)
        j  =idf(i,2)
        IF (ire.eq.0 .or.j.eq.0) cycle
        rate(ire)=rate(ire)*c(idr(ire,j))
      ENDDO

!************************************************************
! compute the production and loss expression for each species
! for the loss term, stoi. coef. is 1, except for the species
! given in idf where it is 2.
! using the index map now

!$OMP PARALLEL DO private(i, p_lpmap, j)
      DO i=1,numsp
        p_lpmap => lpmap(i)
        DO j=1,p_lpmap%nloss
          xfl(i) = xfl(i) + rate(p_lpmap%idl(j))
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO private(i, p_lpmap, j)
      DO i=1,numsp
        p_lpmap => lpmap(i)
        DO j=1,p_lpmap%nprod
          xfr(i) = xfr(i) + rate(p_lpmap%idp(j))*p_lpmap%stpd(j)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!=====END loop over 1,nr

!************************************************************
! correct the loss term for self reaction (relevant for inorg rxns only)
      DO i=1,nf
        ire=idf(i,1)
        j  =idf(i,2)
        IF (ire.eq.0 .or.j.eq.0) cycle
        xfl(idr(ire,j))=xfl(idr(ire,j)) + rate(ire)
      ENDDO

!$OMP PARALLEL DO private(i)
      DO i=1,numsp
!************************************************************
! correct the loss rate from the above expression
!************************************************************
! add the deposition and dilution term to the loss rate
        xfl(i) = xfl(i) / c(i) + rdep(i) + rdil

!************************************************************
! add the emission and exchange term to the production rate
        xfr(i) = xfr(i) + rem(i) + rex(i)
      ENDDO
!$OMP END PARALLEL DO


! save fixed (forced) concs IF applicable
      DO i=1,ncons
        IF (conid(i) .GT. 0) THEN
          cfix(i) = c(conid(i))
        ENDIF
      ENDDO

      IF (noxfix == 1) THEN
        CALL noxconstraint(c, sumnox, idno, idno2, idno3)
      ENDIF

!************************************************************
! compute the new concentration
!$OMP PARALLEL DO private(i)
      DO i=1, numsp
        c(i) = max(0.,(tot(i)+gdt*xfr(i))/(1.+gdt*xfl(i)) )
      ENDDO
!$OMP END PARALLEL DO

!! reinstate fixed (forced) concs
      DO i=1,ncons
        IF (conid(i) .GT. 0) THEN
          c(conid(i)) = cfix(i)
        ENDIF
      ENDDO

      IF (noxfix == 1) THEN
        CALL noxconstraint(c, sumnox, idno, idno2, idno3)
      ENDIF

!************************************************************
! END
      END
