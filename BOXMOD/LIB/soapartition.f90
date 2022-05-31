!**********************************************************************
! This subroutine compute gas/aerosol partitioning assuming an ideal  *
! organic phase.                                                      *
!                                                                     *
! INPUT :                                                             *
!    - nsat : number of species for which psat to be calculated       *
!    - idsat(i) : index of species i                                  * 
!    - psat(i) : saturation vapor pressure                            *
!    - chrsp                                                          *
!    - numsp                                                          *
!                                                                     *
! INPUT/OUTPUT :                                                      *
!    - cgas(i)   : concentration of the species (i) in the gas        *
!                  phase (in molecule.cm-3)                           *
!    - caer(i)   : concentration of the species (i) in the aerosol    *
!                  phase (in molecule.cm-3)                           *
!**********************************************************************
      SUBROUTINE soapartition(temp,nsat,idsat,cgas,caer,psat,cnv)
      USE akparameter_module
      IMPLICIT NONE
      
      INTEGER  nsat,idsat(mxsat)
      REAL     psat(mxsat)
      REAL     cgas(maxsp)
      REAL     caer(mxsat)
      REAL     temp ! temperature
      REAL     cnv  ! non volatile species existing in the aerosol phase

! local
      REAL     cmin, cmax, cinf, csup, ctotit
      REAL     conv, ctot, ratio
      REAL     multifac
      INTEGER  i,niter, it, isp

! --------------------
! saturation vapor pressure(atm)  & temp (k) are passed in
! --------------------
! compute total concentration for each species (gas+aerosol)
! and provide a first estimate for the total number of molecule in 
! the aerosol phase based on the previous time step 
      ctotit=0.
      DO i=1,nsat
        cgas(idsat(i))=cgas(idsat(i))+caer(i)
        ctotit=ctotit+caer(i)
      ENDDO
      ctotit=ctotit+cnv
      WRITE(39,*) '--------------'
      
! R=0.082 L.atm./(K.mol), Na=6.02E23
      multifac = 6.02214E23 / (1000.*temp*0.0820578) 

! --------------------
! CONVERGENCE  LOOP
! --------------------
! intialize parameters for the iteration loops
      niter = 100     ! maximum number of iteration allowed
      cmin = 1E3+cnv      ! min conc. for the sum of condensed species
      cmax = 1E20+cnv     ! max conc. for the sum of condensed species

      cinf = cmin     ! min. conc obtained during the iteration process
      csup = cmax     ! min. conc obtained during the iteration process
      IF (ctotit.LE.cmin) THEN
        ctotit = cmin
      ENDIF

! start iterations
971   DO it = 1,niter
        ctot = 0.
        ratio = 0.

        DO isp = 1,nsat  
          caer(isp) = cgas(idsat(isp))/(1.+(psat(isp)*multifac/ctotit))
          ratio = ratio + cgas(idsat(isp))/(ctotit+(psat(isp)*multifac))
          ctot = ctot + caer(isp)
        ENDDO
        ctot = ctot + cnv
        ratio = ratio + (cnv/ctotit)
        WRITE(39,'(i3,3(1pe16.6))') it,ctotit,ctot, ratio

! check convergence. Convergence reached if total aerosol conc. is 
! modified by the last iteration by less than 1/1000 
        conv = (ctotit-ctot)/ctotit
        IF (abs(conv).lt.1E-3) GOTO 975

! ctot overestimated, real concentration must be between cinf and ctot
        IF (ratio.lt.1.) THEN
          IF (ctotit.le.1000) THEN
            ctot = 0.
            GOTO 975
          ENDIF 
          cinf = cinf
          csup = ctotit
          ctotit = (cinf*csup)**0.5
! ctot underestimated, real concentration must be between ctot and csup
        ELSE
          cinf = ctotit
          csup = csup
          ctotit = (cinf*csup)**0.5
        ENDIF
      ENDDO

! --------------------
! CONVERGENCE  FAILED
! --------------------
! If that point is reached, then convergence failed within the maximum
! number of iteration prescribed (niter). 2 cases can be considered :
! 1) If csup equal Cmax after the iteration process, then Cmax might 
! be too low. 2) the number of iteration (niter) is too low. Stop in 
! any cases
      IF (csup.eq.cmax) THEN
        WRITE(6,*) '--error--, convergence failed with csup=cmax'
        WRITE(6,*) 'change la valeur de cmax'
        STOP
      ELSE
        WRITE(6,*) '--error--, convergence failed it=niter'
        WRITE(6,*) 'change la valeur de niter'
        STOP
      ENDIF

! -----------------------
! CONVERGENCE  SUCCEEDED
! -----------------------
! If that point is reached, then convergence succeeded.
975   CONTINUE

! If a threshold of 1000 molecule.cm-3 condensed is not reached
! then reset aeorols concentration to 0. If above threshold, then
! compute the partioning of every species and return.
      IF (ctot.eq.0) THEN
        DO isp = 1,nsat
          caer(isp) = 0.
        ENDDO
      ELSE
        DO isp = 1,nsat
          caer(isp) = cgas(idsat(isp)) / (1+(psat(isp)*multifac/ctot))
          cgas(idsat(isp)) = cgas(idsat(isp)) - caer(isp)    
        ENDDO            
      ENDIF

! END
      RETURN
      END
