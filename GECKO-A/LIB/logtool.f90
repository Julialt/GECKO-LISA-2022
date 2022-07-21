MODULE logtool
IMPLICIT NONE
CONTAINS
! ======================================================================
! Purpose: write operating conditions and flags for memo to a log file.
! ======================================================================
SUBROUTINE wrtlog()
  USE keyparameter, ONLY: logu
  USE keyflag
  USE keyuser
  IMPLICIT NONE
  
! write mechanism parameters to stdout, for logging
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' GECKO VERSION=v0.0.0 '
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' maxgen = ', maxgen
  WRITE(logu,*) ' default temperature: ', dT
  WRITE(logu,*) ' default 3rd body M (molec cm-3): ', dM
  WRITE(logu,*) '          '

! -----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    ----  MECHANISM  ---- '
  WRITE(logu,*) ' ---------------------------- '

  IF   (g2pfg) THEN; WRITE(logu,*) ' gas/particle mass transfer: included'
  ELSE             ; WRITE(logu,*) ' gas/particle mass transfer: not included'
  ENDIF

  IF   (g2wfg) THEN; WRITE(logu,*) ' gas/wall mass transfer: included'
  ELSE             ; WRITE(logu,*) ' gas/wall mass transfer: not included'
  ENDIF

  IF (dhffg) THEN; WRITE(logu,*) ' DHF isomerisation: yes'
  ELSE           ; WRITE(logu,*) ' DHF isomerisation: no'
  ENDIF
  WRITE(logu,*) '          '

! ----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    --------- SAR -------- '
  WRITE(logu,*) ' ---------------------------- '
  IF     (pvap_sar==1) THEN; WRITE(logu,*) ' vapor pressure scheme=M&Y '
  ELSEIF (pvap_sar==2) THEN; WRITE(logu,*) ' vapor pressure scheme=Nannoolal '
  ELSEIF (pvap_sar==3) THEN; WRITE(logu,*) ' vapor pressure scheme=SIMPOL-1 '
  ENDIF  

  IF     (kohaddfg==1) THEN; WRITE(logu,*) ' OH addition SAR =Peeters 1997 '
  ELSEIF (kohaddfg==2) THEN; WRITE(logu,*) ' OH addition SAR =Ziemann 2009 '
  ELSEIF (kohaddfg==3) THEN; WRITE(logu,*) ' OH addition SAR =Jenkin 2018  '
  ENDIF
  IF (rx_ro2oh) WRITE(logu,*) ' RO2 + OH reactions activated'

  IF     (kisomfg==1) THEN; WRITE(logu,*) ' Alkoxy isomerisation: Atkinson 2007 '
  ELSEIF (kisomfg==2) THEN; WRITE(logu,*) ' Alkoxy isomerisation: Vereecken 2009 '
  ENDIF
  IF     (kdissfg==1) THEN; WRITE(logu,*) ' Alkoxy dissociation: Atkinson 2007 '
  ELSEIF (kdissfg==2) THEN; WRITE(logu,*) ' Alkoxy dissociation: Vereecken 2009 '
  ENDIF

! ----

  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) '    ----- Reductions ----- '
  WRITE(logu,*) ' ---------------------------- '
  WRITE(logu,*) ' critical vapor pressure=', critvp
  IF (isomerfg) THEN; WRITE(logu,*) ' isomerisation allowed '
  ELSE                 ; WRITE(logu,*) ' no isomerisation '
  ENDIF
  WRITE(logu,*) ' high-NOx flag:',highnox
  WRITE(logu,*) ' All classes of RO2 in RO2+RO2: ', multiclass

END SUBROUTINE wrtlog

END MODULE logtool
