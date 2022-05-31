      subroutine setup_GECKO

!-----------------------------------------------------------------------------
!     subroutine to do initial input handling for GECKO
!-----------------------------------------------------------------------------
!----------------------------------------------------------------------------
      USE flags_module,ONLY: printphoto_fg

      IMPLICIT NONE

!====================================================================
! END OF DECLARATIONS
!====================================================================

      PRINT*,"----------------setup_GECKO----------------"

! -------------------------------------
! INITIALIZE DEFAULTS for some values & flags
      CALL init_defaults

! READ FLAGS for simulation from indat*.key
      CALL readkeyflags

! OPEN FILES to read inputs and log run progress
      CALL open_io_files 
  
! READ INPUTS: DICTIONARY & MECHANISM INPUTS & INDEX MAPPING
      CALL get_cheminp

! READ INPUTS: ENVIRONMENTAL INPUTS & CONSTRAINTS
      CALL get_envinp

! -------------------------------------
! SET UP OUTPUTS:
      IF (printphoto_fg) CALL setup_printphoto

      CALL open_op_files ! must be done AFTER reading inputs

! -------------------------------
      END SUBROUTINE setup_GECKO
! -------------------------------

