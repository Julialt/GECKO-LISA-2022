      SUBROUTINE open_io_files
!-----------------------------------------------------------
! open NetCDF input & logfile output files required for GECKO-A box model
!-----------------------------------------------------------

      USE flags_module,ONLY: iofmt_fg,prevflag
      USE io_units_module,ONLY: lout
      USE NetCDF_vars_module,ONLY: ncid_in,ncid_prev

      IMPLICIT NONE

!-----------------------------------------------------------
!--open NetCDF input (link) file
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        PRINT*,"opening NetCDF input file"

        CALL open_ncfile_readonly('indat.nc',ncid_in)

        if (prevflag .eq. 1) then
          CALL open_ncfile_readonly('prevdat.nc',ncid_prev)
        endif

      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

!---------------------------------------------------------------------
! log file (ascii)
      PRINT*,"opening ascii log files"

      OPEN(lout,file='outdat.out',form='formatted',status='unknown')

      RETURN

!---------------------------------------------------------------------
      END SUBROUTINE open_io_files
!---------------------------------------------------------------------
