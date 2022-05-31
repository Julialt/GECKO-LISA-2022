      SUBROUTINE cleanup

! Subroutine to deallocate variables, close files, and 
! perform any other end-of simulation housekeeping tasks

      USE flags_module,ONLY: iofmt_fg,tracer_fg,iscape
      USE io_units_module,only: lout,lpbl,lsoa,ljval,lro2, &
                                lppf,lppa,lpff,lpaa
      USE printphoto_module
      USE NetCDF_vars_module,ONLY: ncid_out
      USE forcing_params_module,ONLY: nbox
      IMPLICIT NONE

!--------------------------------------------------------
!==Deallocate Variables==

      IF(ALLOCATED(cnjv)) DEALLOCATE(cnjv)
      IF(ALLOCATED(idprintphoto)) DEALLOCATE(idprintphoto)
      IF(ALLOCATED(photorates)) DEALLOCATE(photorates)
      IF(ALLOCATED(photoreac_names)) DEALLOCATE(photoreac_names)
      IF(ALLOCATED(cnjv)) DEALLOCATE(cnjv)

!==CLOSE ALL OUTPUT FILES==
!--Ascii
      CLOSE(lout)
      CLOSE(lpbl)
      CLOSE(lsoa)
      CLOSE(lro2)
      CLOSE(ljval)

      CLOSE(38)

      IF (tracer_fg.EQ.1) CLOSE(30)

      IF (iscape.eq.1) THEN
        CLOSE(20)
        CLOSE(67)
        CLOSE(71)
        CLOSE(72)
        CLOSE(73)
        CLOSE(74)
        CLOSE(75)

        IF (nbox.GT.1) THEN
          CLOSE(81)
          CLOSE(82)
          CLOSE(83)
          CLOSE(84)
          CLOSE(85)
        ENDIF

      ENDIF

!--Binary
      IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN
        CLOSE(lppf)
        CLOSE(lppa)
        IF (nbox.GT.1) THEN
          CLOSE(lpff)
          CLOSE(lpaa)
        ENDIF
      ENDIF ! (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2)

!--NetCDF
! (no need to close ncid_in or ncid_prev since they are readonly)
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        CALL close_ncfile(ncid_out)
      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

!--------------------------------------------------------
      END SUBROUTINE cleanup
!--------------------------------------------------------
