      SUBROUTINE write_op_hdrs

      USE flags_module,ONLY: iofmt_fg,soa_fg!,printphoto_fg
      USE io_units_module,ONLY: lppf,lppa,lpff,lpaa!,ljval,lsoa,lro2
      !USE printphoto_module,ONLY: n_printphoto,photoreac_names
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main

      IMPLICIT NONE

!---------------------------------------------------
! BINARY OUTPUT HEADERS
      IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

! gas phase header info 
        CALL wrtlppf(lppf, numsp, chrsp)
        IF (nbox == 2) THEN
           CALL wrtlppf(lpff, numsp, chrsp)
        ENDIF

! particle phase header info
        IF (soa_fg.EQ.1) THEN

          CALL wrtlppa(lppa,nsat,ndim,chrsp,wmol,idsat)
          IF (nbox == 2) THEN
            CALL wrtlppa(lpaa, nsat,ndim, chrsp, wmol, idsat)
          ENDIF

        ELSE IF (soa_fg.EQ.2) THEN

          CALL wrtlppa(lppa,nsat,ndim,chrsp,wmol,idgsat)
          IF (nbox == 2) THEN
            CALL wrtlppa(lpaa,nsat,ndim,chrsp,wmol,idgsat)
          ENDIF

        ENDIF ! (soa_fg.EQ.1)

      ENDIF !(iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2)

!---------------------------------------------------
      END SUBROUTINE write_op_hdrs
!---------------------------------------------------
