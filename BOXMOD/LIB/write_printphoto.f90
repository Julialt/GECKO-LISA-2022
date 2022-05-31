      SUBROUTINE write_printphoto

      USE io_units_module,ONLY: ljval
      USE time_mgmt_module,ONLY: tout,itout
      USE printphoto_module
      USE NetCDF_vars_module,ONLY: ncid_out

      USE flags_module,ONLY: iofmt_fg
      IMPLICIT NONE

      INTEGER :: k

!------------------------------------------------------------------
! output photolysis rates (ascii)
!        IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

          WRITE(ljval,'('//cnjv//'(ES10.3,a1,1x),ES10.3)') &
                        tout,',', &
                        (photorates(k),',', &
                        k=1,n_printphoto-1), &
                        photorates(n_printphoto)
!        ENDIF !(iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

! output photolysis rates (NetCDF)
        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

          CALL eznc_put_1Dreal_into2D(ncid_out,"photorates", &
                                         photorates(1:n_printphoto), &
                                                    1,n_printphoto, &
                                                    itout)
        ENDIF !(iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

! -------------------------------------
      END SUBROUTINE write_printphoto
! -------------------------------------

