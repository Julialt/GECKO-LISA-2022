      SUBROUTINE soa_dyn(ctotaerbox,maerbox)

      USE akparameter_module
      USE flags_module,ONLY: dyn_fg,wall_fg
      USE io_units_module,ONLY: lout
      USE NetCDF_vars_module,ONLY: ncid_prev,ncid_out
      USE solver_params_module,ONLY: dtmin,dtmax,numit
      USE forcing_params_module,ONLY: temp
      USE module_data_gecko_main,ONLY: ratfac,winfac,weqfac, &
                         numain,numaou,numwin,numwou, &
                         idain,idaou,idwin,idwou,aoucf,woucf, &
                         temp,cnv,sumc,ibox,gamm,Mp,Rpo,psat, &
                         wmol,imtr,idrestoi,Rp,qfor 

      IMPLICIT NONE

! LOCAL VARIABLES
      REAL    :: maerbox       ! = maer(ibox)
      REAL    :: ctotaerbox    ! = ctotaer(ibox)

! -------------------------------------
      !PRINT*,"starting dynamic soa"

      IF (dyn_fg.EQ.0) THEN
        ratfac=6.E-3  ! simple multiplicative factor gas-> aero rate
        CALL mtrat_ba(ratfac,winfac,weqfac, &
                       numain,numaou,numwin,numwou, &
                       idain,idaou,idwin,idwou,aoucf,woucf, &
                       temp,ctotaerbox, &
                       qfor, wall_fg)

      ELSE IF (dyn_fg.EQ.1) THEN
        imtr=0      ! 0 mass transfer approach, 
                    ! 1 diffusion-limited, 
                    ! 2 collision-limited
        ratfac=1.   ! simple multiplicative factor gas-> aero rate

        CALL mtrat(ratfac,winfac,weqfac, &
                   numain,numaou,numwin,numwou, &
                   idain,idaou,idwin,idwou, &
                   temp,ctotaerbox,cnv,sumc(ibox), &
                   gamm,Mp,Rpo,maerbox,psat, &
                   wmol,imtr,idrestoi, &
                   Rp,qfor)
      ENDIF

      !PRINT*,"after dynamic soa"

! -------------------------------------
      END SUBROUTINE soa_dyn
! -------------------------------------

