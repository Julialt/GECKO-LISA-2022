      SUBROUTINE wall_Huang2018
        use module_chamber_tools, only: A_V_ratio, vol
        USE forcing_params_module,ONLY: temp
        USE module_data_gecko_main,ONLY: numwin,numwou,idwin,idwou, &
                                       idrestoi,wmol,psat,arrhcf,qfor, &
                                       idpdstoi
        USE fundamental_consts_module,ONLY: Ratm, PI

        implicit none
        real, parameter :: R = 8.314 ! [kg m2 s-2 K-1 mol-1] gas constant
        real, parameter :: Mwall = 200. ! [g mol-1] average molecular weight of teflon wall
        real, parameter :: Dg = 5e-6 ! [m2 s-1] gas diffusivity

        real :: Cw ! [mg m-3] equivalent wall concentraiton
        real :: gamma_inf ! [] activity coefficient in teflon
        real :: Vl ! [m s-1] wall deposition velocity
        real :: aw ! [] wall accomodation coefficient
        real :: cstar ! [ug m-3] vapor saturation concentration
        real :: ke ! [s-1] eddy diffusivity coefficient
        real :: w ! [m s-1] mean molecular velocity


        real :: kf ! [s-1] forward rate
        real :: kb ! [s-1] backward rate

        integer :: ispec, ire, i

! calculate values needed for all
        ke = 0.004 + 10**(-2.25)*vol**(0.74)
        Cw = 10.8*A_V_ratio
! calculate gas -> wall
       do i=1,numwin
         ire = idwin(i)
         ispec = idrestoi(ire, 1)
         cstar = psat(ispec) & 
              * wmol(ispec) &
              /  (Ratm*temp) * 1.E+9
         aw = 10**(-2.744)*cstar**(-0.6566)
         Vl = 1 / (PI / (2 * sqrt(ke*Dg)) + 4 / (aw * Dg))
         qfor(ire) = A_V_ratio * Vl
       enddo

! wall -> gas
      do i = 1, numwou
        ire = idwou(i)
        ispec = idpdstoi(ire, 1)
        cstar = psat(ispec) &
            * wmol(ispec) &
            / (Ratm*temp) * 1.E+9
        gamma_inf = 10**3.299 * cstar ** (-0/6407)
        aw = 10**(-2.744)*cstar**(-0.6566)
        Vl = 1 / (PI / 2 * sqrt(ke*Dg)) + 4 / (aw * Dg)
        qfor(ire) = A_V_ratio*Vl * gamma_inf * cstar * Mwall / &
                    (1e3 * Cw * wmol(ispec))
      enddo
      END SUBROUTINE

