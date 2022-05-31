module inorganic_aer_module
  implicit none
  type aer_population
! constrained
    real :: organic_mass ! ug/m3(air)
    real :: density ! kg/m3(aer)
    real :: pH
    real :: total_sulfates, total_nitrates ! ug/m3(air)
    real :: kappa ! hygroscopicity
    real :: N ! number concentration 1/m3(air)

! calculated
    real :: dry_mass ! ug/m3(air)
    real :: inorg_mass ! ug/m3(air)
    real :: saero ! m2(aer)/m3(air)
    real :: dry_radius ! cm
    real :: wet_radius ! cm
    real :: hp, hso4m, no3m, so4mm ! ions concentrations mol/L(w)
    real :: water_mass ! ug/m3(air)
    real :: vdepaer ! cm/s deposition velocity of aerosol particles

  end type
  real, parameter :: sulfate_pka1 = -3, sulfate_pka2 = 1.99
  real, parameter :: pi = 3.1415926536
  real, parameter :: water_density = 997. ! kg/m3
  real, parameter :: wmol_nitrate = 62., wmol_sulfate = 96. ! in g/mol
  real, parameter :: g = 9.81 ! m/s2
  real, parameter :: n = 1.849e-5 ! air dynamic viscosity at 25C in kg/m/s


  contains

  subroutine init_aer_population(aer, &
                            density, kappa, & ! fixed properties
                            N, pH, total_sulfates, total_nitrates, &! variable constraints
                            organic_mass, &! from model chemistry
                            HR) ! environmental constraint

    type(aer_population), intent(out)  ::  aer
    real, intent(in) :: density, kappa, HR
    real, intent(in) :: N, pH, total_sulfates, total_nitrates
    real, intent(in) :: organic_mass

    aer%density = density
    aer%kappa = kappa
    aer%N = N
    aer%pH = pH
    aer%total_sulfates = total_sulfates
    aer%total_nitrates = total_nitrates
    aer%organic_mass = organic_mass

    ! calculate other properties
    call update_mass(aer)
    call update_water_mass(aer, HR)
    call update_radii(aer)
    call update_saero(aer)
    call update_ions_from_ph(aer)

  end subroutine

  subroutine update_mass(aer)
    type(aer_population), intent(inout)  ::  aer

    aer%inorg_mass = aer%total_sulfates + aer%total_nitrates
    aer%dry_mass = aer%organic_mass + aer%inorg_mass

  end subroutine

  subroutine update_water_mass(aer, HR)
    type(aer_population), intent(inout)  ::  aer
    real, intent(in)  :: HR

    aer%water_mass = aer%dry_mass*aer%kappa*water_density/  &
      (aer%density*(100/HR - 1))

  end subroutine

  subroutine update_radii(aer)
    type(aer_population), intent(inout)  ::  aer

    real :: effective_density
    aer%dry_radius = 3*(aer%dry_mass*1e-9) ! dry mass converted to kg
    aer%dry_radius = (aer%dry_radius/(4*pi*aer%N*aer%density))**(1.0/3.0) ! in m

    aer%wet_radius = 3*(aer%water_mass*1e-9) ! wet mass converted to kg
    effective_density = (aer%water_mass*water_density + aer%dry_mass*aer%density)/ &
                         (aer%water_mass + aer%dry_mass) ! weighted average of aerosol and water density
    aer%wet_radius = (aer%wet_radius/(4*pi*aer%N*effective_density))**(1.0/3.0) ! in m

    if (aer%wet_radius < aer%dry_radius) then
      aer%wet_radius = aer%dry_radius
    endif

    aer%dry_radius = 100.*aer%dry_radius ! convert to cm
    aer%wet_radius = 100.*aer%wet_radius ! convert to cm

  end subroutine


  subroutine update_saero(aer)
    type(aer_population), intent(inout)  ::  aer

    aer%saero = aer%N*4*pi*(aer%wet_radius)*(aer%wet_radius)
  end subroutine

  subroutine update_ions_from_ph(aer)
    type(aer_population), intent(inout)  ::  aer

    real :: sulf_conc, nitrate_conc ! in mol/L
    ! convert ug/m3(air) to mol/L
    nitrate_conc = aer%total_nitrates*water_density/ &
                   (wmol_nitrate*aer%water_mass)
    sulf_conc =  aer%total_sulfates*water_density/ &
                 (wmol_sulfate*aer%water_mass)
    aer%no3m = nitrate_conc


    aer%hp = 10**(-aer%ph)
    aer%hso4m = sulf_conc/(1+10**(sulfate_pka2)/aer%hp)
    aer%so4mm = sulf_conc - aer%hso4m

  end subroutine

  subroutine update_deposition(aer, temp, ResAerdep, pres)
  ! see Zhang et al. (2011), Atmos. Environ. 35 (2001) 549-560
    type(aer_population), intent(inout) :: aer
    real, intent(in) :: temp, ResAerdep, pres ! K, s / m, bar


    real :: Vg ! gravitational settling m/s
    real :: C ! correction factor for small particles
    real :: dp ! particles diameter in m
    real :: lambda ! mean free path of air molecules

    dp = aer%wet_radius * 2. / 100. ! convert radius in cm to diameter in m

    lambda = 1.14e-5 * temp / (pres * 1e5)
    C = 1.0 + 2.0 * lambda * (1.257 + 0.4 * exp(-0.55 * dp / lambda)) / dp

    Vg = aer%density * g * C * dp * dp / (18 * n)

    aer%vdepaer = Vg + 1.0 / (ResAerdep) ! m/s
    ! convert cm/s
    aer%vdepaer = aer%vdepaer * 100


  end subroutine update_deposition
end module inorganic_aer_module
