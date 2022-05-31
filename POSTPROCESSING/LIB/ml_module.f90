module machine_learning
! this module produces output needed by machine learning
use conc
use sorting, only: find_species_index
use io
use bins
contains
subroutine output_environ()

    real, dimension(:,:), allocatable      :: environ_params
    integer, parameter                     :: ind_time = 1, ind_rh = 2, ind_temp = 3, ind_pres = 4, &
    ind_sza = 5, ind_ho = 6, ind_ho2 = 7, ind_no = 8, ind_no2 = 9, ind_coa = 10, ind_o3 = 11, ind_co = 12, &
    ind_ch4 = 13
    character*(header_length)              :: header(ind_ch4)
    integer                                :: itime,iho,iho2,ino,ino2,io3,ico,ich4

    CHARACTER*(10)    :: date,time1
    CALL date_and_time(date,time1)

    allocate(environ_params(ndat,ind_ch4))

    header(ind_time) = "Time [s]"
    header(ind_rh)   = "RH [%]"
    header(ind_temp) = "Temperature [K]"
    header(ind_pres) = "Pressure [hPa]"
    header(ind_sza)  = "Zenith ANgle [deg]"
    header(ind_ho)   = "OH [molec cm-3]"
    header(ind_ho2)  = "HO2 [molec cm-3]"
    header(ind_no)   = "NO [molec cm-3]"
    header(ind_no2)  = "NO2 [molec cm-3]"
    header(ind_coa)  = "COA [molec cm-3]"
    header(ind_o3)   = "O3 [molec cm-3]"
    header(ind_co)   = "CO [molec cm-3]"
    header(ind_ch4)  = "CH4 [molec cm-3]"

    iho  = find_species_index("GHO")
    iho2 = find_species_index("GHO2")
    ino  = find_species_index("GNO")
    ino2 = find_species_index("GNO2")
    io3  = find_species_index("GO3")
    ico  = find_species_index("GCO")
    ich4 = find_species_index("GCH4") 

    do itime = 1, ndat
      environ_params(itime, ind_time) = time(itime)
      environ_params(itime, ind_rh)   = rh(itime)
      environ_params(itime, ind_temp) = temperature(itime)
      environ_params(itime, ind_pres) = pressure(itime)
      environ_params(itime, ind_sza)  = sza(itime)
      environ_params(itime, ind_coa)  = seed(itime)
      environ_params(itime, ind_ho)   = concentrations(itime,iho)
      environ_params(itime, ind_ho2)  = concentrations(itime,iho2)
      environ_params(itime, ind_no)   = concentrations(itime,ino)
      environ_params(itime, ind_no2)  = concentrations(itime,ino2)
      environ_params(itime, ind_o3)   = concentrations(itime,io3)
      environ_params(itime, ind_co)   = concentrations(itime,ico)
      environ_params(itime, ind_ch4)  = concentrations(itime,ich4)
    enddo

    filename = "ml_environ_params"
    call write_2d_array(environ_params, header,filename) 

    deallocate(environ_params)
end subroutine output_environ


subroutine output_bins()
  real, dimension(:,:),  allocatable :: output_array
  character*(header_length), dimension(:), allocatable :: header
  integer, parameter :: ind_time=1,ind_lowerlim=2,ind_upperlim=3,ind_aer=4,ind_gas_aer=5, &
    ind_molw=6,ind_precu=7,ind_heff=8,ind_koh=9,ind_oc=10,ind_hc=11,ind_nc=12
  logical, dimension(:), allocatable :: spec_mask    
  type(bins_class) :: pvap_bins
  real :: minpvap, maxpvap, molec, mass, av_nc
  real, parameter :: ugfac = 1.660578e-12
  integer :: itime, ibin, ispec, iline, nbins

  minpvap = real(floor(minval(log10(species%pvap%pvap_atm_298),mask = species%pvap%pvap_atm_298 > 0.)), kind = 8)  

  maxpvap = real(ceiling(maxval(log10(species%pvap%pvap_atm_298),mask = species%pvap%pvap_atm_298 > 0.)), kind = 8)  
  call pvap_bins%initialize(minpvap, maxpvap, 1.)
  nbins = pvap_bins%nbins

  print *, nbins, minpvap, maxpvap
  allocate(output_array(ndat*nbins, ind_nc), header(ind_nc), spec_mask(numsp))

  header(ind_time) = "Time [s]"
  header(ind_lowerlim) = "Pvap bin lower limit [log10(atm)]"
  header(ind_upperlim) = "Pvap bin upper limit [log10(atm)]"
  header(ind_aer) = "Total aerosol mass [ug m-3]"
  header(ind_gas_aer) = "Total gas + aerosol mass [ug m-3]"
  header(ind_molw) = "Aerosol Molar Weight [g mol-1]"
  header(ind_precu) = "Precursors mass [ug m-3]"
  header(ind_heff) = "Heff [M atm-1]"
  header(ind_koh) = "koh [cm3 molec-1 s-1]"
  header(ind_oc) = "O/C"
  header(ind_hc) = "H/C"
  header(ind_nc) = "N/C"
  
  output_array = 0.
  iline = 0
  do itime =1,ndat
    do ibin = 1, nbins
      iline = iline + 1
      output_array(iline, ind_time) = time(itime)
      output_array(iline, ind_lowerlim) = pvap_bins%lowerlim(ibin)
      output_array(iline, ind_upperlim) = pvap_bins%upperlim(ibin)
    enddo
  enddo

  iline = 0
  do itime = 1, ndat
    do ispec = 1, numsp
      if (species(ispec)%pvap%pvap_atm_298 .eq. -9999.) cycle
      ibin = pvap_bins%find_index(log10(species(ispec)%pvap%pvap_atm_298))
      iline = nbins*(itime-1) + ibin

      molec = concentrations(itime,ispec)
      mass = molec*species(ispec)%molw*ugfac

      if (species(ispec)%phase== AEROSOL) then
        output_array(iline, ind_aer) = output_array(iline, ind_aer) + mass
        output_array(iline, ind_gas_aer) = output_array(iline, ind_gas_aer) + mass
      elseif (species(ispec)%phase== GAS) then
        output_array(iline, ind_gas_aer) = output_array(iline, ind_gas_aer) + mass
      elseif (species(ispec)%phase== PRECURSOR) then
        output_array(iline, ind_precu) = output_array(iline, ind_precu) + mass
      endif
    enddo
  enddo

  do itime = 1, ndat
    do ibin = 1, nbins
      iline = nbins*(itime-1) + ibin
      spec_mask = .false.
      where (species%pvap%pvap_atm_298 .ne. -9999.)
        spec_mask = species%phase == AEROSOL .and. &
                  log10(species%pvap%pvap_atm_298) > pvap_bins%lowerlim(ibin) .and. &
                  log10(species%pvap%pvap_atm_298) <= pvap_bins%upperlim(ibin)
      end where
      output_array(iline, ind_molw) = sum(concentrations(itime,:)*species%molw, &
                   mask = spec_mask) / sum(concentrations(itime,:), mask = spec_mask)
      output_array(iline, ind_heff) = sum(concentrations(itime,:)*species%henry%Keff, &
                   mask = spec_mask.and.species%henry%Keff /= -9999.) &
                 / sum(concentrations(itime,:), &
                   mask = spec_mask.and.species%henry%Keff /= -9999.)
      output_array(iline, ind_koh) = sum(concentrations(itime,:)*species%koh%k_298, &
                   mask = spec_mask.and.species%koh%k_298 /= -9999.) &
                 / sum(concentrations(itime,:), &
                   mask = spec_mask.and.species%koh%k_298 /= -9999.)

      av_nc = sum(concentrations(itime,:)*species%nc, &
                   mask = spec_mask)/sum(concentrations(itime,:), mask = spec_mask)
      output_array(iline, ind_oc) = sum(concentrations(itime,:)*species%no, mask = spec_mask) &
                                 / (sum(concentrations(itime,:), mask = spec_mask) * av_nc)
      output_array(iline, ind_hc) = sum(concentrations(itime,:)*species%nh, mask = spec_mask) &
                                 / (sum(concentrations(itime,:), mask = spec_mask) * av_nc)
      output_array(iline, ind_nc) = sum(concentrations(itime,:)*species%nn, mask = spec_mask) &
                                 / (sum(concentrations(itime,:), mask = spec_mask) * av_nc)
    enddo
  enddo

  filename = "ml_bins"
  call write_2d_array(output_array, header, filename)

  deallocate(header, output_array, spec_mask)

end subroutine output_bins


end module machine_learning
