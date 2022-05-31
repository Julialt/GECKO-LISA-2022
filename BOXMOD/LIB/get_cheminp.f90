      SUBROUTINE get_cheminp

      USE sorting
      USE flags_module,ONLY: iofmt_fg,dimer_fg,soa_fg,pvap_fg,wall_fg
      USE io_units_module,ONLY: lout
      USE NetCDF_vars_module,ONLY: ncid_in !,ncid_out
      USE module_data_gecko_main

      IMPLICIT NONE

!----------------------------------------------------------------------
! initializing the package - get the chemical scheme
!----------------------------------------------------------------------
      WRITE(6,*) 'reading the chemical scheme ...'

      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        PRINT*,"ncid_in = ",ncid_in
!        PRINT*,"ncid_out = ",ncid_out

        CALL spakinit9_ncdf(ncid_in, &
            numsp,numre,num_n,num_m,numfo,numhv, &
            numcvar,numextra,numo2,nummeo2,numreacro2, &
            numiso,idiso, &
            mx12stoi,nauxpar,numstoi,idrestoi,idpdstoi, &
            id_m,idfo, &
            idhv,idcvar,idextra,ido2,idmeo2,idreacro2, &
            idreacdimer,numreacdimer, &
            nself,idselfreac, &
            nt1chromo,chromo1cf,id1chromo, &
            ntmedchromo,chromomedcf,nummedchromo,idmedchromo, &
            chromotopcf,numtopchromo,idtopchromo,             &
            restoicf,pdstoicf,arrhcf,focf,                    &
            hvcf,hvfact,cvarcf,extracf,isocf,                 &
            numain, numaou, numwin, numwou,                   &
            idain,idaou,idwin,idwou, &
            !aoucf,woucf,wincf, &
            woucf,wincf, &
            wmol,chrsp)

      ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
        CALL spakinit9_bin(numsp,numre,num_n,num_m,numfo,numhv,       &
                           numcvar,numextra,numo2,nummeo2,numreacro2, &
                           numiso,idiso, mx12stoi,nauxpar,numstoi,    &
                           idrestoi,idpdstoi, id_m,idfo, idhv,idcvar, &
                           idextra,ido2,idmeo2,idreacro2, idreacdimer,&
                           numreacdimer, nself,idselfreac, nt1chromo, &
                           chromo1cf,id1chromo, ntmedchromo,          &
                           chromomedcf,nummedchromo,idmedchromo,      &
                           chromotopcf,numtopchromo,idtopchromo,      &
                           restoicf,pdstoicf,arrhcf,focf, hvcf,hvfact,&
                           cvarcf,extracf,isocf,                      &
                           numain, numaou, numwin, numwou,            &
                           idain,idaou,idwin,idwou,                   &
                           aoucf,woucf,wincf, wmol,chrsp) !,             &
!                           tralphain,tralphaout,                      &
!                           trdeltahin,trdeltahout,                    &
!                           trhenryin,trhenryout,                      & 
!                           numtrin,numtrout,                          &
!                           idtrin,idtrout,                            &
!                           hydcf,numhyd,idhyd,                        &
!                           ka,numacid,idacid,                         &
!                           ohaq,numohaq,idohaq,                       &
!                           idreacro2_aq,nrpero_aq)
      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

      PRINT*,'after init9'

      ! sort chrsp
      WRITE(6,*) 'sorting chrsp'
      call  sort_species(chrsp)
      WRITE(6,*) 'done'

      DO i=1,nt1chromo
        write(6,*) i,chromo1cf(i),id1chromo(i)
      ENDDO
      write(6,*) '------- med-------'
      DO i=1,ntmedchromo
        write(6,*) i,chromomedcf(i),nummedchromo(i)
      ENDDO
      write(6,*) '------- top -------'
!      DO i=1,mtopchromo
!        write(6,*) i,chromotopcf(i),numtopchromo(i)
!      ENDDO

!----------------------------------------------------------------------
! read data for RO2 counting species
      WRITE(6,*) 'reading the RO2 database ...'

      IF (iofmt_fg.NE.1) THEN
        CALL readro2sof2_ncdf(ncid_in,chrsp,numsp, &
                              nclro2,numchemro2,idchemro2,cro2)
      ELSE ! i.e. bindary (iofmt = 1)
        CALL readro2sof2_bin(chrsp,numsp, &
                             nclro2,numchemro2,idchemro2,cro2)
      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

      WRITE(6,*) 'done RO2 read...'

! read data for dimer counting species
      IF (dimer_fg.EQ.1) THEN
        WRITE(6,*) 'reading the dimer database ...'
        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
! still to create ncdf version
         CALL readdimer(chrsp,numsp, &
                ncldimer, numchemdimer,idchemdimer,cdimer)
        ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
          CALL readdimer(chrsp,numsp, &
                ncldimer, numchemdimer,idchemdimer,cdimer)
        ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
      ENDIF
! initialize # of dimers in NetCDF output
!      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
!        CALL switch_ncfile_to_define_mode(ncid_out)
!        CALL eznc_def_0Dint(ncid_out,"ndim")
!        CALL eznc_def_localatt(ncid_out,"ndim","title", &
!                    "actual number of dimers")
!      ENDIF

!----------------------------------------------------------------------
! read temperature-dependent coefficients (if necessary)
! currently no independent ncdf version exists
      IF (numcvar.gt.0) THEN
        WRITE(6,*) 'reading the stoi. coeff. as a function of T ...'
        CALL readcoeff3(lout,numsp,chrsp,numcvar,cvarcf, &
                        valcoe,ntype,numcoe,nopc,ndatopc,nposopc)
      ENDIF
!----------------------------------------------------------------------
! If particles/walls invoked:
      IF (soa_fg.NE.0.OR.wall_fg.NE.0) THEN

! Read data for Psat evaluation 
        IF (pvap_fg.EQ.1) THEN
          WRITE(6,*) 'reading Myrdal data ...'
          CALL readpvap(chrsp,numsp,nsat,namsat,Tb,HBN,tau,idsat)
        ELSE IF (pvap_fg.EQ.2) THEN
          WRITE(6,*) 'reading Nannoolal data ...'
          IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
            CALL readpvap_nan_ncdf(ncid_in,chrsp, &
                                   numsp,nsat,namsat,  &
                                   Tb,dB,idsat,satid)
          ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
            CALL readpvap_nan_bin(chrsp,        &
                                  numsp,nsat,namsat, &
                                  Tb,dB,idsat)
          ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        ELSE IF (pvap_fg.EQ.3) THEN
          WRITE(6,*) 'reading Simpol data ...'
          IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
            CALL readpvap_sim_ncdf(ncid_in,          &
                                   numsp,nsat,namsat,      &
                                   bk,simpgroup,idsat,satid)
          ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
            CALL readpvap_sim_bin(chrsp,             &
                              numsp,nsat,namsat,      &
                              bk,simpgroup,idsat)
          ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        ENDIF ! (pvap_fg.EQ._)

! Read diffusion volume data
        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
          CALL read_dvsp_ncdf(ncid_in)
        ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
          CALL read_dvsp_bin
        ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

      ENDIF ! (soa_fg.NE.0)

!----------------------------------------------------------------------
! read id for gas, particles, wall and dimers for dynamical mass
! transport
! FOR FUTURE: GET IDASAT,IDWSAT,IDDSAT,SATID INTO LINK FILE
      !IF (soa_fg.eq.1) THEN
        ! need code here to USE AKSPNUM to find idsat. 
        ! satid is not required for soa_fg=1
      !ELSE
      IF (soa_fg.eq.2) THEN
        WRITE(6,*) 'reading idgaw info ...'
        CALL read_idgaw(chrsp,numsp,ndim,nsat,namsat, &
                        idgsat,idasat,idwsat,iddsat,satid)
      ENDIF ! (soa_fg.EQ.2)

      WRITE(6,*) 'end of mechanism data read'

!----------------------------------------------------------------------
! create index map from species to reactions
      CALL map_indices

! find species i.d. for common reactants
      CALL map_sp_ids

!==================================
      END SUBROUTINE get_cheminp
!==================================
