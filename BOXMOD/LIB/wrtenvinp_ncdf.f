      SUBROUTINE wrtenvinp_ncdf(ncid,numsp,
     2                    tstart,tstop,ntstep,ntprint,
     3                    nbox,rtol,atol,dtmin,cbg,cnv,
     3                    gamm, Mp, Rpo,
     3                    emi_spec, cons_spec,
     3                    surface_emi,
     3                    njf,jftim,jfval,szafix,szaval,
     3                    jo2,jh2o,f185,f254,
     4                    nhd,htop,boxt,boxh,nmx,mixt,mixv,vs,
     4                    ndil,diltim,dilval,
     5                    tempm,tempa,temptm,ntk,tktim,tkval,
     6                    rhm,rha,rhtm,nrh,rhtim,rhval,npr,prtim,prval,
     6                    waterfix,water,sumcfix,sumc,noxfix,
     7                    windm,winda,windtm,nws,wstim,wsval,
     8                    saero,nseed,tseed,cseed,nsd,surft,psurf,
     9                   sla,slo,tz,iy,im,id,iscape,iseas,nskip,maskval)

!==================================================================
! PURPOSE: write environmental input variables 
!          into NetCDF-format box model output file
!          Variables already defined in defenvinp_ncdf.f
!          Refer to spreadkey6 for origin of most parameters.
! AUTHOR: Julia Lee-Taylor, NCAR, 123Jan 2018
!==================================================================

      USE akparameter_module
      USE forcing_params_module,ONLY: prconst,presfix
      USE module_data_gecko_main,ONLY: chrsp,conc,isop_fac,mterp_fac
      USE steadystate_vars_module

      INTEGER ncid,numsp
      REAL    maskval
* species background concentrations
      REAL     cbg(maxsp)
! array of background concs, congruent with REAC
      REAL,ALLOCATABLE :: initcbg(:)
* constrained concentrations
      TYPE(species_data), intent(out) :: cons_spec(mbox, maxconst)
* species emissions
      TYPE(species_data) :: emi_spec(maxem)
* surface emissions
      TYPE(surface_data) :: surface_emi(msur)
* solver parameters
      REAL     rtol, atol, dtmin
* parameters for the simulation (starting time, date, coordinate, ...)
      INTEGER  nbox,ntstep,ntprint
      REAL     tstart,tstop,tlen
      REAL     sla,slo,tz
      INTEGER  iy,im,id
* forcing parameters (mixing height, temperature, humidity, ...)
      INTEGER  nhd,nmx,ndil,ntk,nrh,nws,npr
      REAL     htop
      REAL     boxt(mhd),boxh(mhd)
      REAL     mixt(mhd),mixv(mbox,mhd)
      REAL     vs
      REAL     diltim(mtim),dilval(mtim)
      REAL     tempm(mbox),tempa(mbox),temptm(mbox)
      REAL     tktim(mtim),tkval(mbox,mtim)
      REAL     rhm(mbox),rha(mbox),rhtm(mbox)
      REAL     rhtim(mtim),rhval(mbox,mtim)
      REAL     prtim(mtim),prval(mbox,mtim)
      REAL     windm,winda,windtm
      REAL     wstim(mtim),wsval(mtim)
      REAL,DIMENSION(mbox):: water,sumc
      INTEGER  waterfix,sumcfix, noxfix(mbox)
      INTEGER,DIMENSION(non_zero_species_number) :: idinit
* photolysis parameters for OFR simulations
      REAL     jo2,jh2o,f185,f254
* parameter of the ground surface
      INTEGER  nsd
      REAL     psurf(mhd,msur),surft(mhd)
* parameter for ground deposition
      INTEGER  iseas
* thermodynamic
      INTEGER  iscape,nskip
* aerosol surface aera
      REAL     saero
* concentration of nonvolatile cseed aerosol
      INTEGER  nseed
      REAL     cseed(mtim),tseed(mtim)
      REAL    ::  cnv, gamm, Mp, Rpo
* photolysis adjustment factor
      INTEGER  njf
      REAL     jftim(mtim),jfval(mtim)
* flag and value for fixed SZA option
      INTEGER  szafix
      REAL     szaval
* local
      INTEGER i,j, ibox, isurf
* local interpretation of TYPE species_data
      INTEGER nemis       ! # emitted spp
      INTEGER ntem ! # emission times for a species
      REAL,DIMENSION(maxinput) :: emtim,emval ! emissions tables for a species

      INTEGER ncons ! # constrained spp 
      INTEGER ntcons
      REAL, DIMENSION(maxinput) :: constim, consval

!----------------------------------------------------------

! KEYFILE PARAMETERS
      CALL eznc_put_0Dreal(ncid,"tstart",tstart)
      CALL eznc_put_0Dreal(ncid,"tstop",tstop)
  
      tlen = (tstop-tstart)/FLOAT(ntstep)
      CALL eznc_put_0Dreal(ncid,"tlen",tlen)

!      CALL eznc_put_0Dint(ncid,"ntstep",ntstep)
!      CALL eznc_put_0Dint(ncid,"ntprint",ntprint)

      CALL eznc_put_0Dint(ncid,"nbox",nbox)
      CALL eznc_put_0Dreal(ncid,"rtol",rtol)
      CALL eznc_put_0Dreal(ncid,"atol",atol)
      CALL eznc_put_0Dreal(ncid,"dtmin",dtmin)
      CALL eznc_put_0Dreal(ncid,"vs",vs)

      CALL eznc_put_1Dreal(ncid,"cbg",cbg,1,numsp)

      CALL eznc_put_0Dreal(ncid,"cnv",cnv)
      CALL eznc_put_0Dreal(ncid,"Mp",Mp)
      CALL eznc_put_0Dreal(ncid,"Rpo",Rpo)
      CALL eznc_put_0Dreal(ncid,"gamm",gamm)

!-----------------------------------------------------------------
! REAC species initializations defined as in printsteadystate
! cbg is already being written congruent with chrsp:
! here it is again, congruent with REAC

      !! now doing this in defenvinp.f90
      !non_zero_species_number = COUNT(conc(:,1) > 1.)
      CALL eznc_put_0Dint(ncid,"ninit",non_zero_species_number)

      ALLOCATE(initspecnames(non_zero_species_number),
     &          initspecconc(non_zero_species_number,nbox),
     &               initcbg(non_zero_species_number))

      j = 0
      DO i=1,maxsp
        IF (conc(i,1) > 1.) THEN
          j = j+1
          idinit(j) = i
          initspecnames(j) = chrsp(i)
          initcbg(j) = cbg(i)

          PRINT*,"REAC ",idinit(j),initspecnames(j)

          DO ibox = 1,nbox
            initspecconc(j,ibox) = conc(i,ibox)
          ENDDO
        ENDIF
      ENDDO

      CALL eznc_put_1Dint(ncid,"idinit",
     &                          idinit(1:non_zero_species_number),
     &                                 1,non_zero_species_number)

      CALL eznc_put_1Dchar(ncid,"initnam",
     &            initspecnames(1:non_zero_species_number),
     &                   maxlsp,1,non_zero_species_number)

      CALL eznc_put_2Dreal(ncid,"initconc",
     &                 initspecconc(1:non_zero_species_number,1:nbox),
     &                              1,non_zero_species_number,1,nbox)

      CALL eznc_put_1Dreal(ncid,"initcbg",
     &                           initcbg(1:non_zero_species_number),
     &                                   1,non_zero_species_number)

! Deallocations
      IF(ALLOCATED(initspecnames)) DEALLOCATE(initspecnames)
      IF(ALLOCATED(initspecconc))  DEALLOCATE(initspecconc)
      IF(ALLOCATED(initspecconc))  DEALLOCATE(initcbg)

!-----------------------------------------------------------------
!     &              emi_spec, cons_spec,
! nemis and nemit are defined as dimensions (# spp and max # times)
      nemis = COUNT(emi_spec(:)%activefg)
      CALL eznc_put_0Dint(ncid,"nemis",nemis)

      IF(nemis.GT.0)THEN
        CALL eznc_put_1Dint(ncid,"idemis",pack(emi_spec%index,
     &                    mask = emi_spec%activefg),1,nemis)
        CALL eznc_put_1Dchar(ncid,"eminam",pack(emi_spec(:)%name,
     &                    mask = emi_spec%activefg),maxlsp,1,nemis)
        CALL eznc_put_1Dint(ncid, "ntem", emi_spec(:)%npoints, 1, nemis)
        do i = 1, maxem
          if (emi_spec(i)%activefg) then
            ntem = emi_spec(i)%npoints
            emtim(1:ntem) = emi_spec(i)%table(1:ntem, 1)
            emval(1:ntem) = emi_spec(i)%table(1:ntem, 2)
            CALL eznc_put_2Dreal(ncid, "emtim", emtim(1:ntem),
     &                           1, ntem, i, i)
            CALL eznc_put_2Dreal(ncid, "emval", emval(1:ntem),
     &                           1, ntem, i, i)
     
          endif        
        enddo
      ENDIF !(nemis.GT.0)THEN

!-----------------------------------------------------------------
!ncons is defined as a variable (see above note for nemis)
      
      do ibox=1, nbox
        ncons = COUNT(cons_spec(ibox,:)%activefg)
        CALL eznc_put_0Dint_into1D(ncid,"ncons",ncons,ibox)
        
        if (ncons .gt. 0) then
          
          CALL eznc_put_2Dint(ncid,"idcons",
     &         pack(cons_spec(ibox,:)%index,
     &               mask = cons_spec(ibox,:)%activefg),
     &         1, ncons, ibox, ibox)
          CALL eznc_put_2Dchar(ncid,"consnam",
     &         pack(cons_spec(ibox,:)%name,
     &              mask = cons_spec(ibox,:)%activefg),
     &         maxlsp, 1, ncons, ibox, ibox)
          do i = 1, maxconst
            ntcons = cons_spec(ibox,i)%npoints
            constim(1:ntcons) = cons_spec(ibox, i)%table(1:ntcons, 1)
            consval(1:ntcons) = cons_spec(ibox, i)%table(1:ntcons, 2)
            CALL eznc_put_0Dint_into2D(ncid,"ntcons",ntcons,i,ibox)
            CALL eznc_put_3Dreal(ncid,"constim",constim(1:ntcons),
     &                           1,ntcons,i,i,ibox,ibox)
            CALL eznc_put_3Dreal(ncid,"consval",consval(1:ntcons),
     &                           1,ntcons,i,i,ibox,ibox)
          enddo
        endif
      enddo
      
!-----------------------------------------------------------------
! emissions linked to surface nature
      do isurf = 1, msur
        nemis = surface_emi(isurf)%nemis
        CALL eznc_put_0Dint_into1D(ncid,"surf_nemis",
     &                                  surface_emi(isurf)%nemis,
     &                                  isurf)
        
        if (nemis .gt. 0) then
        
          CALL eznc_put_2Dint(ncid,"surf_idemis",
     &         pack(surface_emi(isurf)%emission(:)%index,
     &              mask = surface_emi(isurf)%emission(:)%activefg),
     &         1,nemis,isurf,isurf)
          CALL eznc_put_2Dchar(ncid,"surf_eminam",
     &         pack(surface_emi(isurf)%emission(:)%name,
     &              mask = surface_emi(isurf)%emission(:)%activefg),
     &         maxlsp,1,nemis,isurf,isurf)
          do i = 1, maxem
            if (surface_emi(isurf)%emission(i)%activefg) then
              ntem = surface_emi(isurf)%emission(i)%npoints
              emtim(1:ntem) = surface_emi(isurf)%emission(i)%
     &                                       table(1:ntem,1)
              emval(1:ntem) = surface_emi(isurf)%emission(i)%
     &                                       table(1:ntem,2)
              CALL eznc_put_0Dint_into2D(ncid,"surf_ntem",
     &                                   ntem,i,isurf)
              CALL eznc_put_3Dreal(ncid,"surf_emtim",emtim(1:ntem),
     &                             1,ntem,i,i,isurf,isurf)
              CALL eznc_put_3Dreal(ncid,"surf_emval",emval(1:ntem),
     &                             1,ntem,i,i,isurf,isurf)
            endif
          enddo
        endif      
      enddo


!-----------------------------------------------------------------
!     &              njf,jftim,jfval,szafix,szaval,
      CALL eznc_put_0Dint(ncid,"njf",njf)
      IF (njf.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"jftim",
     &                                   jftim(1:njf),1,njf)
        CALL eznc_put_1Dreal(ncid,"jfval",
     &                                   jfval(1:njf),1,njf)
      ENDIF
      CALL eznc_put_0Dint(ncid,"szafix",szafix)
      IF (szafix.EQ.1) THEN
        CALL eznc_put_0Dreal(ncid,"szaval",szaval)
      ENDIF
!-----------------------------------------------------------------
!     &              jo2,jh2o,f185,f254,
      IF (jo2.GT.0) CALL eznc_put_0Dreal(ncid,"jo2",jo2)
      IF (jh2o.GT.0) CALL eznc_put_0Dreal(ncid,"jh2o",jh2o)
      IF (f185.GT.0) CALL eznc_put_0Dreal(ncid,"f185",f185)
      IF (f254.GT.0) CALL eznc_put_0Dreal(ncid,"f254",f254)


!-----------------------------------------------------------------
!     &              nhd,htop,boxt,boxh,nmx,mixt,mixv,
      CALL eznc_put_0Dreal(ncid,"htop",htop)

      CALL eznc_put_0Dint(ncid,"nhd",nhd)
      IF (nhd.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"boxt",
     &                             boxt(1:nhd),1,nhd)
        CALL eznc_put_1Dreal(ncid,"boxh",
     &                             boxh(1:nhd),1,nhd)
      ENDIF

! table
      CALL eznc_put_0Dint(ncid,"nmx",nmx)
      IF (nmx.NE.0) THEN
        CALL eznc_put_1Dreal(ncid,"mixt",
     &                             mixt(1:nmx),1,nmx)
        CALL eznc_put_2Dreal(ncid,"mixv",
     &                             mixv(1:nbox,1:nmx),1,nbox,1,nmx)
      ENDIF

      CALL eznc_put_0Dint(ncid,"ndil",ndil)
      IF (ndil.NE.0) THEN
        CALL eznc_put_1Dreal(ncid,"diltim",
     &                             diltim(1:ndil),1,ndil)
        CALL eznc_put_1Dreal(ncid,"dilval",
     &                             dilval(1:ndil),1,ndil)
      ENDIF

!-----------------------------------------------------------------
!     5              tempm,tempa,temptm,ntk,tktim,tkval,
! sine
      IF (tempm(1).GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"tempm",
     &                                   tempm(1:nbox),1,nbox)
        CALL eznc_put_1Dreal(ncid,"tempa",
     &                                   tempa(1:nbox),1,nbox)
        CALL eznc_put_1Dreal(ncid,"temptm",
     &                                   temptm(1:nbox),1,nbox)

!        CALL eznc_put_1Dreal_into2D(ncid,"TEMP",
!     &                                   tempm(1:nbox),1,1,1,nbox)
!        CALL eznc_put_1Dreal_into2D(ncid,"TEMP",
!     &                                   tempa(1:nbox),2,2,1,nbox)
!        CALL eznc_put_1Dreal_into2D(ncid,"TEMP",
!     &                                   temptm(1:nbox),3,3,1,nbox)
      ENDIF

! table
      CALL eznc_put_0Dint(ncid,"ntk",ntk)
      IF (ntk.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"tktim",
     &                                   tktim(1:ntk),1,ntk)
        CALL eznc_put_2Dreal(ncid,"tkval",
     &                                   tkval(1:nbox,1:ntk),
     &                                         1,nbox,1,ntk)
      ENDIF

!-----------------------------------------------------------------
!     6              rhm,rha,rhtm,nrh,rhtim,rhval,
! sine
      IF (rhm(1).GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"rhm",rhm(1:nbox),1,nbox)
        CALL eznc_put_1Dreal(ncid,"rha",rha(1:nbox),1,nbox)
        CALL eznc_put_1Dreal(ncid,"rhtm",rhtm(1:nbox),1,nbox)

!        CALL eznc_put_1Dreal_into2D(ncid,"RH",rhm(1:nbox),1,1,1,nbox)
!        CALL eznc_put_1Dreal_into2D(ncid,"RH",rha(1:nbox),2,2,1,nbox)
!        CALL eznc_put_1Dreal_into2D(ncid,"RH",rhtm(1:nbox),3,3,1,nbox)
      ENDIF

! table
      CALL eznc_put_0Dint(ncid,"nrh",nrh)
      IF (nrh.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"rhtim",
     &                                   rhtim(1:nrh),1,nrh)
        CALL eznc_put_2Dreal(ncid,"rhval",
     &                                   rhval(1:nbox,1:nrh),
     &                                         1,nbox,1,nrh)
      ENDIF

      CALL eznc_put_0Dint(ncid,"npr",npr)
      IF (npr.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"prtim",
     &                                   prtim(1:npr),1,npr)
        CALL eznc_put_2Dreal(ncid,"prval",
     &                                   prval(1:nbox,1:npr),
     &                                         1,nbox,1,npr)
      ENDIF


!-----------------------------------------------------------------
!     6              waterfix,water,sumcfix,sumc,noxfix,
      CALL eznc_put_0Dint(ncid,"waterfix",waterfix)
      IF (waterfix.EQ.1) THEN
        CALL eznc_put_1Dreal(ncid,"water",water(1:nbox),1,nbox)
      ENDIF
      CALL eznc_put_0Dint(ncid,"sumcfix",sumcfix)
      IF (sumcfix.EQ.1) THEN
        CALL eznc_put_1Dreal(ncid,"sumc",sumc(1:nbox),1,nbox)
      ENDIF
      CALL eznc_put_0Dint(ncid,"presfix",presfix)
      IF (presfix.EQ.1) THEN
        CALL eznc_put_1Dreal(ncid,"prconst",prconst(1:nbox),1,nbox)
      ENDIF
      CALL eznc_put_1Dint(ncid,"noxfix",noxfix(1:nbox),1,nbox)

!-----------------------------------------------------------------
! sine
      IF (windm.GT.0) THEN
        CALL eznc_put_0Dreal(ncid,"windm",windm)
        CALL eznc_put_0Dreal(ncid,"winda",winda)
        CALL eznc_put_0Dreal(ncid,"windtm",windtm)
      ENDIF

! table
      CALL eznc_put_0Dint(ncid,"nws",nws)
      IF (nws.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"wstim", wstim(1:nws),1,nws)
        CALL eznc_put_1Dreal(ncid,"wsval", wsval(1:nws),1,nws)
      ENDIF

!-----------------------------------------------------------------
      CALL eznc_put_0Dint(ncid,"nseed",nseed)
      IF (nseed.GT.0) THEN
        CALL eznc_put_1Dreal(ncid,"tseed",tseed(1:nseed),1,nseed)
        CALL eznc_put_1Dreal(ncid,"cseed",cseed(1:nseed),1,nseed)
      ENDIF

      CALL eznc_put_0Dint(ncid,"nsd",nsd)
      IF (nsd.GT.0) THEN

        CALL eznc_put_1Dreal(ncid,"surft",surft(1:nsd),1,nsd)
        CALL eznc_put_2Dreal(ncid,"psurf",psurf(1:nsd,1:msur),
     &                                          1,nsd,1,msur)
      ENDIF

      CALL eznc_put_0Dreal(ncid,"isop_fac",isop_fac)
      CALL eznc_put_0Dreal(ncid,"mterp_fac",mterp_fac)

!-----------------------------------------------------------------
      IF(sla.GT.0)THEN
      CALL eznc_put_0Dreal(ncid,"sla",sla)
      CALL eznc_put_0Dreal(ncid,"slo",slo)
      CALL eznc_put_0Dreal(ncid,"tz",tz)
      CALL eznc_put_0Dint(ncid,"iy",iy)
      CALL eznc_put_0Dint(ncid,"im",im)
      CALL eznc_put_0Dint(ncid,"id",id)
      ENDIF

!-----------------------------------------------------------------
      CALL eznc_put_0Dint(ncid,"iscape",iscape)
      CALL eznc_put_0Dint(ncid,"iseas",iseas)
! nskip as input by user is ONE LESS than nskip as used by the model
      CALL eznc_put_0Dint(ncid,"nskip",nskip-1)

!==end write keyfile inputs

!-----------------------------------------------------------------
!==WRITE MASK THRESHOLD FOR "SMALL" OUTPUT VALUES
      CALL eznc_put_0Dreal(ncid,"maskval",maskval)

      END SUBROUTINE wrtenvinp_ncdf
