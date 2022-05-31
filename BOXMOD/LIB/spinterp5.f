************************************************************************
* Cette routine a pour objet le calcul des constantes photolytiques    *
* en fonction de l'heure. Les constantes sont calculees par            *
* interpolation des constantes tabulees.                               *
* INPUT :                                                              *
*   time        : temps (mudulo 24*3600)                               *
*   sla         : latitude                                             *
*   slo         : longitude                                            *
*   tz          : decalage horaire par rapport au meridien de          *
*                 greenwich en TU                                      *
*   iy          : annee                                                *
*   im          : mois                                                 *
*   id          : jour                                                 *
*   numtet      : nombre de d'angle zenithaux tabules                  *
*   xang(j)     : valeur des angles zenithaux tabules                  *
*   ratpho(i,j) : frequence de photolyse de la reaction i              *
*                 a l'angle xang(j)                                    *
*   npos_cf     : numero de la reaction ayant le drapeau /HV/ 1        *
*                 (correspond inperativement a O3 -> 2 OH              *
*   coefpho(i,*): valeur des coef. de l'interpolation des              *
*                 frequences de photolyse pour la reaction i           *
*                 (l'interpolation est un polynome de degre 2          *
*                 sur chaque intervalle de frequence de photolyse      *
*                 soit 3*numtet coefficient                            *
* OUTPUT :                                                             *
*    les constantes de photolyse sont directement affectees comme      *
*    le premier parametre de la loi d'arrhenius                        *
************************************************************************
      SUBROUTINE interp5(
     1                   numhv,idhv,hvfact,
     7                   nt1chromo,id1chromo,
     7                   ntmedchromo,nummedchromo,idmedchromo,
     7                   numtopchromo,idtopchromo,
     3                   time,
     4                   sla,slo,tz,iy,im,id,
     5                   szafix,szaval,
     5                   numtet,xang,
     3                   rat1pho,coef1pho,ratmedpho,coefmedpho,
     4                   rattoppho,coeftoppho,
     6                   arrhcf, xd, cbox)
      USE akparameter_module
      USE flags_module,ONLY: iofmt_fg, jall_fg, OFR_fg, lightsonoff_fg
      USE io_units_module,ONLY: ljall
      USE time_mgmt_module,ONLY:daytime_fg,iskip,nskip,itout,tout,tstart
      USE netcdf_vars_module,ONLY: ncid_out
      USE forcing_params_module,ONLY: jfac
      USE module_data_gecko_main,ONLY:chromo1cf,chromomedcf,
     &                    chromotopcf,nchrom,idchrom,small
      USE module_chamber_tools,ONLY: run_chamber_lights

      IMPLICIT NONE

* INPUT
      INTEGER numhv,numtet,idhv(maxhv)
      REAL    hvfact(maxhv)
      REAL    time 
      REAL    sla,slo,tz
      INTEGER iy,im,id
      REAL    xang(maxang),ratpho(mchromo,maxang),coefpho(mchromo,nlo)
c      INTEGER  ntchromo,numchromo(mchromo)
c      INTEGER  idchromo(mchromo,mspchromo)

      INTEGER  numtopchromo(mtopchromo)
      INTEGER  idtopchromo(mtopchromo,msptopchromo)

      INTEGER  nummedchromo(mmedchromo)
      INTEGER  idmedchromo(mmedchromo,mspmedchromo)
      INTEGER  ntmedchromo

      INTEGER  nt1chromo,id1chromo(mchromo)

      REAL     rat1pho(mchromo,maxang),coef1pho(mchromo,nlo)
      REAL     ratmedpho(mmedchromo,maxang),coefmedpho(mmedchromo,nlo)
      REAL     rattoppho(mtopchromo,maxang),coeftoppho(mtopchromo,nlo)

      INTEGER  szafix
      REAL     szaval

      REAL,DIMENSION(maxsp) :: cbox ! = conc(1:numsp,ibox)

* OUTPUT
      REAL    arrhcf(maxre,3)
      REAL, intent(out) :: xd ! computed or fixed solar zenital angle

* LOCAL
      INTEGER iti,imn
      REAL    xmn,ti
      REAL    xc
      REAL    cor  ! correcting factor applied to all J values
      REAL    cck(nlo),yy(maxang),v(5),ratact(maxhv)
      !REAL    small
      REAL    lnsmall, k_temp
      INTEGER i,j,k,ire
      REAL    valj(nchrom),jchrom(nchrom)

!=========================================
! flag activates ASCII output of instantaneous j-values by chromophore
      IF(jall_fg.EQ.1)THEN

! write jids once
!! NOW DONE IN open_op_files.f90 and setup_ncdf_op.f90
        IF(time.EQ.tstart)THEN  ! itout = 1 is at time 0
          IF(iofmt_fg.NE.0) THEN
            WRITE(ljall,*)"0",(idchrom(i),i=1,nchrom)
          ENDIF ! iofmt (binary)
          IF(iofmt_fg.NE.1) THEN
            CALL eznc_put_1Dint(ncid_out,"idchrom",
     &                                    idchrom(1:nchrom),
     &                                            1,nchrom)
          ENDIF ! iofmt (NetCDF)
        ENDIF ! iskip = 0 
      ENDIF ! jall_fg
!=========================================
* set minimum value
      !small=1E-32
      lnsmall=log(small)

! set scaling factor (dangerous hardwire)
c      cor=0.252  ! scaling factor for blacklight in Camergie mellon chamber
       cor=1.0  ! scaling factor for no corrections

* change time format for subroutine solar
      iti=int(time/3600.)
      xmn=time/60.
      imn=int(xmn-(iti*60))
      ti = iti*100.0 +imn

* compute the zenithal angle
      xc = 0.0
      CALL solar(sla,slo,tz,iy,im,id,ti,xc,5)
      xd = 90.-xc
** JMLT EDIT FOR FIXED SZA **
      IF(szafix.EQ.1) xd = szaval

      IF (lightsonoff_fg .gt. 0) THEN
        CALL run_chamber_lights(xd) 
      ENDIF
** END EDIT **
      v(1) = xd

* check if it is night or day (no photolysis at night)
      IF (xd.GT.90.0) THEN
        DO i=1,numhv
          ire=idhv(i)
          arrhcf(ire,1)=lnsmall
          arrhcf(ire,2)=0.
          arrhcf(ire,3)=0.
        ENDDO
        daytime_fg = .FALSE.
        RETURN
      ENDIF
      daytime_fg = .TRUE.

* -----------------------------------
* loop over the photolytic reactions
* -----------------------------------

!---------------------------------------------
**** LOOP OVER THE MOST USED LABEL (TOP TABLES)
        !print*,"1st"
      DO i=1,mtopchromo
        k=i

* set the polynomial coef. in J value for reaction i
        DO j=1,nlo
          cck(j)=coeftoppho(i,j)
        ENDDO
        DO j=1,numtet
          yy(j)=rattoppho(i,j)
        ENDDO

* compute J value at the given zenithal angle (v(1))
        CALL splnb(numtet,xang,yy,cck,V)
        !print*,v
        valj(k) = ABS(v(2))
        DO j=1,numtopchromo(i)
          ratact(idtopchromo(i,j)) = valj(k)
        ENDDO
      ENDDO

!---------------------------------------------
***** LOOP OVER THE REGULARLY USED LABEL (MED TABLES)
        !print*,"2nd"
      DO i=1,ntmedchromo
        k=mtopchromo+i

* set the polynomial coef. in J value for reaction i
        DO j=1,nlo
          cck(j)=coefmedpho(i,j)
        ENDDO
        DO j=1,numtet
          yy(j)=ratmedpho(i,j)
        ENDDO

* compute J value at the given zenithal angle (v(1))
        CALL splnb(numtet,xang,yy,cck,V)
        !print*,v
        valj(k) = ABS(v(2))
        DO j=1,nummedchromo(i)
          ratact(idmedchromo(i,j)) = valj(k)
        ENDDO
      ENDDO

!---------------------------------------------
**** LOOP OVER USED ONLY ONCE LABEL (1 TABLES)
        !print*,"3rd"
      DO i=1,nt1chromo
        k=mtopchromo+ntmedchromo+i

* set the polynomial coef. in J value for reaction i
        DO j=1,nlo
          cck(j)=coef1pho(i,j)
        ENDDO
        DO j=1,numtet
          yy(j)=rat1pho(i,j)
        ENDDO

* compute J value at the given zenithal angle (v(1))
        CALL splnb(numtet,xang,yy,cck,V)
        !print*,v
        valj(k) = ABS(v(2))
        ratact(id1chromo(i)) = valj(k)
      ENDDO

!---------------------------------------------------
! revise photolysis rates if in OFR mode
! affects valj and ratact

      IF (OFR_fg.EQ.1) CALL calc_phot_OFR(cbox,valj,ratact)
        
* ---------------------------------
* set J value to the arrhenius parameter of the reaction
* ---------------------------------
!debug
!      rewind(94)
!      write(94,*)  time, numhv, cor
!      DO i=1,numhv
!        write(94,*) i, idhv(i), ratact(i), hvfact(i)
!      ENDDO
!end debug      
      
      DO i=1,numhv

        ire=idhv(i)
        k_temp = ratact(i)*hvfact(i)*cor*jfac

! DEBUG
!        IF(ratact(i).GT.1e-20) PRINT*,i, k_temp

        if (k_temp .gt. 0.0) then
          !k_temp = k_temp * jfac
          arrhcf(ire,1)=log(k_temp)
        else
          arrhcf(ire,1)=lnsmall
        endif
        arrhcf(ire,2)=0.
        arrhcf(ire,3)=0.
      ENDDO

* ---------------------------------
* set J value to the arrhenius parameter of the reaction
* ---------------------------------
!DEBUG
!      rewind(94)
!      write(94,*)  time, numhv, cor
!      DO i=1,numhv
!        write(94,*) i, idhv(i), ratact(i), hvfact(i)
!      ENDDO
!END DEBUG      
      
      DO i=1,numhv

        ire=idhv(i)
        k_temp = ratact(i)*hvfact(i)*cor*jfac

! DEBUG
!        IF(ratact(i).GT.1e-20) PRINT*,i, k_temp
!END DEBUG      

        if (k_temp .gt. 0.0) then
          arrhcf(ire,1)=log(k_temp)
        else
          arrhcf(ire,1)=lnsmall
        endif
        arrhcf(ire,2)=0.
        arrhcf(ire,3)=0.
      ENDDO

!=========================================
! write j-values each output time
      IF(jall_fg.EQ.1)THEN

        IF(iskip.EQ.nskip)THEN

          jchrom = valj*cor*jfac

          IF(iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN
            WRITE(ljall,*)tout,(jchrom(i),i=1,nchrom)
          ENDIF ! iofmt (binary)

          IF(iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

           CALL eznc_put_1Dreal_into2D(ncid_out,"jchrom",
     &                                           jchrom(1:nchrom), 
     &                                                  1,nchrom, itout)
          ENDIF ! iofmt (NetCDF)

        ENDIF ! iskip = 0 or nskip

      ENDIF ! jall_fg
!=========================================

      END
