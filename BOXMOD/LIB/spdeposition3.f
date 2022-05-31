*******************************************************************
* Cette routine a pour objet de calculer les vitesses de depot    *
* sec des especes en phase gazeuse                                *
*                                                                 *
* INPUT :                                                         *
*   ndepspe  : nombre d'espece qui se depose                      *
*   time        : temps (mudulo 24*3600)                          *
*   temp        : temperature                                     *
*   rik(i)   : resistance stomatale minimum pour la vapeur d'eau  *
*              pour la surface (i)                                *
*   rlu0k(i) : resistance des feuilles de surface dans la canope  *
*              superieure ou resistance cuticulaire des feuilles  *
*              pour la surface (i)                                *
*   rack(i)  : resistance de transfert (depend de la hauteur et   *
*              la densite de la canope) pour la surface (i)       *
*   rgsSk(i) : resistance de piegage au sol pour SO2              *
*              pour la surface (i)                                *
*   rgsOk(i) : resistance de piegage au sol pour O3               *
*              pour la surface (i)                                *
*   rclSk(i) : resistance des composants de la partie inferieure  *
*              de la canopee (ecorce, brindilles...) pour SO2     *
*              pour la surface (i)                                *
*   rclOk(i) : resistance des composants de la partie inferieure  *
*              de la canopee (ecorce, brindilles...) pour O3      *
*              pour la surface (i)                                *
*   z0k(i)   : hauteur de rugosite pour la surface (i)            *
*   sla      : latitude                                           *
*   slo      : longitude                                          *
*   tz       : decalage horaire par rapport au meridien de        *
*              greenwich en TU                                    *
*   iy       : annee                                              *
*   im       : mois                                               *
*   id       : jour                                               *
*   depnamspe(i) : nom de la ieme espece qui se depose            *
*   depdatspe(i,j) : parametre j de la ieme espece                *
*                    j=1 => D_H2O/D_x                             *
*                    j=2 => constante de Henry                    *
*                    j=3 => facteur de reactivite                 *
*   iddepspe(i)    : numero d'identification de l'espece dans le  *
*                    mecanisme chimique                           *
*   windm          : vitesse moyenne du vent (m/s)                *
*   winda          : amplitude de variation de la vitesse du vent *
*   windtm         : temps pour lequel la vitesse est maximum     *
*   => nws   : number of tabular wind speeds                      *
*   => wstim : time for tabular wind speeds (s)                   *
*   => wsval : tabular wind speeds (m.s-1)                        *
*                                                                 *
* OUTPUT :                                                        *
*   Vd(i)          : vitesse de depot de l-espece de numero i     *
*   ResAerdep      : sum of resistances to be used to calculate   *
*                    aer dep later
*******************************************************************
      SUBROUTINE deposition3(time, temp, xsurf, Rik, Rlu0k,
     &                      Rack, RgsSk, RgsOk, RclSk, RclOk, z0k,
     &                      ndepspe,depnamspe, depdatspe,iddepspe,
     &                      windm,winda,windtm,nws,wstim,wsval,
     &                      iseas,
     &                      sla,slo,tz,iy,im,id,Vd, ResAerdep)

      USE akparameter_module
      USE module_data_gecko_main,ONLY: small

      IMPLICIT NONE

* INPUT
      CHARACTER(maxlsp) depnamspe(mxdep)
      INTEGER   iseas, nws
      REAL      time, temp,windm,winda,windtm
      REAL      wstim(mtim),wsval(mtim)
      REAL      xsurf(msur)
      INTEGER   ndepspe,iddepspe(mxdep)
      REAL      depdatspe(mxdep,3)
      REAL      Rik(msur),Rlu0k(msur),Rack(msur),RgsSk(msur),
     &          RgsOk(msur),RclSk(msur),RclOk(msur),z0k(msur)
      REAL      sla,slo,tz
      INTEGER   iy,im,id

* OUTPUT
      REAL      Vd(mxdep), ResAerdep

* LOCAL
      INTEGER   i,j,iti,imn
      REAL      pi,cvk,T0,dwat,xmu,gslope,zs,tempfac
      REAL      fconv,b,wind,xmn,ti,xc,zen,w
      REAL      Ri,Rlu0,Rac,RgsS,RgsO,RclS,RclO,z0,rt,rs,rdc
      REAL      Rsm,Rlu,Rcl,Rgs
      REAL      V,Vdtyp(mxdep,msur), yy(mxdep)
      INTEGER   iter,lout
      REAL      ra,rb,rc,rtot
      REAL      ustarcor,ustar,xx,Z_L,HV,phiH,phiM,ZZ,err,DH2O
      REAL      xmin,xmax,ustark
      !REAL      small
      REAL      drdt
      REAL      logzsz0


*********************************
* set some value                *
*********************************
      !small = TINY(1.0)
* pi
      pi=3.141592

* Von Karman constant
      cvk=0.4

* reference temperature (273 K)
      T0=273.

* diffusion coefficient for H2O at 273 K (m2/s)
      Dwat=0.23E-4

* air viscosity (m2/s)
      xmu=1.461E-5

* slope at ground (angle)
      gslope=0.

* reference height (height at which the wind is given) (m)
      zs=10.

* temperature factor
      tempfac=100.

*****************************************
* initialise                            *
*****************************************
      Vd = 0.
      Vdtyp = 0.

******************************************
* compute some time dependent parameters *
******************************************

* interpolate wind speed (from input table) (m/s)
      IF(nws.GT.0)THEN
        IF(wstim(nws).lt.time)THEN
          WRITE (lout,*) '--error--, in deposition3 -wind speed- '
          WRITE (lout,*) 'upper limit for time not found'
          STOP
        ENDIF
        DO i=1,nws-1
          IF (wstim(i+1).GT.time) THEN
            drdt=(wsval(i+1)-wsval(i))/(wstim(i+1)-wstim(i))
            wind=wsval(i) + drdt*(time-wstim(i))
            EXIT
          ENDIF
        ENDDO
      ELSE
* OR compute wind speed (sinusoidal expression) (m/s)
        fconv=7.2722E-5
        b=(pi/2.)-fconv*windtm
        wind=windm + winda*SIN((fconv*time)+b)
      ENDIF

* compute zenith angle (zen) and solar rayonnement (W/m2)
      iti=int(time/3600.)
      xmn=time/60.
      imn=int(xmn-(iti*60))
      ti = iti*100.0 +imn
      xc = 0.0
      CALL solar(sla,slo,tz,iy,im,id,ti,xc,5)
      zen = 90.-xc
      W=1370.*(0.25+0.50*0.6)*cos(pi*zen/180.)
      IF (zen.GE.90.) W=0.

****************************************
* loop over the different ground type  *
****************************************

      DO 100 i=1,msur

* go to next surface if surface 'i' is set to 0
        IF (xsurf(i).EQ.0.) GOTO 100

* set the surface resistance and rougthness for the given surface
        Ri=Rik(i)
        Rlu0=Rlu0k(i)
        Rac=Rack(i)
        RgsS=RgsSk(i)
        RgsO=RgsOk(i)
        RclS=RclSk(i)
        RclO=RclOk(i)
        z0=z0k(i)
        if (z0 .eq. 0) GOTO 100

* make some correction for temperature dependence
        rt=1000.*EXP(-(temp-273.)-4.)
        Rlu0=Rlu0+rt
        RclS=RclS+rt
        RclO=RclO+rt
        RgsS=RgsS+rt
        RgsO=RgsO+rt

* set maximum and minimum value
        IF (Ri.GE.9998.)   Ri=1E5
        IF (Rlu0.GE.9998.) Rlu0=1E5
        IF (Rac.GE.9998.)  Rac=1E5
        IF (RgsS.GE.9998.) RgsS=1E5
        IF (RgsO.GE.9998.) RgsO=1E5
        IF (RclS.GE.9998.) RclS=1E5
        IF (RclO.GE.9998.) RclO=1E5
        IF (rac.LT.1.)     Rac=1.
        IF (rgsS.LT.1.)    RgsS=1.
        Rs=Ri

        IF (ri.GE.9998.) GOTO 42
        IF ( ((temp-273.).LE.0.) .OR. ((temp-273.).GE.40.) ) GOTO 41
        tempfac=400./( (temp-273.)*(40.-(temp-273.)) )
41      Rs=Ri*(1.+(200./(W+0.1))**2.)*tempfac
42      Rdc=100.*( 1+1000./(W+10.) )/(1.+1000.*gslope)


* compute aerodynamic resistance
* -------------------------------

* compute "flux de chaleur virtuel"
* iseas=1 : hiver    ////  iseas=2 : summer
        IF (iseas.eq.1) THEN
          xmin=200.
          xmax=50.
        ENDIF

        IF (iseas.eq.2) THEN
          xmin=50.
          xmax=250.
        ENDIF

        IF (zen.GE.90.) THEN
          HV=-xmin
        else
          HV=xmax*COS(pi*zen/180.)-xmin
        ENDIF

* vitesse de friction (without correction at first)
        logzsz0=log(max(zs/z0, small))
        ustar=wind*cvk/logzsz0

* make iteration to evaluate correction function
        ustarcor=ustar
        iter=0

16      xx=ustarcor
* correction for the wind

* JMLT 2020: avoid FPEs by pre-filtering HV
        IF(-0.032*HV.LE.-1)THEN
          Z_L=-0.99
          ZZ=log(-Z_L)
          phiM=exp(0.032+0.448*ZZ-0.132*ZZ*ZZ)
          GOTO 14
        ELSEIF(-0.032*HV.GE.1)THEN
          Z_L=0.99
          phiM=-5.*Z_L
          GOTO 14
        ELSE

* original equation &  post-filters
          Z_L=-0.032*HV/max((ustarcor**3.)*temp,small)

          IF (Z_L.LT.0.) GOTO 13
          IF (Z_L.GE.1.) Z_L=0.99
          phiM=-5.*Z_L
          GOTO 14

13        IF (Z_L.le.-1.) Z_L=-0.99
          ZZ=log(-Z_L)
          phiM=exp(0.032+0.448*ZZ-0.132*ZZ*ZZ)

        ENDIF

14      CONTINUE

* estimation of the gap with the previous estimate of ustar
        ustarcor=wind*cvk/(logzsz0-phiM)
        err=abs((ustarcor-xx)/max(xx,small))
        iter=iter+1
        IF (iter.GT.10) GOTO 15
        IF (err.GT.0.01) GOTO 16

15      CONTINUE

* end iteration to evaluate correction function
        phiH=phiM
        IF (Z_L.LT.0) phiH=exp(0.598+0.39*ZZ-0.09*ZZ*ZZ)

* set max value for Ra
        ustark=max(ustarcor*cvk,small)

        IF(ustark.LE.small)THEN
          Ra = 9999.
        ELSE
          Ra=(logzsz0-phiH)/ustark
        ENDIF

* loop over the different species
        DO 200 j=1, ndepspe

* compute laminar resistance
* --------------------------

          DH2O=Dwat*(temp/T0)**1.75

          IF(ustark.LE.small)THEN
            Rb = 9999.
          ELSE
            Rb=(2./max(ustark,small))*((xmu/DH2O)**(2./3.))*
     $          (depdatspe(j,1)**(2./3.))
          ENDIF

* compute surface resistance
* --------------------------

* compute O3 surface resistance
          IF (depnamspe(j)(1:3).EQ.'O3 ') THEN
            Rc=1./(Rs*depdatspe(j,1)) + 1./rlu0 + 1./(Rdc+RclO) +
     &         1./(Rac+RgsO)
            Rc=1./Rc
            GOTO 984
          ENDIF

* compute SO2 surface resistance
          IF (depnamspe(j)(1:4).EQ.'SO2 ') THEN
            Rc=1./(Rs*depdatspe(j,1)) + 1./rlu0 + 1./(Rdc+RclS) +
     &         1./(Rac+RgsS)
            Rc=1./Rc
            GOTO 984
          ENDIF

* compute surface resistance for other species
          Rsm=Rs*depdatspe(j,1) +
     &        1./( depdatspe(j,2)/3000. + 100.*depdatspe(j,3) )
          Rlu=rlu0/( depdatspe(j,2)*1E-5 + depdatspe(j,3) )
          Rcl=1./( (depdatspe(j,2)/1E5/RclS) + (depdatspe(j,3)/RclO) )
          Rgs=1./( (depdatspe(j,2)/1E5/RgsS) + (depdatspe(j,3)/RgsO) )
          Rc=1./( (1./Rsm) + (1./Rlu) + (1./(Rdc+Rcl)) +
     &          (1./(Rac+Rgs)) )

984       CONTINUE
* test minimum and maximum value
          IF (Rc.GT.9999.) Rc=9999.
          IF (Rc.LT.10.) Rc=10.

* compute deposition velocity
* ----------------------------

          Rtot=Ra+Rb+Rc

c          IF (i.eq.4) THEN
c            IF (depnamspe(j)(1:4).EQ.'NO2 ') THEN
c              write(53,'(6(E11.3))')
c     &              time/3600.,ra,rb,rc,rtot,(1./rtot)*100.
c            ENDIF
c            IF (depnamspe(j)(1:5).EQ.'HONO ') THEN
c              write(54,'(6(E11.3))')
c     &              time/3600.,ra,rb,rc,rtot,(1./rtot)*100.
c            ENDIF
c            IF (depnamspe(j)(1:5).EQ.'HNO3 ') THEN
c              write(55,'(6(E11.3))')
c     &              time/3600.,ra,rb,rc,rtot,(1./rtot)*100.
c            ENDIF
c          ENDIF

          V=1./Rtot
* convert metres to centimeters
          Vdtyp(j,i)=V*100.

200     CONTINUE
100   CONTINUE


**************************************************
* compute deposition velocities for the scenario *
**************************************************
      DO i=1,msur
        DO j=1,ndepspe
          Vd(j)=Vd(j)+xsurf(i)*Vdtyp(j,i)
        ENDDO
      ENDDO

      ResAerdep = (Ra + Rs)

c      vd(8)=vd(8)*2.

c      DO j=1,ndepspe
c        yy(j)=0.
c        DO i=1,msur
c          yy(j)=yy(j) + xsurf(i)*Vdtyp(j,i)
c        ENDDO
c      ENDDO
c      write(53,'(20(E10.3))') time/3600.,(yy(j),j=1,12)


      END
