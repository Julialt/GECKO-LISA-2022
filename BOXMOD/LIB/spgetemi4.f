*******************************************************************
* Cette routine a pour objet de calculer les flux d'emission      *
* relatifs aux differents types de surface.                       *
* surface 1 => surface urbaine                                    *
* surface 2 => surface terre cultive                              *
* surface 3 => surface foret de feuillu                           *
* surface 4 => surface foret de connifere                         *
*                                                                 *
* INPUT :                                                         *
*   lout   : numero de fichier pour ecriture des erreurs          *
*   temp   : temperature                                          *
*   time   : temps (modulo 24h)                                   *
*   iscape : drapeau pour calcul thermo
*                                                                 *
* EMISSION ....                                                   *
*   nenott    : nombre d'espece emise dont le facteur d'emission  *
*               ne depend pas du temps (SO2,CH4,Ca,K,CL,Na,Mg)    *
*   idenott(i): identification number de la ieme espece dans le   *
*               schema chimique                                   *
*   enott(i,j): facteur d'emission de de la ieme espece pour      *
*               le jieme environnement                            *
*   idnh3     : identification number de NH3                      *
*   enh3(j)   : facteur d'emission de NH3 pour le jieme           *
*               environnement                                     *
*   idisop    : identification number de l'isoprene               *
*   idapin    : identification number de l'alpha-pinene           *
*   idbpin    : identification number du beta-pinene              *
*   eisop(j)  : facteur d'emission standart de l'isoprene pour    *
*               le jieme environnement                            *
*   eapin(j)  : facteur d'emission standart de l'a-pinene pour    *
*               le jieme environnement                            *
*   ebpin(j)  : facteur d'emission standart du b-pinene pour      *
*               le jieme environnement                            *
*   ecour     : facteur d'emission du CO pour environnement       *
*               urbain                                            *
*   enoxur    : facteur d'emission des NOx pour environnement     *
*               urbain                                            *
*   senoxur(i): speciation des NOx pour environnement             *
*               urbain (i=1=> NO, i=2=> NO2, i=3=> HONO)          *
*   nshcur    : nombre d'hydrocarbure dans la speciation des HC   *
*               pour le scenario urbain                           *
*   idhcur(i) : identification number de la ieme espece de la     *
*               la speciation des HC urbain                       *
*   ehcur     : facteur d'emission des HC                         *
*   ehcdatur(i,j) : speciation pour le ieme HC (j=1) et nombre de *
*                   carbone (j=2)                                 *
*   cscoef3(i,j) : parametre pour interpolation des facteurs      *
*                 d'emission des activites urbaines               *
*                                                                 *
* Pour le calcul du PAR                                           *
*   sla         : latitude                                        *
*   slo         : longitude                                       *
*   tz          : decalage horaire par rapport au meridien de     *
*                 greenwich en TU                                 *
*   iy          : annee                                           *
*   im          : mois                                            *
*   id          : jour                                            *
*                                                                 *
* OUTPUT :                                                        *
*   eflux(i) : flux d'emission pour l'espece i                    *
*******************************************************************
      SUBROUTINE getemi4(chrsp,numsp,temp,time,xsurf,iscape,
     1                   sla,slo,tz,iy,im,id,szafix,szaval,
     2                   ecour,enoxur,ehcur,senoxur,nshcur,
     3                   idehcur,ehcdatur,cscoef3,
     4                   nenott,idenott,enott,
     5                   idnh3,enh3,
     6                   idisop,idapin,idbpin,eisop,eapin,ebpin,
     7                   eflux)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      CHARACTER(maxlsp) chrsp(maxsp)
      INTEGER  numsp,iscape
      REAL     temp,time
      REAL     sla,slo,tz
      INTEGER  iy,im,id
      INTEGER  nshcur,nenott,idnh3,idisop,idapin,idbpin
      INTEGER  idehcur(maxem),idenott(maxem)
      REAL     enoxur,ecour,ehcur
      REAL     ehcdatur(maxem,2),senoxur(3)
      REAL     enott(maxem,msur)
      REAL     enh3(msur),eisop(msur),eapin(msur),ebpin(msur)
      REAL     xsurf(msur)
      REAL     cscoef3(4,25)
      INTEGER  szafix
      REAL     szaval

* OUTPUT
      REAL     eflux(maxsp)

* LOCAL
      CHARACTER(maxlsp) namesp
      REAL     coefbnh3(24), par
      INTEGER  i,idco,idno,idno2,idhono,nt,itime
      REAL     xtime,pcoe,em,f_isotemp,f_isolight,f_tertemp

* facteur d'emission pour NH3 biogenique
      data coefbnh3 /0.65, 0.65, 0.65, 0.65, 0.65, 0.65,
     &               0.75, 0.90, 1.05, 1.15, 1.25, 1.40,
     &               1.50, 1.60, 1.60, 1.45, 1.35, 1.25,
     &               1.10, 1.00, 0.85, 0.70, 0.70, 0.70/

* initialise
      DO i=1,maxsp
        eflux(i)=0.
      ENDDO

* get some ID number that are not provided by the input
      namesp='GCO '
      CALL akspnum(namesp,chrsp,numsp,idco)
      IF (idco.eq.0) THEN
        WRITE(6,*) '--error--, in getemi3'
        WRITE(6,*) '           species CO not found'
        STOP
      ENDIF

      namesp='GNO '
      CALL akspnum(namesp,chrsp,numsp,idno)
      IF (idno.eq.0) THEN
        WRITE(6,*) '--error--, in getemi3'
        WRITE(6,*) '           species NO not found'
        STOP
      ENDIF

      namesp='GNO2 '
      CALL akspnum(namesp,chrsp,numsp,idno2)
      IF (idno2.eq.0) THEN
        WRITE(6,*) '--error--, in getemi3'
        WRITE(6,*) '           species NO2 not found'
        STOP
      ENDIF

      namesp='GHNO2 '
      CALL akspnum(namesp,chrsp,numsp,idhono)
      IF (idco.eq.0) THEN
        WRITE(6,*) '--error--, in getemi3'
        WRITE(6,*) '           species HONO not found'
        STOP
      ENDIF

* compute PAR
      CALL get_par(time,par,sla,slo,tz,iy,im,id,szafix,szaval)

*********************************************************
* COMPUTE EMISSION FOR THE URBAN ENVIRONNEMENT (FLAG=1) *
*********************************************************

      IF (xsurf(1).gt.0.) THEN

* calcul des emissions ayant un facteur d'emission INDEPENDANT du temps
* =====================================================================

* species in the enott table (SO2,CH4,Ca, ...)
        DO i=1,nenott
          eflux(idenott(i))=eflux(idenott(i)) + (enott(i,1)*xsurf(1))
        ENDDO

* NH3 emission (only if iscape set to 1)
        IF (iscape.eq.1) THEN
          eflux(idnh3)=eflux(idnh3) + (enh3(1)*xsurf(1))
        ENDIF

* calcul des emissions ayant un facteur d'emission DEPENDANT du temps
* ===================================================================

* calcul du coefficient horaire
        itime=int(time/3600.)
        xtime=time/3600.-real(itime)
        itime=itime+1
        pcoe=cscoef3(1,itime) + cscoef3(2,itime)*xtime +
     &        0.5*cscoef3(3,itime)*(xtime**2.) +
     &       (1./6.)*cscoef3(4,itime)*(xtime**3.)
        pcoe=pcoe*0.24

* emission des NOx
        eflux(idno)=eflux(idno)+(enoxur*senoxur(1)*pcoe*xsurf(1))
        eflux(idno2)=eflux(idno2)+(enoxur*senoxur(2)*pcoe*xsurf(1))
        eflux(idhono)=eflux(idhono)+(enoxur*senoxur(3)*pcoe*xsurf(1))

* emission de CO
        eflux(idco)=eflux(idco)+(ecour*pcoe*xsurf(1))

* emission des COV
        DO i=1,nshcur
          eflux(idehcur(i))= eflux(idehcur(i)) +
     &             (ehcur*pcoe*xsurf(1)*(ehcdatur(i,1)/ehcdatur(i,2)))
        ENDDO
c        DO i=1,nshcur
c          eflux(idehcur(i))= eflux(idehcur(i)) +
c     &    (ehcur*pcoe*xsurf(1)*6.02E23*(ehcdatur(i,1)/ehcdatur(i,2)))
c        ENDDO
      ENDIF

*********************************************************
* COMPUTE EMISSION FOR THE CULTURE (FLAG=2)             *
*********************************************************

      IF (xsurf(2).gt.0.) THEN

* calcul des emissions ayant un facteur d'emission INDEPENDANT du temps
* =====================================================================

* species in the enott table (SO2,CH4,Ca, ...)
        DO i=1,nenott
          eflux(idenott(i))=eflux(idenott(i))+ (enott(i,2)*xsurf(2))
        ENDDO

* calcul des emissions ayant un facteur d'emission DEPENDANT du temps
* ===================================================================

* NH3 emission (only if iscape set to 1)
        IF (iscape.eq.1) THEN
          nt=int(time/3600. + 1)
          eflux(idnh3)=eflux(idnh3) + (enh3(2)*coefbnh3(nt))*xsurf(2)
        ENDIF

* NO emission
        em=xsurf(2)*0.90*EXP(0.071*(0.67*(temp-273.))+8.8)
        eflux(idno)=eflux(idno)+em

* isoprene emission
        f_isotemp=exp((9.5E4*(temp-303.))/(8.314*303.*temp))/
     &            (1+exp((2.3E5*(temp-314.))/(8.314*303.*temp)))
        f_isolight=(0.0027*1.066*par)/((1.+7.29E-6*par**2)**0.5)
        em=xsurf(2)*eisop(2)*f_isotemp*f_isolight
        eflux(idisop)=eflux(idisop)+em

* alpha pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(2)*eapin(2)*f_tertemp
        eflux(idapin)=eflux(idapin)+em

* beta pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(2)*ebpin(2)*f_tertemp
        eflux(idbpin)=eflux(idbpin) + em

      ENDIF

*********************************************************
* COMPUTE EMISSION FOR THE LEAF FOREST (FLAG=3)         *
*********************************************************

      IF (xsurf(3).gt.0.) THEN

* calcul des emissions ayant un facteur d'emission INDEPENDANT du temps
* =====================================================================

* species in the enott table (SO2,CH4,Ca, ...)
        DO i=1,nenott
          eflux(idenott(i))=eflux(idenott(i))+ (enott(i,3)*xsurf(3))
        ENDDO

* calcul des emissions ayant un facteur d'emission DEPENDANT du temps
* ===================================================================

* NH3 emission (only if iscape set to 1)
        IF (iscape.eq.1) THEN
          nt=int(time/3600. + 1)
          eflux(idnh3)=eflux(idnh3) + (enh3(3)*coefbnh3(nt))*xsurf(3)
        ENDIF

* NO emission
        em=xsurf(3)*0.07*EXP(0.071*(0.84*(temp-273.))+3.6)
        eflux(idno)=eflux(idno)+em

* isoprene emission
        f_isotemp=exp((9.5E4*(temp-303.))/(8.314*303.*temp))/
     &            (1+exp((2.3E5*(temp-314.))/(8.314*303.*temp)))
        f_isolight=(0.0027*1.066*par)/((1.+7.29E-6*par**2)**0.5)
        em=xsurf(3)*eisop(3)*f_isotemp*f_isolight
        eflux(idisop)=eflux(idisop)+em

* alpha pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(3)*eapin(3)*f_tertemp
        eflux(idapin)=eflux(idapin)+em

* beta pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(3)*ebpin(3)*f_tertemp
        eflux(idbpin)=eflux(idbpin) + em

      ENDIF

*********************************************************
* COMPUTE EMISSION FOR THE CONIFERE FOREST (FLAG=4)     *
*********************************************************

      IF (xsurf(4).gt.0.) THEN

* calcul des emissions ayant un facteur d'emission INDEPENDANT du temps
* =====================================================================

* species in the enott table (SO2,CH4,Ca, ...)
        DO i=1,nenott
          eflux(idenott(i))=eflux(idenott(i))+ (enott(i,4)*xsurf(4))
        ENDDO

* calcul des emissions ayant un facteur d'emission DEPENDANT du temps
* ===================================================================

* NH3 emission (only if iscape set to 1)
        IF (iscape.eq.1) THEN
          nt=int(time/3600. + 1)
          eflux(idnh3)=eflux(idnh3) + (enh3(4)*coefbnh3(nt))*xsurf(4)
        ENDIF

* NO emission
        em=xsurf(4)*0.07*EXP(0.071*(0.84*(temp-273.))+3.6)
        eflux(idno)=eflux(idno)+em

* isoprene emission
        f_isotemp=exp((9.5E4*(temp-303.))/(8.314*303.*temp))/
     &            (1+exp((2.3E5*(temp-314.))/(8.314*303.*temp)))
        f_isolight=(0.0027*1.066*par)/((1.+7.29E-6*par**2)**0.5)
        em=xsurf(4)*eisop(4)*f_isotemp*f_isolight
        eflux(idisop)=eflux(idisop)+em

* alpha pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(4)*eapin(4)*f_tertemp
        eflux(idapin)=eflux(idapin)+em

* beta pinene
        f_tertemp=EXP(0.07*((temp-273.)-30.))
        em=xsurf(4)*ebpin(4)*f_tertemp
        eflux(idbpin)=eflux(idbpin) + em

      ENDIF

      END

*******************************************************************
*******************************************************************

*******************************************************************
* cette routine a pour objet de calculer le PAR                   *
*                                                                 *
* INPUT :                                                         *
*   time   : temps (modulo 24h)                                   *
*   sla    : latitude                                             *
*   slo    : longitude                                            *
*   tz     : decalage horaire par rapport au meridien de          *
*            greenwich en TU                                      *
*   iy     : annee                                                *
*   im     : mois                                                 *
*   id     : jour                                                 *
*                                                                 *
* OUTPUT :                                                         *
*   par
*******************************************************************
      SUBROUTINE get_par(time,par,sla,slo,tz,iy,im,id,szafix,szaval)
      IMPLICIT NONE

* INPUT
      REAL     time
      REAL     sla,slo,tz
      INTEGER  iy,im,id
      INTEGER  szafix
      REAL     szaval

* OUTPUT
      REAL     par

* LOCAL
      INTEGER  iti,imn
      REAL     ti,xmn,xc,xd,So,a,c,d,an,sw

* time format for subroutine solar
      iti=int(time/3600.)
      xmn=time/60.
      imn=int(xmn-(iti*60))
      ti = iti*100.0 +imn

* compute zenithal angle
      xc = 0.0
      CALL solar(sla,slo,tz,iy,im,id,ti,xc,5)
      xd = 90.-xc
** JMLT EDIT FOR FIXED SZA **
      IF(szafix.EQ.1) xd = szaval
** END EDIT **

* Ref: J.A. Foley JGR 99 (D10) pp 20773-20783, 1994
* par : fournie en Ein.m-2.s-1 = moles.m-2.s-1
* conversion em micromoles.m-2.s-1
* nul si c'est la nuit, xd superieur a 90

      IF (xd.gt.85.0) THEN
        par=0.
*        rsol=0.
      ELSE
        So=1370.
        a=0.15
        c=0.25
        d=0.50
        an=0.6

        sw=So*(1-a)*(c+d*an)*COS(3.1415927*xd/180.)
        par=(sw*1E6)/(2.*2.18E5)
*        rsol=So*(c+d*an)*COS(3.1415927*xd/180.)
      ENDIF

      END


