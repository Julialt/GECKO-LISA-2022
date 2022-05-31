*******************************************************************
* Cette routine a pour objet de lire l'ensemble des donnees       *
* relative au different type de surface.                          *
* surface 1 => surface urbaine                                    *
* surface 2 => surface terre cultive                              *
* surface 3 => surface foret de feuillu                           *
* surface 4 => surface foret de connifere                         *
*                                                                 *
* INPUT :                                                         *
*   lout   : numero de fichier pour ecriture des erreurs          *
*   iscape : drapeau si le calcul thermo est demande              *
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
*               urbain (i=1=> NO, i=2=> NO2, i=3=> HONO)                      *
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
* RESISTANCE DE SURFACE .....                                     *
* OUTPUT :                                                        *
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
*******************************************************************
      SUBROUTINE datsurf(lout,chrsp,numsp,iscape,iseas,
     1                   ecour,enoxur,ehcur,senoxur,nshcur,
     2                   idehcur,ehcdatur,cscoef3,
     3                   nenott,idenott,enott,
     4                   idnh3,enh3,
     5                   idisop,idapin,idbpin,eisop,eapin,ebpin,
     6                   rik,rlu0k,rack,rgsSk,rgsOk,rclSk,rclOk,z0k)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lout,numsp,iscape,iseas
      CHARACTER(maxlsp) :: chrsp(maxsp)

* OUTPUT
      INTEGER  nshcur, nenott, idnh3, idisop, idapin,idbpin
      INTEGER  idehcur(maxem),idenott(maxem)
      REAL     ecour,enoxur,ehcur
      REAL     ehcdatur(maxem,2),senoxur(3)
      REAL     enott(maxem,msur)
      REAL     enh3(msur),eisop(msur),eapin(msur),ebpin(msur)
      REAL     cscoef3(4,25)
      REAL     rik(msur),rlu0k(msur),rack(msur),rgsSk(msur)
      REAL     rgsOk(msur),rclSk(msur),rclOk(msur),z0k(msur)

* local
      CHARACTER(76) line
      CHARACTER(4)  keyword
      CHARACTER(maxlsp) ehcnamcu(maxem)
      CHARACTER(maxlsp) ehcnamff(maxem)
      CHARACTER(maxlsp) ehcnamfc(maxem)
      CHARACTER(maxlsp) ehcnamur(maxem)
      CHARACTER(maxlsp) namesp
      CHARACTER(200) filename
      INTEGER  nchar
      REAL     ehcdatcu(maxem,2),idehccu(maxem)
      REAL     ehcdatff(maxem,2),idehcff(maxem)
      REAL     ehcdatfc(maxem,2),idehcfc(maxem)
      INTEGER  i,j,nsur, i_val, ierr, llen, isp, lenstr
      INTEGER  nshccu, nshcff,nshcfc
      REAL     r_val
      REAL     eso2ur,ech4ur,enh3ur,ecaur,ekur,eclur,enaur,emgur
      REAL     ech4cu,ehccu,enh3cu,ecacu,ekcu,eclcu,enacu,emgcu
      REAL     ech4ff,ehcff,enh3ff,ecaff,ekff,eclff,enaff,emgff
      REAL     ech4fc,ehcfc,enh3fc,ecafc,ekfc,eclfc,enafc,emgfc

      nchar=76

* initialise the output data
      enoxur=0.
      ecour=0.
      ehcur=0.
      senoxur(1)=0.
      senoxur(2)=0.
      senoxur(3)=0.
      nshcur=0
      nenott=0
      idnh3=0
      idisop=0
      idapin=0
      idbpin=0
      DO i=1,maxem
        idehcur(i)=0
        ehcdatur(i,1)=0.
        ehcdatur(i,2)=0.
        idenott(i)=0
      ENDDO
      DO i=1,msur
        enh3(i)=0.
        eisop(i)=0.
        eapin(i)=0.
        ebpin(i)=0.
        rik(i)=0.
        rlu0k(i)=0.
        rack(i)=0.
        rgsSk(i)=0.
        rgsOk(i)=0.
        rclSk(i)=0.
        rclOk(i)=0.
        z0k(i)=0.
        DO j=1,maxem
          enott(j,i)=0.
        ENDDO
      ENDDO

!     RETURN

**************************************************
* read data for urban surface (FLAG=1)           *
**************************************************

      write(6,*) '      ... reading urban data'

* initialise
      eso2ur=0.
      enoxur=0.
      ecour=0.
      ech4ur=0.
      ehcur=0.
      enh3ur=0.
      ecaur=0.
      ekur=0.
      eclur=0.
      enaur=0.
      emgur=0.
      senoxur(1)=0.
      senoxur(2)=0.
      senoxur(3)=0.
      nshcur=0
      DO i=1,maxem
        ehcnamur(i)=' '
        idehcur(i)=0
        ehcdatur(i,1)=0.
        ehcdatur(i,2)=0.
      ENDDO

* open the file
      nsur=1
      OPEN (12,file='urbain.sur',status='OLD')

* return label to read new keyword
90    CONTINUE
      READ (12,'(a4,(a))') keyword,line
c      write(6,'(a4,(a))') keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'/') GOTO 90

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 100

* emission de SO2
      ELSE IF (keyword(1:4).EQ.'ESO2') THEN
        eso2ur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword ESO2'
          STOP
        ENDIF
        GOTO 90
* emission de NOx
      ELSE IF (keyword(1:4).EQ.'ENOX') THEN
        enoxur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword ENOX'
          STOP
        ENDIF
        GOTO 90
* emission de CO
      ELSE IF (keyword(1:4).EQ.'EMCO') THEN
        ecour=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EMCO'
          STOP
        ENDIF
        GOTO 90
* emission de CH4
      ELSE IF (keyword(1:4).EQ.'ECH4') THEN
        ech4ur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword ECH4'
          STOP
        ENDIF
        GOTO 90
* emission de non methane HC
      ELSE IF (keyword(1:4).EQ.'EMHC') THEN
        ehcur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EMHC'
          STOP
        ENDIF
        GOTO 90
* emission de NH3
      ELSE IF (keyword(1:4).EQ.'ENH3') THEN
        enh3ur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword ENH3'
          STOP
        ENDIF
        GOTO 90
* emission de calcium particulaire
      ELSE IF (keyword(1:4).EQ.'EPCA') THEN
        ecaur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EPCA'
          STOP
        ENDIF
        GOTO 90
* emission de potassium particulaire
      ELSE IF (keyword(1:4).EQ.'EMPK') THEN
        ekur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EMPK'
          STOP
        ENDIF
        GOTO 90
* emission de chlorure particuliare
      ELSE IF (keyword(1:4).EQ.'EPCL') THEN
        eclur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EPCL'
          STOP
        ENDIF
        GOTO 90
* emission de sodium particulaire
      ELSE IF (keyword(1:4).EQ.'EPNA') THEN
        enaur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EPNA'
          STOP
        ENDIF
        GOTO 90
* emission de magnesium particulaire
      ELSE IF (keyword(1:4).EQ.'EPMG') THEN
        emgur=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword EPMG'
          STOP
        ENDIF
        GOTO 90
* speciation des NOx a l'emission
      ELSE IF (keyword(1:4).EQ.'SNOX') THEN
        senoxur(1)=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword SNOX'
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        senoxur(2)=r_val(line,1,nchar,2,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword SNOX'
          WRITE(lout,*)'            second value not read correctly'
          STOP
        ENDIF
        senoxur(3)=r_val(line,1,nchar,3,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword SNOX'
          WRITE(lout,*)'            third value not read correctly'
          STOP
        ENDIF
        GOTO 90
* speciation des HC a l'emission
      ELSE IF (keyword(1:4).EQ.'SPHC') THEN
        nshcur=i_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)'            keyword SPHC'
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        DO i=1,nshcur
          READ(12,'(a)') line
          ehcnamur(i)=line(1:16)
          CALL akspnum(ehcnamur(i),chrsp,numsp,isp)
          idehcur(i)=isp
          IF (isp.eq.0)then
            WRITE(lout,*)' --error--  while reading urbain.sur'
            WRITE(lout,*)' species unkown=',idehcur(i)
            STOP
          ENDIF
          ehcdatur(i,1)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)'            species name=',ehcnamur(i)
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          ehcdatur(i,2)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)'            species name=',ehcnamur(i)
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
        ENDDO
        GOTO 90

* surface resistance : season 1
      ELSE IF (keyword(1:4).EQ.'RES1') THEN
        IF (iseas.EQ.1) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 90

* surface resistance : season 2
      ELSE IF (keyword(1:4).EQ.'RES2') THEN
        IF (iseas.EQ.2) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in urbain.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 90

* roughness height
      ELSE IF (keyword(1:4).EQ.'HRUG') THEN
        z0k(nsur)=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in urbain.sur :'
          WRITE(lout,*)' line= ',line
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        GOTO 90

* keyword unknown : stop
      ELSE
        WRITE(lout,*)' --error--  while reading in urbain.sur'
        WRITE(lout,*)'            keyword unkown: ',keyword
        STOP
      ENDIF

100   CONTINUE
      CLOSE (12)

C * open the file that gives the emission coefficient for the urban environment
C        CALL GETENV("HOME",line)
C        llen=lenstr(line,nchar)
C        IF (llen.lt.1) THEN
C           WRITE(lout,*) 'no environment variable HOME'
C           WRITE(lout,*) 'emission coefficient for the urban scenario'
C           WRITE(lout,*) 'can not be read'
C           STOP
C        ENDIF
C        filename=line(1:llen)//"/V1/INPUT/EMISSION/emi_coef.dat"
C        llen=lenstr(filename,200)
C        open(12,file=filename(1:llen),form='formatted',status='old')

C * read the coefficient
C       DO j=1,25
C          DO i=1,4
C            READ(12,*) cscoef3(i,j)
C          ENDDO
C       ENDDO
C       CLOSE (12)

**************************************************
* read data for culture surface (FLAG=2)         *
**************************************************

      write(6,*) '      ... reading crops data'

* initialise
      ech4cu=0.
      ehccu=0.
      enh3cu=0.
      ecacu=0.
      ekcu=0.
      eclcu=0.
      enacu=0.
      emgcu=0.
      nshccu=0
      DO i=1,maxem
        ehcnamcu(i)=' '
        idehccu(i)=0
        ehcdatcu(i,1)=0.
        ehcdatcu(i,2)=0.
      ENDDO

* open the file
      nsur=2
      OPEN (12,file='culture.sur',status='OLD')

* return label to read new keyword
91    CONTINUE
      READ (12,'(a4,(a))') keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'/') GOTO 91

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 101

* emission de CH4
      ELSE IF (keyword(1:4).EQ.'ECH4') THEN
        ech4cu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword ECH4'
          STOP
        ENDIF
        GOTO 91
* emission de non methane HC
      ELSE IF (keyword(1:4).EQ.'EBHC') THEN
        ehccu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EBHC'
          STOP
        ENDIF
        GOTO 91
* emission de NH3
      ELSE IF (keyword(1:4).EQ.'ENH3') THEN
        enh3cu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword ENH3'
          STOP
        ENDIF
        GOTO 91
* emission de calcium particulaire
      ELSE IF (keyword(1:4).EQ.'EPCA') THEN
        ecacu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EPCA'
          STOP
        ENDIF
        GOTO 91
* emission de potassium particulaire
      ELSE IF (keyword(1:4).EQ.'EMPK') THEN
        ekcu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EMPK'
          STOP
        ENDIF
        GOTO 91
* emission de chlorure particuliare
      ELSE IF (keyword(1:4).EQ.'EPCL') THEN
        eclcu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EPCL'
          STOP
        ENDIF
        GOTO 91
* emission de sodium particulaire
      ELSE IF (keyword(1:4).EQ.'EPNA') THEN
        enacu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EPNA'
          STOP
        ENDIF
        GOTO 91
* emission de magnesium particulaire
      ELSE IF (keyword(1:4).EQ.'EPMG') THEN
        emgcu=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword EPMG'
          STOP
        ENDIF
        GOTO 91
* speciation des HC a l'emission
      ELSE IF (keyword(1:4).EQ.'SBHC') THEN
        nshccu=i_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)'            keyword SBHC'
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        DO i=1,nshccu
          READ(12,'(a)') line
          ehcnamcu(i)=line(1:16)
          CALL akspnum(ehcnamcu(i),chrsp,numsp,isp)
          idehccu(i)=isp
          IF (isp.eq.0)then
            WRITE(lout,*)' --error--  while reading culture.sur'
            WRITE(lout,*)' species unkown=',idehccu(i)
            STOP
          ENDIF
          ehcdatcu(i,1)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)'            species name=',ehcnamcu(i)
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          ehcdatcu(i,2)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)'            species name=',ehcnamcu(i)
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
        ENDDO
        GOTO 91

* surface resistance : season 1
      ELSE IF (keyword(1:4).EQ.'RES1') THEN
        IF (iseas.EQ.1) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 91

* surface resistance : season 2
      ELSE IF (keyword(1:4).EQ.'RES2') THEN
        IF (iseas.EQ.2) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in culture.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 91

* roughness height
      ELSE IF (keyword(1:4).EQ.'HRUG') THEN
        z0k(nsur)=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in culture.sur :'
          WRITE(lout,*)' line= ',line
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        GOTO 91

* keyword unknown : stop
      ELSE
        WRITE(lout,*)' --error--  while reading in culture.sur'
        WRITE(lout,*)'            keyword unkown: ',keyword
        STOP
      ENDIF

101   CONTINUE
      CLOSE (12)

**************************************************
* read data for leaf forest surface (FLAG=3)     *
**************************************************

      write(6,*) '      ... reading leaf forest data'

* initialise
      ech4ff=0.
      ehcff=0.
      enh3ff=0.
      ecaff=0.
      ekff=0.
      eclff=0.
      enaff=0.
      emgff=0.
      nshcff=0
      DO i=1,maxem
        ehcnamff(i)=' '
        idehcff(i)=0
        ehcdatff(i,1)=0.
        ehcdatff(i,2)=0.
      ENDDO

* open the file
      nsur=3
      OPEN (12,file='foretfeu.sur',status='OLD')

* return label to read new keyword
92    CONTINUE
      READ (12,'(a4,(a))') keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'/') GOTO 92

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 102

* emission de CH4
      ELSE IF (keyword(1:4).EQ.'ECH4') THEN
        ech4ff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword ECH4'
          STOP
        ENDIF
        GOTO 92
* emission de non methane HC
      ELSE IF (keyword(1:4).EQ.'EBHC') THEN
        ehcff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EMHC'
          STOP
        ENDIF
        GOTO 92
* emission de NH3
      ELSE IF (keyword(1:4).EQ.'ENH3') THEN
        enh3ff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword ENH3'
          STOP
        ENDIF
        GOTO 92
* emission de calcium particulaire
      ELSE IF (keyword(1:4).EQ.'EPCA') THEN
        ecaff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EPCA'
          STOP
        ENDIF
        GOTO 92
* emission de potassium particulaire
      ELSE IF (keyword(1:4).EQ.'EMPK') THEN
        ekff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EMPK'
          STOP
        ENDIF
        GOTO 92
* emission de chlorure particuliare
      ELSE IF (keyword(1:4).EQ.'EPCL') THEN
        eclff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EPCL'
          STOP
        ENDIF
        GOTO 92
* emission de sodium particulaire
      ELSE IF (keyword(1:4).EQ.'EPNA') THEN
        enaff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EPNA'
          STOP
        ENDIF
        GOTO 92
* emission de magnesium particulaire
      ELSE IF (keyword(1:4).EQ.'EPMG') THEN
        emgff=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword EPMG'
          STOP
        ENDIF
        GOTO 92
* speciation des HC a l'emission
      ELSE IF (keyword(1:4).EQ.'SBHC') THEN
        nshcff=i_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)'            keyword SBHC'
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        DO i=1,nshcff
          READ(12,'(a)') line
          ehcnamff(i)=line(1:16)
          CALL akspnum(ehcnamff(i),chrsp,numsp,isp)
          idehcff(i)=isp
          IF (isp.eq.0)then
            WRITE(lout,*)' --error--  while reading foretfeu.sur'
            WRITE(lout,*)' species unkown=',idehcff(i)
            STOP
          ENDIF
          ehcdatff(i,1)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)'            species name=',ehcnamff(i)
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          ehcdatff(i,2)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)'            species name=',ehcnamff(i)
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
        ENDDO
        GOTO 92

* surface resistance : season 1
      ELSE IF (keyword(1:4).EQ.'RES1') THEN
        IF (iseas.EQ.1) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 92

* surface resistance : season 2
      ELSE IF (keyword(1:4).EQ.'RES2') THEN
        IF (iseas.EQ.2) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 92

* roughness height
      ELSE IF (keyword(1:4).EQ.'HRUG') THEN
        z0k(nsur)=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretfeu.sur :'
          WRITE(lout,*)' line= ',line
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        GOTO 92

* keyword unknown : stop
      ELSE
        WRITE(lout,*)' --error--  while reading in foretfeu.sur'
        WRITE(lout,*)'            keyword unkown: ',keyword
        STOP
      ENDIF

102   CONTINUE
      CLOSE (12)

***************************************************
* read data for connifere forest surface (FLAG=4) *
***************************************************

      write(6,*) '      ... reading connifere forest data'

* initialise
      ech4fc=0.
      ehcfc=0.
      enh3fc=0.
      ecafc=0.
      ekfc=0.
      eclfc=0.
      enafc=0.
      emgfc=0.
      nshcfc=0
      DO i=1,maxem
        ehcnamfc(i)=' '
        idehcfc(i)=0
        ehcdatfc(i,1)=0.
        ehcdatfc(i,2)=0.
      ENDDO

* open the file
      nsur=4
      OPEN (12,file='foretcon.sur',status='OLD')

* return label to read new keyword
93    CONTINUE
      READ (12,'(a4,(a))') keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'/') GOTO 93

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 103

* emission de CH4
      ELSE IF (keyword(1:4).EQ.'ECH4') THEN
        ech4fc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword ECH4'
          STOP
        ENDIF
        GOTO 93
* emission de non methane HC
      ELSE IF (keyword(1:4).EQ.'EBHC') THEN
        ehcfc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EMHC'
          STOP
        ENDIF
        GOTO 93
* emission de NH3
      ELSE IF (keyword(1:4).EQ.'ENH3') THEN
        enh3fc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword ENH3'
          STOP
        ENDIF
        GOTO 93
* emission de calcium particulaire
      ELSE IF (keyword(1:4).EQ.'EPCA') THEN
        ecafc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EPCA'
          STOP
        ENDIF
        GOTO 93
* emission de potassium particulaire
      ELSE IF (keyword(1:4).EQ.'EMPK') THEN
        ekfc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EMPK'
          STOP
        ENDIF
        GOTO 93
* emission de chlorure particuliare
      ELSE IF (keyword(1:4).EQ.'EPCL') THEN
        eclfc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EPCL'
          STOP
        ENDIF
        GOTO 93
* emission de sodium particulaire
      ELSE IF (keyword(1:4).EQ.'EPNA') THEN
        enafc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EPNA'
          STOP
        ENDIF
        GOTO 93
* emission de magnesium particulaire
      ELSE IF (keyword(1:4).EQ.'EPMG') THEN
        emgfc=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword EPMG'
          STOP
        ENDIF
        GOTO 93
* speciation des HC a l'emission
      ELSE IF (keyword(1:4).EQ.'SBHC') THEN
        nshcfc=i_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)'            keyword SBHC'
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        DO i=1,nshcfc
          READ(12,'(a)') line
          ehcnamfc(i)=line(1:16)
          CALL akspnum(ehcnamfc(i),chrsp,numsp,isp)
          idehcfc(i)=isp
          IF (isp.eq.0)then
            WRITE(lout,*)' --error--  while reading foretcon.sur'
            WRITE(lout,*)' species unkown=',idehcfc(i)
            STOP
          ENDIF
          ehcdatfc(i,1)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)'            species name=',ehcnamfc(i)
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          ehcdatfc(i,2)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0)then
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)'            species name=',ehcnamfc(i)
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
        ENDDO
        GOTO 93

* surface resistance : season 1
      ELSE IF (keyword(1:4).EQ.'RES1') THEN
        IF (iseas.EQ.1) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 93

* surface resistance : season 2
      ELSE IF (keyword(1:4).EQ.'RES2') THEN
        IF (iseas.EQ.2) THEN
          rik(nsur)=r_val(line,1,nchar,1,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            first value not read correctly'
            STOP
          ENDIF
          rlu0k(nsur)=r_val(line,1,nchar,2,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            second value not read correctly'
            STOP
          ENDIF
          rack(nsur)=r_val(line,1,nchar,3,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            third value not read correctly'
            STOP
          ENDIF
          rgsSk(nsur)=r_val(line,1,nchar,4,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            forth value not read correctly'
            STOP
          ENDIF
          rgsOk(nsur)=r_val(line,1,nchar,5,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            fiveth value not read correctly'
            STOP
          ENDIF
          rclSk(nsur)=r_val(line,1,nchar,6,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            sixth value not read correctly'
            STOP
          ENDIF
          rclOk(nsur)=r_val(line,1,nchar,7,ierr)
          IF(ierr.ne.0) THEN
            WRITE(lout,*)' --error--  while reading in foretcon.sur :'
            WRITE(lout,*)' line= ',line
            WRITE(lout,*)'            seventh value not read correctly'
            STOP
          ENDIF
        ENDIF
        GOTO 93

* roughness height
      ELSE IF (keyword(1:4).EQ.'HRUG') THEN
        z0k(nsur)=r_val(line,1,nchar,1,ierr)
        IF(ierr.ne.0)then
          WRITE(lout,*)' --error--  while reading in foretcon.sur :'
          WRITE(lout,*)' line= ',line
          WRITE(lout,*)'            first value not read correctly'
          STOP
        ENDIF
        GOTO 93

* keyword unknown : stop
      ELSE
        WRITE(lout,*)' --error--  while reading in foretcon.sur'
        WRITE(lout,*)'            keyword unkown: ',keyword
        STOP
      ENDIF

103   CONTINUE
      CLOSE (12)

**************************************************
* end of data to be read                         *
**************************************************

* put together emission that are of similar nature

* The species below have emission factors that do not depend on time
* This species are stored using directly the ID number of the species
* in the mechanism
      nenott=0
      DO i=1,maxem
        idenott(i)=0
        DO j=1,msur
          enott(i,j)=0.
        ENDDO
      ENDDO

C * SO2 emission
C       namesp='GSO2 '
C       CALL akspnum(namesp,chrsp,numsp,isp)
C       nenott=nenott+1
C       idenott(nenott)=isp
C       IF (isp.eq.0)then
C          WRITE(lout,*)' --error--  in subroutine datsuf'
C          WRITE(lout,*)' species unkown=SO2'
C          STOP
C       ENDIF
C       enott(nenott,1)=eso2ur

C * CH4 emission
C       namesp='GCH4 '
C       CALL akspnum(namesp,chrsp,numsp,isp)
C       nenott=nenott+1
C       idenott(nenott)=isp
C       IF (isp.eq.0)then
C          WRITE(lout,*)' --error--  in subroutine datsuf'
C          WRITE(lout,*)' species unkown=CH4'
C          STOP
C       ENDIF
C       enott(nenott,1)=ech4ur
C       enott(nenott,2)=ech4cu
C       enott(nenott,3)=ech4ff
C       enott(nenott,4)=ech4fc

C * calcium particulaire
C       IF (iscape.eq.1) THEN
C         namesp='ACA'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         nenott=nenott+1
C         idenott(nenott)=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=ACA'
C            STOP
C         ENDIF
C         enott(nenott,1)=ecaur
C         enott(nenott,2)=ecacu
C         enott(nenott,3)=ecaff
C         enott(nenott,4)=ecafc
C       ENDIF

C * potassium particulaire
C       IF (iscape.eq.1) THEN
C         namesp='AK'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         nenott=nenott+1
C         idenott(nenott)=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=AK'
C            STOP
C         ENDIF
C         enott(nenott,1)=ekur
C         enott(nenott,2)=ekcu
C         enott(nenott,3)=ekff
C         enott(nenott,4)=ekfc
C       ENDIF

C * chlorure particulaire
C       IF (iscape.eq.1) THEN
C         namesp='AHCL'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         nenott=nenott+1
C         idenott(nenott)=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=AHCL'
C            STOP
C         ENDIF
C         enott(nenott,1)=eclur
C         enott(nenott,2)=eclcu
C         enott(nenott,3)=eclff
C         enott(nenott,4)=eclfc
C       ENDIF

C * emission de sodium particulaire
C       IF (iscape.eq.1) THEN
C         namesp='ANA'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         nenott=nenott+1
C         idenott(nenott)=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=ANA'
C            STOP
C         ENDIF
C         enott(nenott,1)=enaur
C         enott(nenott,2)=enacu
C         enott(nenott,3)=enaff
C         enott(nenott,4)=enafc
C       ENDIF

C * emission de magnesium particulaire
C       IF (iscape.eq.1) THEN
C         namesp='AMG'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         nenott=nenott+1
C         idenott(nenott)=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=AMG'
C            STOP
C         ENDIF
C         enott(nenott,1)=emgur
C         enott(nenott,2)=emgcu
C         enott(nenott,3)=emgff
C         enott(nenott,4)=emgfc
C       ENDIF

C * put together emission of NH3
C       IF (iscape.eq.1) THEN
C         namesp='NH3'
C         CALL akspnum(namesp,chrsp,numsp,isp)
C         idnh3=isp
C         IF (isp.eq.0)then
C            WRITE(lout,*)' --error--  in subroutine datsuf'
C            WRITE(lout,*)' species unkown=NH3'
C            STOP
C         ENDIF
C         enh3(1)=enh3ur
C         enh3(2)=enh3cu
C         enh3(3)=enh3ff
C         enh3(4)=enh3fc
C       ENDIF

C * standart emission of isoprene, alpha-pinene and beta-pinene
C * data are converted from microg m-2 h-1 to molecule cm-2 s-1
C c      namesp='ISO '
C c      CALL akspnum(namesp,chrsp,numsp,idisop)
C c      IF (idisop.eq.0) THEN
C c         WRITE(lout,*)' --error--  in subroutine datsuf'
C c         WRITE(lout,*)' species unkown=C5H8'
C c         STOP
C c      ENDIF
C       namesp='API '
C       CALL akspnum(namesp,chrsp,numsp,idapin)
C       IF (idapin.eq.0) THEN
C          WRITE(lout,*)' --error--  in subroutine datsuf'
C          WRITE(lout,*)' species unkown=APINENE'
C          STOP
C       ENDIF
C       namesp='GBPIN '
C       CALL akspnum(namesp,chrsp,numsp,idbpin)
C       IF (idbpin.eq.0) THEN
C          WRITE(lout,*)' --error--  in subroutine datsuf'
C          WRITE(lout,*)' species unkown=BPINENE'
C          STOP
C       ENDIF

C * culture
C       DO i=1,nshccu
C         IF (idehccu(i).eq.idisop) THEN
C           eisop(2)=ehccu*ehcdatcu(i,1)
C           eisop(2)=(eisop(2)*6.02E23)/(3600*1E4*1E6*ehcdatcu(i,2))
C           idehccu(i)=0
C         ENDIF
C         IF (idehccu(i).eq.idapin) THEN
C           eapin(2)=ehccu*ehcdatcu(i,1)
C           eapin(2)=(eapin(2)*6.02E23)/(3600*1E4*1E6*ehcdatcu(i,2))
C           idehccu(i)=0
C         ENDIF
C         IF (idehccu(i).eq.idbpin) THEN
C           ebpin(2)=ehccu*ehcdatcu(i,1)
C           ebpin(2)=(ebpin(2)*6.02E23)/(3600*1E4*1E6*ehcdatcu(i,2))
C           idehccu(i)=0
C         ENDIF
C       ENDDO
C * check that all species are known
C       DO i=1,nshccu
C         IF (idehccu(i).ne.0) THEN
C           WRITE(lout,*)' --error--  in subroutine datsuf'
C           WRITE(lout,*)' species with ID number=',idehccu(i)
C           WRITE(lout,*)' can not be taken into account in'
C           WRITE(lout,*)' the -culture- environment'
C           WRITE(lout,*)' please change the programm accordingly'
C           STOP
C         ENDIF
C       ENDDO

C * foret leaf
C       DO i=1,nshcff
C         IF (idehcff(i).eq.idisop) THEN
C           eisop(3)=ehcff*ehcdatff(i,1)
C           eisop(3)=(eisop(3)*6.02E23)/(3600*1E4*1E6*ehcdatff(i,2))
C           idehcff(i)=0
C         ENDIF
C         IF (idehcff(i).eq.idapin) THEN
C           eapin(3)=ehcff*ehcdatff(i,1)
C           eapin(3)=(eapin(3)*6.02E23)/(3600*1E4*1E6*ehcdatff(i,2))
C           idehcff(i)=0
C         ENDIF
C         IF (idehcff(i).eq.idbpin) THEN
C           ebpin(3)=ehcff*ehcdatff(i,1)
C           ebpin(3)=(ebpin(3)*6.02E23)/(3600*1E4*1E6*ehcdatff(i,2))
C           idehcff(i)=0
C         ENDIF
C       ENDDO
C * check that all species are known
C       DO i=1,nshcff
C         IF (idehcff(i).ne.0) THEN
C           WRITE(lout,*)' --error--  in subroutine datsuf'
C           WRITE(lout,*)' species with ID number=',idehcff(i)
C           WRITE(lout,*)' can not be taken into account in'
C           WRITE(lout,*)' the -foret feuillu- environment'
C           WRITE(lout,*)' please change the programm accordingly'
C           STOP
C         ENDIF
C       ENDDO

C * foret conifere
C       DO i=1,nshcfc
C         IF (idehcfc(i).eq.idisop) THEN
C           eisop(4)=ehcfc*ehcdatfc(i,1)
C           eisop(4)=(eisop(4)*6.02E23)/(3600*1E4*1E6*ehcdatfc(i,2))
C           idehcfc(i)=0
C         ENDIF
C         IF (idehcfc(i).eq.idapin) THEN
C           eapin(4)=ehcfc*ehcdatfc(i,1)
C           eapin(4)=(eapin(4)*6.02E23)/(3600*1E4*1E6*ehcdatfc(i,2))
C           idehcfc(i)=0
C         ENDIF
C         IF (idehcfc(i).eq.idbpin) THEN
C           ebpin(4)=ehcfc*ehcdatfc(i,1)
C           ebpin(4)=(ebpin(4)*6.02E23)/(3600*1E4*1E6*ehcdatfc(i,2))
C           idehcfc(i)=0
C         ENDIF
C       ENDDO
C * check that all species are known
C       DO i=1,nshcfc
C         IF (idehcfc(i).ne.0) THEN
C           WRITE(lout,*)' --error--  in subroutine datsuf'
C           WRITE(lout,*)' species with ID number=',idehcfc(i)
C           WRITE(lout,*)' can not be taken into account in'
C           WRITE(lout,*)' the -foret conifere- environment'
C           WRITE(lout,*)' please change the programm accordingly'
C           STOP
C         ENDIF
C       ENDDO

      END

