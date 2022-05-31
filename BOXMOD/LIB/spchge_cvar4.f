***********************************************************************
* Cette routine a pour objet de changer les coef. stoechio. qui       * 
* dependent de la temperature                                         *  
* INPUT :                                                             *
*   lout   : numero de fichier pour ecriture des erreurs              *
*   t      : temperature                                              *
* (in the common block)                                               *
*   valcoe(i,j,k) : valeur des coef. stoech. a la ieme temperature,   *
*                   pour la kieme reaction cvar et pour la jieme      *
*                   espece                                            *
*   numcoe(k)     : nbre de coef. stoe. tabule pour la kieme reaction *
*   nopc(k)       : nbre d'espece compteur de la kieme reaction qui   *
*                   doivent etre reevalue a partir des coef. stoe.    *
*                   associes aux operateurs                           *
*   ndatopc(k,j)  : nbre de donnees pour reevaluer le jieme compteur  *
*                   de la kieme reaction                              *
*   nposopc(k,j,l): position (l) des operateurs pour la reevaluation  *
*                   du jieme compteur de la kieme reaction            *
*                                                                     *
* OUTPUT : (dans akcommon.h)                                          *
*   stoicf(ire,j,2) : coefficient stoechiometrique de la reaction     *
*                     ire, jieme espece                               * 
***********************************************************************
      SUBROUTINE chge_cvar4(maxre,maxcvar,mxright,
     1                      mopc,mpos,nset,maxcoe,
     2                      lout,numcvar,idcvar,numstoi,
     3                      numcoe,nopc,nposopc,ndatopc,valcoe,
     4                      temp,pdstoicf)
      IMPLICIT NONE

*INPUT
      INTEGER  maxre,maxcvar,mxleft,mxright,mopc,mpos,nset,maxcoe
      INTEGER  lout,numcvar,idcvar(maxcvar)
      INTEGER  numstoi(maxre,2)
      REAL     temp
      INTEGER  numcoe(maxcvar), nopc(maxcvar)
      INTEGER  ndatopc(maxcvar,mopc), nposopc(maxcvar,mopc,mpos)
      REAL     valcoe(nset,maxcoe,maxcvar)

* OUTPUT
      REAL     pdstoicf(maxre,mxright)   

* LOCAL
      INTEGER  ii,ire,jj,j,k,kk,nposc,nposo,npos1,npos2

c      write(6,*) 'inside cvar'
c      WRITE(6,*) ' toto --- nopc(1)=',nopc(1)

* test si  la temperature est dans le bon intervalle
      IF (temp.LT.260.0.or.temp.GT.320.0) THEN
        WRITE(lout,'(a55)')
     &  '--error--,temp. must be greater than 260K or less 320K'      
        STOP
      ENDIF

* loop over cvar reaction
      DO 10 ii=1,numcvar
      ire=idcvar(ii)

* change the stoechiometric coefficient
      IF (temp.lt.280.0) THEN
         DO 20 jj=1,numcoe(ii) 
            pdstoicf(ire,jj)=valcoe(1,jj,ii)+
     &      ((valcoe(2,jj,ii)-valcoe(1,jj,ii))/20.)*(temp-260.)
20       CONTINUE

      ELSE IF (temp.lt.300.0) THEN
         DO 30 jj=1,numcoe(ii) 
            pdstoicf(ire,jj)=valcoe(2,jj,ii)+
     &      ((valcoe(3,jj,ii)-valcoe(2,jj,ii))/20.)*(temp-280.) 
30       CONTINUE

      ELSE IF (temp.lt.320.0) THEN
         DO 40 jj=1,numcoe(ii) 
            pdstoicf(ire,jj)=valcoe(3,jj,ii)+
     &      ((valcoe(4,jj,ii)-valcoe(3,jj,ii))/20.)*(temp-300.) 
40       CONTINUE
      ENDIF
      
* overwrite the stoechiometric coefficient for counting species
      IF (nopc(ii).ne.0.) THEN
        DO j=1,nopc(ii)
          nposc=nposopc(ii,j,1)
          pdstoicf(ire,nposc)=0.
          DO k=2,ndatopc(ii,j)+1
            nposo=nposopc(ii,j,k)
            pdstoicf(ire,nposc)=
     &             pdstoicf(ire,nposc)+pdstoicf(ire,nposo)
          ENDDO
        ENDDO
      ENDIF

* change the stoicf for parameter 3
c      npos1=1
c      npos2=1
c      DO 60 kk=1,numstoi(ire,3)
c      IF (idstoi(ire,kk,3).eq.idstoi(ire,npos1,1).AND.
c     &    idstoi(ire,kk,3).eq.idstoi(ire,npos2,2)) THEN
c         stoicf(ire,kk,3)=stoicf(ire,npos2,2)-stoicf(ire,npos1,1)
c         npos1=npos1+1
c         npos2=npos2+1
c      ELSE IF (idstoi(ire,kk,3).eq.idstoi(ire,npos2,2)) THEN
c         stoicf(ire,kk,3)=stoicf(ire,npos2,2)
c         npos2=npos2+1
c      ELSE IF (idstoi(ire,kk,3).eq.idstoi(ire,npos1,1)) THEN
c         stoicf(ire,kk,3)=-stoicf(ire,npos1,1)
c         npos1=npos1+1
c      ELSE 
c         WRITE (lout,*) '--error-- in chge_cvar'
c         STOP
c      ENDIF
c60    CONTINUE

10    CONTINUE
      END 

