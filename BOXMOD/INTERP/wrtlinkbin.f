      SUBROUTINE wrtlinkbin(llink,numsp,numre,num_n,num_m,
     &            numfo,numhv,numcvar,numextra,numo2,nummeo2,
     &            numain,numaou,numwin,numwou, nauxpar,
     &            numstoi,numiso,mx12stoi,mx1stoi,mx2stoi,
     &            chrsp,id_n,id_m,
     &            idfo,idhv,idcvar,idextra,ido2,idmeo2,
     &            idiso,idain,idaou,idwin,idwou,
     &            idrestoi,restoicf,idpdstoi,pdstoicf,
     &            arrhcf,focf,hvcf,hvfact,cvarcf,extracf, 
!     &            isocf,aincf,aoucf,woucf,wincf,
     &            isocf,aoucf,woucf,wincf,
     &            nrpero,idreacro2,nrdimer,idreacdimer,wmol)

      USE akparameter_module
      IMPLICIT NONE

      INTEGER llink
      INTEGER numsp, numre, num_n, numo2, nummeo2,numiso
      INTEGER num_m, numfo, numhv, numcvar, numextra
      INTEGER mx12stoi, mx1stoi, mx2stoi, mx3stoi
      INTEGER numain, numaou, numwin, numwou
      INTEGER nauxpar(maxaux)
      INTEGER numstoi(maxre,2)
      INTEGER idrestoi(maxre,mxleft)
      INTEGER idpdstoi(maxre,mxright)
      INTEGER id_n(maxre), id_m(max_m)
      INTEGER idfo(maxfo,3)
      INTEGER idhv(maxhv), idcvar(maxcvar)
      INTEGER ido2(maxo2)
      INTEGER idiso(maxiso)
      INTEGER idmeo2(mxrpero)
      INTEGER idextra(maxextra)
      INTEGER idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      !INTEGER nrpero(mxrpero)
      INTEGER nrpero(maxro2)
      INTEGER idreacro2(mxrpero,maxro2)
      !INTEGER nrdimer(mxrdimer)
      INTEGER nrdimer(maxdimer)
      INTEGER idreacdimer(mxrdimer,maxdimer)

      REAL focf(maxaux+3,maxfo)
      REAL extracf(maxaux,maxextra)
! NB: must use compilation option real-8
!      DOUBLE PRECISION extracf(maxaux,maxextra)
      REAL isocf(maxaux,maxiso)
      REAL restoicf(maxre,mxleft)
      REAL pdstoicf(maxre,mxright)
      REAL arrhcf(maxre,3)
      REAL hvcf(maxhv), hvfact(maxhv), cvarcf(maxcvar)
      REAL wmol(maxsp)
      REAL aoucf(2,maxt),woucf(3,maxt),wincf(3,maxt)

      CHARACTER*(maxlsp)  chrsp(maxsp)

      INTEGER i,ire,k


* ----------------------
* CREATING LINK-FILE
* ----------------------

      WRITE(6,*) 'in subroutine and writing the output ...'

      OPEN(llink,file='outdat.akli',form='UNFORMATTED')

      WRITE(llink)maxlsp,numsp,numre,num_n,
     &            num_m,numfo,numhv,numcvar,numextra,numo2,nummeo2,
     &            numain,numaou,numwin,numwou,
     &            numiso,mx12stoi,mx1stoi,mx2stoi,
     &            maxaux,
     &            (nauxpar(i),i=1,maxaux)

      WRITE(llink)(chrsp(i),i=1,numsp)

      WRITE(llink)(id_n(i),i=1,num_n),(id_m(i),i=1,num_m),
     &            ((idfo(i,k),k=1,3),i=1,numfo),
     &            (idhv(i),i=1,numhv),(idcvar(i),i=1,numcvar),
     &            (idextra(i),i=1,numextra)

      WRITE(llink) (ido2(i),i=1,numo2)

      WRITE(llink) (idmeo2(i),i=1,nummeo2)
      WRITE(llink) (idiso(i),i=1,numiso)

      WRITE(llink) (idain(i),i=1,numain)
      WRITE(llink) (idaou(i),i=1,numaou)
      WRITE(llink) (idwin(i),i=1,numwin)
      WRITE(llink) (idwou(i),i=1,numwou)

      WRITE(llink)((arrhcf(ire,k),k=1,3),ire=1,numre)
      WRITE(llink)((numstoi(ire,k),k=1,2),ire=1,numre)

      WRITE(llink)((idrestoi(ire,i),restoicf(ire,i),i=1,numstoi(ire,1)),
     &             (idpdstoi(ire,i),pdstoicf(ire,i),i=1,numstoi(ire,2)),
     &              ire=1,numre)

      WRITE(llink)((focf(i,k),i=1,maxaux+3),k=1,numfo)

      WRITE(llink)(hvcf(k),k=1,numhv)
      WRITE(llink)(hvfact(k),k=1,numhv)

      WRITE(llink)(cvarcf(k),k=1,numcvar)

      WRITE(llink)((extracf(i,k),i=1,maxaux),k=1,numextra)
      WRITE(llink)((isocf(i,k),i=1,maxaux),k=1,numiso)

!      WRITE(llink)((aincf(i,k),i=1,3),k=1,numain)
      WRITE(llink)((aoucf(i,k),i=1,2),k=1,numaou)
      WRITE(llink)((woucf(i,k),i=1,3),k=1,numwou)
      WRITE(llink)((wincf(i,k),i=1,3),k=1,numwin)

      WRITE(llink) (nrpero(k),k=1,maxro2)
      WRITE(llink) ((idreacro2(i,k),i=1,nrpero(k)),k=1,maxro2)

      WRITE(llink) (nrdimer(k),k=1,maxdimer)
      WRITE(llink) ((idreacdimer(i,k),i=1,nrdimer(k)),k=1,maxdimer)

      WRITE(llink)(wmol(i),i=1,numsp)

      CLOSE(llink)

* -------------------------------
      RETURN

      END SUBROUTINE wrtlinkbin

