* INCA : subroutine chkaux
*
* Purpose : 
* check the auxiliary information provided to a given reaction.
*
*
* INPUT :
*  maxaux    : maximum number of data that can be given between slashes
*              (and maximum number of keyword available)
*  maxfo     : maximum number of fall off reactions
*  maxre     : maximum number of reactions
*  maxextra  : maximum number of extra reactions
*  maxcvar   : maximum number of cvar reactions
*  max_m     : maximum number of TBODY (formerly '+M') reactions
*  maxhv     : maximum number of hv reactions
*  numre     : total number of reactions
*  xauxcf(i,j): jth data for the ith keyword 
*
* INPUT/OUTPUT
*  numfo     : number of fo reactions
*  numcvar   : number of cvar reactions
*  num_m     : number of TBODY ('+M') reactions
*  numextra  : number of EXTRA reactions
*  numhv     : number of HV reactions
*  numo2     : number of O2 reaction
*  lout      : unit number for the output file
*  idfo(i)   : reaction number of the ith fall off reaction
*  idhv(i)   : reaction number of the ith hv reaction
*  idextra(i): reaction number of the ith extra reaction
*  idcvar(i) : reaction number of the ith cvar reaction
*  id_m(i)   : reaction number of the ith TBODY ('+M') reaction
*  ido2(i)   : reaction number of
*  cvarcf(i) : label for the ith cvar reaction
*  extracf(j,i): jth data for the ith extra reaction
*  isocf(j,i): jth data for the ith isomerisation reaction
*  focf(j,i) : jth data for the ith fo reaction
*  hvcf(i)   : label for the ith hv reaction
*  hvfact(i) : coefficient for the ith hv reaction
*  locheck(i): logical check for the ith reaction
*  lo_m      : logical for TBODY ('+M') reaction
*  lofo      : logical for fall off reaction
*  lohv      : logical for hv reaction
*  loaux     : logical for auxiliary information (if true=>data read)
*  locvar    : logical for CVAR reaction
*  loextra,  : logica for extra reaction
*  lo_iso    : logical for izomerisation reaction (with Vereecken SAR)
*  itype(i)  : type for the ith reaction (1=M,3=FO,4=O2,5=PERO,6=HV,
*              7=EXTRA,8=CVAR,9=AOU,10=AIN,11=WOU,12=WIN,13=ISOM,14=ISOM)
* nrpero(i)  : total number of reaction involving PEROi
* idreacro2(j,i) : number of jth reaction involving PEROi
* nrdimer(i)  : total number of reaction involving DIM_i
* idreacdimer(j,i) : number of jth reaction involving DIM_i
*
*
* OUTPUT
*  lostop    : logical to stop program if error
**********************************************************************
      SUBROUTINE chkaux (
     1   maxaux, maxfo, maxre, 
     2   maxextra, maxcvar,
     3   max_m, maxhv, maxo2, mxrpero, maxro2,maxiso,
     3   maxt,mxrdimer,maxdimer,
     4   numre, numfo, numcvar, num_m,
     5   numextra, numhv, numo2, nummeo2,numiso,
     5   numain, numaou, numwin, numwou,
     6   lout,
     7   idfo, idhv,
     9   idextra, idcvar, id_m, ido2, idmeo2,idiso,
     9   idain,idaou,idwin,idwou,
     9   cvarcf, extracf,isocf,
     8   focf, hvcf, hvfact,
     8   woucf,wincf,
     &   nrpero,idreacro2,
     &   nrdimer,idreacdimer,
     7   lostop, locheck, lo_m, lofo,lo_iso,
     6   lohv,loaux,locvar,loextra,lo_o2,lo_meo2,lopero,ipero,
     6   lodimer,idimer,
     7   loain,loaou,lowin,lowou,
     l   itype, xauxcf)
      IMPLICIT NONE

* input
      INTEGER maxaux, maxfo, maxre
      INTEGER maxextra
      INTEGER maxcvar, max_m, maxhv, maxo2, mxrpero, maxro2,maxiso
      INTEGER mxrdimer,maxdimer
      INTEGER maxt
      INTEGER numre, numfo, numcvar
      INTEGER num_m, numextra, numhv, numo2, nummeo2,numiso
      INTEGER numain, numaou, numwin, numwou
      INTEGER lout
      INTEGER idfo(maxfo,3)
      INTEGER idhv(maxhv)
      INTEGER ido2(maxo2)
      INTEGER idiso(maxiso)
      INTEGER idmeo2(mxrpero)
      INTEGER idextra(maxextra) 
      INTEGER idcvar(maxcvar)
      INTEGER id_m(max_m)
      INTEGER idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      REAL xauxcf(0:maxaux,maxaux)
!      DOUBLE PRECISION xauxcf(0:maxaux,maxaux)

* input/output
      REAL cvarcf(maxcvar)
      REAL extracf(maxaux,maxextra)
! USE COMPILATION OPTION real-8
!      DOUBLE PRECISION extracf(maxaux,maxextra)
      REAL isocf(maxaux,maxiso)
      REAL focf(maxaux+3,maxfo)
      REAL hvcf(maxhv)
      REAL hvfact(maxhv)
      REAL woucf(1,maxt),wincf(1,maxt)
      LOGICAL loaux,lofo,lohv,locvar,loextra,lo_o2,lopero,lo_meo2
      LOGICAL lo_iso,lodimer
      INTEGER ipero,idimer
      LOGICAL loain,loaou,lowin,lowou
      LOGICAL lo_m, lostop, locheck(maxre)
      INTEGER itype(maxre)
      !INTEGER nrpero(mxrpero),nrdimer(mxrdimer)
      INTEGER nrpero(maxro2),nrdimer(maxdimer)
      INTEGER idreacro2(mxrpero,maxro2),idreacdimer(mxrdimer,maxdimer)

* local
      INTEGER i

* ===================================================================
* If auxiliary information was not provided (loaux=false) then the
* reaction cannot be a falloff, hv, cvar, extra reaction or 
* isomerisation. Case
* allowed is only a regular TBODY ('+M') reaction or the reaction with oxygen.
* If simple thermal reaction then exit 8000. Note that this point is 
* not checked.
* ===================================================================

      IF (.not.loaux) THEN
        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a fall-off reaction must have',
     &               '  auxiliary information'
          GOTO 9000        

        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a "cvar" reaction must have',
     &                 '  auxiliary information'
          GOTO 9000        

        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a "extra" reaction must have',
     &                 ' auxiliary information'
          GOTO 9000        

        ELSE IF (lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a "hv" reaction must have',
     &                 ' auxiliary information'
          GOTO 9000        

        ELSE IF (lo_iso) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a ISOM reaction must have',
     &                 ' auxiliary information'
          GOTO 9000        

        ELSE IF (lo_m) THEN
          IF (num_m.ge.max_m) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations',
     &                   ' with "M" is greater than max_m=',max_m
            GOTO 9000
          ELSE
            itype(numre)=1
            num_m=num_m+1
            id_m(num_m)=numre
            GOTO 8000
          ENDIF

        ELSE IF (lo_o2) THEN
          IF (numo2.ge.maxo2) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations with',
     &                   '  "OXYGEN" is greater than maxo2=',maxo2
            GOTO 9000
          ELSE
            itype(numre)=4
            numo2=numO2+1
            ido2(numo2)=numre
            GOTO 8000
          ENDIF

        ELSE IF (lo_meo2) THEN
          IF (nummeo2.ge.mxrpero) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations with',
     &                   '  "MEPERO" is greater than mxrpero=',mxrpero
            GOTO 9000
          ELSE
            itype(numre)=5
            nummeo2=nummeo2+1
            idmeo2(nummeo2)=numre
            GOTO 8000
          ENDIF

        ELSE IF (lopero) THEN
          IF ( (ipero.lt.1).or.(ipero.gt.9) ) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   ID number with pero not set'
            GOTO 9000
          ENDIF
          IF (nrpero(ipero).gt.mxrpero) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations with',
     &                  'PERO',ipero,' is greater than mxpero=',mxrpero
            GOTO 9000
          ELSE
            nrpero(ipero)=nrpero(ipero)+1
            idreacro2(nrpero(ipero),ipero)=numre
            itype(numre)=5
            GOTO 8000
          ENDIF

        ELSE IF (lodimer) THEN
          IF ( (idimer.lt.1).or.(idimer.gt.4) ) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   ID number with dimer not set'
            GOTO 9000
          ENDIF
          IF (nrdimer(idimer).gt.mxrdimer) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations with',
     &             'DIM_',idimer,' is greater than mxrdimer=',mxrdimer
            GOTO 9000
          ELSE
            nrdimer(idimer)=nrdimer(idimer)+1
            idreacdimer(nrdimer(idimer),idimer)=numre
            itype(numre)=14
            GOTO 8000
          ENDIF


! mass transfer to wall or aerosols 
! remember the "type" table : 
!  itype(i)  : type for the ith reaction (1=M,3=FO,4=O2,5=PERO,6=HV,
!              7=EXTRA,8=CVAR,9=AOU,10=AIN,11=WOU,12=WIN,13=ISOM,14=DIMER)
        ELSE IF (loain) THEN
          IF (numain.ge.maxt) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   number of equations with',
     &                   '  "AIN" is greater than maxt=',maxt
            GOTO 9000
          ELSE
            itype(numre)=10
            numain=numain+1
            idain(numain)=numre
            GOTO 8000
          ENDIF

        ELSE IF (lowin) THEN
            WRITE(lout,*)
            WRITE(lout,*)'  --error--   a "WIN" reaction must have',
     &                   ' auxiliary information'
            GOTO 9000

        ELSE IF (lowou) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   a "WOU" reaction must have',
     &                 ' auxiliary information'
          GOTO 9000        
         ELSE
           GOTO 8000
         ENDIF
       ENDIF

* ===================================================================
* following statement are read only if auxiliary infomation was 
* provided for the reaction. The program check the consistency
* of the information provided and update the table.
* 
* REMEMBER : link between keywords and number in xauxcf
*     1 = FALLOFF    |   2 = not used
*     3 = not used   |   4 = not used
*     5 = HV         |   6 = EXTRA
*     7 = CVAR       |   8 = AOU
*     9 = WOU        |  10 = WIN
*    11 = ISOM 
!*    11 = AIN        |  12 = ISOM
* ===================================================================

* -------------------
* PRELIMINARY CHECKS
* -------------------

* check that a reaction was written first (before any aux. info.)
      IF (numre.EQ.0) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   auxiliary information',
     &               ' must be after reaction information.'
        STOP
      ENDIF

* keyword LOW (2), SRI (3) and LT (4) cannot be used anymore
      IF (xauxcf(0,2).NE.0.0) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   keyword LOW (2) cannot be used'
        STOP
      ENDIF
      IF (xauxcf(0,3).NE.0.0) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   keyword SRI (3) cannot be used'
        STOP
      ENDIF
      IF (xauxcf(0,4).NE.0.0) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   keyword LT (4) cannot be used'
        STOP
      ENDIF

* At least one of the xauxcf(0,*) must be a non-zero value
      IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &    xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &    xauxcf(0,7).EQ.0.0.AND.xauxcf(0,8).EQ.0.0.AND.
     &    xauxcf(0,9).EQ.0.0.AND.xauxcf(0,10).EQ.0.0.AND.
     &    xauxcf(0,11).EQ.0.0)THEN
!     &    xauxcf(0,11).EQ.0.0.AND.xauxcf(0,12).EQ.0.0)THEN
        WRITE(lout,*)
        WRITE(lout,*)'   --error--     third body species can only',
     &               ' be used within fall-off reactions'
        WRITE(lout,*)'                 or in conection with the',
     &               ' keyword "extra"'
        STOP
      ENDIF

* -----------------------------
* FALL OFF REACTION (formerly LOW, TROE)
* -----------------------------

      IF (xauxcf(0,1).EQ.1.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &          xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &          xauxcf(0,7).EQ.0.0)THEN

* check that the only keyword in the reaction is FALLOFF (previously '(+M)').
* Number of fo reactions must not exceed maxfo. If any problem => goto 9000
        IF (.not.lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   key-word "FALLOFF"',
     &                 '        can only be used in a fall-off reaction'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "hv" can not be used',
     &                 '                in a fall-off reaction'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "cvar" can not be used',
     &                 '                in a fall-off reaction'
          GOTO 9000
        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "extra" can not be used',
     &                 '               in a fall-off reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   an additional "m" can not',
     &                 '               be used in a fall-off reaction'
          GOTO 9000
        ELSE IF (numfo.ge.maxfo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of fall-off equations',
     &                 '               is greater than maxfo=',maxfo
          GOTO 9000

* load the data in the tables holding the info for fall off reactions
* and jump to 8000 (no problem in auxiliary info)
        ELSE
          itype(numre)=3
          numfo=numfo+1
          idfo(numfo,1)=numre
          idfo(numfo,2)=2
          DO i=1,3
            focf(i,numfo)=xauxcf(i,1)
          ENDDO
          DO i=1,maxaux
            focf(i+3,numfo)=xauxcf(i,2)
          ENDDO
          GOTO 8000
        ENDIF

* -----------------------------
* PHOTOLYTIC REACTION (HV)
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.1.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &         xauxcf(0,7).EQ.0.0) THEN

* check that the only keyword found in the reaction is HV. Number
* of hv reaction must not exceed maxhv. If any problem => goto 9000
        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "HV" can not be used within',
     &                 '                a fall-off reaction'
          GOTO 9000
        ELSE IF (.not.lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "HV"',
     &                 '               in order to include "HV" in the'
          WRITE(lout,*)'               auxiliary information'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "CVAR" and "HV" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "EXTRA" and "HV" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "M" and "HV" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (numhv.ge.maxhv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of "HV" equations',
     &                 '               is greater than maxhv=',maxhv
          GOTO 9000

* load the data in the tables holding the info for hv reactions
* and jump to 8000 (no problem in auxiliary info)
        ELSE
!          WRITE(6,*) numre
          itype(numre)=6
          numhv=numhv+1
          idhv(numhv)=numre
          hvcf(numhv)=xauxcf(1,5)
          hvfact(numhv)=xauxcf(2,5)
          GOTO 8000
        ENDIF

* -----------------------------
* EXTRA REACTION
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,3).EQ.0.0.AND.xauxcf(0,4).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.1.0.AND.
     &         xauxcf(0,7).EQ.0.0) THEN

* check that the only keyword found in the reaction is EXTRA. Number
* of EXTRA reaction must not exceed maxextra. If any problem => goto 9000
        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "extra" can not be used',
     &                 ' in a fall-off reaction'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "hv" and "extra" can not',
     &                 ' be used within the same reaction'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "cvar" and "extra" can not',
     &                 ' be used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "m" and "extra" can not',
     &                 ' be used within the same reaction'
          GOTO 9000
        ELSE IF (.not.loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "extra"',
     &                 ' in order to include "extra"'
          WRITE(lout,*)'               in the auxiliary information'
          GOTO 9000
        ELSE IF (numextra.ge.maxextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of equations',
     &                 ' with "extra" is greater than maxextra=',
     &                  maxextra
          GOTO 9000

* load the data in the tables holding the info for extra reactions
* and jump to 8000 (no problem in auxiliary info)
        ELSE
          itype(numre)=7
          numextra=numextra+1
          idextra(numextra)=numre
          DO i=1,maxaux
              extracf(i,numextra)=xauxcf(i,6)
          ENDDO
          GOTO 8000
        ENDIF

* -----------------------------
* CVAR REACTION
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &         xauxcf(0,7).EQ.1.0) THEN

* check that the only keyword found in the reaction is CVAR. Number
* of CVAR reaction must not exceed maxextra. If any problem => goto 9000
        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "cvar" can not be used',
     &                 ' within a fall-off reaction'
          GOTO 9000
        ELSE IF (.not.locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "cvar"',
     &                 ' in order to include "cvar" in the'
          WRITE(lout,*)'               auxiliary information'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "hv" and "cvar" can not be',
     &                 ' used in the same reaction'
          GOTO 9000
        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "extra" and "cvar" can not',
     &                 ' be used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "m" and "cvar" can not be',
     &                 ' used in the same reaction'
          GOTO 9000
        ELSE IF (numcvar.ge.maxcvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of "cvar"',
     &                 ' equations is greater than maxcvar=',maxcvar
          GOTO 9000

* load the data in the tables holding the info for cvar reactions
* and jump to 8000 (no problem in auxiliary info)
        ELSE
          itype(numre)=8
          numcvar=numcvar+1
          idcvar(numcvar)=numre
          cvarcf(numcvar)=xauxcf(1,7)
          GOTO 8000
        ENDIF


!* -----------------------------
!* AIN REACTION (xauxcf=11,itype=10)
! mass transfer to wall or aerosols 
* NO LONGER HAS COEFFICIENTS
* -----------------------------
* AOU REACTION (xauxcf = 8,itype = 9)
! mass transfer to wall or aerosols 
* NO LONGER HAS COEFFICIENTS
* -----------------------------

* -----------------------------
* WOU REACTION (xauxcf = 9, itype = 11)
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &         xauxcf(0,7).EQ.0.0.AND.xauxcf(0,8).EQ.0.0.AND.
     &         xauxcf(0,9).EQ.1.0.AND.xauxcf(0,10).EQ.0.0) THEN

        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "WOU" can not be used within',
     &                 '                a fall-off reaction'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)'   --error--   "WOU" can not be used within',
     &                 '                a fall-off reaction'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "CVAR" and "WOU" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "EXTRA" and "WOU" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "M" and "WOU" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (loaou) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "WOU" and "AOU" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (.not.lowou) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "WOU"',
     &                 '               in order to include "WOU" in the'
          WRITE(lout,*)'               auxiliary information'
          GOTO 9000
        ELSE IF (numwou.ge.maxt) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of "WOU" equations',
     &                 '               is greater than maxt=',maxt
          GOTO 9000

* load the data in the tables holding the info for hv reactions
* and jump to 8000 (no problem in auxiliary info)
! remember the "type" table : 
!  itype(i)  : type for the ith reaction (1=M,3=FO,4=O2,5=PERO,6=HV,
!              7=EXTRA,8=CVAR,9=AOU,10=AIN,11=WOU,12=WIN,13=ISOM)
        ELSE
!          WRITE(6,*) numre
          itype(numre)=11
          numwou=numwou+1
          idwou(numwou)=numre
          woucf(1,numwou)=xauxcf(1,9)
          GOTO 8000
        ENDIF

* -----------------------------
* WIN REACTION (xauxcf = 10, itype = 12)
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &         xauxcf(0,7).EQ.0.0.AND.xauxcf(0,8).EQ.0.0.AND.
     &         xauxcf(0,9).EQ.0.0.AND.xauxcf(0,10).EQ.1.0) THEN

        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "WIN" can not be used within',
     &                 '                a fall-off reaction'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)'   --error--   "WIN" can not be used within',
     &                 '                a fall-off reaction'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "CVAR" and "WIN" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (loextra) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "EXTRA" and "WIN" can not be',
     &                 '                used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "M" and "WIN" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (loaou) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "WIN" and "AOU" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (lowou) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "WIN" and "WOU" can not be used',
     &                 '                in the same reaction'
          GOTO 9000
        ELSE IF (.not.lowin) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "WIN"',
     &                 '               in order to include "WIN" in the'
          WRITE(lout,*)'               auxiliary information'
          GOTO 9000
        ELSE IF (numwin.ge.maxt) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of "WIN" equations',
     &                 '               is greater than maxt=',maxt
          GOTO 9000

* load the data in the tables holding the info for hv reactions
* and jump to 8000 (no problem in auxiliary info)
! remember the "type" table : 
!  itype(i)  : type for the ith reaction (1=M,3=FO,4=O2,5=PERO,6=HV,
!              7=EXTRA,8=CVAR,9=AOU,10=AIN,11=WOU,12=WIN)
        ELSE
!          WRITE(6,*) numre
          itype(numre)=12
          numwin=numwin+1
          idwin(numwin)=numre
          wincf(1,numwin)=xauxcf(1,10)
          GOTO 8000
        ENDIF

* -----------------------------
* ISOMERISATION REACTION
* -----------------------------

      ELSE IF (xauxcf(0,1).EQ.0.0.AND.xauxcf(0,2).EQ.0.0.AND.
     &         xauxcf(0,3).EQ.0.0.AND.xauxcf(0,4).EQ.0.0.AND.
     &         xauxcf(0,5).EQ.0.0.AND.xauxcf(0,6).EQ.0.0.AND.
     &         xauxcf(0,7).EQ.0.0.AND.xauxcf(0,8).EQ.0.0.AND.
     &         xauxcf(0,9).EQ.0.0.AND.xauxcf(0,10).EQ.0.0.AND.
     &         xauxcf(0,11).EQ.1.0) THEN
!     &         xauxcf(0,11).EQ.0.0.AND.xauxcf(0,12).EQ.1.0) THEN

* check that the only keyword found in the reaction is ISOM. Number
* of ISOM reaction must not exceed maxiso. If any problem => goto
* 9000
        IF (lofo) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "isom" can not be used',
     &                 ' in a fall-off reaction'
          GOTO 9000
        ELSE IF (lohv) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "hv" and "isom" can not',
     &                 ' be used within the same reaction'
          GOTO 9000
        ELSE IF (locvar) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "cvar" and "isom" can not',
     &                 ' be used in the same reaction'
          GOTO 9000
        ELSE IF (lo_m) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   "m" and "isom" can not',
     &                 ' be used within the same reaction'
          GOTO 9000
        ELSE IF (.not.lo_iso) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   reaction must have a "ISOM"',
     &                 ' in order to include "ISOM"'
          WRITE(lout,*)'               in the auxiliary information'
          GOTO 9000
        ELSE IF (numiso.ge.maxiso) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   number of equations',
     &                 ' with "ISOM" is greater than maxiso=',
     &                  maxiso
          GOTO 9000
* load the data in the tables holding the info for isom reactions
* and jump to 8000 (no problem in auxiliary info)
        ELSE
          itype(numre)=13
          numiso=numiso+1
          idiso(numiso)=numre
          DO i=1,maxaux
            !isocf(i,numiso)=xauxcf(i,12)
            isocf(i,numiso)=xauxcf(i,11)
          ENDDO
          GOTO 8000
        ENDIF

* -----------------------------
* UNIDENTIFIED REACTION
* -----------------------------

      ELSE
        WRITE(lout,*)
        WRITE(lout,*)'   --error--   an incorrect combination',
     &               '               of key-words is used within the'
        WRITE(lout,*)'               auxiliary information'
        GOTO 9000
      ENDIF

* ===================================================================
* final statements => initialize the logical for the next reaction.
* If 8000 is executed, then no problem was found in the provided 
* auxiliary information. If 9000 is executed => error occured
* ===================================================================

8000  lo_m=.false.
      lofo=.false.
      lohv=.false.
      loaux=.false.
      locvar=.false.
      loextra=.false.
      lo_o2=.false.
      lo_meo2=.false.
      lopero=.false.
      lodimer=.false.
      loain=.false.
      loaou=.false.
      lowin=.false.
      lowou=.false.
      lo_iso=.false.
      RETURN

9000  WRITE(lout,*)
      lo_m=.false.
      lofo=.false.
      lohv=.false.
      loaux=.false.
      locvar=.false.
      loextra=.false.
      lo_o2=.false.
      lo_meo2=.false.
      lopero=.false.
      lodimer=.false.
      loain=.false.
      loaou=.false.
      lowin=.false.
      lowou=.false.
      lo_iso=.false.
      locheck(numre)=.true.
      lostop=.true.

      END
