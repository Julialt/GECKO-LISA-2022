! INCA : subroutine getreac
!
! Purpose : 
! read a given reaction to find reactants, products, associated
! stoichimetric coefficient, keyword and arrhenius parameter.
!
! INPUT :
!  line      : the line to be read
!  lenline   : maximum length of the line
!  lout      : unit number for the output file
!  maxre     : maximum number of reactions
!  maxsp     : maximum number of species
!  maxstoi   : maximum number of species in a given reaction
!  numre     : total number of reactions
!  numsp     : number of species in the chemical scheme
!  chrsp(i)  : name of the species i
!  sortsp(j) : sorted list of the species
!  isortlk(j): link between the sorted list and the unsorted list
!
! INPUT/OUTPUT
!  numstoi(i,j) : number of species reacting/produced at side j 
!                 for reaction i.
!  idstoi(i,k,j) : id number of the kth species involved in reaction i
!                  at side j.
!  stoicf(i,k,j) : stoi coef of the kth species involved in reaction i
!                  at side j.
!  mx1stoi   : maximum number of reactant found in a given reaction
!             (full chemical scheme considered)
!  mx2stoi   : maximum number of product found in a given reaction
!             (full chemical scheme considered)
!  mx3stoi   : maximum number of distinct species found in a given 
!              reaction (full chemical scheme considered)
!  mx12stoi  : maximum between mx1stoi and mx2stoi
!  lostop    : logical to stop program if error
!  locheck(i): logical check fot ith reaction
!  lo_m      : logical for +M reaction
!  lofo      : logical for fall off reaction
!  lohv      : logical for hv reaction
!  locvar    : logical for CVAR reaction
!  loextra,  : logica for extra reaction
!  lo_o2     : logical for reaction with O2
!  lo_iso    : logical for izomerisation reaction (with Vereecken SAR)
!  itype(i)  : type gor the ith reaction (1=M,3=FO,4=O2, 6=HV,8=CVAR,
!             7=EXTRA)
!  arrhcf(i,j): jth arrhenius coef. for reaction i
!
! CALLED SUBROUTINE   : errline
! CALLED FUNCTION     : r_val
!
!*********************************************************************
      SUBROUTINE getreac (
     1  line, lenline, lout, 
     3  numre, numsp, numstoi, 
     4  idrestoi, idpdstoi, restoicf, pdstoicf, chrsp,
     5  mx12stoi, mx1stoi, mx2stoi, mx3stoi,
     6  lostop, locheck, lo_m, lofo,lo_o2, lo_meo2,lo_iso,
     7  lohv, locvar, loextra, lopero, ipero,lodimer,idimer,
     7  loain,loaou,lowin,lowou,
     8  itype, 
     9  arrhcf
     1  )

      USE sorting, ONLY : srtid_chrsp, find_species_index
      USE akparameter_module
      IMPLICIT NONE
      INCLUDE 'general.h'

! input
      INTEGER lenline, lout
      INTEGER numre, numsp
      CHARACTER*(*) line
      CHARACTER*(*) chrsp(maxsp)

! input/output
      INTEGER mx12stoi,mx1stoi,mx2stoi,mx3stoi
      INTEGER itype(maxre)
      INTEGER numstoi(maxre,2)
      INTEGER idrestoi(maxre,mxleft)
      INTEGER idpdstoi(maxre,mxright)
      REAL arrhcf(maxre,3)
      REAL restoicf(maxre,mxleft)
      REAL pdstoicf(maxre,mxright)
      LOGICAL lostop,locheck(maxre)
      LOGICAL lo_m, lofo, lohv, locvar, loextra, lo_o2, lo_meo2
      LOGICAL loain,loaou,lowin,lowou,lo_iso
      LOGICAL lopero
      INTEGER ipero,idpero
      LOGICAL lodimer
      INTEGER idimer,iddimer


! local
      INTEGER i, j, k, ii, jj, iii, iloc
      INTEGER iterm, ib, ierr, iside
      INTEGER numdata, idlopar, lenmin
      INTEGER ipart
      REAL xcoeff
      LOGICAL lomemo, lonumber, lopar, lonothing
      CHARACTER*(maxlsp) tempsp

      REAL r_val
      INTEGER search
      INTEGER iholdsp(mxright,3), cpiholdsp(mxright,3)
      INTEGER ncoef
      REAL    coefsp(mxright,3), cpcoefsp(mxright,3)
      INTEGER minsp, nspe(3), ipos(2)

      !print*,"-----------getreac----------"
! initialize
      lo_m =.false.
      lofo =.false.
      lohv =.false.
      locvar =.false.
      lomemo =.false.
      lonumber =.false.
      lopar =.false.
      loextra =.false.
      lonothing =.false.
      lo_o2=.false.
      lo_iso=.false.
      lo_meo2=.false.
      loain=.false.
      loaou=.false.
      lowin=.false.
      lowou=.false.
      iterm=lenline
      ib=0
      ierr=0
      iside=1
      numdata=0
      idlopar=0
      ncoef=0
      ipero=-1
      idimer=-1
      iholdsp=0
      cpiholdsp=0
      coefsp=0.
      cpcoefsp=0.

! add reaction 
      numre=numre+1
      itype(numre)=0
      IF (numre.gt.maxre) THEN
        WRITE(lout,*)
        WRITE(lout,*) line
        WRITE(lout,*)'   --error--   number of reactions can',
     &               ' not be greater than maxre=',maxre
        WRITE(lout,*)
        lostop=.true.
        RETURN
      ENDIF

! --------------------------
! READ ARRHENIUS PARAMETER
! --------------------------

! read the line from the end and put the cursor (iterm) 3 block of 
! data before end of line
 
      DO i=lenline,1,-1
        IF (line(i:i).eq.' ') THEN
          IF (lomemo) THEN
            lomemo=.false.
            iterm=i
            numdata=numdata+1
            IF (numdata.eq.3) GOTO 3250
          ENDIF
        ELSE
          lomemo=.true.
        ENDIF
      ENDDO

! If the following statement are executed, expected data for arrhenius
! parameterization was not found => error 

      WRITE(lout,*)
      WRITE(lout,*) line
      WRITE(lout,*)'  --error--  not all reaction coefficients found'
      WRITE(lout,*)
      locheck(numre)=.true.
      lostop=.true.
      RETURN

! read the arrhenius coefficient

3250  CONTINUE
      DO i=1,3
        arrhcf(numre,i)=r_val(line,iterm,lenline,i,ierr)
        IF (ierr.ne.0) THEN
          WRITE(lout,*)
          WRITE(lout,*) line
          WRITE(lout,*)' --error--  while reading reaction coeff. ',i
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
        ENDIF
      ENDDO

! ----------------------------------------------
! READ THE REACTION (SPECIES AND STOECHIOMETRIE)
! ----------------------------------------------

      ii=0
      xcoeff=+1.

! 3300 is the re-entry point.

3300  CONTINUE

      ii=ii+1
      !print*,'3300',ii, line(ii:ii)

      IF (ii.ge.iterm) THEN
        WRITE(lout,*) line
        CALL errline(lout,line,1,iterm,ii)
        WRITE(lout,*)'   --error--   reaction not properly',
     &               ' completed'
        WRITE(lout,*)
        locheck(numre)=.true.
        lostop=.true.
        RETURN
      ENDIF

! find what must be read at the pointer position (ii).
! if a species is read => must start with 'A-Z', 'a-z', or '(<[{' 
      IF (line(ii:ii).eq.' ') THEN
        GOTO 3300

      ELSE IF (line(ii:ii).eq.'+' ) THEN
        IF (.not.lonumber) xcoeff=+1. 
        GOTO 3300

! JMLT 2021: dashes are ALLOWED in MECHGEN names
      ELSE IF (line(ii:ii).eq.'-'.AND.lco.NE.8 ) THEN
        IF (.not.lonumber) xcoeff=-1.
        GOTO 3300

      ELSE IF (line(ii:ii).eq.'.'.or .
     &        (line(ii:ii).ge.'0'.and.line(ii:ii).le.'9') ) THEN
        IF (.not.lonumber) THEN
          lonumber=.true.
          ib=ii
        ENDIF
        GOTO 3300

      ELSE IF( (line(ii:ii).ge.'A'.and.line(ii:ii).le.'Z') .or.
     &         (line(ii:ii).ge.'a'.and.line(ii:ii).le.'z') .or.
     &          line(ii:ii).eq.'('.or .line(ii:ii).eq.'<'  .or.
     &          line(ii:ii).eq.'['.or .line(ii:ii).eq.'{' ) THEN
        GOTO 3400

      ELSE
        CALL errline(lout,line,1,iterm,ii)
        WRITE(lout,*)'   --error1--   unexpected character',
     &               ' in reaction equation'
        WRITE(lout,*)
        locheck(numre)=.true.
        lostop=.true.
        RETURN
      ENDIF

! read a "block" of data (species or keyword + corresponding stoi. coe.
! ====================================================================

! find delimiter (end of species or keyword)
3400  CONTINUE
      tempsp=' '
      DO jj=ii+1,iterm
        !print*,'3400', jj, line(jj:jj)
        IF ( (line(jj:jj).eq.' ') .or.
     &       (line(jj:jj).eq.'+') .or.
     &       (line(jj:jj).eq.'-'.AND.lco.NE.8) .or.
     &       (line(jj:jj+1).eq.'(+') .or.
     &       (line(jj:jj+1).eq.'=>') ) THEN
          IF (jj-ii.gt.maxlsp) THEN
            WRITE(lout,*) '--error--, separator not found '
            WRITE(lout,*) 'species may exceed maxlsp in line: '
            WRITE(lout,*) line
            WRITE(lout,*) 'at character chain : ',line(ii:jj)
            lostop=.true.
            RETURN
          ENDIF
          tempsp=line(ii:jj-1)
          GOTO 3410
        ENDIF
      ENDDO
3410  CONTINUE
      !print*,'3410: tempsp = ', tempsp(1:maxlsp)

! check if keyword. Length of the keyword is used to identify it in
! the next section (should be changed in a futur version of the program) 
      lenmin=0
      idpero=-1
      iddimer=-1
      IF (tempsp(1:4).eq.'PERO') THEN 
        idpero=ICHAR(tempsp(5:5)) - ICHAR('0')
        IF ( (idpero.lt.0) .or. (idpero.gt.9) ) THEN
          WRITE(lout,*) '--error--, keyword PERO not associated with '
          WRITE(lout,*) 'a number (0-9 expected) at line :  '
          WRITE(lout,*) line
          WRITE(lout,*) 'at character chain : ',line(ii:jj)
          lostop=.true.
          RETURN
        ENDIF
        lenmin=5
      ELSE IF (tempsp(1:4).eq.'DIM_') THEN
        iddimer=ICHAR(tempsp(5:5)) - ICHAR('0')
        IF ( (iddimer.lt.0) .or. (iddimer.gt.9) ) THEN
          WRITE(lout,*) '--error--, keyword DIM_ not associated with '
          WRITE(lout,*) 'a number (0-9 expected) at line :  '
          WRITE(lout,*) line
          WRITE(lout,*) 'at character chain : ',line(ii:jj)
          lostop=.true.
          RETURN
        ENDIF
        lenmin=5
      ELSE IF (tempsp(1:6).eq.'MEPERO') THEN
         lenmin=6
      ELSE IF (tempsp(1:3).eq.'HV ') THEN
         lenmin=2
      ELSE IF (tempsp(1:7).eq.'NOTHING') THEN
         lenmin=7
      ELSE IF (tempsp(1:6).eq.'OXYGEN') THEN
         lenmin=6
      ELSE IF (tempsp(1:6).eq.'EXTRA ') THEN 
         lenmin=5
      ELSE IF (tempsp(1:5).eq.'CVAR ') THEN
         lenmin=4
      ELSE IF (tempsp(1:5).eq.'ISOM ') THEN
         lenmin=4
      ELSE IF (tempsp(1:6).eq.'TBODY ') THEN
         lenmin=5
      ELSE IF (tempsp(1:8).eq.'FALLOFF ') THEN
         lenmin=7
      ELSE IF (tempsp(1:4).eq.'AIN ') THEN
         lenmin=3
      ELSE IF (tempsp(1:4).eq.'AOU ') THEN
         lenmin=3
      ELSE IF (tempsp(1:4).eq.'WIN ') THEN
         lenmin=3
      ELSE IF (tempsp(1:4).eq.'WOU ') THEN
         lenmin=3
      ENDIF

! If keyword found => set corresponding logical parameter
      IF (lenmin.ne.0) THEN
        IF (lonumber) THEN
          CALL errline(lout,line,1,iterm,ib)
          WRITE(lout,*)'   --error--   "M", "HV","CVAR", "NOTHING"',
     &                 ' or "EXTRA" can not be preceded'
          WRITE(lout,*)'               by a stoichiometric coefficient'
          WRITE(lout,*)
          lonumber=.false.
          locheck(numre)=.true.
          lostop=.true.
        ENDIF

! lenmin=7 can be due to keyword 'NOTHING' or 'FALLOFF' => check which
        IF (lenmin.eq.7) THEN
          IF (INDEX(line,"FALLOFF").NE.0) THEN
            lofo=.true.
          ELSE
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   the key-word NOTHING can',
     &                     ' not be used as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE IF (iside.eq.1) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   the key-word NOTHING can',
     &                     ' not be used on the left side'
              WRITE(lout,*)'               of the reaction'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE
              lonothing=.true.
            ENDIF
          ENDIF

! lenmin=6 can be due to keyword 'MEPERO' or 'OXYGEN' => check which
! keyword is ok. temporary version, must be changed in the futur.
        ELSE IF (lenmin.eq.6) THEN 
          IF (tempsp(1:6).eq.'MEPERO') THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   MEPERO can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE
              lo_meo2=.true.
            ENDIF
          ELSE
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   OXYGEN can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE
              lo_o2=.true.
            ENDIF
          ENDIF

! lenmin=5 can be due to keyword 'EXTRA', 'PEROx' or 'TBODY' => check which
        ELSE IF (lenmin.eq.5) THEN
          IF ((idpero.eq.-1).AND.(iddimer.eq.-1)) THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   EXTRA can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE IF (INDEX(line,'EXTRA').NE.0) THEN
              loextra=.true.
            ENDIF
          ELSE IF (iddimer.eq.-1) THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   PERO can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
              RETURN
            ENDIF
            IF (iside.eq.2) THEN
              CALL errline(lout,line,1,iterm,ib)
               WRITE(lout,*)'   --error--   the key-word PERO can',
     &                      ' not be used on the left side'
              WRITE(lout,*)'               of the reaction'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
              RETURN
            ENDIF
            lopero=.true.
            ipero=idpero
          ELSE IF (idpero.eq.-1) THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   DIM_ can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
              RETURN
            ENDIF
            IF (iside.eq.2) THEN
              CALL errline(lout,line,1,iterm,ib)
               WRITE(lout,*)'   --error--   the key-word DIM_ can',
     &                      ' not be used on the left side'
              WRITE(lout,*)'               of the reaction'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
              RETURN
            ENDIF
            lodimer=.true.
            idimer=iddimer
          ENDIF

        ELSE IF (lenmin.eq.4) THEN
          IF (tempsp(1:4).eq.'CVAR') THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   CVAR can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE
              locvar=.true.
            ENDIF
          ELSE IF (tempsp(1:4).eq.'ISOM') THEN
            IF (lopar) THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   ISOM can not be',
     &                     ' use as a "fall-off" species'
              WRITE(lout,*)
              locheck(numre)=.true.
              lostop=.true.
            ELSE
              lo_iso=.true.
            ENDIF
          ENDIF

        ELSE IF (lenmin.eq.2) THEN
          IF (lopar) THEN
            CALL errline(lout,line,1,iterm,ib)
            WRITE(lout,*)'   --error--   HV can not be',
     &                   ' use as a "fall-off" species'
            WRITE(lout,*)
            locheck(numre)=.true.
            lostop=.true.
          ELSE
            lohv=.true.
          ENDIF

        ELSE IF (lenmin.eq.3) THEN
          IF (lopar) THEN
            CALL errline(lout,line,1,iterm,ib)
            WRITE(lout,*)'   --error--   transfer can not be',
     &                   ' use as a "fall-off" species'
            WRITE(lout,*)
            locheck(numre)=.true.
            lostop=.true.
          ELSE
           IF (tempsp(1:3).eq.'AIN') loain=.true.
           IF (tempsp(1:3).eq.'AOU') loaou=.true.
           IF (tempsp(1:3).eq.'WIN') lowin=.true.
           IF (tempsp(1:3).eq.'WOU') lowou=.true.
          ENDIF


        ELSE IF (lenmin.eq.1) THEN
          IF (.not.lopar) lo_m=.true.
        ENDIF
        iii=ii+lenmin-1

! else check if species is known. If not found => error and return.
! Get stoi. coef. and store the data to the tables
      ELSE
        iloc = find_species_index(tempsp,chrsp)
        if (iloc > 4456016) then
          !print*, tempsp, chrsp(iloc)
          stop
        endif
        IF (iloc.le.0) THEN
          WRITE(lout,*) line
          WRITE(lout,*) '--error--, species not found:"',tempsp,'"'
          WRITE(lout,*) 'iloc=',iloc
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ENDIF

        k=INDEX(chrsp(iloc),' ')
        IF (k.eq.0) k=maxlsp+1 
        iii=ii+k-2
        jj=iloc

        IF (lopar) THEN
          CALL errline(lout,line,1,iterm,ii)
          WRITE(lout,*)'   --error--   only "M" can',
     &                 ' be used as the "fall-off" species'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ENDIF

! read stoichiometric coefficient
        IF (lonumber) THEN
          lonumber=.false.
          xcoeff=r_val(line,ib,ii-1,1,ierr)*xcoeff
          IF (ierr.ne.0) THEN
            CALL errline(lout,line,1,iterm,ib)
            WRITE(lout,*)' --error--   while reading stoichio. coeff.'
            WRITE(lout,*)
            locheck(numre)=.true.
            lostop=.true.
          ENDIF
        ENDIF

        ncoef=ncoef+1
        IF (ncoef.gt.mxright) THEN
          WRITE(lout, *) '--error--, species in react. exceed mxright'
          locheck(numre)=.true.
          lostop=.true.
          WRITE(lout, *) line
          WRITE(lout, *) 'iside = ', iside
          DO j=1,mxright
            WRITE(lout,*) 'j=',j,' ',chrsp(iholdsp(j,iside))
          ENDDO
          RETURN
        ENDIF
        iholdsp(ncoef,iside)=jj
        coefsp(ncoef,iside)=xcoeff

      ENDIF

! set pointer to the new position
      ii=iii

! get the next species in the reaction line 
! ==========================================
! if a parenthesis is open (that does not belong to the name of 
! the species), then search for the closing parenthesis. If not found
! then exit
      IF (lopar) THEN
        lopar=.false.
3470    CONTINUE
        ii=ii+1
        !print*,'3470',ii, line(ii:ii)
        IF (ii.gt.iterm) THEN
          CALL errline(lout,line,1,iterm,ii)
          WRITE(lout,*)'   --error--   reaction not properly completed'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ELSE IF (line(ii:ii).eq.' ') THEN
          GOTO 3470
        !ELSE IF (line(ii:ii).ne.')') THEN
        ELSE IF (line(ii:ii+1).ne.'F )') THEN
          CALL errline(lout,line,1,iterm,ii)
          WRITE(lout,*)'   --error--   reaction not properly completed'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ENDIF
      ENDIF

! At that point, a species and the corresponding stoi coef has been
! found => find the next delimiter (i.e. "+","-","(","=>") and go
! back to 3300 to read the next species. If end of the reaction line
! is reached (iterm), then goto 3550 for final check
3500  CONTINUE
      ii=ii+1
      !print*,'3500', ii, line(ii:ii)
      
      IF (ii.gt.iterm) THEN
        GOTO 3550
      ELSE IF (line(ii:ii).eq.' ') THEN
        GOTO 3500
      ELSE IF (line(ii:ii).eq.'(') THEN
        IF (idlopar.eq.iside) THEN
          CALL errline(lout,line,1,iterm,ii)
          WRITE(lout,*)'   --error--   too many "fall-off" indications',
     &                 ' in the reaction'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ELSE
          idlopar=iside
        ENDIF
        lopar=.true.
        lofo=.true.
      !print*,'3500', ii, line(ii:ii)," lofo"
        GOTO 3500
      ELSE IF (line(ii:ii).eq.'+') THEN
        xcoeff=+1.
        GOTO 3300
      ELSE IF (line(ii:ii).eq.'-'.AND.lco.NE.8) THEN
        xcoeff=-1.
        GOTO 3300
      ELSE IF (line(ii:ii+1).eq.'=>') THEN
        xcoeff=+1.
        IF(iside.eq.2)THEN
          CALL errline(lout,line,1,iterm,ii)
          WRITE(lout,*)'   --error--   more than one => in reaction'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
          RETURN
        ENDIF
        ii=ii+1
        iside=2
        ncoef=0
        GOTO 3300
       ELSE
        CALL errline(lout,line,1,iterm,ii)
          PRINT*,'   --error2--   unexpected character',
     &               ' in reaction equation'
          PRINT*,' ii,iterm = ',ii,iterm
          STOP !DEBUG

        WRITE(lout,*)'   --error2--   unexpected character',
     &               ' in reaction equation'
        WRITE(lout,*)' ii,iterm = ',ii,iterm
        WRITE(lout,*)
        locheck(numre)=.true.
        lostop=.true.
        RETURN
      ENDIF

! -----------------------------------------------
! CHECK THAT REACTION IS COMPLETE AND STORE DATA
! -----------------------------------------------

3550  CONTINUE
      !print*,'3550'


! npse(i) is the number of species in the reaction at side i
      DO j=1,3
        nspe(j)=0
      ENDDO

! search for duplicate species (in the same side only)
      DO j=1,2
        DO i=1,mxright
          IF (iholdsp(i,j).ne.0) THEN
            nspe(j)=nspe(j)+1
            DO k=i+1,mxright
              IF (iholdsp(i,j).eq.iholdsp(k,j) ) THEN
                coefsp(i,j) = coefsp(i,j)+coefsp(k,j)
                iholdsp(k,j) = 0
                coefsp(k,j) = 0.
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO

! sort the species according to their declaration order in the input
! chemical scheme => store sorted tables in cpiholdsp and cpcoefsp
      DO j=1,2
        IF  (nspe(j).gt.0) THEN
          DO k=1,nspe(j)
            minsp=numsp+10
            iloc=0
            DO i=1,mxright
              IF (iholdsp(i,j).ne.0) THEN
                IF (iholdsp(i,j).lt.minsp) THEN
                  minsp=iholdsp(i,j)
                  iloc=i
                ENDIF
              ENDIF
            ENDDO
            IF (minsp.eq.0) STOP 'minsp=0'
            IF (iloc.eq.0) STOP 'iloc=0'
            IF (minsp.gt.numsp) STOP 'minsp=numsp+1'
            cpiholdsp(k,j)=iholdsp(iloc,j)
            cpcoefsp(k,j)=coefsp(iloc,j)
            iholdsp(iloc,j)=0
          ENDDO
        ENDIF
      ENDDO

! put the values in the output tables
!      DO j=1,2
!        ipart=0
!        IF  (nspe(j).gt.0) THEN
!          DO i=1,nspe(j)
!              idstoi(numre,i,j)=cpiholdsp(i,j)
!              stoicf(numre,i,j)=cpcoefsp(i,j)
!          ENDDO
!        ENDIF
!        numstoi(numre,j)=nspe(j)
!        mx12stoi=max(numstoi(numre,j),mx12stoi)
!      ENDDO

! put the values in the output tables - reactant side.
! Up to now, number of reactant could be lesser or equal mxright. 
! Check now that the number of reactant does in fact not exceed mxleft
      IF (nspe(1).gt.mxleft) THEN
        WRITE(lout,*) ' '
        WRITE(lout,*) line
        WRITE(lout,*)'  --error--  number of reactant exceed mxleft'
        WRITE(lout,*)
        locheck(numre)=.true.
        lostop=.true.
        RETURN
      ENDIF
      IF  (nspe(1).gt.0) THEN
        DO i=1,nspe(1)
           idrestoi(numre,i)=cpiholdsp(i,1)
           restoicf(numre,i)=cpcoefsp(i,1)
        ENDDO
      ENDIF
      numstoi(numre,1)=nspe(1)
      mx12stoi=max(numstoi(numre,1),mx12stoi)

! put the values in the output tables - product side
      IF (nspe(2).gt.0) THEN
        DO i=1,nspe(2)
           idpdstoi(numre,i)=cpiholdsp(i,2)
           pdstoicf(numre,i)=cpcoefsp(i,2)
        ENDDO
      ENDIF
      numstoi(numre,2)=nspe(2)
      mx12stoi=max(numstoi(numre,2),mx12stoi)

      mx1stoi=max(numstoi(numre,1),mx1stoi)
      mx2stoi=max(numstoi(numre,2),mx2stoi)

! check if reactants were given
      IF (numstoi(numre,1).eq.0) THEN
        CALL errline(lout,line,1,iterm,0)
        WRITE(lout,*)'   --error--   reaction equation has no',
     &               ' reactants'
        WRITE(lout,*)
        locheck(numre)=.true.
        lostop=.true.
      ENDIF

! check if products were given 
      IF (lonothing) THEN
        IF (numstoi(numre,2).ne.0) THEN
          CALL errline(lout,line,1,iterm,0)
          WRITE(lout,*)'   --error--   reaction can not have any',
     &                 ' products if key-word NOTHING is used'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
        ENDIF
      ELSE
        IF (numstoi(numre,2).eq.0) THEN
          CALL errline(lout,line,1,iterm,0)
          WRITE(lout,*)'  --error--  reaction equation has no products.'
          WRITE(lout,*)'             Use keyword NOTHING if needed'
          WRITE(lout,*)
          locheck(numre)=.true.
          lostop=.true.
        ENDIF
      ENDIF

! COMPUTE DATA FOR SIDE 3 (side3=side2-side1)

! copy tables
!      ipart=0
!      DO j=1,2
!        DO i=1,nspe(j)
!          iholdsp(i,j)=cpiholdsp(i,j)
!          coefsp(i,j)=cpcoefsp(i,j)
!        ENDDO
!      ENDDO

! sort species accoding to declaration in the input file using
! both side 1 and 2
!      DO k=1,nspe(1)+nspe(2)
!        minsp=numsp+10
!        DO j=1,2
!          DO i=1,nspe(j)
!              IF (iholdsp(i,j).ne.0) THEN
!                IF (iholdsp(i,j).lt.minsp) THEN
!                  minsp=iholdsp(i,j)
!                ENDIF
!              ENDIF
!          ENDDO
!        ENDDO
!
!        DO j=1,2
!          ipos(j)=0
!          DO i=1,nspe(j)
!            IF (iholdsp(i,j).eq.minsp) THEN
!              ipos(j)=i
!            ENDIF
!          ENDDO
!        ENDDO

!        IF (ipos(1)+ipos(2).gt.0) THEN
!          ipart=ipart+1
!          IF (ipos(1).gt.0) THEN
!            idstoi(numre,ipart,3) = cpiholdsp(ipos(1),1)
!            stoicf(numre,ipart,3) = -cpcoefsp(ipos(1),1)
!            iholdsp(ipos(1),1)=0
!          ENDIF
!          IF (ipos(2).gt.0) THEN
!            idstoi(numre,ipart,3) = cpiholdsp(ipos(2),2)
!            stoicf(numre,ipart,3) = 
!     &                      stoicf(numre,ipart,3)+cpcoefsp(ipos(2),2)
!            iholdsp(ipos(2),2)=0
!          ENDIF
!        ENDIF
!      ENDDO

! store data for side 3
!      numstoi(numre,3)=ipart
!      mx3stoi=max(numstoi(numre,3),mx3stoi)

      END
