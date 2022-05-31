      SUBROUTINE akkrat6(chrsp,numsp,numre,num_m,numfo,numextra,numo2,
     &                   nummeo2,numiso,idiso,
     &                   nclro2,numreacro2,idrestoi,id_m,idfo,idextra,
     &                   ncldimer,numreacdimer,idreacdimer,cdimer,
     &                   ido2,idmeo2,idreacro2,
     &                   arrhcf,focf,extracf,isocf,wmol,
     &                   vd,cro2,cmeo2,
     &                   temp,sumc,ibox,water,height,saero,qfor,
     &                   inorg_aer)

      USE flags_module,ONLY: soa_fg,wall_fg,vbs_fg
      USE akparameter_module
      USE inorganic_aer_module,ONLY : aer_population
      USE module_data_gecko_main,ONLY: small,idain,idaou,idwin,idwou,
     &                                     numain,numaou,numwin,numwou
      USE vbs_module, ONLY: init_vbs_rates
      !$ use OMP_LIB

      IMPLICIT NONE
!-------------------------------------------------------------
* INPUT

      CHARACTER(maxlsp) chrsp(maxsp)
      INTEGER  ibox
      INTEGER  numsp, numre
      INTEGER  num_m, numfo, numextra
      INTEGER  nclro2,ncldimer
      INTEGER  idrestoi(maxre,mxleft)
      INTEGER  id_m(max_m)
      INTEGER  idfo(maxfo,3)
      INTEGER  idextra(maxextra)
      INTEGER  numo2,ido2(maxo2)
      INTEGER  numiso,idiso(maxiso)
      INTEGER  nummeo2,idmeo2(mxrpero)
      INTEGER  numreacro2(maxro2),idreacro2(mxrpero,maxro2)
      INTEGER  numreacdimer(maxdimer),idreacdimer(mxrdimer,maxdimer)
      
      REAL     arrhcf(maxre,3)
      REAL     focf(maxaux+3,maxfo)
      REAL     extracf(maxaux,maxextra)
      REAL     temp,height,saero
! USE REAL UNDER COMPILATION OPTION real-8
!      DOUBLE PRECISION     extracf(maxaux,maxextra)
      REAL     isocf(maxaux,maxiso)
      REAL     wmol(maxsp)
      REAL     cro2(maxro2)
      REAL     cdimer(maxdimer)
      REAL     Vd(maxsp)
      REAL     cmeo2
      REAL,DIMENSION(mbox) :: sumc,water

      TYPE(aer_population), INTENT(IN) :: inorg_aer(mbox)

* OUTPUT
      REAL     qfor(maxre)

* LOCAL
      INTEGER  ire, i,j
      REAL     tln,tinv,sumcln,lno2,lncro2(maxro2),lncdimer(maxdimer)
!      REAL     small
      REAL     cmeo2ln,totcro2

! PRINT*,"initialize"
      !small= TINY(1.0)
!      small = 1e-32

! PRINT*,"initialize 1"
      qfor = 0.0

! PRINT*,"initialize 2" 
      tln=log(max(temp, small))
      tinv=1.0/temp

! PRINT*,"initialize 3"
      sumcln=log(max(sumc(ibox), small))
      lno2=log(max(sumc(ibox)*0.2, small*0.2))

! PRINT*,"initialize 4"
      cmeo2ln=log(max(cmeo2, small))

      lncro2 = -700.
      DO i=1, nclro2
        IF(cro2(i) .GT. 0) THEN
          lncro2(i) = log(cro2(i))
        ENDIF
      ENDDO

! PRINT*,"initialize 4"
      lncdimer(1:ncldimer) = log(max(cdimer(1:ncldimer), small))

! PRINT*,"compute rate constant"
      DO ire=1,numre
        qfor(ire)=arrhcf(ire,1)+arrhcf(ire,2)*tln-arrhcf(ire,3)*tinv
      ENDDO

! PRINT*,"reaction with third body M"
      DO i=1,num_m
        ire=id_m(i)
        qfor(ire)=qfor(ire)+sumcln
      ENDDO

! PRINT*,"reaction with O2"
      DO i=1,numo2
        ire=ido2(i)
        qfor(ire)=qfor(ire)+lno2
      ENDDO

! PRINT*,"reaction with CH3O2"
      DO i=1,nummeo2
        ire=idmeo2(i)
        qfor(ire)=qfor(ire)+cmeo2ln
      ENDDO

! PRINT*,"reaction with RO2"
      DO i=1,nclro2
        DO j=1,numreacro2(i)
          ire=idreacro2(j,i)
          qfor(ire)=qfor(ire)+lncro2(i)
        ENDDO
      ENDDO

! PRINT*,"reaction with DIM_x"
      DO i=1,ncldimer
        DO j=1,numreacdimer(i)
          ire=idreacdimer(j,i)
          qfor(ire)=qfor(ire)+lncdimer(i)
        ENDDO
      ENDDO

! PRINT*,"Isomerisation reaction"
      DO i=1,numiso
          ire=idiso(i)
          qfor(ire)=qfor(ire)+log(max(isocf(1,i)*temp**4 +
     &              isocf(2,i)*temp**3 + isocf(3,i)*temp**2 +
     &              isocf(4,i)*temp + isocf(5,i), small))
      ENDDO

! PRINT*,"fall off reaction"
      DO i=1,numfo
        ire=idfo(i,1)
!JMLT: MECHGEN, DO NOT CALL akkfo IF NO-INORGS
        IF(EXP(arrhcf(ire,1)).GT.small) THEN ! JMLT MECHGEN
        !PRINT*,"ire = ",ire, EXP(arrhcf(ire,1)),arrhcf(ire,2:3)
        CALL akkfo(maxre,maxfo,maxaux,
     1             ire,i,focf,temp,sumc(ibox),arrhcf,qfor(ire))
        ENDIF ! JMLT for MECHGEN
      ENDDO

! PRINT*,"extra reaction"
      if (vbs_fg .gt. 0) then
        call init_vbs_rates(temp)
      endif

      DO i=1,numextra
        ire=idextra(i)
!JMLT: MECHGEN, DO NOT CALL akkextra4 IF NO-INORGS
        IF(EXP(arrhcf(ire,1)).GT.small) THEN ! JMLT MECHGEN
        !PRINT*,"ire = ",ire, extracf(1,i)
        CALL akkextra4
     &           (chrsp,numsp,idrestoi,extracf,wmol,Vd,
     &            ire,i,temp,sumc,ibox,water,height,saero,qfor(ire),
     &            inorg_aer)
        ENDIF ! JMLT for MECHGEN
      ENDDO

      DO ire=1,numre
        qfor(ire)=exp(qfor(ire))
      ENDDO

! zero out soa rates unless using dynamic aerosol
      IF(soa_fg.NE.2)THEN
!$OMP PARALLEL DO private(i,ire)
        DO i = 1,numain
          ire = idain(i)
          qfor(ire) = 0.
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(i,ire)
        DO i = 1,numaou
          ire = idaou(i)
          qfor(ire) = 0.
        ENDDO
!$OMP END PARALLEL DO
      ENDIF
! zero out wall exchange rates if not required 
      IF(wall_fg.EQ.0)THEN
!$OMP PARALLEL DO private(i,ire)
        DO i = 1,numwin
          ire = idwin(i)
          qfor(ire) = 0.
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO private(i,ire)
        DO i = 1,numwou
          ire = idwou(i)
          qfor(ire) = 0.
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      END
