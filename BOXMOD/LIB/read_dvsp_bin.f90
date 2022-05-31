      SUBROUTINE read_dvsp_bin
! read diffusion volumes and
! redistribute from mxsat into maxsp space 

      USE flags_module,ONLY: wall_fg
      USE io_units_module,ONLY: lread,lout
      USE akparameter_module,ONLY: mxsat,maxsp,maxlsp
      USE module_data_gecko_main,ONLY: namdif,chrsp,numsp, &
                                       dvsp,difvol,difid

      IMPLICIT NONE

      INTEGER:: i,isp,j,j1,k,k2
      CHARACTER(70) :: line
      CHARACTER(20) :: filnam
      CHARACTER(maxlsp) :: tmpname

! -------------------------------------
! retrieve diffusion volume
      filnam="difv.dat"
      OPEN(lread,FILE=filnam,STATUS='OLD')

      k2=2
      IF(wall_fg.eq.1) k2=3

!! read dicnams, data
      DO i=1,mxsat
        READ(lread,'(a)') line
        IF(line(1:3).EQ.'END')EXIT
        READ(line,*,err=9) namdif(i), difvol(i)

!! OLD CODE, PRODUCES ERRORS
!! call akspnum to assign species index
!! Only assignes FIRST index encountered!
!!        CALL akspnum(namdif(i),chrsp,numsp,isp)
!!        dvsp(isp)=difvol(i) 

! NEW CODE
! set up indexing as per satid in wrtlinkncdf: 
! difid = dif index of chrsp species
        j1 = 1
        isp = j1
! END NEW INITIALIZATION

        DO k=1,k2
          tmpname = ''
          if (k == 1) then
            tmpname(1:1) = 'G'
          else if (k == 2) then
            tmpname(1:1) = 'A'
          else if (k==3) then
            tmpname(1:1) = 'W'
          endif
          tmpname(2:maxlsp) = namdif(i)(2:maxlsp) 
          call akspnum(tmpname(1:maxlsp), chrsp, numsp, isp)
          dvsp(isp)=difvol(i)
          difid(i) = i 
        ENDDO
! set up indexing as per satid in wrtlinkncdf: 
! difid = dif index of chrsp species
!        j1 = 1
!        DO k=1,k2
!          DO j=j1,maxsp
!            IF (namdif(i)(2:maxlsp).eq.chrsp(j)(2:maxlsp)) THEN
!              isp=j
!              j1=j+1
!              difid(j)=i
!              dvsp(isp)=difvol(i) 
!              EXIT
!            ELSE
!              isp=0
!            ENDIF
!          ENDDO
!        ENDDO
        IF(isp.EQ.0)THEN
          WRITE(6,*) '--error--, sat species unidentified in readdif'
          WRITE(6,'(a)') namdif(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF

        CYCLE ! bypass error handling if read successful

9       WRITE(lout,*) '-error--, while reading difv.dat'
        WRITE(lout,*) '        , at line : ',i
        STOP
      ENDDO

      CLOSE(lread)

! -------------------------------------
      END SUBROUTINE read_dvsp_bin
! -------------------------------------

