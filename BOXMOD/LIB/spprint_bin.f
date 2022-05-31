*
* ===========================================================================
*
      SUBROUTINE wrtlout(chrsp,numsp,lout,time,c,temp,sumc,ibox)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      CHARACTER(maxlsp) chrsp(maxsp)
      INTEGER  numsp,lout,ibox
      REAL     time,c(maxsp),temp,sumc

* LOCAL
      INTEGER npl,nprint,j,imin,imax,i
      REAL    r
      CHARACTER(1)::cbox
!---------------------------------------

      npl=3
      nprint=numsp/npl
      r=mod(numsp,npl)
      IF (r.ne.0.0) nprint=nprint+1

* write to output-file
      WRITE(lout,*) ' '

      WRITE(cbox,'(i1)') ibox
      WRITE(lout,*)'Box //',cbox

      WRITE(lout,*) ' '
      WRITE(lout,1000) time,sumc,temp
      DO j=1,nprint
        imin=npl*(j-1)+1
        imax=imin+npl-1
        IF(imax.gt.numsp)imax=numsp
        WRITE(lout,2000)(c(i),chrsp(i),i=imin,imax)
      ENDDO

* format
1000  FORMAT(/,' t(s)=',1pe11.4,'        sumc=',1pe11.4,
     &       '      T(K)=',1pe11.4)
2000  FORMAT(1x,4(1pe11.4,1x,a15))

      END
*
* =========================================================================
*      
      SUBROUTINE wrtlppf(lpp, numsp, chrsp)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lpp, numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* LOCAL
      CHARACTER(80) nameothers(3)
      INTEGER  i
!---------------------------------------

c      nameothers(1)='TIME'
c      nameothers(2)='PRESSURE'
c      nameothers(3)='TEMPERATURE'
c
c* write to ppf-file
c      WRITE(lpp)numsp,3,maxlsp,0,0
c      WRITE(lpp)(nameothers(i)(1:maxlsp),i=1,3),
c     &           (chrsp(i),i=1,numsp)

      nameothers(1)='TIME'
      nameothers(2)='TEMPERATURE'
      nameothers(3)='HUMIDITY'

* write to ppf-file
      WRITE(lpp)numsp,3,maxlsp,0,0
      WRITE(lpp)(nameothers(i)(1:maxlsp),i=1,3),
     &           (chrsp(i),i=1,numsp)

      END
*
* =========================================================================
*
      SUBROUTINE wrtlppa(lpa, nsat,ndim, chrsp, wmol, idsat)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lpa, nsat, idsat(mxsat),ndim
      CHARACTER(maxlsp) chrsp(maxsp)
      REAL     wmol(maxsp)
* LOCAL
      CHARACTER(80) nameothers(3)
      INTEGER  i
!---------------------------------------

      nameothers(1)='TIME'
      nameothers(2)='TEMPERATURE'
      nameothers(3)='HUMIDITY'

* write to ppa file
      WRITE(lpa)nsat,ndim,2,maxlsp,0,0
      WRITE(lpa)(nameothers(i)(1:maxlsp),i=1,2),
     &           (chrsp(idsat(i)),i=1,nsat),
     &           (wmol(idsat(i)),i=1,nsat)

      END
