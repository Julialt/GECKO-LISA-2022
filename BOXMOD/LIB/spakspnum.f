      SUBROUTINE akspnum(line,chrsp,numsp,nsp)

      USE sorting, ONLY: find_species_index

      USE akparameter_module
      IMPLICIT NONE

* INPUT
      CHARACTER(maxlsp) line
      CHARACTER(maxlsp),INTENT(IN) :: chrsp(maxsp)
      INTEGER  numsp
      INTEGER  stridx(maxsp) ! sorted index of chrsp

* OUTPUT
      INTEGER  nsp

* LOCAL
      CHARACTER(maxlsp) ichrsp
      INTEGER  i, ibgn, iendp1, iend, isplen, isp


      ichrsp=' '
      nsp=0

      DO 100 i=1,len(line)
        ibgn=i
        IF (line(i:i).NE.' ') GOTO 200
100   CONTINUE
      RETURN
200   CONTINUE

      DO 300 i=ibgn+1,len(line)
        iendp1=i
        IF (line(i:i).EQ.' ') GOTO 400
300   CONTINUE

      iendp1=len(line)+1
400   CONTINUE



      iend=iendp1-1
      isplen=iendp1-ibgn
      IF (isplen.GT.maxlsp) RETURN
      ichrsp=line(ibgn:iend)
      
      nsp = find_species_index(ichrsp, chrsp)

      END

