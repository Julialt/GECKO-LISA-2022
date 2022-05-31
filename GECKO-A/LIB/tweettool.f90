MODULE tweettool
IMPLICIT NONE
CONTAINS

! ======================================================================
SUBROUTINE rdtweet()
  USE keyparameter, ONLY: tfu1
  USE references, ONLY: mxtweet,mxlcod,mxltweet,ntweet,code,tweet
  USE sortstring              ! to sort the tweet codes
  USE searching, ONLY: srh5   ! to seach in a sorted list
  IMPLICIT NONE

  CHARACTER(LEN=100) :: filename
  INTEGER       :: i,n1,ierr,ilin,lcode,cblank
  CHARACTER(LEN=420) :: line                   ! string to be read
  CHARACTER(LEN=mxlcod+10) :: leftline         ! code side string
  CHARACTER(LEN=mxltweet+5):: rightline        ! tweet side string
  CHARACTER(LEN=mxlcod)    :: tpcode(mxtweet)  ! temporary copy (unsorted) of code
  CHARACTER(LEN=mxltweet)  :: tptweet(mxtweet) ! temporary copy (unsorted) of code
  
  ntweet=0 ; tpcode(:)=' ' ; tptweet(:)=' '
  
! open the file
  filename='../DATA/tweetdic.dat'
  OPEN(tfu1,FILE=filename, FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF (ierr/=0) THEN
    WRITE(6,*) '--error--, while trying to open: ',TRIM(filename)
    STOP "in rdtweet" 
  ENDIF

! Read the file
! ---------------
  ilin=0
  rdloop: DO
    leftline=' ' ; rightline=' '
    ilin = ilin+1
    READ(tfu1,'(a)',IOSTAT=ierr) line
    IF (ierr/=0) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'keyword "END" missing ?'
      STOP "in rdtweet" 
    ENDIF
    IF (line(1:1)=='!') CYCLE rdloop
    IF (line(1:3)=='END') EXIT rdloop
    
! read line
    line=ADJUSTL(line) ; n1=INDEX(line,':') 
    IF (n1<=1) THEN 
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Expect ":" separator in line: ',TRIM(line)
      STOP "in rdtweet" 
    ENDIF

! read the code
    leftline=line(1:n1-1) ; rightline=ADJUSTL(line(n1+1:))
    lcode=LEN_TRIM(leftline)              ! check the length of the code
    IF (lcode > mxlcod) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Length of the code exceed mxlcod in: ',TRIM(line)
      STOP "in rdtweet" 
    ENDIF
    cblank=INDEX(leftline(1:lcode),' ')   ! check for " " char
    IF (cblank /=0 ) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Unexpected " " char in code: ',TRIM(line)
      STOP "in rdtweet" 
    ENDIF
    ntweet=ntweet+1
    IF (ntweet > mxtweet) THEN             ! check # of tweets
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'too many tweets (check table size). mxtweet= ',mxtweet
      STOP "in rdtweet" 
    ENDIF
    tpcode(ntweet)=leftline(1:mxlcod)

! read the tweet
    IF (LEN_TRIM(leftline)>mxltweet) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*) 'at line number: ',ilin
      WRITE(6,*) 'Length of the code exceed mxltweet in: ',TRIM(line)
      STOP "in rdtweet" 
    ENDIF
    tptweet(ntweet)=rightline(1:mxltweet)
  ENDDO rdloop
  CLOSE(tfu1)

! load sorted tables tables
! -------------------------

! sort the codes
  code(:)=tpcode(:)
  CALL sort_string(code(1:ntweet)) 

! check for duplicate 
  ierr=0
  DO i=1,ntweet-1
    IF (code(i)==code(i+1)) THEN
      WRITE(6,*) '--error--, while reading file: ',TRIM(filename)
      WRITE(6,*)  'Following code identified 2 times: ',TRIM(code(i))
      ierr=1
    ENDIF
  ENDDO
  IF (ierr/=0) STOP "in rdtweet" 
  
! sort the tables according to formula
  DO i=1,ntweet
    ilin=srh5(tpcode(i),code,ntweet)
    IF (ilin <= 0) THEN
      WRITE(6,*) '--error--, while sorting tweet in: ',TRIM(filename)
      WRITE(6,*) 'tweet "lost" after sorting the list: ', TRIM(tpcode(i))
      STOP "in rdtweet"
    ENDIF
    tweet(ilin)=tptweet(i)
  ENDDO

END SUBROUTINE rdtweet ! -----------------------------------------------

CHARACTER*(mxltweet) FUNCTION fullref(shortcode)
 USE references
 IMPLICIT NONE
 CHARACTER*(mxlcod) :: shortcode
 INTEGER            :: i
 
 fullref=' '
 DO i=1,mxtweet
   IF (code(i)==shortcode) THEN
     fullref=tweet(i)
     RETURN
   ENDIF 
 ENDDO
END FUNCTION

END MODULE tweettool
