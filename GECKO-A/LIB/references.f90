MODULE references
  IMPLICIT NONE
  INTEGER,PARAMETER :: mxtweet=1000 ! max # of comment/references 
  INTEGER,PARAMETER :: mxlcod=10    ! max length of comment's code in database
  INTEGER,PARAMETER :: mxltweet=400 ! max length of a comment/reference
  
  INTEGER,SAVE           :: ntweet         ! # of available "tweet" (references)  
  CHARACTER(LEN=mxlcod)  :: code(mxtweet)  ! Code for each "tweet" 
  CHARACTER(LEN=mxltweet):: tweet(mxtweet) ! Full text correxponding to the tweet 

END MODULE references
