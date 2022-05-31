!------------------------------------------
      PROGRAM read_rates_auto
!------------------------------------------
! program to find instantaneous reaction rates for all reactions
! involving target species "findspec"
! compile & run IN THE RESULTS DIRECTORY using command sequence:
! gfortran `nc-config --flibs` `nc-config --fflags` read_rate_auto.f
! ./a.out
!------------------------------------------
      USE NETCDF

      IMPLICIT NONE

      INTEGER :: ncid
      INTEGER :: ntout,nbox,maxre,maxsp
      INTEGER :: numre,numsp
      INTEGER :: mxleft,mxright
      INTEGER :: maxro2, mxrpero
      INTEGER :: i,i1,i2,j,k,p1,p2,isp,ire,itime,ialert,ipero

      REAL,     ALLOCATABLE :: time(:)
      INTEGER,  ALLOCATABLE :: numstoi(:,:)
      INTEGER,  ALLOCATABLE :: idrestoi(:,:)
      INTEGER,  ALLOCATABLE :: idpdstoi(:,:)
      INTEGER,  ALLOCATABLE :: idreacro2(:,:)
      INTEGER,  ALLOCATABLE :: idmeo2(:)
      REAL,     ALLOCATABLE :: restoicf(:,:)
      REAL,     ALLOCATABLE :: pdstoicf(:,:)
      REAL,     ALLOCATABLE :: reacrate(:,:,:)
      REAL,     ALLOCATABLE :: aggrate(:)
      REAL,     ALLOCATABLE :: netrate(:)

      CHARACTER(LEN=1) :: cpero
      CHARACTER(LEN=9) :: filename

!----------------------------------------
! USER EDIT SECTION
!----------------------------------------
! choose correct species length for model
      INTEGER,PARAMETER :: mxlsp =8 ! GECKO
!      INTEGER,PARAMETER :: mxlsp =10 ! MG
      CHARACTER(LEN=mxlsp),ALLOCATABLE :: chrsp(:)
! list of species to assess
      INTEGER,PARAMETER :: nfind = 10
      CHARACTER(LEN=mxlsp-2),PARAMETER :: findspec(nfind) = 
     & (/"C03000","203000","203001","O03001","O03000","P03000",
     &   "D03000","CH3O2 ","3K3000","3K2000"/)
!     & (/"RO2-0001"/)
!     & (/"RO2-0001","I-C3-OH ","N-C3-OH ","PPN     ","PROPALD ",
!     &   "CH3O2   ","RO3-0001","RO3-0002"/)
      CHARACTER(LEN=mxlsp-2) :: teststr
! last time for which you want output
!      REAL,PARAMETER :: reftime = 19800
!----------------------------------------
      filename = "outdat.nc"
!----------------------------------------

      CALL open_ncfile_readonly(filename,ncid)

!---- read dimensions ----
      CALL eznc_get_dimension(ncid,"maxre",maxre)
      CALL eznc_get_dimension(ncid,"maxsp",maxsp)
      CALL eznc_get_dimension(ncid,"mxleft",mxleft)
      CALL eznc_get_dimension(ncid,"mxright",mxright)
      CALL eznc_get_dimension(ncid,"maxro2",maxro2)
      CALL eznc_get_dimension(ncid,"mxrpero",mxrpero)
      CALL eznc_get_dimension(ncid,"nbox",nbox)
      CALL eznc_get_dimension(ncid,"ntout",ntout)

      CALL eznc_get_0Dint(ncid,"numre",numre)
      CALL eznc_get_0Dint(ncid,"numsp",numsp)

!---- allocate variable arrays ----
      ALLOCATE(time(ntout),numstoi(maxre,mxleft),
     $                     idrestoi(maxre,mxleft),
     &                     idpdstoi(maxre,mxright),
     &                     idreacro2(mxrpero,maxro2),
     &                     idmeo2(mxrpero),
     &                     restoicf(maxre,mxleft),
     &                     pdstoicf(maxre,mxright),
     &                     reacrate(maxre,nbox,ntout),
     &                     aggrate(ntout),
     &                     netrate(ntout))
      ALLOCATE(character(len=mxlsp) :: chrsp(maxsp))

!---- read variable arrays ----
      CALL eznc_get_1Dreal(ncid,"time",ntout,  
     $                          time(1:ntout),1,ntout)

      CALL eznc_get_2Dint(ncid,"numstoi",maxre,mxleft,
     $                    numstoi(1:numre,1:mxleft),    
     $                    1,numre,1,mxleft)

      CALL eznc_get_2Dint(ncid,"idrestoi",maxre,mxleft,
     $                    idrestoi(1:numre,1:mxleft),    
     $                    1,numre,1,mxleft)

      CALL eznc_get_2Dint(ncid,"idpdstoi",maxre,mxright,
     $                    idpdstoi(1:numre,1:mxright),    
     $                    1,numre,1,mxright)

      CALL eznc_get_2Dint(ncid,"idreacro2",mxrpero,maxro2,
     $                    idreacro2(1:mxrpero,1:maxro2),    
     $                    1,mxrpero,1,maxro2)

      CALL eznc_get_1Dint(ncid,"idmeo2",mxrpero,
     $                    idmeo2(1:mxrpero),    
     $                    1,mxrpero)

      CALL eznc_get_2Dreal(ncid,"restoicf",maxre,mxleft,
     $                    restoicf(1:numre,1:mxleft),    
     $                    1,numre,1,mxleft)

      CALL eznc_get_2Dreal(ncid,"pdstoicf",maxre,mxright,
     $                    pdstoicf(1:numre,1:mxright),    
     $                    1,numre,1,mxright)

      CALL eznc_get_3Dreal(ncid,"reacrate",maxre,nbox,ntout, 
     $                    reacrate(1:numre,1:nbox,1:ntout),    
     $                    1,numre,1,nbox,1,ntout)

      CALL eznc_get_1Dchar(ncid,"chrsp",mxlsp,maxsp,
     $                    chrsp,    1,numsp)

!---- find index of "findspec" ----
      DO k = 1,nfind

! loop (test) all species
        DO i = 1,numsp
          teststr = findspec(k)
          IF(INDEX(chrsp(i),teststr(1:LEN_TRIM(teststr))).GT.0)THEN
            isp = i
            PRINT*,"species found: ",chrsp(isp),isp
            EXIT
          ENDIF
        ENDDO ! numsp

! write species header
        WRITE(55,*)"time(s), rates_for ",findspec(k)(:),time(:)
        netrate(:) = 0.
        aggrate(:) = 0.

! check that reaction (1) does NOT include findspec(k)
        DO j = 1,numstoi(1,1)
          IF (idrestoi(1,j).EQ.isp) THEN
            PRINT*,"ERROR: REACTION 1 CONTAINS TARGET SPECIES"
            STOP
          ENDIF
        ENDDO
        
! loop all subsequent reactions
        !DO i = 2,numre+1
        DO i = 118,127
          i1 = i-1

! check if it's a PERO reaction
          ialert = 0
          ipero = 0
          DO p1 = 1,mxrpero
            IF(i.EQ.idmeo2(p1)) THEN
                ialert = 1
                ipero = 0
                EXIT
            ENDIF
            DO p2 = 1,maxro2
              IF(i.EQ.idreacro2(p1,p2)) THEN
                ialert = 1
                ipero = p2
                EXIT
              ENDIF
            ENDDO
          ENDDO

! if you're out of reactions, write rate (if you have one) and exit
          IF(i1.EQ.numre) THEN
            IF(aggrate(ntout).NE.0)THEN
              IF (numstoi(i2,1).EQ.1)THEN
                WRITE(55,*)i2,"loss ",
     & chrsp(idrestoi(i2,1))(1:LEN_TRIM(chrsp(idrestoi(i2,1))))//" ",
     &                                  aggrate(:)
              ELSE
!! bimolecular loss reaction
                WRITE(55,*)i2,"loss ",
     & chrsp(idrestoi(i2,1))(1:LEN_TRIM(chrsp(idrestoi(i2,1))))//"+"//
     & chrsp(idrestoi(i2,2))(1:LEN_TRIM(chrsp(idrestoi(i2,2))))//" ",
     &                                  aggrate(:)
              ENDIF ! numstoi(i1,1).EQ.1
            ENDIF ! aggrate(ntout).NE.0
            GOTO 99
          ENDIF

! NET rates are for all reactions of species "findspec(k)"
! AGGREGATE rates are for reactions of "findspec(k) with same reactants
! if new reactant side, WRITE PREVIOUS aggregated totals
! loop reactants
          DO j = 1,numstoi(i,1)
            IF (idrestoi(i,j).EQ.isp) THEN
              i2 = i ! index for last found loss reaction

! check if we're aggregating or starting fresh
              IF((j.EQ.1.AND.idrestoi(i,2).NE.idrestoi(i1,2)).OR.
     &           (j.EQ.2.AND.idrestoi(i,1).NE.idrestoi(i1,1)).OR.
     &           (ialert.EQ.1))THEN
! start fresh
! write PREVIOUS aggregated rate if you have one
                IF(aggrate(ntout).NE.0)THEN
!! unimolecular loss reaction
                  IF (numstoi(i,1).EQ.1.AND.ialert.EQ.0)THEN
                   WRITE(55,*)i1,"loss ",
     & chrsp(idrestoi(i1,1))(1:LEN_TRIM(chrsp(idrestoi(i1,1))))//" ",
     &                                  aggrate(:)
                  ELSEIF (ialert.EQ.1)THEN
!! bimolecular loss reaction
                    WRITE(cpero,'(i1)')ipero
                    IF (ipero.EQ.0)THEN
                    WRITE(55,*)i1,"loss ",
     & chrsp(idrestoi(i1,1))(1:LEN_TRIM(chrsp(idrestoi(i1,1))))//"+"//
     & "MEPERO"//" ", aggrate(:)
                    ELSE
                    WRITE(55,*)i1,"loss ",
     & chrsp(idrestoi(i1,1))(1:LEN_TRIM(chrsp(idrestoi(i1,1))))//"+"//
     & "PERO"//cpero(1:1)//" ", aggrate(:)
                    ENDIF
                  ELSE
!! bimolecular loss reaction
                    WRITE(55,*)i1,"loss ",
     & chrsp(idrestoi(i1,1))(1:LEN_TRIM(chrsp(idrestoi(i1,1))))//"+"//
     & chrsp(idrestoi(i1,2))(1:LEN_TRIM(chrsp(idrestoi(i1,2))))//" ",
     &                                  aggrate(:)
                  ENDIF ! numstoi(i1,1).EQ.1
                ENDIF ! aggrate(ntout).NE.0
! reset aggrate
                aggrate(:)=0
              ENDIF

              aggrate(:) = aggrate(:)-reacrate(i,1,:)*restoicf(i,j)
              netrate(:) = netrate(:)-reacrate(i,1,:)*restoicf(i,j) 
              
            ENDIF ! idrestoi(i,j).EQ.isp
          ENDDO ! j=numstoi(i,1)

! loop products
          DO j = 1,numstoi(i,2)
            IF (idpdstoi(i,j).EQ.isp) THEN
              netrate(:)= netrate(:)+reacrate(i,1,:)*pdstoicf(i,j) 

!! unimolecular production reaction
              IF (numstoi(i,1).EQ.1)THEN
                WRITE(55,*)i,"prod ",
     & chrsp(idrestoi(i,1))(1:LEN_TRIM(chrsp(idrestoi(i,1))))//" ",
     &                      reacrate(i,1,:)*pdstoicf(i,j)
              ELSE
!! bimolecular production reaction
                WRITE(55,*)i,"prod ",
     & chrsp(idrestoi(i,1))(1:LEN_TRIM(chrsp(idrestoi(i,1))))//"+"//
     & chrsp(idrestoi(i,2))(1:LEN_TRIM(chrsp(idrestoi(i,2))))//" ",
     &                      reacrate(i,1,:)*pdstoicf(i,j)
              ENDIF ! numstoi(i,1).EQ.1
            ENDIF ! idpdstoi(i,j).EQ.isp
          ENDDO ! j=numstoi(i,2)

        ENDDO ! i=numre

99      CONTINUE

! write net rates for species
        WRITE(55,*)i,"net_rate ",findspec(k),netrate(:)
        WRITE(55,*)" "
      ENDDO ! nfind

      DEALLOCATE(time,reacrate,netrate,aggrate,chrsp,numstoi,
     &                               idrestoi,idpdstoi,
     &                               idmeo2,idreacro2,
     &                               pdstoicf,restoicf)
      
      END PROGRAM

!==============================================================!
! PURPOSE: COLLECTION OF ROUTINES TO HANDLE NETCDF-            !
!          FORMAT INPUT & OUTPUT                               !
!          FOR GECKO-A BOX MODEL                               !
! Author: Julia Lee-Taylor, NCAR, January 2018                 !
! Updates:                                                     !
! 2019-April-11 : removed unused variables                     !
!==============================================================!
! ASSUMPTION: variable dimensions, size, type are known        !
!             If this is not the case, the info can be         ! 
!             found using 2 successive calls, e.g.             !
!      status = nf90_inq_varid(ncid, "errnam", VarId)          ! 
!      status = nf90_inquire_variable(ncid,VarId,              !  
!     $         xtype,numDims,dimids,numAtts)                  ! 
!                                                              !
! SUBROUTINE LIST:                                             !
! Subroutine names generally begin with "eznc_" (Easy NetCDF)  !
! to distinguish them from NetCDF system calls                 !
! Error handling:
!      SUBROUTINE eznc_handle_err(status,text)
!
! File handling (open/close):
!      SUBROUTINE open_ncfile_new(filename,ncid)
!      SUBROUTINE open_ncfile_readonly(filename,ncid)
!      SUBROUTINE open_ncfile_readwrite(filename,ncid)
!      SUBROUTINE close_ncfile(ncid)
!
! File handling (sync data to avoid loss in case run terminates
!                unexpectedly):
!      SUBROUTINE sync_ncfile(ncid)
!
! File handling ("define" variable space then "data" for write/read):
!      SUBROUTINE switch_ncfile_to_define_mode(ncid)
!      SUBROUTINE switch_ncfile_to_data_mode(ncid)
!
! Attribute handling ("global" for file metadata)
!      SUBROUTINE eznc_def_sysglobatt(ncid,AttName,Command)
!      SUBROUTINE eznc_def_globalatt(ncid,AttName,AttTxt)
!      SUBROUTINE eznc_get_globalatt(ncid,AttName,AttTxt)
!      SUBROUTINE eznc_del_globalatt(ncid,AttName)
!      SUBROUTINE eznc_copy_globalatt(ncid_in,AttName,ncid_out)
!      SUBROUTINE eznc_rename_globalatt(ncid,OldName,NewName)
!
! Attribute handling ("local" for variable metadata)
!      SUBROUTINE eznc_def_localatt(ncid,VarName,AttName,AttTxt)
!      SUBROUTINE eznc_get_localatt(ncid,VarName,AttName,AttTxt)
!
! Dimension handling:
!      SUBROUTINE eznc_def_dim(ncid,DimName,dimval)
!      SUBROUTINE eznc_get_dimid(ncid,DimName,dimid)
!      SUBROUTINE eznc_get_dimension(ncid,DimName,dimval)
!
! Variable handling - scalars (define,put,get = declare,write,read):
!      SUBROUTINE eznc_def_0Dint(ncid,VarName)
!      SUBROUTINE eznc_put_0Dint(ncid,VarName,value)
!      SUBROUTINE eznc_get_0Dint(ncid,VarName,value)
!      SUBROUTINE eznc_def_0Dreal(ncid,VarName)
!      SUBROUTINE eznc_put_0Dreal(ncid,VarName,value)
!      SUBROUTINE eznc_get_0Dreal(ncid,VarName,value)
!
! Variable handling - arrays (define,put,get = declare,write,read):
!      SUBROUTINE eznc_def_1Dint(ncid,VarName,DimName)
!      SUBROUTINE eznc_put_1Dint(ncid,VarName,invals,start1,end1)
!      SUBROUTINE eznc_get_1Dint(ncid,VarName,dimval,outvals,
!     $                          start1,end1)
!      SUBROUTINE eznc_def_1Dreal(ncid,VarName,DimName)
!      SUBROUTINE eznc_put_1Dreal(ncid,VarName,invals,start1,end1)
!      SUBROUTINE eznc_get_1Dreal(ncid,VarName,dimval,outvals,
!     $                           start1,end1)
!      SUBROUTINE eznc_def_1Dchar(ncid,VarName,DimName1,DimName2)
!      SUBROUTINE eznc_put_1Dchar(ncid,VarName,values,
!     $                                        nchar,start1,end1)
!      SUBROUTINE eznc_get_1Dchar(ncid,VarName,nchar,dimval,
!     $                                values,start1,end1)
!      SUBROUTINE eznc_def_2Dint(ncid,VarName,DimName1,DimName2)
!      SUBROUTINE eznc_put_2Dint(ncid,VarName,invals,
!     $                          start1,end1,start2,end2)
!      SUBROUTINE eznc_get_2Dint(ncid,VarName,dim1,dim2,outvals,
!     $                          start1,end1,start2,end2)
!      SUBROUTINE eznc_def_2Dreal(ncid,VarName,DimName1,DimName2)
!      SUBROUTINE eznc_put_2Dreal(ncid,VarName,invals,
!     $                           start1,end1,start2,end2)
!      SUBROUTINE eznc_get_2Dreal(ncid,VarName,dim1,dim2,outvals,
!     $                           start1,end1,start2,end2)
!      SUBROUTINE eznc_def_2Ddbl(ncid,VarName,DimName1,DimName2)
!      SUBROUTINE eznc_def_2Dchar(ncid,VarName,DimName1,DimName2,DimName3)
!      SUBROUTINE eznc_put_2Dchar(ncid,VarName,values,
!     $                                nchar,start1,end1,start2,end2)
!      SUBROUTINE eznc_get_2Dchar(ncid,VarName,nchar,dimval,
!     $                                values,start1,end1,start2,end2)
!      SUBROUTINE eznc_def_3Dint(ncid,VarName,DimName1,DimName2,DimName3)
!      SUBROUTINE eznc_put_3Dint(ncid,VarName,invals,
!     $                          start1,end1,start2,end2,start3,end3)
!      SUBROUTINE eznc_get_3Dint(ncid,VarName,dim1,dim2,dim3,outvals,
!     $                         start1,end1,start2,end2,start3,end3)
!      SUBROUTINE eznc_def_3Dreal(ncid,VarName,DimName1,DimName2,DimName3)
!      SUBROUTINE eznc_put_3Dreal(ncid,VarName,invals,
!     $                           start1,end1,start2,end2,start3,end3)
!      SUBROUTINE eznc_get_3Dreal(ncid,VarName,dim1,dim2,dim3,outvals,
!     $                           start1,end1,start2,end2,start3,end3)
!
! NetCDF Fortan manual may be found at:                        !
! https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90/index.html#Top
!-------------------------------------------------------!

      SUBROUTINE eznc_handle_err(status,text)

! Return a NetCDF error code for failed NetCDF operation,
!   and terminate program.
! A dictionary of codes is available at:
! https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/NetCDF_002d3-Error-Codes.html#NetCDF_002d3-Error-Codes

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER status
      CHARACTER(*) text

! INTERNAL VARIABLES
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(16) errnam
      CHARACTER(60) report

      status = -(status)
      errnam = "              "    
      report = "                "

      SELECT CASE (status)
        CASE (0)
          errnam = "NC_NOERR        "    
          report = "No Error      "
          RETURN
        CASE (2)
          errnam = "xx_BADFILENAME  "    
          report = "No file found with that name      "
        CASE (33)
          errnam = "NC_EBADID       "
          report = "Not a netcdf id"
        CASE (34)
          errnam = "NC_ENFILE       "
          report = "Too many netcdfs open"
        CASE (35)
          errnam = "NC_EEXIST       "
          report = "netcdf file already exists (NC_NOCLOBBER)"
        CASE (36)
          errnam = "NC_EINVAL       "
          report = "Invalid Argument"
        CASE (37)
          errnam = "NC_EPERM        "
          report = "Write to read only"
        CASE (38)
          errnam = "NC_ENOTINDEFINE "
          report = "Operation not allowed in data mode"
        CASE (39)
          errnam = "NC_EINDEFINE    "
          report = "Operation not allowed in define mode"
        CASE (40)
          errnam = "NC_EINVALCOORDS "
          report = "Index exceeds dimension bound"
        CASE (41)
          errnam = "NC_EMAXDIMS     "
          report = "NC_MAX_DIMS exceeded"
        CASE (42)
          errnam = "NC_ENAMEINUSE   "
          report = "String match to name in use"
        CASE (43)
          errnam = "NC_ENOTATT      "
          report = "Attribute not found"
        CASE (44)
          errnam = "NC_EMAXATTS     "
          report = "Not a netcdf data type"
        CASE (46)           
          errnam = "NC_EBADDIM      "
          report = "Invalid dimension id or errnam"
        CASE (47)           
          errnam = "NC_EUNLIMPOS    "
          report = "NC_UNLIMITED in the wrong index"
        CASE (48)           
          errnam = "NC_EMAXVARS     "
          report = "NC_MAX_VARS exceeded"
        CASE (49)           
          errnam = "NC_ENOTVAR      "
          report = "Variable not found"
        CASE (50)           
          errnam = "NC_EGLOBAL      "
          report = "Action prohibited on NC_GLOBAL varid"
        CASE (51)           
          errnam = "NC_ENOTNC       "
          report = "Not a netcdf file"
        CASE (52)           
          errnam = "NC_ESTS         "
          report = "In Fortran, string too short"
        CASE (53)           
          errnam = "NC_EMAXNAME     "
          report = "NC_MAX_NAME exceeded"
        CASE (54)           
          errnam = "NC_EUNLIMIT     "
          report = "NC_UNLIMITED size already in use"
        CASE (55)           
          errnam = "NC_ENORECVARS   "
          report = "nc_rec op when there are no record vars"
        CASE (56)           
          errnam = "NC_ECHAR        "
          report = "Attempt to convert between text & numbers"
        CASE (57)           
          errnam = "NC_EEDGE        "
          report = "Edge+start exceeds dimension bound"
        CASE (58)           
          errnam = "NC_ESTRIDE      "
          report = "Illegal stride"
        CASE (59)           
          errnam = "NC_EBADNAME     "
          report ="Attribute or variable errnam contains illegal chars"
        CASE (60)           
          errnam = "NC_ERANGE       "
          report = "Math result not representable"
        CASE (61)           
          errnam = "NC_ENOMEM       "
          report = "Memory allocation (malloc) failure"
        CASE (62)           
          errnam = "NC_EVARSIZE     "
          report="One or more variable sizes violate format constraints"
        CASE (63)           
          errnam = "NC_EDIMSIZE     "
          report = "Invalid dimension size"
        CASE (64)           
          errnam = "NC_ETRUNC       "
          report = "File likely truncated or possibly corrupted"
        CASE (101)
          errnam = "NC_EHDFERR      "
        CASE (102)
          errnam = "NC_ECANTREAD    "
        CASE (103)
          errnam = "NC_ECANTWRITE   "
        CASE (104)
          errnam = "NC_ECANTCREATE  "
        CASE (105)
          errnam = "NC_EFILEMETA    "
        CASE (106)
          errnam = "NC_EDIMMETA     "
        CASE (107)
          errnam = "NC_EATTMETA     "
        CASE (108)
          errnam = "NC_EVARMETA     "
        CASE (109)
          errnam = "NC_ENOCOMPOUND  "
        CASE (110)
          errnam = "NC_EATTEXISTS   "
        CASE (111) 
          errnam = "NC_ENOTNC4      "
          report = "Attempting netcdf-4 operation on netcdf-3 file."
        CASE (112) 
          errnam = "NC_ESTRICTNC3   "
          report = 
     $    "Attempting netcdf-4 operation on strict nc3 netcdf-4 file."
        CASE (113) 
          errnam = "NC_EBADGRPID    "
          report = "Bad group id. Bad!"
        CASE (114) 
          errnam = "NC_EBADTYPEID   "
          report = "Bad type id."
        CASE (115) 
          errnam = "NC_EBADFIELDID  "
          report = "Bad field id."
        CASE (116)
          errnam = "NC_EUNKNAME     "
        CASE DEFAULT
          status = -(status)
          errnam = "NC_UNDEFINED_ERR"
          report = "Undefined error"
      END SELECT

      WRITE(lout,*) "NetCDF error: ",text,status,errnam," ",report
      PRINT*, "NetCDF error: ",text,status
      PRINT*, errnam," : ",report
      STOP

      END SUBROUTINE eznc_handle_err

* -----------------------------------------------------
      SUBROUTINE open_ncfile_new(filename,ncid)

! create a new netcdf file (replace any existing file)
! file is automatically created in DEFINE mode
!------------------------------------------------------
! cmode : OPTIONS for dealing with an existing file of same name:
! NF90_CLOBBER = overwrite (replace) existing file
! NF90_NOCLOBBER = return an error and stop IF file already exists
!-----------------
! other options :
! NF90_SHARE = allow different processes to read / write simulataneously
! NF90_64BIT_OFFSET = create 64-bit offest format file (useful for >2GB)
! NF90_HDF5 = output is in netCDF-4/HDF5 format
!-----------------

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      CHARACTER(*) filename
! OUTPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      print*,"open_ncfile_new : ",filename,ncid

      text   ="NF90_CREATE, "//filename
      status = NF90_CREATE(path=filename,
     &                     cmode=NF90_CLOBBER,
     &                     ncid=ncid)
!      status = NF90_CREATE(path=filename,
!     &                     cmode=OR(NF90_CLOBBER,NF90_64BIT_OFFSET),
!     &                     ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! file automatically CREATEd in DEFINE mode: no need to switch

      END SUBROUTINE open_ncfile_new

* -------------------------------
      SUBROUTINE open_ncfile_readonly(filename,ncid)

! open an existing netcdf file for reading only (in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      CHARACTER(*) filename
! OUTPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      text   ="NF90_OPEN, "//filename
      status = NF90_OPEN(path=filename,mode=NF90_NOWRITE,ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      print*,"open_ncfile_readonly : ",filename,ncid

! switch to DATA mode just in case it's not already there
      CALL switch_ncfile_to_data_mode(ncid)

      END SUBROUTINE open_ncfile_readonly

* -------------------------------
      SUBROUTINE open_ncfile_readwrite(filename,ncid)

! open an existing netcdf file for reading and writing (in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      CHARACTER(*) filename
! OUTPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      print*,"open_ncfile_readwrite : ",filename,ncid

      text   ="NF90_OPEN, "//filename
      status = NF90_OPEN(path=filename,mode=NF90_WRITE,ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! switch to DEFINE mode just in case it's not already there
      IF (status.EQ.-39) RETURN  ! already in DEFINE mode, no action needed
      CALL switch_ncfile_to_define_mode(ncid)

      END SUBROUTINE open_ncfile_readwrite

* -------------------------------
      SUBROUTINE open_ncfile_writedata(filename,ncid)

! open an existing netcdf file for writing data (in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      CHARACTER(*) filename
! OUTPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      print*,"open_ncfile_writedata : ",filename,ncid

      text   ="NF90_OPEN, "//filename
      status = NF90_OPEN(path=filename,mode=NF90_WRITE,ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! switch to DATA mode if not already there
      IF (status.EQ.-38) RETURN  ! already in DATA mode, no action needed
      CALL switch_ncfile_to_data_mode(ncid)

      END SUBROUTINE open_ncfile_writedata

* -------------------------------
      SUBROUTINE sync_ncfile(ncid)

! sync data to an open netcdf file 

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

!      print*,"sync_ncfile : ",ncid

      text   ="NF90_SYNC"
      status = NF90_SYNC(ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE sync_ncfile

* -------------------------------
      SUBROUTINE close_ncfile(ncid)

! close an open netcdf file 

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      print*,"close_ncfile : ",ncid

      text   ="NF90_CLOSE"
      status = NF90_CLOSE(ncid=ncid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE close_ncfile

* -------------------------------
      SUBROUTINE switch_ncfile_to_data_mode(ncid)

! switch an open netcdf file to DATA mode, 
! allowing values to be written

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      print*,"switch_ncfile_to_data_mode : ",ncid

      text   ="NF90_ENDDEF"
      status = NF90_ENDDEF(ncid=ncid)
      IF (status.EQ.-38) RETURN  ! already in DATA mode, no action needed
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE switch_ncfile_to_data_mode

* -------------------------------
      SUBROUTINE switch_ncfile_to_define_mode(ncid)

! switch an open netcdf file to DEFINE mode, 
! allowing array dimensions or new variables to be specifed

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(25) text

      !print*,"switch_ncfile_to_define_mode : ",ncid

      text   ="NF90_REDEF"
      status = NF90_REDEF(ncid=ncid)
      IF (status.EQ.-39) RETURN  ! already in DEFINE mode, no action needed
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE switch_ncfile_to_define_mode

* -------------------------------
      SUBROUTINE eznc_def_globalatt(ncid,AttName,AttTxt)

! define global attribute using specified text ("AttTxt")
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) AttName
      CHARACTER(*) AttTxt
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(40) text

      text   ="NF90_PUT_ATT "//AttName
      status = NF90_PUT_ATT(ncid,NF90_GLOBAL,AttName,AttTxt)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,AttName,AttTxt)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_globalatt

* -------------------------------
      SUBROUTINE eznc_def_sysglobatt(ncid,AttName,Command)

! define global attribute using results of a system call
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) Command
      !CHARACTER(*) SysCall
      CHARACTER(*) AttName
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(60) AttTxt
      CHARACTER(60) text
      CHARACTER(60) SysCall

      SysCall = 'echo '//Command//' > dummy.99'
      CALL SYSTEM(SysCall)

      OPEN(99,file='dummy.99')
        READ(99,'(a60)') AttTxt
        PRINT*,"AttTxt = ",AttTxt
      CLOSE(99)

      text   ="NF90_PUT_ATT  "//AttName
      status = NF90_PUT_ATT(ncid,NF90_GLOBAL,AttName,AttTxt)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_PUT_ATT(ncid,NF90_GLOBAL,AttName,AttTxt)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_sysglobatt

* -------------------------------
      SUBROUTINE eznc_copy_globalatt(ncid_in,AttName,ncid_out)

! copy global attribute from one dataset to another
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid_in,ncid_out
      CHARACTER(*) AttName
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(40) text

! need a checker to make sure variable exists...

      text   ="NF90_COPY_ATT "//AttName
      status = NF90_COPY_ATT(ncid_in,NF90_GLOBAL,AttName,
     &                       ncid_out,NF90_GLOBAL)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid_out)
        status = NF90_COPY_ATT(ncid_in,NF90_GLOBAL,AttName,
     &                         ncid_out,NF90_GLOBAL)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_copy_globalatt

* -------------------------------
      SUBROUTINE eznc_rename_globalatt(ncid,OldName,NewName)

! rename existing global attribute 
! note that dataset must be in DEFINE mode to do this, unless
! len(Name_old) /< len(Name_new)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) OldName
      CHARACTER(*) NewName
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(40) text

10    CONTINUE

      text   ="NF90_RENAME_ATT "//NewName
      status = NF90_RENAME_ATT(ncid,NF90_GLOBAL,OldName,NewName)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        GOTO 10
! if "new" attribute name is already in use, delete it
      ELSE IF (status.EQ.-42) THEN  
        status = NF90_DEL_ATT(ncid,NF90_GLOBAL,NewName)
        GOTO 10
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_rename_globalatt

* -------------------------------
      SUBROUTINE eznc_del_globalatt(ncid,AttName)

! delete an existing global attribute to allow writing of a new one
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) AttName
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(40) text

      text   ="NF90_DEL_ATT "//AttName
      status = NF90_DEL_ATT(ncid,NF90_GLOBAL,AttName)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_DEL_ATT(ncid,NF90_GLOBAL,AttName)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_del_globalatt

* -------------------------------
      SUBROUTINE eznc_def_localatt(ncid,VarName,AttName,AttTxt)

! define variable-specific attribute using specified text ("AttTxt")
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
      CHARACTER(*) AttName
      CHARACTER(*) AttTxt
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(40) text

! find variable ID from its name
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! assuming variable exists, apply attribute to it
      text   ="NF90_PUT_ATT "//AttName
      status = NF90_PUT_ATT(ncid,varid,AttName,AttTxt)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_PUT_ATT(ncid,varid,AttName,AttTxt)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_localatt

* -------------------------------
      SUBROUTINE eznc_def_dim(ncid,DimName,dim_in)

! define a dimension, given its name and value
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER,INTENT(IN) :: dim_in
      CHARACTER(*) DimName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER dimid
      INTEGER dimval
      CHARACTER(25) text

      IF(dim_in.EQ.0) THEN
         dimval = 1
      ELSE
         dimval = dim_in
      ENDIF
      print*,"eznc_def_dim:",DimName,DimVal

! write dimension
      text   ="NF90_DEF_DIM "//DimName
      status = NF90_DEF_DIM(ncid,DimName,dimval,dimid)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_DEF_DIM(ncid,DimName,dimval,dimid)
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_dim

* -------------------------------
      SUBROUTINE eznc_def_dim_unlim(ncid,DimName)

! define an UNLIMITED dimension, given its name
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) DimName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER dimid
      CHARACTER(25) text

      print*,"eznc_def_dim_unlim:",DimName

! write dimension
      text   ="NF90_DEF_DIM "//DimName
      status = NF90_DEF_DIM(ncid,DimName,NF90_UNLIMITED,dimid)

! if file is in DATA mode, switch to DEFINE mode and try again
      IF (status.EQ.-38) THEN  
        CALL switch_ncfile_to_define_mode(ncid)
        status = NF90_DEF_DIM(ncid,DimName,NF90_UNLIMITED,dimid)
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_dim_unlim

* -------------------------------
      SUBROUTINE eznc_def_0Dint(ncid,VarName)

! define a scalar integer, given its name
! dataset must be in DEFINE mode
! write the value later, using subroutine ncput_0Dint 

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(25) text

      !print*,"eznc_def_0Dint:",VarName

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_INT,varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_INT,varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_0Dint

* -------------------------------
      SUBROUTINE eznc_def_0Dreal(ncid,VarName)

! define a scalar real, given its name
! dataset must be in DEFINE mode
! write the value later, using subroutine ncput_0Dreal

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(25) text

      !print*,"eznc_def_0Dreal:",VarName

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_FLOAT,varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_FLOAT,varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_0Dreal

* -------------------------------
      SUBROUTINE eznc_def_1Dint(ncid,VarName,DimName)

! define a 1-D integer array, given its name and dimension name 
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_1Dint

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
      CHARACTER(*) DimName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid
      CHARACTER(25) text

      !print*,"eznc_def_1Dint:",VarName,"(",DimName,")"

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimid of named dimension
        text   ="NF90_INQ_DIMID "//DimName
        status = NF90_INQ_DIMID(ncid,DimName,dimid)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_INT,(/dimid/),varid)

! > if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_INT,(/dimid/),varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_1Dint

* -------------------------------
      SUBROUTINE eznc_def_1Dreal(ncid,VarName,DimName)

! define a 1-D array of reals
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_1Dreal

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
      CHARACTER(*) DimName
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid
      CHARACTER(25) text

      !print*,"eznc_def_1Dreal:",VarName,"(",DimName,")"

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimid of named dimension
        text   ="NF90_INQ_DIMID "//DimName
        status = NF90_INQ_DIMID(ncid,DimName,dimid)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_FLOAT,(/dimid/),varid)

! > if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_FLOAT,(/dimid/),varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_1Dreal

* -------------------------------
      SUBROUTINE eznc_def_1Dchar(ncid,VarName,DimName1,DimName2)

! define a 1-D array (size Dim2) of character strings with length Dim1
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2
      CHARACTER(25) text

      !print*,"eznc_def_1Dchar:",VarName

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_CHAR,
     $                       (/dimid1,dimid2/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_CHAR,
     $                         (/dimid1,dimid2/),varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_1Dchar

* -------------------------------
      SUBROUTINE eznc_def_2Dint(ncid,VarName,DimName1,DimName2)

! define a 2-D array of integers, given name and dimension names
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_2Dint

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2
      CHARACTER(25) text


      !print*,"eznc_def_2Dint:",VarName,"(",DimName1,",",DimName2,")"

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_INT,
     &                       (/dimid1,dimid2/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_INT,
     &                         (/dimid1,dimid2/),varid)
        ENDIF
      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_2Dint

* -------------------------------
      SUBROUTINE eznc_def_2Dreal(ncid,VarName,DimName1,DimName2)

! define a 2-D array of reals, given name and dimension names
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_2Dreal

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2
      CHARACTER(25) text

      !print*,"eznc_def_2Dreal:",VarName,"(",DimName1,",",DimName2,")"

! check if variable already defined. 
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_REAL,
     &                       (/dimid1,dimid2/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_REAL,
     &                         (/dimid1,dimid2/),varid)
        ENDIF
      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_2Dreal

* -------------------------------
      SUBROUTINE eznc_def_2Ddbl(ncid,VarName,DimName1,DimName2)

! USE SPARINGLY !
! define a 2-D array of DP values, given name and dimension names
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_2Ddbl

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2
      CHARACTER(25) text

      !print*,"eznc_def_2Ddbl:",VarName

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_DOUBLE,
     &                       (/dimid1,dimid2/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_DOUBLE,
     &                         (/dimid1,dimid2/),varid)
        ENDIF
      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_2Ddbl

* -------------------------------
      SUBROUTINE eznc_def_2Dchar(ncid,VarName,
     &                           DimName1,DimName2,DimName3)

! define a 2-D array (size Dim2,Dim3) of character strings with length Dim1
! (dataset must be in DEFINE mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2,DimName3
! OUTPUT VARIABLES
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2,dimid3
      CHARACTER(25) text

      !print*,"eznc_def_2Dchar:",VarName

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName3
        status = NF90_INQ_DIMID(ncid,DimName3,dimid3)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_CHAR,
     $                       (/dimid1,dimid2,dimid3/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_CHAR,
     $                         (/dimid1,dimid2,dimid3/),varid)
        ENDIF

      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_2Dchar

* -------------------------------
      SUBROUTINE eznc_def_3Dint(ncid,VarName,
     &                          DimName1,DimName2,DimName3)

! define a 3-D array of integers, given name and dimension names
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_2Dint

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2,DimName3
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2,dimid3
      CHARACTER(25) text


      !print*,"eznc_def_3Dint:",VarName,"(",DimName1,",",DimName2,
!     &                                     DimName3,")"

! check if variable already defined. If yes, return
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName3
        status = NF90_INQ_DIMID(ncid,DimName3,dimid3)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_INT,
     &                       (/dimid1,dimid2,dimid3/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_INT,
     &                         (/dimid1,dimid2,dimid3/),varid)
        ENDIF
      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_3Dint

* -------------------------------
      SUBROUTINE eznc_def_3Dreal(ncid,VarName,
     &                           DimName1,DimName2,DimName3)

! define a 3-D array of reals, given name and dimension names
! dataset must be in DEFINE mode
! write the values later, using subroutine eznc_put_2Dreal

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName,DimName1,DimName2,DimName3
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER dimid1,dimid2,dimid3
      CHARACTER(25) text

!      print*,"eznc_def_3Dreal:",VarName,"(",DimName1,",",DimName2,",",
!     &                                      DimName3,")"

! check if variable already defined. 
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)

! if variable not found, define it
      IF (status.EQ.-49) THEN 

! > find dimids of named dimensions
        text   ="NF90_INQ_DIMID "//DimName1
        status = NF90_INQ_DIMID(ncid,DimName1,dimid1)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName2
        status = NF90_INQ_DIMID(ncid,DimName2,dimid2)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

        text   ="NF90_INQ_DIMID "//DimName3
        status = NF90_INQ_DIMID(ncid,DimName3,dimid3)
        IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! > define variable
        text   ="NF90_DEF_VAR "//VarName
        status = NF90_DEF_VAR(ncid,VarName,NF90_REAL,
     &                       (/dimid1,dimid2,dimid3/),varid)

! if file is in DATA mode, switch to DEFINE mode and try again
        IF (status.EQ.-38) THEN  
          CALL switch_ncfile_to_define_mode(ncid)
          status = NF90_DEF_VAR(ncid,VarName,NF90_REAL,
     &                         (/dimid1,dimid2,dimid3/),varid)
        ENDIF
      ENDIF

! Handle any other inquiry errors
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_def_3Dreal

* -------------------------------
      SUBROUTINE eznc_put_0Dint(ncid,VarName,value)

! write a value to a previously-defined scalar integer
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER value
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(25) text

      !print*,"eznc_put_0Dint:",VarName,value

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,value)

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,value)
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_0Dint

* -------------------------------
      SUBROUTINE eznc_put_0Dreal(ncid,VarName,value)

! write a value to a previously-defined scalar real number
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      REAL    value
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(25) text

      !print*,"eznc_put_0Dreal:",VarName,value

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,value)

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,value)
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_0Dreal

* -------------------------------
      SUBROUTINE eznc_put_1Dint(ncid,VarName,invals,start1,end1)

! write values to a previously-defined 1-D integer array 
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,end1,invals(end1-start1+1)
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      CHARACTER(25) text

      !print*,"eznc_put_1Dint:",VarName

      count1=end1-start1+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1/),
     $                      count=(/count1/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals(start1:end1),
     $                        start=(/start1/),
     $                        count=(/count1/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dint

* -------------------------------
      SUBROUTINE eznc_put_1Dreal(ncid,VarName,invals,start1,end1)

! write values to a previously-defined 1-D real array 
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,end1
      REAL,DIMENSION(end1-start1+1) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      CHARACTER(25) text

      !print*,"eznc_put_1Dreal:",VarName !,invals

      count1=end1-start1+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals(start1:end1),
     $                      start=(/start1/),
     $                      count=(/count1/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals(start1:end1),
     $                        start=(/start1/),
     $                        count=(/count1/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dreal

* -------------------------------
      SUBROUTINE eznc_put_1Dchar(ncid,VarName,values,nchar,start1,end1)

! write values to a previously-defined 1-D array of string variables
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER nchar,start1,end1
      CHARACTER(nchar) values(end1-start1+1)
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      CHARACTER(25) text

      !print*,"eznc_put_1Dchar:",VarName
      count1=end1-start1+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,values,
     $                      start=(/1,start1/),
     $                      count=(/nchar,count1/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,values,
     $                        start=(/1,start1/),
     $                        count=(/nchar,count1/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dchar

* -------------------------------
      SUBROUTINE eznc_put_2Dint(ncid,VarName,invals,
     $                       start1,end1,start2,end2)

! write values to a previously-defined 2-D integer array
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,start2,end1,end2
      INTEGER,DIMENSION((end1-start1+1),(start2-end2+1)) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      CHARACTER(25) text

      !print*,"eznc_put_2Dint:",VarName

      count1=end1-start1+1
      count2=end2-start2+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals,
     $                        start=(/start1,start2/),
     $                        count=(/count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_2Dint

* -------------------------------
      SUBROUTINE eznc_put_2Dreal(ncid,VarName,invals,
     $                       start1,end1,start2,end2)

! write values to a previously-defined 2-D real array
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,start2,end1,end2
      REAL,DIMENSION((end1-start1+1),(start2-end2+1)) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      CHARACTER(25) text

      !print*,"eznc_put_2Dreal:",VarName

      count1=end1-start1+1
      count2=end2-start2+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals,
     $                        start=(/start1,start2/),
     $                        count=(/count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_2Dreal

* -------------------------------
      SUBROUTINE eznc_put_2Dchar(ncid,VarName,values,nchar,
     &                           start1,end1,start2,end2)

! write values to a previously-defined 1-D array of string variables
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER nchar,start1,end1,start2,end2
      CHARACTER(nchar) values(end1-start1+1,end2-start2+1)
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      CHARACTER(25) text

      !print*,"eznc_put_2Dchar:",VarName
      count1=end1-start1+1
      count2=end2-start2+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,values,
     $                        start=(/1,start1,start2/),
     $                      count=(/nchar,count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,values,
     $                        start=(/1,start1,start2/),
     $                      count=(/nchar,count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_2Dchar

* -------------------------------
      SUBROUTINE eznc_put_3Dint(ncid,VarName,invals,
     $                       start1,end1,start2,end2,start3,end3)

! write values to a previously-defined 3-D integer array
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,start2,start3,end1,end2,end3
      INTEGER,DIMENSION((end1-start1+1),(start2-end2+1),
     &                  (start3-end3+1)) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2,count3
      CHARACTER(25) text

      !print*,"eznc_put_3Dint:",VarName

      count1=end1-start1+1
      count2=end2-start2+1
      count3=end3-start3+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1,start2,start3/),
     $                      count=(/count1,count2,count3/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals,
     $                        start=(/start1,start2,start3/),
     $                        count=(/count1,count2,count3/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_3Dint

* -------------------------------
      SUBROUTINE eznc_put_3Dreal(ncid,VarName,invals,
     $                       start1,end1,start2,end2,start3,end3)

! write values to a previously-defined 3-D real array
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,start2,start3,end1,end2,end3
      REAL,DIMENSION((end1-start1+1),(start2-end2+1),(start3-end3+1)) 
     &                        :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2,count3
      CHARACTER(25) text

      !print*,"eznc_put_3Dreal:",VarName

      count1=end1-start1+1
      count2=end2-start2+1
      count3=end3-start3+1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1,start2,start3/),
     $                      count=(/count1,count2,count3/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,invals,
     $                      start=(/start1,start2,start3/),
     $                      count=(/count1,count2,count3/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_3Dreal

* -------------------------------
* -------------------------------
!      SUBROUTINE eznc_put_2Ddbl(........
!
!! CAUTION !! DO NOT USE - DEFINE AS NF90_DOUBLE BUT PUT AS "REAL"
!             UNDER COMPILATION OPTION real-8
!
* -------------------------------
      SUBROUTINE eznc_put_0Dreal_into1D(ncid,VarName,inval,start1)

! write a scalar value to a previously-defined 1-D real array 
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1
      REAL :: inval
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
!      REAL,DIMENSION(start1) :: ncvals
      CHARACTER(25) text

      !print*,"eznc_put_0Dreal_into1D:",VarName,inval

! convert scalar into 1-D array entry
!      ncvals(start1)=inval
      count1=1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,(/inval/),
     $                      start=(/start1/),
     $                      count=(/count1/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,(/inval/),
     $                        start=(/start1/),
     $                        count=(/count1/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_0Dreal_into1D

* -------------------------------
      SUBROUTINE eznc_put_0Dreal_into2D(ncid,VarName,inval,
     &                                       start1,start2)

! write a scalar value to a previously-defined 2-D real array 
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,start2
      REAL :: inval
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      CHARACTER(25) text

      !print*,"eznc_put_0Dreal_into2D:",VarName,inval

! convert scalar into 2-D array entry
      count1=1
      count2=1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,(/inval/),
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,(/inval/),
     $                        start=(/start1,start2/),
     $                        count=(/count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_0Dreal_into2D

* -------------------------------
      SUBROUTINE eznc_put_1Dreal_into2D(ncid,VarName,invals,
     $                               start1,end1,start2)

! write 1D REAL data slab into a previously-defined 2-D real array
! as ncvals(islab,start:end)
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,end1,start2
      REAL,DIMENSION(end1-start1+1) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      !REAL,DIMENSION(end1-start1+1,start2) :: ncvals
      REAL,DIMENSION(end1,1) :: ncvals
      CHARACTER(25) text

      !print*,"eznc_put_1Dreal_into2D:",VarName
      !print*,"invals",invals
      !print*,"start1,end1,start2",start1,end1,start2

! convert 1-D array into 2-D array slab of dimensions (count1,1)
      ncvals(start1:end1,1) = invals

      count1=end1-start1+1
      count2=1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,
     $                      ncvals(start1:end1,1),
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,
     $                        ncvals(start1:end1,1),
     $                        start=(/start1,start2/),
     $                        count=(/count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dreal_into2D

* -------------------------------
      SUBROUTINE eznc_put_1Dreal_into2D_mask(ncid,VarName,invals,
     $                               start1,end1,start2,maskval)

! write 1D REAL data slab into a previously-defined 2-D real array
! as ncvals(islab,start:end)
! (dataset must be in DATA mode)
! Writes ZERO ("0.") for values under threshold "maskval"

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,end1,start2
      REAL,DIMENSION(end1-start1+1) :: invals
      REAL    maskval
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      REAL,DIMENSION(end1,1) :: ncvals
      CHARACTER(25) text

      !print*,"eznc_put_1Dreal_into2D_mask:",VarName,varid,maskval
      !print*,invals

! zero out undesirable values from 1-D input array
      WHERE (invals.LT.maskval) invals=0.

! convert 1-D array into 2-D array slab of dimensions (count1,1)
      ncvals(start1:end1,1) = invals

      count1=end1-start1+1
      count2=1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid,
     $                      ncvals(start1:end1,1),
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,
     $                        ncvals(start1:end1,1),
     $                        start=(/start1,start2/),
     $                        count=(/count1,count2/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dreal_into2D_mask

* -------------------------------
      SUBROUTINE eznc_put_1Dreal_into3D(ncid,VarName,invals,
     $                               start1,end1,start2,start3)

! write 1D REAL data slab into a previously-defined s32-D real array
! as ncvals(islab,start:end)
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      INTEGER start1,end1,start2,start3
      REAL,DIMENSION(end1-start1+1) :: invals
      CHARACTER(*) VarName
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2,count3
      REAL,DIMENSION(end1,1,1) :: ncvals
      CHARACTER(25) text

      !print*,"eznc_put_1Dreal_into3D:",VarName
      !print*,invals

! convert 1-D array into 3-D array slab of dimensions (count1,1,1)
      ncvals(start1:end1,1,1) = invals

      count1=end1-start1+1
      count2=1
      count3=1

! find variable ID
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! write variable
      text   ="NF90_PUT_VAR "//VarName
      status = NF90_PUT_VAR(ncid,varid, 
     $                  ncvals(start1:end1,1,1),
     $                  start=(/start1,start2,start3/),
     $                  count=(/count1,count2,count3/))

! if file is in DEFINE mode, switch to DATA mode and try again
      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_PUT_VAR(ncid,varid,
     $                  ncvals(start1:end1,1,1),
     $                  start=(/start1,start2,start3/),
     $                  count=(/count1,count2,count3/))
      ENDIF

      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_put_1Dreal_into3D
* -------------------------------
      SUBROUTINE eznc_get_globalatt(ncid,AttName,AttTxt)

! retrieve text of global attribute using attribute name

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) AttName
      CHARACTER(*) AttTxt
! INTERNAL VARIABLES
      INTEGER status
      CHARACTER(40) text

      !print*,"eznc_get_globalatt : ",AttName,ncid

! find global attribute
      text   ="NF90_GET_ATT "//Attname
      status = NF90_GET_ATT(ncid,NF90_GLOBAL,AttName,AttTxt)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_get_globalatt

* -------------------------------
      SUBROUTINE eznc_get_localatt(ncid,VarName,AttName,AttTxt)

! retrieve attribute text using variable and attribute names

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
      CHARACTER(*) AttName
      CHARACTER(*) AttTxt
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(40) text

      !print*,"eznc_get_localatt : ",VarName," ",AttName,ncid

! find variable ID from its name
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

! assuming variable exists, find its attribute
      text   ="NF90_GET_ATT "//AttName
      status = NF90_GET_ATT(ncid,varid,AttName,AttTxt)
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      END SUBROUTINE eznc_get_localatt

* -------------------------------
      SUBROUTINE eznc_get_dimension(ncid,DimName,dimval)

! read the value of a dimension whose name is already known
! (dataset can be in either DEFINE or DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER :: ncid
      CHARACTER(*) :: DimName
! OUTPUT VARIABLES
      INTEGER :: dimval
! INTERNAL VARIABLES
      INTEGER :: status
      INTEGER :: dimid
      CHARACTER(40) :: text

      text   ="NF90_INQ_DIMID "//DimName
      status = NF90_INQ_DIMID(ncid,DimName,dimid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      text   ="NF90_INQUIRE_DIMENSION "//DimName
      status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=dimval)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      !print*,"eznc_get_dimension: ",DimName,dimval

      END SUBROUTINE eznc_get_dimension

* -------------------------------
      SUBROUTINE eznc_get_dimid(ncid,DimName,dimid)

! read the I.D (only) of a dimension whose name is already known
! (dataset can be in either DEFINE or DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER :: ncid
      CHARACTER(*) :: DimName
! OUTPUT VARIABLES
      INTEGER :: dimid
! INTERNAL VARIABLES
      INTEGER :: status
      CHARACTER(40) :: text


      text   ="NF90_INQ_DIMID "//DimName
      status = NF90_INQ_DIMID(ncid,DimName,dimid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      !print*,"eznc_get_dimid: ", DimName,": dimid=",dimid

      END SUBROUTINE eznc_get_dimid

* -------------------------------
      SUBROUTINE eznc_get_0Dint(ncid,VarName,value)

! read the value of a scalar integer whose name is already known
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      INTEGER value
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(40) text

      !print*,"eznc_get_0Dint: ",VarName

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,value)

      IF (status.EQ.-39) THEN  
        CALL switch_ncfile_to_data_mode(ncid)
        status = NF90_GET_VAR(ncid,varid,value)
      ENDIF

      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      !print*,"VarName=",VarName,": varid=",varid,": Value=",value
      !print*,"eznc_get_0Dint: ",VarName,value

      END SUBROUTINE eznc_get_0Dint

* -------------------------------
      SUBROUTINE eznc_get_0Dreal(ncid,VarName,value)

! read the value of a scalar integer whose name is already known
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      REAL value
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      CHARACTER(40) text

      !print*,"eznc_get_0Dreal: ",VarName

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,value)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      !print*,"VarName=",VarName,": varid=",varid,": Value=",value
      !print*,"eznc_get_0Dreal: ",VarName,value

      END SUBROUTINE eznc_get_0Dreal

* -------------------------------
      SUBROUTINE eznc_get_1Dint(ncid,VarName,dimval,outvals,
     $                     start1,end1)
! read a 1-D integer array whose name & dimension are already known
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,dimval,start1,end1
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      INTEGER,DIMENSION(end1-start1+1) :: outvals
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      INTEGER,DIMENSION(dimval) :: values
      CHARACTER(40) text

      !print*,"eznc_get_1Dint: ",VarName

      values(:) = 0
      outvals(:) = 0
      count1=end1-start1+1

      IF(count1.LT.1) RETURN

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,values(start1:end1),
     $                      start=(/start1/),
     $                      count=(/count1/))
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      outvals=values(start1:end1)
      !print*,"VarName=",VarName,": varid=",varid,": values=",outvals
      !print*,outvals

      END SUBROUTINE eznc_get_1Dint

* -------------------------------
      SUBROUTINE eznc_get_1Dreal(ncid,VarName,dimval,outvals,
     $                           start1,end1)

! read a 1-D real array whose name & dimension are already known
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,dimval,start1,end1
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      REAL,DIMENSION(end1-start1+1) :: outvals
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      REAL,DIMENSION(dimval) :: values
      CHARACTER(40) text

      !print*,"eznc_get_1Dreal: ",VarName

      values(:) = 0.
      outvals(:) = 0.
      count1=end1-start1+1

      IF(count1.LT.1) RETURN

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,values(start1:end1),
     $                      start=(/start1/),
     $                      count=(/count1/))
      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)

      outvals=values(start1:end1)

      !print*,"VarName=",VarName,": varid=",varid,": values=",outvals
      !print*,outvals

      END SUBROUTINE eznc_get_1Dreal

* -------------------------------
      SUBROUTINE eznc_get_1Dchar(ncid,VarName,nchar,dimval,
     $                                values,start1,end1)

! read a 1-D array of character strings, 
! whose name & dimensions are already known
! NOTE THAT THIS SUBROUTINE EXPECTS TO START AT LIST ITEM (1)
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,nchar,dimval,start1,end1
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      CHARACTER(nchar) values(dimval)
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1
      INTEGER,PARAMETER :: firstchar = 1
      CHARACTER(40) text

      !print*,"eznc_get_1Dchar: ",VarName

      values(:) = ""
      count1=end1-start1+1

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,values,
     $                      start=(/firstchar,start1/),
     $                      count=(/nchar,count1/))
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)


      !print*,"VarName=",VarName,": varid=",varid,": values=",values
      !print*,values

      END SUBROUTINE eznc_get_1Dchar

* -------------------------------
      SUBROUTINE eznc_get_2Dint(ncid,VarName,dim1,dim2,outvals,
     $                     start1,end1,start2,end2)

! read a 2-D integer array whose name & dimensions are already known
! dimension values are used in variable declaration
! THIS FORM USED FOR RETRIEVING SUBSECTION OF STORED ARRAY
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,dim1,dim2
      INTEGER start1,start2
      INTEGER end1,end2
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      INTEGER,DIMENSION((end1-start1+1),(end2-start2+1)) :: outvals
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      INTEGER,DIMENSION(dim1,dim2) :: values
      CHARACTER(40) text
      INTEGER, PARAMETER ::  stride1=1, stride2=1

      !print*,"eznc_get_2Dint: ",VarName

      values(:,:) = 0
      outvals(:,:) = 0
      count1=end1-start1+1
      count2=end2-start2+1

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      STATUS = NF90_INQ_VARID(ncid,VarName,varid)
      IF (STATUS/=NF90_NOERR) call eznc_handle_err(STATUS,text)

      PRINT*, 'eznc_get_2Dint:', varname ,count1,count2

      IF(count1.LT.1.OR.count2.LT.1) RETURN

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,
     $                      values(start1:end1,start2:end2),
     $                      start=(/start1,start2/),
     $                      count=(/count1,count2/),
     $                      stride=(/stride1,stride2/))
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)
      outvals(1:count1,1:count2)=values(start1:end1,start2:end2)

      END SUBROUTINE eznc_get_2Dint

* -------------------------------
      SUBROUTINE eznc_get_2Dreal(ncid,VarName,dim1,dim2,outvals,
     $                     start1,end1,start2,end2)

! read a 2-D real array whose name & dimensions are already known
! THIS FORM USED FOR RETRIEVING SUBSECTION OF STORED ARRAY
! "start" and "end" refer to the indices of the STORED array, 
! not the subsection.
!!!! Remember that the dimensions are quoted IN REVERSE ORDER 
!    in the NetCDF file !!!
! (dataset must be in DATA mode)

      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,dim1,dim2
      INTEGER start1,start2
      INTEGER end1,end2
      CHARACTER(*) VarName
! OUTPUT VARIABLES
      REAL,DIMENSION((end1-start1+1),(end2-start2+1)) :: outvals
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1,count2
      REAL,DIMENSION(dim1,dim2) :: values
      CHARACTER(40) text
      INTEGER tend, tcount, tstart
      INTEGER, PARAMETER :: stride1=1, stride2=1


      values(:,:) = 0
      outvals(:,:) = 0
      count1=end1-start1+1
      count2=end2-start2+1

      PRINT*,"eznc_get_2Dreal:",VarName,count1,count2

      IF(count1.LT.1.OR.count2.LT.1) RETURN

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,
     $                    values(start1:end1,start2:end2),
     $                    start=(/start1,start2/),
     $                    count=(/count1,count2/))
!     $                    stride=(/stride1,stride2/))
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      outvals(1:count1,1:count2)=values(start1:end1,start2:end2)

      !print*,"VarName=",VarName,": varid=",varid,": values=",outvals
      !print*,"VarName=",VarName,": starts=",start1,start2,": values="
      !print*,outvals

      END SUBROUTINE eznc_get_2Dreal

* ----------------------------
      SUBROUTINE eznc_get_3Dreal(ncid, VarName, dim1, dim2, dim3, 
     &          outvals,start1, end1, start2, end2, start3, end3)
! read a 3-D real array whose name & dimensions are already known
! THIS FORM USED FOR RETRIEVING SUBSECTION OF STORED ARRAY
! "start" and "end" refer to the indices of the STORED array, 
! not the subsection.
!!!! Remember that the dimensions are quoted IN REVERSE ORDER 
!    in the NetCDF file !!!
! (dataset must be in DATA mode)
      USE netcdf
      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER ncid,dim1,dim2, dim3
      INTEGER start1,start2, start3
      INTEGER end1,end2, end3
      CHARACTER(*) VarName

! OUTPUT VARIABLES
      REAL,DIMENSION(end1-start1+1,end2-start2+1,end3-start3+1) ::
     &               outvals
! INTERNAL VARIABLES
      INTEGER status
      INTEGER varid
      INTEGER count1, count2, count3
      CHARACTER(40) text
      
      print*,"eznc_get_3Dreal: ",VarName
      count1=end1-start1+1
      count2=end2-start2+1
      count3=end3-start3+1
      print*,start1,start2,start3
      print*,count1,count2,count3
      print*,dim1,dim2,dim3
      
      IF(count1.LT.1.OR.count2.LT.1.or.count3.LT.1) RETURN

      outvals(:,:,:) = 0

! find ID of variable
      text   ="NF90_INQ_VARID "//VarName
      status = NF90_INQ_VARID(ncid,VarName,varid)
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

! read variable
      text   ="NF90_GET_VAR "//VarName
      status = NF90_GET_VAR(ncid,varid,
     $                      outvals,
     $                      start=(/start1,start2,start3/),
     $                      count=(/count1,count2,count3/))
      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)

      
      END SUBROUTINE eznc_get_3Dreal

* -------------------------------
!-------------------------------------------------------!
!     BUFFER: USEFUL CODE SNIPPETS
!             USEFUL LINES FOR DEBUGGING
!-------------------------------------------------------!
* -------------------------------
!!-----find dimension ID------!
!      text   ="NF90_INQ_DIMID "//DimName
!      status = NF90_INQ_DIMID(ncid,DimName,dimid)
!      IF (status/=NF90_NOERR) CALL eznc_handle_err(status,text)
* -------------------------------
* -------------------------------
! TROUBLESHOOTING: 
!!-----investigate dimensions of known variable---------!
! VARIABLES FOR TROUBLESHOOTING
!      INTEGER :: numdims,dimsize,i
!      INTEGER,PARAMETER :: maxdims = 4, maxchars = 20
!      INTEGER,DIMENSION(maxdims) :: dimids
!      CHARACTER(maxchars) DimName
!
!      text   ="NF90_INQUIRE_VARIABLE "//VarName
!      status = NF90_INQUIRE_VARIABLE(ncid,varid,ndims=numDims,natts=i)
!      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)
!      text   ="NF90_INQUIRE_VARIABLE "//VarName
!      status = NF90_INQUIRE_VARIABLE(ncid,varid,dimids=dimids(:numDims))
!      IF (status/=NF90_NOERR) call eznc_handle_err(status,text)
!      DO i = 1,numdims
!        text   ="NF90_INQUIRE_DIMENSION "
!        status = NF90_INQUIRE_DIMENSION(ncid,dimids(i),
!     &                                  name=DimName,len=DimSize)
!        !print*,"dimension ",i,"=",dimname,dimsize
!      ENDDO

!-------------------------------------------------------!
