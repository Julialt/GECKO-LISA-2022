! The aim of this module is to contain all things related
! to constraining and forcing parameters
! such as concentrations, emissions, environmental values
!CMV! 20180822: start with inorganic aerosol constraints pH, sulf, nitrates, kappa
!           before generalizing and centralizing other things

MODULE constraints_module
USE akparameter_module
implicit none
  type constrained_thing
    real :: val ! current_value, supposed to be updated for the current timestep
    integer  :: npoints
    real, dimension(maxinput) :: time, values
  end type constrained_thing
  
  contains
  
  subroutine update_constraint(const_thing, cur_time)
    type(constrained_thing), intent(inout) :: const_thing
    real, intent(in) :: cur_time
    
    const_thing%val = interpolate_value(const_thing, cur_time)
  end subroutine
  
  function interpolate_value(const_thing, cur_time) result(cur_val)
    type(constrained_thing), intent(in) :: const_thing
    real, intent(in) :: cur_time
    real  :: cur_val
    
    !real :: t1, t2, v1, v2, m, p
    logical :: found
    integer :: j
    ! linearly interpolate between the two closest points in time
    found = .false.
    DO j=1, const_thing%npoints-1
      IF (const_thing%time(j).le.cur_time .and. const_thing%time(j+1).gt.cur_time) THEN
        found = .true.
        cur_val = linear_interp(cur_time, &
                                const_thing%time(j), const_thing%time(j+1), &
                                const_thing%values(j), const_thing%values(j+1))
        exit
      ENDIF
    ENDDO
    if (.not. found) then
      write(6,*) 'could not interpolate value in module constraint_module'
      stop
    endif
        
  end function
  
  function linear_interp(t, t1, t2, v1, v2) result(v)
    real, intent(in) :: t, t1, t2, v1, v2
    real :: v
    real :: m, p
    
    m  = (v2-v1)/(t2-t1)
    p  = v2-m*t2
    v = m*t + p
  end function
!----------------------------------------------------
! function to estimate SIMPOL Cstar from Nannoolal values
! currently based on standard products of 6-to-12-alkanols
! Update as more info becomes available
! JMLT, march 4, 2019
!----------------------------------------------------
  function nan_to_sim(Cstar)  result(tCstar)
    REAL,INTENT(IN) :: Cstar
    REAL :: tCstar

    tCstar = MIN(Cstar,(Cstar+209.4)/1.81)

  end function

END MODULE constraints_module
