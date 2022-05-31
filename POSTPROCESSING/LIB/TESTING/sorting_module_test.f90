MODULE sorting

! module for sorting the species arrray based on code names
! routine found on http://fortranwiki.org/fortran/show/qsort_inline
  USE DICO, only : species, numsp, species_detailed_info
  USE CONC, only : concentrations, rdat
  USE RATE, only : idrestoi,idpdstoi
  USE USER_INPUT, only : flag_rates
  IMPLICIT NONE
  integer, parameter   :: QSORT_THRESHOLD = 16
  integer   :: array_size
  contains
  
  function find_species_index(key)  result(pos)
    character*(*), intent(in)  :: key
    integer :: range, start, finish, mid, pos
    
! binary search in ordered species%code(1:numsp) for code == key
    start = 1
    finish = numsp
    range = finish - start
    mid = (start+finish)/2
    
    do while (trim(species(mid)%code) /= trim(key) .and. range > 0 )
      if (trim(key) > trim(species(mid)%code)) then
        start = mid+1
      else
        finish = mid -1
      endif
      range = finish - start
      mid = (start+finish)/2
      if (mid == 0) exit
    enddo
    if (mid > 0) then
      if(trim(species(mid)%code) /= trim(key)) then
        pos = 0
      else
        pos = mid
      endif   
    else
      pos = 0
    endif
  end function 

  
  subroutine init()
    array_size = numsp
  end subroutine
  
  logical function less_than(a,b)
    integer, intent(in) :: a, b
    less_than = (species(a)%code < species(b)%code)

    
  end function
  
  subroutine swap(a,b)
    integer, intent(in) :: a, b
    type(species_detailed_info)  :: hold_spec
    real, dimension(rdat)        :: hold_conc
    
    hold_spec  = species(a)
    species(a) = species(b)
    species(b) = hold_spec
    
    hold_conc           = concentrations(:,a)
    concentrations(:,a) = concentrations(:,b)
    concentrations(:,b) = hold_conc

    if (flag_rates) then
      WHERE(idrestoi.eq.a)
          idrestoi = b
      ELSEWHERE(idrestoi.eq.b)
          idrestoi = a
      ENDWHERE

      WHERE(idpdstoi.eq.a)
          idpdstoi = b
      ELSEWHERE(idpdstoi.eq.b)
          idpdstoi = a
      ENDWHERE
    endif

    return
    
  end subroutine
  
  subroutine rshift(left, right)
    integer, intent(in) :: left, right
    type(species_detailed_info) :: hold_spec
    real, dimension(rdat)       :: hold_conc
    integer  ::  i
    hold_spec = species(right)
    hold_conc = concentrations(:, right)

    do i=right, left+1, -1
      species(i) = species(i-1)
      concentrations(:,i) = concentrations(:,i-1)

      if (flag_rates) then
        WHERE(idrestoi.eq.i)
            idrestoi = i-1
        ELSEWHERE(idrestoi.eq.i-1)
            idrestoi = i
        ENDWHERE

        WHERE(idpdstoi.eq.i)
            idpdstoi = i-1
        ELSEWHERE(idpdstoi.eq.i-1)
            idpdstoi = i
        ENDWHERE
      endif

    enddo
    species(left) = hold_spec
    concentrations(:,left) = hold_conc
    
  end subroutine rshift
    
  subroutine sort_species()
   
!======================================================================
! Fast in-line QSORT+INSERTION SORT for Fortran.
! Author: Joseph M. Krahn
! FILE: qsort_inline.inc
! PURPOSE:
! Generate a custom array sort procedure for a specific type,
! without the comparison-callback overhead of a generic sort procedure.
! This is essentially the same as an in-line optimization, which generally
! is not feasible for a library-based generic sort procedure.
!
! This implementation is as generic as possible, while avoiding the need
! for a code pre-processor. The success of this approach assumes that
! internal procedures are always in-lined with optimized Fortran compilation.
!
! USAGE:
! This file contains the sort subroutine body. You must supply
! an integer parameter QSORT_THRESHOLD, and internal procedures:
!    subroutine INIT()
!    logical function LESS_THAN(a,b)
!    subroutine SWAP(a,b)
!    subroutine RSHIFT(left,right)
!
! Descriptions:
!
! SUBROUTINE INIT()
!   Any user initialization code. This is needed because executable
!   statements cannot precede this code, which begins with declarations.
!   In many cases, this is just an empty procedure.
!
! LOGICAL FUNCTION LESS_THAN(a,b)
!   Return TRUE if array member 'a' is less than array member 'b'.
!   Only a TRUE value causes a change in sort order. This minimizes data
!   manipulation, and maintains the original order for values that are
!   equivalent by the sort comparison. It also avoids the need to
!   distinguish equality from greater-than.
!
! SUBROUTINE SWAP(A,B)
!   Swap array members 'a' and 'b'
!
! SUBROUTINE RSHIFT(LEFT,RIGHT)
!   Perform a circular shift of array members LEFT through RIGHT,
!   shifting the element at RIGHT back to the position at LEFT.
!
! QSORT_THRESHOLD:
!   The QSORT is used down to the QSORT_THRESHOLD size sorted blocks.
!   Then insertion sort is used for the remainder, because it is faster
!   for small sort ranges. The optimal size is not critical. Most of
!   the benefit is in blocks of 8 or less, and values of 16 to 128
!   are generally about equal speed. However, the optimal value
!   depends a lot on the hardware and the data being sorted, so this
!   is left as a tunable parameter for cases where ther is an
!   effect on performance.
!
!---------------------------------------------------------------------
! NOTES:
! The procedure uses a optimized combination of QSORT and INSERTION
! sorting. The algorithm is based on code used in GLIBC. 
! A stack is used in place of recursive calls. The stack size must
! be at least as big as the number of bits in the largest array index.
!
! Sorting vectors of a multidimensional allocatable array can be
! VERY slow. In this case, or with large derived types, it is better
! to sort a simple derived type of key/index pairs, then reorder
! tha actual data using the sorted indices.
!
!---------------------------------------------------------------------
  integer :: stack_top, right_size, left_size
  integer :: mid, left, right, low, high

! A stack of 32 can handle the entire extent of a 32-bit
! index, so this value is fixed. If you have 64-bit indexed
! arrays, which might contain more than 2^32 elements, this
! should be set to 64.
  integer, parameter :: QSORT_STACK_SIZE = 32
  
  type qsort_stack; integer :: low, high; end type
  type(qsort_stack) :: stack(QSORT_STACK_SIZE)

  
  call init()

  if (array_size > QSORT_THRESHOLD) then
    low = 1
    high = array_size
    stack_top = 0

    QSORT_LOOP: &
    do
      mid = (low + high)/2
      if (LESS_THAN (mid, low)) then
        call SWAP(mid,low)
      end if
      if (LESS_THAN (high, mid)) then
        call SWAP(high,mid)
        if (LESS_THAN (mid, low)) then
          call SWAP(mid,low)
        end if
      end if
      left  = low + 1
      right = high - 1

      COLLAPSE_WALLS: &
      do
        do while (LESS_THAN (left, mid))
          left=left+1
        end do
        do while (LESS_THAN (mid, right))
          right=right-1
        end do
        if (left < right) then
          call SWAP(left,right)
          if (mid == left) then
            mid = right
          else if (mid == right) then
            mid = left
          end if
          left=left+1
          right=right-1
        else
          if (left == right) then
            left=left+1
            right=right-1
          end if
          exit COLLAPSE_WALLS
        end if
      end do COLLAPSE_WALLS

! Set up indices for the next iteration.
! Determine left and right partition sizes.
! Defer partitions smaller than the QSORT_THRESHOLD.
! If both partitions are significant,
! push the larger one onto the stack.
      right_size = right - low
      left_size = high - left
      if (right_size <= QSORT_THRESHOLD) then
        if (left_size <= QSORT_THRESHOLD) then
          ! Ignore both small partitions: Pop a partition or exit.
          if (stack_top<1) exit QSORT_LOOP
          low=stack(stack_top)%low; high=stack(stack_top)%high
          stack_top=stack_top-1
        else
          ! Ignore small left partition.
          low = left
        end if
      else if (left_size <= QSORT_THRESHOLD) then
        ! Ignore small right partition.
        high = right
      else if (right_size > left_size) then
        ! Push larger left partition indices.
        stack_top=stack_top+1
        stack(stack_top)=qsort_stack(low,right)
        low = left
      else
        ! Push larger right partition indices.
        stack_top=stack_top+1
        stack(stack_top)=qsort_stack(left,high)
        high = right
      end if
    end do QSORT_LOOP
  end if ! (array_size > QSORT_THRESHOLD)

! Sort the remaining small partitions using insertion sort,
! which should be faster for partitions smaller than the
! appropriate QSORT_THRESHOLD.

! First, find smallest element in first QSORT_THRESHOLD and
! place it at the array's beginning. This places a lower
! bound 'guard' position, and speeds up the inner loop
! below, because it will not need a lower-bound test.
  low = 1
  high = array_size

! left is the MIN_LOC index here:
  left=low
  do right = low+1, min(low+QSORT_THRESHOLD,high)
    if (LESS_THAN(right,left)) left=right
  end do
  if (left/=low) call SWAP(left,low)

! Insertion sort, from left to right.
! (assuming that the left is the lowest numbered index)
  INSERTION_SORT: &
  do right = low+2,high
    left=right-1
    if (LESS_THAN(right,left)) then
      do while (LESS_THAN(right,left-1))
        left=left-1
      end do
      call RSHIFT(left,right)
    end if
  end do INSERTION_SORT
    
    
  end subroutine
!--------------------------------------------------------------  

END MODULE
