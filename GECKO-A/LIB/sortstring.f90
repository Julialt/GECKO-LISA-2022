! the subroutine sort_string(string) return the array "string" provided
! as input as a sorted array (alphabetical order).
MODULE sortstring
IMPLICIT NONE
INTEGER, PARAMETER :: QSORT_THRESHOLD = 16

CONTAINS

LOGICAL FUNCTION less_than(a,b,string)
  INTEGER, INTENT(IN) :: a, b
  CHARACTER(LEN=*), INTENT(IN) :: string(:)
  less_than = (string(a) < string(b))
END FUNCTION
  
SUBROUTINE swap(a,b,string)
  INTEGER, INTENT(IN) :: a, b
  CHARACTER(LEN=*),INTENT(INOUT) :: string(:)
  CHARACTER(LEN=LEN(string(1)))  :: hold_string
  hold_string = string(a)
  string(a) = string(b)
  string(b) = hold_string
  RETURN
END SUBROUTINE

SUBROUTINE r_shift(left,right,string)
  INTEGER, INTENT(IN) :: left, right
  CHARACTER(LEN=*),INTENT(INOUT) :: string(:)
  INTEGER  ::  i
  CHARACTER(LEN=LEN(string(1)))  :: hold_string
  hold_string = string(right)
  DO i=right, left+1, -1
    string(i) = string(i-1)
  ENDDO
  string(left) = hold_string
END SUBROUTINE r_shift


SUBROUTINE sort_string(string)
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
  CHARACTER(LEN=*),INTENT(INOUT) :: string(:)
  INTEGER :: array_size
  INTEGER :: stack_top, right_size, left_size
  INTEGER :: mid, left, right, low, high

! A stack of 32 can handle the entire extent of a 32-bit
! index, so this value is fixed. If you have 64-bit indexed
! arrays, which might contain more thant 2^32 elements, this
! should be set to 64.
  integer, parameter :: QSORT_STACK_SIZE = 32
  
  type qsort_stack; integer :: low, high; end type
  type(qsort_stack) :: stack(QSORT_STACK_SIZE)

  array_size=SIZE(string,1)
!  call init()

  if (array_size > QSORT_THRESHOLD) then
    low = 1
    high = array_size
    stack_top = 0

    QSORT_LOOP: &
    do
      mid = (low + high)/2
!      if (LESS_THAN (mid, low)) then
      if (LESS_THAN (mid, low,string)) then
!        call SWAP(mid,low)
        call SWAP(mid,low,string)
      end if
!      if (LESS_THAN (high, mid)) then
      if (LESS_THAN (high, mid,string)) then
!        call SWAP(high,mid)
        call SWAP(high,mid,string)
!        if (LESS_THAN (mid, low)) then
        if (LESS_THAN (mid, low,string)) then
!          call SWAP(mid,low)
          call SWAP(mid,low,string)
        end if
      end if
      left  = low + 1
      right = high - 1

      COLLAPSE_WALLS: &
      do
!        do while (LESS_THAN (left, mid))
        do while (LESS_THAN (left, mid,string))
          left=left+1
        end do
!        do while (LESS_THAN (mid, right))
        do while (LESS_THAN (mid, right, string))
          right=right-1
        end do
        if (left < right) then
!          call SWAP(left,right)
          call SWAP(left,right,string)
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

! Sort the remaining small partitions using insertion sort, which should 
! be faster for partitions smaller than the appropriate QSORT_THRESHOLD.

! First, find smallest element in first QSORT_THRESHOLD and place it at
! the array's beginning. This places a lower bound 'guard' position, and
! speeds up the inner loop below, because it will not need a lower-bound
! test.
  low = 1
  high = array_size

! left is the MIN_LOC index here:
  left=low
  do right = low+1, min(low+QSORT_THRESHOLD,high)
!    if (LESS_THAN(right,left)) left=right
    if (LESS_THAN(right,left,string)) left=right
  end do
!  if (left/=low) call SWAP(left,low)
  if (left/=low) call SWAP(left,low,string)

! Insertion sort, from left to right.
! (assuming that the left is the lowest numbered index)
  insertion_sort: &
  do right = low+2,high
    left=right-1
!    if (LESS_THAN(right,left)) then
    if (LESS_THAN(right,left,string)) then
!      do while (LESS_THAN(right,left-1))
      do while (LESS_THAN(right,left-1,string))
        left=left-1
      end do
!      call RSHIFT(left,right)
      call R_SHIFT(left,right,string)
    end if
  end do insertion_sort
    
    
END SUBROUTINE

END MODULE sortstring
