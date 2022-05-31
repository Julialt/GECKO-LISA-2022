MODULE BINS
  IMPLICIT NONE
  private
  
  type, public :: bins_class
    private
    real(kind=8), public :: min_val, max_val, resolution
    real(kind=8), dimension(:), allocatable, public :: lowerlim, upperlim
    integer, public :: nbins
    contains
      procedure, public :: initialize => initialize
      procedure, public :: find_index => find_bin
      procedure, public :: destroy => destroy
      procedure :: calc_nbins => calc_nbins
      procedure :: fill_bin_boundaries => fill_bin_boundaries
  end type bins_class
  
  CONTAINS
  
  subroutine initialize(this, min_val, max_val, resolution)
    class(bins_class), intent(inout) :: this
    real(kind=8), intent(in) :: min_val, max_val, resolution
    this%min_val = min_val
    this%max_val = max_val
    this%resolution = resolution
    call this%calc_nbins()    
    call this%fill_bin_boundaries()
    
  end subroutine
  
  subroutine destroy(this)
    class(bins_class), intent(inout) :: this
    
    if (allocated(this%lowerlim)) deallocate(this%lowerlim)
    if (allocated(this%upperlim)) deallocate(this%upperlim)
  end subroutine
  
  function find_bin(this, value) result(ibin)
    class(bins_class), intent(in) :: this
    real(kind=8), intent(in) :: value
    integer :: ibin
    integer :: i
    ibin = -1
    do i=this%nbins,1,-1
      if (this%lowerlim(i) <= value) then
        ibin = i
        exit
      endif        
    enddo
    
    if (ibin < 1 .or. ibin > this%nbins) then
      write(6,*) "error in finding bin index for value:", value
      write(6,*) "ibin=",ibin
      stop
    endif
    return    
    
  end function find_bin
  
  subroutine calc_nbins(this)
    class(bins_class), intent(inout) :: this
    
    this%nbins = ceiling((this%max_val - this%min_val)/this%resolution)
    
  end subroutine calc_nbins
  
  subroutine fill_bin_boundaries(this)
    class(bins_class), intent(inout) :: this
    integer :: ibin
    
    if (allocated(this%lowerlim)) then
      deallocate(this%lowerlim)
    endif
    if (allocated(this%upperlim)) then
      deallocate(this%upperlim)
    endif
    allocate(this%lowerlim(this%nbins), this%upperlim(this%nbins))
    
    do ibin = 1, this%nbins
      this%lowerlim(ibin) = this%min_val + (ibin-1)*this%resolution
      this%upperlim(ibin) = this%min_val + ibin*this%resolution
    enddo
  
  end subroutine

END MODULE BINS
