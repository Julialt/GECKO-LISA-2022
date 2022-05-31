module dynamic_filter
  USE CONC, ONLY: numsp
  USE DICO, ONLY: species_detailed_info, species, func_codes, function_types
  USE PARAMETERS, ONLY: mfilt
  
  
  IMPLICIT NONE
  
  contains
  
  subroutine calc_anyphase_filters(phase, dyn_filter, species_mask)
    integer, intent(in) :: phase
    character*(*)   :: dyn_filter
  
    logical, dimension(:), intent(out)  :: species_mask 
    
    ! first filter phase
    
    species_mask = (species%phase == phase) .and. create_spec_mask(dyn_filter)
        
  end subroutine

  function create_spec_mask(text_line) result(spec_mask)
    character*(*), intent(in)  :: text_line
    logical, dimension(:), allocatable  :: spec_mask
    character*(mfilt)  :: atom_desc
    character*(mfilt)  :: min_func_desc
    character*(mfilt)  :: max_func_desc
    integer            :: i,j
    
    i = index(text_line, ";")
    j = index(text_line, "!")
    
    if (i == 0 .or. j == 0) then
      print *,'error, expected ; and ! in spec filter'
      print *, text_line
      stop 
    endif
    
    atom_desc = text_line(1:i-1)
    min_func_desc = text_line(i+1:j-1)
    max_func_desc = text_line(j+1:len_trim(text_line))
  
    allocate(spec_mask(numsp))
    
    spec_mask = .true.
    if (len_trim(atom_desc) .gt. 0) then
      spec_mask = spec_mask .and. atoms_mask(atom_desc)
    endif
    if (len_trim(min_func_desc) .gt. 0) then
      spec_mask = spec_mask .and. min_func_mask(min_func_desc)
    endif
    if (len_trim(max_func_desc) .gt. 0) then
      spec_mask = spec_mask .and. max_func_mask(max_func_desc)
    endif
        
  end function create_spec_mask
  
  function get_natoms(atom,atom_desc) result(n)
    character*(1), intent(in) :: atom
    character*(*), intent(in) :: atom_desc
    
    integer :: n
    character*(4) :: fmt
    integer       :: i, j, ind, nlen
    
    nlen = 0
    ! -1 = no atom constraint
    n = -1
    i = index(atom_desc, atom)
    if (i > 0) then
      ind = i
      ! try to read numbers until not possible
      do
        ind = ind + 1
!        if (ind > len(atom_desc) + 1) then
!          write(6,*) atom_desc
!          write(6,*) 'syntax error in dynamic filter'
!          write(6,*) ' cannot count number of atom in constraint for ', atom
!          stop
!        endif
        j = 0
        if (ind <= len(atom_desc)) then
          j = index("0123456789", atom_desc(ind:ind))
        endif
        if (j .le. 0 .or. ind .ge. len(atom_desc)) exit
      enddo  
      
      nlen = ind - i -1
      if (nlen <= 0) then
        write(6,*) atom_desc
        write(6,*) 'syntax error in dynamic filter'
        write(6,*) ' constraint on number of atoms does not work for ', atom
        stop
      endif
      write(fmt,'(A2,I1,A1)') '(i',nlen,')'
      read(atom_desc(i+1:ind-1), fmt) n
    endif
    
  
  end function 
  
  function atoms_mask(atom_desc) result(spec_mask)
    character*(*), intent(in)           :: atom_desc
    logical, dimension(:), allocatable  :: spec_mask
    integer :: nc, no, nh, nn !, i, j, ind
    
    ! -1 = no constraint
    nc = get_natoms('C', atom_desc)
    no = get_natoms('O', atom_desc)
    nh = get_natoms('H', atom_desc)
    nn = get_natoms('N', atom_desc)
       
    
    allocate(spec_mask(numsp))
    spec_mask = .true.
    
    if (nc > -1) then
      spec_mask = spec_mask .and. species%nc == nc
    endif
    if (no > -1) then
      spec_mask = spec_mask .and. species%no == no
    endif
    if (nh > -1) then
      spec_mask = spec_mask .and. species%nh == nh
    endif
    if (nn > -1) then
      spec_mask = spec_mask .and. species%nn == nn
    endif    
  
  end function atoms_mask
  
  function min_func_mask(func_desc) result(spec_mask)
    character*(*), intent(in)           :: func_desc
    logical, dimension(numsp)  :: spec_mask
    character*(1)         :: f
    integer               :: ifunc, min_nf, nf, ispec
  
    spec_mask = .true.
    
    do ifunc = 1, function_types
      f = func_codes(ifunc)
      min_nf = count(transfer(func_desc, 'a', len(func_desc)) == f)
      if (min_nf <= 0) cycle
      do ispec = 1, numsp
        if (len(species(ispec)%functions) <= 0)  cycle
        nf = count(transfer(species(ispec)%functions, 'a', len(species(ispec)%functions)) == f)
        if (nf < min_nf) then
          spec_mask(ispec) = .false.
        endif
      enddo            
    enddo
      
  end function min_func_mask
  
  function max_func_mask(func_desc) result(spec_mask)
    character*(*), intent(in)           :: func_desc
    logical, dimension(numsp)  :: spec_mask
    character*(1)         :: f
    integer               :: ifunc, max_nf, nf, ispec
  
    spec_mask = .true.
    
    do ifunc = 1, function_types
      f = func_codes(ifunc)
      max_nf = count(transfer(func_desc, 'a', len(func_desc)) == f)
      if (max_nf <= 0) cycle
      do ispec = 1, numsp
        if (len(species(ispec)%functions) <= 0)  cycle
        nf = count(transfer(species(ispec)%functions, 'a', len(species(ispec)%functions)) == f)
        if (nf >= max_nf) then
          spec_mask(ispec) = .false.
        endif
      enddo            
    enddo
  
  end function max_func_mask  

end module dynamic_filter
