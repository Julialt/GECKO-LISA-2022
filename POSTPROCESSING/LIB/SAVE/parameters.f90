MODULE PARAMETERS

! code length of a given spcies
integer, parameter     ::  lco    = 6
!integer, parameter     ::  lco    = 8  ! MECHGEN
! max length of species name
integer, parameter     ::  maxlsp = lco+2
! max number of species
integer, parameter     ::  maxsp  = 9000000
! length of string in dictionary
integer, parameter     ::  ldi    = 146
! formula length of a given species
integer, parameter     ::  lfo    = 120
! functionalities list length
integer, parameter     ::  lfl    = 15
! max length of group
integer, parameter     ::  lgr    = 21
! max number of groups
integer, parameter     :: mca     = 29
! maximum number of data in the result file
integer, parameter     ::  mdat   = 20000
! maximum typical species per requested set
integer, parameter     :: mtpcal  = 1000
! maximum typical species per bin
integer, parameter     :: mbinspe = 20,   mspe = 10000
! maximum number of single species to plot
!integer, parameter     :: mxselspe = 100
integer, parameter     :: mxselspe = 650
! maximum number of precursors (primary species) in a simulation
integer, parameter     :: mxprecu = 100
! max number of carbon by species
integer                :: kmax
! max number of functions by species
integer                :: pmax
! maximum number of data in the pvap/cstar distribution (is min log10(pvap)+1)
integer                :: mpc
! for real size
integer, parameter     :: real_kind = kind(1.0)
! maximum number of nitrates to look at on a molecule
integer, parameter     :: max_nitrates = 10
! max length of dynamic filters
integer, parameter     :: mfilt = 300
! maximum rings allowed, and ring-joining characters
integer, parameter     :: mri = 4
character*(1), parameter :: digit(mri) = (/'1','2','3','4'/)
! maximum # of copies of formula allowed
integer, parameter     :: mco = 99


END MODULE
