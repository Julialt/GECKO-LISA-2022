!**************************************************************
!  HEADER INFORMATION: general.h
!  BROADLY USED PARAMETERS :
!**************************************************************


! dimension for table
! ---------------
      INTEGER mca, mri, mnr, mnp, mcp, mco, mni, mlv, mra, mfn, mkr, mnt
      INTEGER mnic, mbg, mrd, mps, mog, maq
      INTEGER mhyd, mhiso

! maximum carbons allowed
      !PARAMETER (mca=13)
      PARAMETER (mca=29)  ! < for Leeds MCM

! maximum rings allowed, and ring-joining characters
      PARAMETER (mri=4)
      CHARACTER(LEN=1)   digit(mri)
      DATA digit /'1','2','3','4'/

! maximum number of reactions for a given species
      PARAMETER (mnr=30)

! maximum number of products allowed in a given reaction
      PARAMETER (mnp=20)

! maximum number of coproducts allowed for a given species
      PARAMETER (mcp=5)

! maximum # of copies of formula allowed
      PARAMETER (mco=99)

! maximum # of species allowed in the mechanism for a given number of C
!      PARAMETER (mnic=447392)

! maximum # of species allowed in the mechanism
!      PARAMETER (mni=1000000)
      PARAMETER (mni=10000000)

! maximum number of aqueous species allowed in the mechanism
      PARAMETER (maq=1000000)

! maximum number of species in the voc stack
      PARAMETER (mlv=1500000)

! maximum number of species in the radical stack
      PARAMETER (mra=200)

! maximum number of inorganic species and in fixedname.dat and in
! special_dict.dat
      PARAMETER (mfn=3000)

! maximum number of known reactions in a given series (e.g.
! NO3+VOC, OH+VOC or photolytic reaction
      PARAMETER (mkr=110)

! maximum number of benson groups
      parameter (mbg=500)

! maximum number of rate data as input for a given series
      parameter (mrd=500)

! maximum number of "primary" species that can be given as input
      parameter (mps=300)

! maximum number of reactions linked to special species (i.e. not set
! by the generator - "og"=out generator)
      parameter (mog=10000)

! maximum number of hydrate a molecule can have (number of carbonyl)

      parameter (mhyd=10)
! maximum number of hydrate isomer a molecule can have
! for example, for 70 distinct can be made with 4 addition of water
! from 7 possible slots

      parameter (mhiso=700)
! maximum number of species a total species can contain
      parameter (mnt=100)


! dimension for string
! --------------------
      INTEGER lgr, lgb, lco, lfo, lcf, ldi, lfl, llin, lst

! maximum length of group string (previously also used for benson group)
      PARAMETER (lgr=15)
! maximum length of Benson group string
      PARAMETER (lgb=24)
! code length of a given species
      !PARAMETER (lco=14) ! Jenkins
      PARAMETER (lco=6) ! SAPRC
      !PARAMETER (lco=8) ! MECHGEN
! formula length of a given species
!      PARAMETER (lfo=100)
      PARAMETER (lfo=120)
! length of code+formula
!      PARAMETER (lcf=106)
      PARAMETER (lcf=126)
! length of string in the dictionary
!      PARAMETER (ldi=149) ! NCAR version: includes generation #
!      PARAMETER (ldi=178) ! PARIS version: no generation #, does include
!                         ! mass, radflg, #C,#H,#N,#O,4 more integer molecule fields
      PARAMETER (ldi=182) ! NCAR-PARIS version: includes generation #, 
!                         ! mass, radflg, #C,#H,#N,#O,4 more integer molecule fields
! length of a string in the stack (code+formula+i3+i3)
!      PARAMETER (lst=112)
      PARAMETER (lst=132)
! length of functionalities list
      PARAMETER (lfl=15)
! length of a typical input line (e.g. filename, ...)
      PARAMETER (llin=100)
! length of typical subroutine name
      INTEGER,PARAMETER :: lsb=20
! length of typical error message
      INTEGER,PARAMETER :: ler=60

! flag to write debug output from ncutil.f
      LOGICAL,PARAMETER :: debug_nc_fg = .TRUE.
      !LOGICAL,PARAMETER :: debug_nc_fg = .FALSE.
