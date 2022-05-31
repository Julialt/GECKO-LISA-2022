!akparameter.h file for mechanism propane_vd
! max length of species names
      INTEGER,PARAMETER :: maxlsp = 8
! max length of a printed reaction
      INTEGER,PARAMETER :: maxreac_char =  90
! max number of species
      INTEGER,PARAMETER :: maxsp = 1191
! max number of reactants in a reaction
      INTEGER,PARAMETER :: mxleft = 2
! max number of products in a reaction
      INTEGER,PARAMETER :: mxright = 4
! max number of reactions
      INTEGER,PARAMETER :: maxre = 3175
! max number of reactions with "M"
      INTEGER,PARAMETER :: max_m = 5
! max number of fall-off reactions 
      INTEGER,PARAMETER :: maxfo = 23
! max number of reactions with "HV"
      INTEGER,PARAMETER :: maxhv = 326
! max number of reactions with "CVAR"
      INTEGER,PARAMETER :: maxcvar = 1
! max number of reactions with "EXTRA"
      INTEGER,PARAMETER :: maxextra = 6
! max number of reactions with "OXYGEN"
      INTEGER,PARAMETER :: maxo2 = 36
! max number of ISOMERIZATIONS
      INTEGER,PARAMETER :: maxiso = 1
! max number of species undergoing phase equilibrium
      INTEGER,PARAMETER :: maxt = 190
! max number of different types of auxiliary information
      INTEGER,PARAMETER :: maxaux = 12
! max number of different classes of RO2
      INTEGER,PARAMETER :: maxro2 = 9
! max number of PEROx in a class
      INTEGER,PARAMETER :: mxro2cl = 42
! max number of reactions with PEROx
      INTEGER,PARAMETER :: mxrpero = 105
! max number of different classes of dimer
      INTEGER,PARAMETER :: maxdimer = 4
! max number of dimers in a given class
      INTEGER,PARAMETER :: mxrdimer = 1
! max number of variable coefficients in CVAR type reaction
      INTEGER,PARAMETER :: maxcoe = 30
! max number of data set in function of temp. in CVAR type reaction
      INTEGER,PARAMETER :: nset = 4
! max number of angles in "HV" function
      INTEGER,PARAMETER :: maxang = 12
! -------------------------- DATA FOR CHROMOPHORE ---------------
! max # of different types of chromophore
      INTEGER,PARAMETER :: mchromo = 60
! number of chromophores  to be stored in "most used"  chromophore
!     (the "top" tables)
      INTEGER,PARAMETER :: mtopchromo = 10
! max # of species that can be stored in "most used"  chromophore
!     (the "top" tables)
      INTEGER,PARAMETER :: msptopchromo = 326
! number of chromophores  to be stored in "regularly used"  chromophore
!     (the "med" tables)
      INTEGER,PARAMETER :: mmedchromo = 500
! max # of species that can be stored in "regularly used"  chromophore
!     (the "med" tables)
      INTEGER,PARAMETER :: mspmedchromo = 326
! -------------------------- 
! max coefficient for interpolation of the photolytic frequencies
      INTEGER,PARAMETER :: nlo = maxang*3
! max number of boxes in the model
      INTEGER,PARAMETER :: mbox = 2
! max number of data to compute mixing height
      INTEGER,PARAMETER :: mhd = 60
! max number of isurface types
      INTEGER,PARAMETER :: msur = 4
! max number of counting species for which stoe. coff. need to be
!     evaluated in CVAR application
      INTEGER,PARAMETER :: mopc = 4
! max number of positions used to evaluate stoe. coff. for counting
!     species from the operator species
      INTEGER,PARAMETER :: mpos = 6
! max number of self reactions
      INTEGER,PARAMETER :: mself = 21
! -------------------------- 
! max number of species for which Psat is evaluated
      INTEGER,PARAMETER :: mxsat = 192
! max number of species for which Vdep is evaluated
      INTEGER,PARAMETER :: mxdep = 199
! -------------------------- 
! max number of emitted species
      INTEGER,PARAMETER :: mes = 55
      INTEGER,PARAMETER :: maxem = 55
! max number of emitted times
      INTEGER,PARAMETER :: mtim = 25
! max number of spp assessing prod & loss rates
      INTEGER,PARAMETER :: mtr = 100
! max number of species that can be constrained
      INTEGER,PARAMETER :: maxconst = 20
! max number of lines for datapoints for constrained species
      INTEGER,PARAMETER :: maxinput = 1500
! -------------------------- 
! structure for external data used to constrain
! species emissions, concentrations, deposition...
! with possibility for using an input file
      TYPE species_data 
        LOGICAL              :: activefg 
        CHARACTER(LEN=maxlsp):: name 
        CHARACTER(LEN=10)    :: unit 
        INTEGER              :: index,npoints 
        REAL                 :: table(maxinput,2) 
      END TYPE 
!structure for storing surface data (emissions for now
! TODO later: add deposition data...
      TYPE surface_data 
        INTEGER              :: nemis 
        TYPE(species_data)   :: emission(maxem) 
      END TYPE 
