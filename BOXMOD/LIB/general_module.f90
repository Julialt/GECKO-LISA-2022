!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  HEADER INFORMATION: general_module
!  BROADLY USED PARAMETERS :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE general_module
      IMPLICIT NONE

! dimension for table
! ---------------
      INTEGER mca, mri, mnr, mnp, mcp, mco, mni, mlv, mra, mfn, mkr 
      INTEGER mnic, mbg, mrd, mps, mog, mhy, mnitro, mox

! maximum carbons allowed
      !PARAMETER (mca=13)
      PARAMETER (mca=25)  ! < for Leeds MCM
      
! maximum rings allowed, and ring-joining characters
      PARAMETER (mri=4)
      CHARACTER(1)   digit(mri)
      DATA digit /'1','2','3','4'/
      
! maximum number of reactions for a given species
      PARAMETER (mnr=20)
      
! maximum number of products allowed in a given reaction
      PARAMETER (mnp=18)

! maximum number of coproducts allowed for a given species
      PARAMETER (mcp=5)
      
! maximum # of copies of formula allowed
      PARAMETER (mco=99)
      
! maximum # of species allowed in the mechanism for a given number of C
!      PARAMETER (mnic=447392)

! maximum # of species allowed in the mechanism
!      PARAMETER (mni=1000000)
      PARAMETER (mni=10000000)
      
! maximum number of species in the voc stack
      PARAMETER (mlv=1500000)

! maximum number of species in the radical stack
      PARAMETER (mra=200)

! maximum number of inorganic species and in fixedname.dat and in 
! special_dict.dat
      PARAMETER (mfn=1000)
      
! maximum number of known reactions in a given series (e.g.
! NO3+VOC, OH+VOC or photolytic reaction
      PARAMETER (mkr=110)

! maximum number of benson groups
      parameter (mbg=310)     

! maximum number of rate data as input for a given series
      parameter (mrd=500)     

! maximum number of "primary" species that can be given as input 
      parameter (mps=300)     
      
! maximum number of reactions linked to special species (i.e. not set
! by the generator - "og"=out generator)
      parameter (mog=10000)     
      
! dimension for string
! --------------------
      INTEGER lgr, lco, lfo, lcf, ldi, lfl, llin, lst

! maximum length of group string (also used for benson group)
      PARAMETER (lgr=21)
! code length of a given species
      !PARAMETER (lco=14) ! Jenkins
      PARAMETER (lco=6) ! SAPRC
      !PARAMETER (lco=8) ! MECHGEN
! formula length of a given species
      PARAMETER (lfo=120)
! length of code+formula
      PARAMETER (lcf=126)
! length of string in the dictionnary
      PARAMETER (ldi=146)
!      PARAMETER (ldi=100)
! length of a string in the stack (code+formula+i3+i3)
      PARAMETER (lst=112)
! length of functionnalities list
      PARAMETER (lfl=15)
! length of a typical input line (e.g. filename, ...)
      PARAMETER (llin=100) 
     
      END MODULE general_module
