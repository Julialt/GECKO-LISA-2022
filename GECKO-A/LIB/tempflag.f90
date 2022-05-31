MODULE tempflag
  IMPLICIT NONE 

! ignore reaction products (e.g. after a given # of generation)
  INTEGER,SAVE :: iflost  ! flag raised if products are not written 
  INTEGER,SAVE :: xxc     ! # of C atom in the reactant
  INTEGER,SAVE :: xxh     ! # of H atom in the reactant
  INTEGER,SAVE :: xxn     ! # of N atom in the reactant
  INTEGER,SAVE :: xxo     ! # of O atom in the reactant
  INTEGER,SAVE :: xxr     ! # of O atom in the reactant
  INTEGER,SAVE :: xxs     ! # of O atom in the reactant
  INTEGER,SAVE :: xxfl    ! # of O atom in the reactant
  INTEGER,SAVE :: xxbr    ! # of O atom in the reactant
  INTEGER,SAVE :: xxcl    ! # of O atom in the reactant

END MODULE tempflag
