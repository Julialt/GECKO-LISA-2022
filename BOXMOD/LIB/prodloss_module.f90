!structure for storing index map from species to reactions
MODULE prodloss_module
      TYPE spec_reac_map
        INTEGER              :: nloss, nprod
        INTEGER, ALLOCATABLE :: idl(:) !id of loss reactions
        INTEGER, ALLOCATABLE :: idp(:) !id of prod reactions
        REAL, ALLOCATABLE :: stpd(:) ! prod associated stoechiometric coefficient
      END TYPE

END MODULE prodloss_module

