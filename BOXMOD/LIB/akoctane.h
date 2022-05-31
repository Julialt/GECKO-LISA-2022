* max length of species names
      integer maxlsp
      parameter(maxlsp=8)

* max number of species
      integer maxsp
      parameter(maxsp=1300000)

* maximum number of reactants in a reaction
      integer mxleft
      parameter (mxleft=2)

* maximum number of products in a reaction
      integer mxright
      parameter (mxright=6)

* max number of reactions
      integer maxre
      parameter(maxre=7900000)

* max number of reactions with "M"
      integer max_m
      parameter(max_m=10)

* max number of fall-off reactions
      integer maxfo
      parameter(maxfo=30)

* max number of reactions with "HV"
      integer maxhv
      parameter(maxhv=1723000)

* max number of reactions with "CVAR"
      integer maxcvar
      parameter(maxcvar=1)

* max number of reactions with "EXTRA"
      integer maxextra
      parameter(maxextra=10)

* max number of reactions with "OXYGEN"
      integer maxo2
      parameter(maxo2=73000)

* max number of different types of auxiliary information
      integer maxaux
      parameter(maxaux=7)

* max number of different class of RO2
      integer maxro2
      parameter(maxro2=9)

* max number of reaction using the keyword PEROx
      integer mxrpero
      parameter(mxrpero=350000)

* max number of variable coefficient in CVAR type reaction
      integer maxcoe
      parameter(maxcoe=30)

* max of data set in fonction of temp. in CVAR type reaction
      integer nset
      parameter(nset=4)

* max of angle in "HV" function
      integer maxang
      parameter(maxang=20)

* max of different type of chomophore
      integer mchromo
      parameter(mchromo=140)

* max number of species that can be stored in a given chromophore
      integer mspchromo
      parameter(mspchromo=110000)

* max coefficient for interpolation of the photolytic frequencies
      integer nlo
      parameter(nlo=maxang*3)

* max number of boxes in the model
      integer mbox
      parameter(mbox=2)

* max number of data to compute mixing height
      integer mhd
      parameter(mhd=10)

* max number of surface type
      integer msur
      parameter(msur=4)

* max number of counting species for which stoe. coef. need to be 
* evaluated in CVAR application
      integer mopc
      parameter(mopc=4)

* max number of position used to evaluate stoe. coef. for counting 
* species from the operator species 
      integer mpos
      parameter(mpos=6)

* max number of emmitted species
      integer mes
      parameter(mes=500)

* max number of self reaction
      integer mself
      parameter(mself=20)
