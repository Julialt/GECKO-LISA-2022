***************************************************************
*  HEADER INFORMATION: general.h
*  BROADLY USED PARAMETERS :
***************************************************************


* dimension for table
* ---------------
      INTEGER mca, mri, mnr, mnp, mcp, mco, mni, mlv, mra, mfn, mkr 
      INTEGER mbg, mrd, mps, mog

* maximum carbons allowed
      PARAMETER (mca=22)
      !PARAMETER (mca=15)  ! < for Leeds MCM
      
* maximum rings allowed, and ring-joining characters
      PARAMETER (mri=4)
      CHARACTER*(1)   digit(mri)
      DATA digit /'1','2','3','4'/
      
* maximum number of reactions for a given species
      PARAMETER (mnr=20)
      
* maximum number of products allowed in a given reaction
      PARAMETER (mnp=18)

* maximum number of coproducts allowed for a given species
      PARAMETER (mcp=5)
      
* maximum # of copies of formula allowed
      PARAMETER (mco=99)
      
* maximum # of species allowed in the mechanism
      PARAMETER (mni=10000000)
      
* maximum number of species in the voc stack
      PARAMETER (mlv=1500000)

* maximum number of species in the radical stack
      PARAMETER (mra=200)

* maximum number of inorganic species and in fixedname.dat and in 
* special_dict.dat
      PARAMETER (mfn=1000)
      
* maximum number of known reactions in a given series (e.g.
* NO3+VOC, OH+VOC or photolytic reaction
      PARAMETER (mkr=110)

* maximum number of benson groups
      parameter (mbg=310)     

* maximum number of rate data as input for a given series
      parameter (mrd=500)     

* maximum number of "primary" species that can be given as input 
      parameter (mps=300)     
      
* maximum number of reactions linked to special species (i.e. not set
* by the generator - "og"=out generator)
      parameter (mog=10000)     
      
* dimension for string
* --------------------
      INTEGER lgr, lco, lfo, lcf, ldi, lfl, llin, lst

* maximum length of group string (also used for benson group)
      PARAMETER (lgr=21)
* code length of a given species
      !PARAMETER (lco=14) ! Jenkins
      PARAMETER (lco=6) ! SAPRC
* formula length of a given species
      PARAMETER (lfo=100)
* length of code+formula
      PARAMETER (lcf=106)
* length of string in the dictionnary
      PARAMETER (ldi=126)
*      PARAMETER (ldi=100)
* length of a string in the stack (code+formula+i3+i3)
      PARAMETER (lst=112)
* length of functionnalities list
      PARAMETER (lfl=15)
* length of a typical input line (e.g. filename, ...)
      PARAMETER (llin=100) 
     
***************************************************************
*  HEADER INFORMATION: general.h
*  BROADLY USED VARIABLES :
***************************************************************
* boolean to write or not information during the generation of the
* scheme
      INTEGER         wtflag, wtopeflag
* boolean to consider only reaction with NOx for peroxys and acylperoxys
      INTEGER         high_NOxfg, zero_NOXfg
* boolean to consider only 3 types of peroxys for recombinaison
* reactions
      INTEGER         per3_typefg
* number of species in : dict,chmlst,namlst,dbrch
c      INTEGER         nrec(mca),nrectot
      INTEGER         nrec,nrectot

* table of inorganic species potentially produced
      INTEGER         ninorg
      CHARACTER*(ldi) inorglst(mfn)

* special species
      INTEGER         nspsp
      LOGICAL         lospsp(mfn)
      CHARACTER*(ldi) dictsp(mfn)
      

* table of species having a special name
      COMMON          nspsp,lospsp,dictsp

      COMMON          nrectot,nrec
      COMMON          inorglst,ninorg
* flag to allow the writting in various files
      COMMON          wtflag
      COMMON          wtopeflag
* flag to consider only reaction with NOx for peroxys and acylperoxys
      COMMON          high_NOxfg, zero_NOxfg
* flag to consider only 3 types of peroxys 
      COMMON         per3_typefg
      
**********************************      
      
            
* ------------ define organic functionalities -----------------
*
      CHARACTER*(1)  fluorine, sulfur, aromatic
      CHARACTER*(2)  bromine, chlorine, amine, carbonyl, secondary
      CHARACTER*(3)  ketene, aldehyde, nitroso(2), ketone, acyl,
     &             primary, methyl, ether
      CHARACTER*(4)  alkoxy, hydroxy
      CHARACTER*(5)  alkyl_peroxy, hydro_peroxide, nitro, nitrite
      CHARACTER*(6)  carboxylic_acid, nitrate, criegee, t_peroxy, 
     &             acyl_oxy, s_alcohol
      CHARACTER*(7)  acyl_peroxy, peroxy_acid, s_peroxy, p_alcohol,
     &             hot_criegee, peroxy_nitrate
      CHARACTER*(8)  ext_criegee, int_criegee, p_peroxy
      CHARACTER*(9)  pan, cold_criegee
C
      CHARACTER(61)    alfa
      CHARACTER(70)    pri*70
C
      INTEGER      nalfa
C
C
      DATA fluorine        /'F'/
      DATA sulfur          /'S'/
      DATA bromine         /'Br'/
      DATA chlorine        /'Cl'/
      DATA amine           /'NH'/
      DATA carbonyl        /'CO'/
      DATA secondary       /'CH'/
      DATA ketene          /'CdO'/
      DATA aldehyde        /'CHO'/
      DATA ketone          /'CO('/
      DATA acyl            /'CO.'/
      DATA primary         /'CH2'/
      DATA methyl          /'CH3'/
      DATA nitroso(1)      /'NO '/,    nitroso(2)     /'NO)'/
      DATA alkoxy          /'(O.)'/
      DATA hydroxy         /'(OH)'/
      DATA alkyl_peroxy    /'(OO.)'/
      DATA hydro_peroxide  /'(OOH)'/
      DATA nitro           /'(NO2)'/
      DATA nitrite         /'(ONO)'/
      DATA carboxylic_acid /'CO(OH)'/
      DATA nitrate         /'(ONO2)'/
      DATA criegee         /'.(OO.)'/
      DATA t_peroxy        /'C(OO.)'/ 
      DATA acyl_oxy        /'CO(O.)'/
      DATA s_alcohol       /'CH(OH)'/
      DATA acyl_peroxy     /'CO(OO.)'/
      DATA peroxy_acid     /'CO(OOH)'/
      DATA peroxy_nitrate  /'(OONO2)'/
      DATA hot_criegee     /'.(OO.)*'/
      DATA s_peroxy        /'CH(OO.)'/
      DATA p_alcohol       /'CH2(OH)'/
      DATA ext_criegee     /'CH.(OO.)'/
      DATA int_criegee     /'C.(OO.) '/
      DATA p_peroxy        /'CH2(OO.)'/
      DATA pan             /'CO(OONO2)'/
      DATA cold_criegee    /'CH2.(OO.)'/
      DATA ether           /'-O-'/
      DATA aromatic        /'c'/

***********************************************************************

      DATA pri 
     &  /'CHOCHC CH(CH2 CH3O.)OO. OOH ONOF  Br2Br)Cl2Cl)NO2NO)OH)'/
      DATA alfa 
     & /'0123456789ABCDEFGHIJKLMNOPQRSTUWXYZabcdefghijklmnopqrstuvwxyz'/
      DATA nalfa /61/
c      DATA alfa 
c     & /'0123456789ABCDEFGHIJKLMNOPQRSTUWXYZ'/
c      DATA nalfa /35/

********************************************************************
