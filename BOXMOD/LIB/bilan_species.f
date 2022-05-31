****************************************************************************
*                                                                          *
*         compute loss and production rates for a given species            *
*                                                                          *
****************************************************************************
      SUBROUTINE bilan_species(qfor,numre,numstoi,idrestoi,idpdstoi,
     &                         restoicf,pdstoicf,idXX,cbot,chrsp,
     &                         sourceXX,puitsXX,time,nbsource,nbpuits,
     &                         kicovi,sum_source,sum_puits)
      USE akparameter_module
      IMPLICIT NONE
      INCLUDE 'general.h'
   
* input:
      INTEGER  numstoi(maxre,2)
      INTEGER  idrestoi(maxre,mxleft)
      INTEGER  idpdstoi(maxre,mxright)
      REAL     restoicf(maxre,mxleft),pdstoicf(maxre,mxright)
      INTEGER  numre,idXX
      CHARACTER(maxlsp) chrsp(maxsp)
      REAL     cbot(maxsp)
      REAL     time
      REAL     qfor(maxre)

* input/output
      REAL     sourceXX(45),puitsXX(45)
      REAL     kicovi
      INTEGER  nbsource,nbpuits

* internal
      INTEGER  i,j,k
      REAL     sum_source,sum_puits,rate(maxre)
      REAL     sum_source_XX,sum_puits_XX
      
* initialize
      sum_source = 0
      sum_puits = 0
      sum_source_XX = 0
      sum_puits_XX = 0
      rate = 0   
      sourceXX = 0
      puitsXX = 0  
      kicovi = 0   
*******************
      WRITE(50,*) '-----------------------------'
      WRITE(50,*) '------ sources de ',chrsp(idXX),' au temps ',time 

      WRITE(51,*) '-----------------------------'
      WRITE(51,*) '------ puits de ',chrsp(idXX),' au temps ',time 

**************************************************************************************
*        calcul de la somme des sources et des puits de l'espece recherchée          *
**************************************************************************************
      DO 23 i=1,numre
        rate(i)=qfor(i)
        DO k=1,numstoi(i,1)
          rate(i)=rate(i)*cbot(idrestoi(i,k))                              ! calcul de la vitesse associée à la réaction i
          IF (restoicf(i,k).EQ.2) rate(i)=rate(i)*cbot(idrestoi(i,k))
        ENDDO

        DO j=1,numstoi(i,1)
          DO k=1,numstoi(i,2)
            IF ((idrestoi(i,j).EQ.idpdstoi(i,k)).AND.
     &          (idrestoi(i,j).EQ.idXX).AND.
     &          (restoicf(i,j).EQ.pdstoicf(i,k))) GOTO 23                  ! si une espèce est à la fois réactif et produits, on ne la compte pas ds les sources/puits
          ENDDO
        ENDDO

        DO k=1,numstoi(i,2)                                                ! boucle sur le nombre de produits
          IF (idpdstoi(i,k).eq.idXX) THEN
            sum_source=sum_source+rate(i)*pdstoicf(i,k)                    ! si un des produits est l'espèce recherchée, on somme les sources de i * le coeff stochio
          ENDIF
        ENDDO

        IF ((idrestoi(i,1).eq.idXX).OR.(idrestoi(i,2).eq.idXX))THEN
            sum_puits=sum_puits+rate(i)                                    ! si un des réactifs est l'espèce recherchée, on somme les puits de i
            IF ((idrestoi(i,1).eq.idXX).AND.(idrestoi(i,2).eq.idXX)) 
     &        sum_puits=sum_puits+rate(i)                                  ! si A + A -> B le puits de A est multiplié par 2
        ENDIF
23    CONTINUE

*******************************************************************************************
* calcul de la contribution de l'espece recherchée par rapport au total des sources/puits *
*******************************************************************************************
      nbsource=1
      nbpuits=1
      DO 24 i=1,numre
        
        DO j=1,numstoi(i,1)
          DO k=1,numstoi(i,2)
            IF ((idrestoi(i,j).EQ.idpdstoi(i,k)).AND.
     &          (idrestoi(i,j).EQ.idXX).AND.
     &          (restoicf(i,j).EQ.pdstoicf(i,k))) GOTO 24                  ! si une espèce est à la fois réactif et produits, on ne la compte pas ds les sources/puits
          ENDDO
        ENDDO

        IF (i.LT.46) THEN                    !!!!!!!!!!! SCHEMA INORGANIQUE !!!!!!!!!!!!!!!!        
          DO k=1,numstoi(i,2)                                              ! boucle sur le nombre de produits pour chaque réaction
            IF (idpdstoi(i,k).eq.idXX) THEN
c              sourceXX(nbsource)=100*rate(i)*pdstoicf(i,k)/sum_source      ! la source i est exprimée en pourcentage par rapport à la somme des sources de l'espèce
              sourceXX(nbsource)=rate(i)*pdstoicf(i,k)

              WRITE(50,'(i2,2x,i2,2x,f13.6,2x,e10.3,2x,a7,a3,a7,a2,a7,         
     &        a3,a7)')  nbsource,i,100000*rate(i)/sum_source,rate(i),
     &              chrsp(idrestoi(i,1)),' + ', chrsp(idrestoi(i,2)),
     &         '=>',chrsp(idpdstoi(i,1)),' + ', chrsp(idpdstoi(i,2))       ! ce fichier ne sert qu'à identifier quelle réaction correspond à la source "nbsource"

              nbsource=nbsource+1
            ENDIF
          ENDDO

          IF ((idrestoi(i,1).eq.idXX).OR.                                 ! si un des réactifs est l'espece recherchée, la réaction i est un puits
     &      (idrestoi(i,2).eq.idXX)) THEN
c            puitsXX(nbpuits)=100*rate(i)/sum_puits                         ! le puits i est exprimé en pourcentage par rapport à la somme des puits de l'espèce
            puitsXX(nbpuits)=rate(i)
            IF ((idrestoi(i,1).eq.idXX).AND.(idrestoi(i,2).eq.idXX)) 
     &       puitsXX(nbpuits)=puitsXX(nbpuits)*2                         ! si A + A -> B le puits de A est multiplié par 2

            
            WRITE(51,'(i2,2x,i2,2x,f13.6,2x,e10.3,2x,a7,a3,a7,a2,a7,
     &       a3,a7)')  nbpuits,i,100000*rate(i)/sum_puits,rate(i),
     &            chrsp(idrestoi(i,1)),' + ', chrsp(idrestoi(i,2)),
     &       '=>',chrsp(idpdstoi(i,1)),' + ', chrsp(idpdstoi(i,2))         ! ce fichier ne sert qu'à identifier quelle réaction correspond au puits "nbpuits"

            nbpuits=nbpuits+1

          ENDIF

        ELSE                                    !!!!!!!!!!! SCHEMA ORGANIQUE !!!!!!!!!!!!!!!!

          DO k=1,numstoi(i,2)                                               ! boucle sur le nombre de produits pour chaque réaction
            IF (idpdstoi(i,k).eq.idXX) THEN
              sum_source_XX=sum_source_XX+rate(i)*pdstoicf(i,k)             ! la source i est exprimée en pourcentage par rapport à la somme des sources de l'espèce
            ENDIF
          ENDDO
          IF ((idrestoi(i,1).eq.idXX).OR.                             
     &           (idrestoi(i,2).eq.idXX)) THEN                              ! si un des réactifs est l'espece recherchée, la réaction i est un puit
            sum_puits_XX=sum_puits_XX+rate(i)                               ! les puits organiques sont sommés
          ENDIF


          IF ((chrsp(idrestoi(i,1))(2:2).NE.'1').AND.
     &        (chrsp(idrestoi(i,1))(2:2).NE.'2').AND.
     &        (chrsp(idrestoi(i,1))(2:2).NE.'3')) THEN                      ! sommes des kiCOVi
            kicovi=kicovi+qfor(i)*cbot(idrestoi(i,1))
          ENDIF

        ENDIF            
24    CONTINUE


c      sourceXX(nbsource)=100*sum_source_XX/sum_source                       ! la source organique est exprimée en pourcentage par rapport à la somme des sources de l'espèce                  
c      puitsXX(nbpuits)=100*sum_puits_XX/sum_puits                           ! le puits organique est exprimé en pourcentage par rapport à la somme des puits de l'espèce
      sourceXX(nbsource)=sum_source_XX              
      puitsXX(nbpuits)=sum_puits_XX

      RETURN

      END
