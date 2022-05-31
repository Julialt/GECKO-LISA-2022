***********************************************************************
* calcul du pourcentage de surface au temps time                      *
*     => time  : temps local                                          *
*     => mhd   : nombre maximum de donnees surf=f(time)               *
*     => msur   : nombre maximum de type de surface                   *
*     => surft(i) : temps (s) pour lequel la surface est tabulee      *
*     => psurf(i,j) : pourcentage de la surface j au temps i          *
* variables sorties                                                   *
*     => xsurf(j) : pourcentage de surface j au temps time            *
***********************************************************************
      SUBROUTINE timesurf3(lout,time,mhd,msur,nsd,psurf,surft,xsurf)
      IMPLICIT NONE

* input
      INTEGER lout,mhd,msur,nsd
      REAL    time
      REAL    psurf(mhd,msur),surft(mhd)

* output
      REAL    xsurf(msur)

* local
      INTEGER i,j
      REAL    slope
       
* controle que le temps est plus grand ou egal a la
* permiere valeur en temps tabule
      IF (time.lt.surft(1)) THEN
         WRITE (lout,*) '--error--, in timesurf3'
         WRITE (lout,*) 'time lesser than surft(1)'
         stop
      ENDIF

* calcul de la surface au temps time
      DO i=1,nsd-1
        IF (surft(i+1).GT.time) THEN
          DO j=1,msur
            slope=(psurf(i+1,j)-psurf(i,j))/(surft(i+1)-surft(i))
            xsurf(j)=psurf(i,j) + slope*(time-surft(i))
          ENDDO
          GOTO 10
        ENDIF 
      ENDDO
* intervalle en temps non trouve
      WRITE (lout,*) '--error--, in timesurf3'
      WRITE (lout,*) 'upper limit for time not found'
      write(lout,*) time, surft(i+1)
      stop
10    CONTINUE

      END
