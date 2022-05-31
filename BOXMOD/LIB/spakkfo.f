      SUBROUTINE akkfo(maxre,maxfo,maxaux,
     1                 ire,ifo,focf,temp,sumcb,arrhcf,qln)
      IMPLICIT NONE

* INPUT
      INTEGER  maxre,maxfo,maxaux
      INTEGER  ire,ifo
      REAL     temp, sumcb
      REAL     focf(maxaux+3,maxfo)
      REAL     arrhcf(maxre,3)

* OUTPUT
      REAL     qln

* LOCAL
      REAL    effko,rapk,facteur,fcent


* expression fall off type TROE
* effko=constante efficace basse pression
* rapk=effko/Kinfinie
* Pour le coefficient F: 2 parametrisations:
* fcent=constante ---> T*** doit etre egal a 0
* fcent est fonction de la temperature ---> expression type chemkin

      effko = focf(1,ifo) * ((temp/300.)**focf(2,ifo)) *
     &        exp(-focf(3,ifo)/temp)*sumcb
      rapk = effko / (exp(arrhcf(ire,1))*((temp/300.)**arrhcf(ire,2))*
     &                exp(-arrhcf(ire,3)/temp))
      facteur=1./(1.+(log10(rapk))**2.)

      IF (focf(5,ifo).eq.0. .or. focf(6,ifo).eq. 0) THEN
         qln=log(effko/(1.+rapk))+facteur*log(focf(4,ifo))
      ELSE
         fcent=(1.-focf(4,ifo))*exp(-temp/focf(5,ifo))+focf(4,ifo)*
     &         exp(-temp/focf(6,ifo))+exp(-focf(7,ifo)/temp)
         IF (focf(7,ifo).eq.0.) fcent=fcent-1.
         qln=log(effko/(1.+rapk))+facteur*log(fcent)
      ENDIF
      END
