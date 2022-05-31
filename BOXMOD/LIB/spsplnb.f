      SUBROUTINE SPLNB (N,X,Y,c,V)                                   

      USE akparameter_module
c      implicit double precision (a-h,o-z)

      DIMENSION X(maxang), Y(maxang), V(5),c(nlo)                           
      V(5)=2.0                                                        
      LIM=N-1                                                         
C                  --------------------------------------------
C                  DETERMINE IN WHICH INTERVAL THE INDEPENDENT        
C                  VARIABLE,V(1),LIES.                               
C                  --------------------------------------------
      DO 10 I=2,LIM                                                   
      IF (V(1).LT.X(I)) GO TO 20                                      
   10 CONTINUE                                                         
      I=N                                                              
      IF (V(1).GT.X(N)) V(5)=3.0                                      
      GO TO 30                                                       
   20 IF (V(1).LT.X(1)) V(5)=1.0                                     
C                  ---------------------------------------------
C                  Q IS THE SIZE OF THE INTERVAL CONTAINING V(1).     
C                  ---------------------------------------------
C                  Z IS A LINEAR TRANSFORMATION OF THE INTERVAL       
C                  ONTO (0,1) AND IS THE VARIABLE FOR WHICH           
C                  THE COEFFICIENTS WERE COMPUTED BY SPLNA.         
C                  --------------------------------------------
   30 Q=X(I)-X(I-1)                                                   
      Z=(V(1)-X(I-1))/Q                                               
      V(2)=((Z*C(3*I-3)+C(3*I-4))*Z+C(3*I-5))*Z+Y(I-1)              
      V(3)=((3.*Z*C(3*I-3)+2.0*C(3*I-4))*Z+C(3*I-5))/Q              
      V(4)=(6.*Z*C(3*I-3)+2.0*C(3*I-4))/(Q*Q)                       
      RETURN                                                         
      END                                                           
