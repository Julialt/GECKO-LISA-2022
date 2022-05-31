      SUBROUTINE SPLNA (N,X,Y,c,J)                                 

      USE akparameter_module
c      implicit double precision (a-h,o-z)

      DIMENSION X(maxang), Y(maxang), D(2)                        
      DIMENSION c(nlo), W(nlo)                        
c
       d(1)=0.
       d(2)=0.
      !do 5 i=1,nlo
      do i=1,nlo
5      c(i)=0.
      enddo
C     -----------------------------------------------------------
C       OVER THE INTERVAL X(I) TO X(I+1), THE INTERPOLATING          
C       POLYNOMIAL                                                   
c            Y=Y(I)+A(I)*Z+B(I)*Z**2+E(I)*Z**3                    
C            WHERE Z=(X-X(I))/(X(I+1)-X(I))                            
C     IS USED. THE COEFFICIENTS A(I),B(I) AND E(I) ARE COMPUTED
C     BY SPLNA AND STORED IN LOCATIONS C(3*I-2),C(3*I-1) AND    
C     C(3*I) RESPECTIVELY.                                    
C     WHILE WORKING IN THE ITH INTERVAL,THE VARIABLE Q WILL     
C     REPRESENT Q=X(I+1) - X(I), AND Y(I) WILL REPRESENT       
C             Y(I+1)-Y(I)                                              
C     -------------------------------------------------------------
C                                                                      
      Q=X(2)-X(1)                                                       
      YI=Y(2)-Y(1)                                                     
      IF (J.EQ.2) GO TO 10                                             
C     ------------------------------------------------------------
C             IF THE FIRST DERIVATIVE AT THE END POINTS IS GIVEN,      
C             A(1) IS KNOWN, AND THE SECOND EQUATION BECOMES           
C             MERELY B(1)+E(1)=YI - Q*D(1).                           
C     ------------------------------------------------------------
      C(1)=Q*D(1)                                                      
      C(2)=1.0                                                        
      W(2)=YI-C(1)                                                     
      GO TO 20                                                       
C     -------------------------------------------------------------
C             IF THE SECOND DERIVATIVE AT THE END POINTS IS GIVEN     
C             B(1) IS KNOWN, THE SECOND EQUATION BECOMES              
C             A(1)+E(1)=YI-0.5*Q*Q*D(1). DURING THE SOLUTION OF        
C             THE 3N-4 EQUATIONS,A1 WILL BE KEPT IN CELL C(2)          
C             INSTEAD OF C(1) TO RETAIN THE TRIDIAGONAL FORM OF THE    
C             COEFFICIENT MATRIX.                                     
C     -------------------------------------------------------------
   10 C(2)=0.0                                                         
      W(2)=0.5*Q*Q*D(1)                                               
   20 M=N-2                                                            
      IF (M.LE.0) GO TO 40                                             
C     ------------------------------------------------------------
C             UPPER TRIANGULARIZATION OF THE TRIDIAGONAL SYSTEM OF     
C             EQUATIONS FOR THE COEFFICIENT MATRIX FOLLOWS--           
C     ------------------------------------------------------------
      !DO 30 I=1,M                                                     
      DO I=1,M                                                     
        AI=Q                       
        Q=X(I+2)-X(I+1)           
        H=AI/Q 
        C(3*I)=-H/(2.0-C(3*I-1)) 
        W(3*I)=(-YI-W(3*I-1))/(2.0-C(3*I-1))                            
        C(3*I+1)=-H*H/(H-C(3*I))
        W(3*I+1)=(YI-W(3*I))/(H-C(3*I))                                 
        YI=Y(I+2)-Y(I+1)                                               
        C(3*I+2)=1.0/(1.0-C(3*I+1))                                    
   30   W(3*I+2)=(YI-W(3*I+1))/(1.0-C(3*I+1))                          
      ENDDO
C     ------------------------------------------------------------
C             E(N-1) IS DETERMINED DIRECTLY FROM THE LAST EQUATION    
C             OBTAINED ABOVE, AND THE FIRST OR SECOND DERIVATIVE       
C             VALUE GIVEN AT THE END POINT.                           
C     ------------------------------------------------------------
   40 IF (J.EQ.1) GO TO 50                                             
      C(3*N-3)=(Q*Q*D(2)/2.0-W(3*N-4))/(3.0-C(3*N-4))                
      GO TO 60                                                        
   50 C(3*N-3)=(Q*D(2)-YI-W(3*N-4))/(2.0-C(3*N-4))                   
   60 M=3*N-6                                                          
      IF (M.LE.0) GO TO 80                                            
C     -----------------------------------------------------------
C             BACK SOLUTION FOR ALL COEFFICENTS EXCEPT                
C             A(1) AND B(1) FOLLOWS--                                  
C     -----------------------------------------------------------
      !DO 70 II=1,M                                                  
      DO II=1,M                                                  
        I=M-II+3                                                        
   70   C(I)=W(I)-C(I)*C(I+1)                                           
      ENDDO
   80 IF (J.EQ.1) GO TO 90                                            
C     -----------------------------------------------------------
C             IF THE SECOND DERIVATIVE IS GIVEN AT THE END POINTS,  
C             A(1) CAN NOW BE COMPUTED FROM THE KNOWN VALUES OF       
C             B(1) AND E(1). THEN A(1) AND B(1) ARE PUT INTO THEIR    
C             PROPER PLACES IN THE C ARRAY.                            
C     -----------------------------------------------------------
      C(1)=Y(2)-Y(1)-W(2)-C(3)                                         
      C(2)=W(2)                                                        
      RETURN                                                          
   90 C(2)=W(2)-C(3)                                                  
      RETURN                                                          
      END                                                           
