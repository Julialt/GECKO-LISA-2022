      SUBROUTINE SOLAR (SLA,SLO,TZ,IY,IM,ID,TIME,D,NV)                
c      SAVE                                                          
C***                                                                 
C***     SLA...  LATITUDE (DEG)  SOUTH = MINUS                      
C***     SLO...  LONGITUDE (DEG)  EAST = MINUS                       
C***     TZ...  TIME ZONE                                          
C***            ALSO INCLUDES FRACTION IF LOCAL TIME IS NOT           
C***             STANDARD MERIDIAN TIME.  E.G. POONA, INDIA 5.5      
C***     IY..  YEAR                                                  
C***     IM..  MONTH                                                 
C***     ID..  DAY                                                   
C***     TIME.. LOCAL STANDARD TIME IN HOURS AND MINUTES.            
C***            1 30 PM = 1330  ** STANDARD TIME **                   
C***     D..  RETURNED VALUE                                         
C***     NV..  VALUE TO BE RETURNED, SELECTED AS FOLLOWS....         
C***           1...  DECLINATION (DEG.)                            
C***           2...  EQUATION OF TIME ADJUSTMENT (HRS.)               
C***           3...  TRUE SOLAR TIME (HRS.)                            
C***           4...  HOUR ANGLE (DEG.)                                 
C***           5...  SOLAR ELEVATION (DEG.)                           
C***           6...  OPTICAL AIRMASS                               
C***     0 ( NV ( 7.  OTHERWISE, D = 9999.                            
C***                                                                
c      implicit double precision (a-h,o-z)
      DIMENSION MD(11)                                               
      DATA MD/31,29,31,30,31,30,2*31,30,31,30/                        
      DATA A,B,C,SIGA/0.15,3.885,1.253,279.9348/                      
c      RAD=572957.75913E-4                                             
      RAD=572957.76E-4                                             
c      SDEC=39784.988432E-5                                            
      SDEC=39784.988E-5                                            
      RE=1.                                                         
      IF (SLO.LT.0.) RE=-1.                                        
      KZ=TZ                                                           
      TC=(TZ-FLOAT(KZ))*RE                                           
      TZZ=FLOAT(KZ)*RE                                                
      SLB=SLA/RAD                                                    
      K=ID                                                           
      TIMH=TIME/100.                                                  
      I=TIMH                                                        
      TIMLOC=(TIMH-FLOAT(I))/0.6+FLOAT(I)+TC                         
      IMC=IM-1                                                       
      IF (IMC.LT.1) GO TO 20                                          
      !DO 10 I=1,IMC                                                    
      DO I=1,IMC   
        K=K+MD(I) 
   10 ENDDO  
   20 LEAP=1  
      NL=MOD(IY,4)                                                      
      IF (NL.LT.1) LEAP=2                                             
      SMER=TZZ*15.                                                     
      TK=((SMER-SLO)*4.)/60.                                          
      KR=1                                                              
      IF (K.GE.61.AND.LEAP.LT.2) KR=2                                 
      DAD=(TIMLOC+TZZ)/24.                                          
      DAD=DAD+FLOAT(K-KR)                                            
      DF=DAD*360./365.242                                              
      DE=DF/RAD                                                       
      DESIN=SIN(DE)                                                   
      DECOS=COS(DE)                                                    
      DESIN2=SIN(DE*2.0)                                               
      DECOS2=COS(DE*2.0)                                               
      SIG=SIGA+DF+1.914827*DESIN-0.079525*DECOS+0.019938*DESIN2
     &-0.00162*DECOS2                                                          
      SIG=SIG/RAD                                                    
      DECSIN=SDEC*SIN(SIG)                                            
      EFFDEC=ASIN(DECSIN)                                             
      IF (NV.NE.1) GO TO 30                                          
      D=EFFDEC*RAD                                                    
      RETURN                                                          
   30 EQT=0.12357*DESIN-0.004289*DECOS+0.153809*DESIN2+0.060783*
     &DECOS2 
      IF (NV.NE.2) GO TO 40                                           
      D=EQT                                                          
      RETURN                                                        
   40 TST=TK+TIMLOC-EQT                                            
      IF (NV.NE.3) GO TO 50                                           
      D=TST                                                            
      IF (D.LT.0.) D=D+24.                                            
      IF (D.GE.24.) D=D-24.                                                        
      RETURN                                                           
   50 HRANGL=ABS(TST-12.)*15.                                         
      IF (NV.NE.4) GO TO 60                                           
      D=HRANGL                                                         
      RETURN                                                           
   60 HRANGL=HRANGL/RAD                                              
      SOLSIN=DECSIN*SIN(SLB)+COS(EFFDEC)*COS(SLB)*COS(HRANGL)         
      SOLEL=ASIN(SOLSIN)*RAD                                           
      IF (NV.NE.5) GO TO 70                                           
      D=SOLEL
c      write(11,*)'angle',d,90.-d                                                         
      RETURN                                                           
   70 IF (NV.NE.6) GO TO 80                                           
      IF (SOLEL.LE.0.) GO TO 80                                       
      TK=SOLEL+B                                                      
      E=1./TK**C                                                      
      D=1./(A*E+SOLSIN)                                               
      RETURN                                                           
   80 D=9999.                                                          
      RETURN                                                          
      END                                                            
