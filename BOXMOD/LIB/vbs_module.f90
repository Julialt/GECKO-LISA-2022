      MODULE vbs_module 
            
      implicit none
      
      real, save :: no_br, fragbound
      real, parameter :: Aro2ho2 = 2.91e-13
      real, parameter :: Aro2ro2 = 2.5e-13
      real, parameter :: Aro2no  = 2.7e-12
      real, parameter :: Aro2no3 = 2.3e-12
      !real, parameter :: small = 1e-35
      
      
      CONTAINS
      
      SUBROUTINE init_vbs_rates(temp)
        USE module_data_gecko_main, only : conc, idho2,idno, cro2, ibox, idno3, small
        real, intent(in) :: temp
        
        real :: kro2ho2, kro2ro2, kro2no, kro2no3
        real :: ho2_rate, ro2_rate, no_rate,no3_rate
        
        kro2ho2 = Aro2ho2*exp(1300./temp)
        kro2ro2 = Aro2ro2
        kro2no =  Aro2no*exp(360./temp)
        kro2no3 = Aro2no3
        
        ho2_rate = conc(idho2, ibox)*kro2ho2
        ro2_rate = sum(cro2(:))*kro2ro2
        no_rate  = conc(idno, ibox)*kro2no
        no3_rate = conc(idno3, ibox)*kro2no3
        
        no_br = (no_rate + no3_rate)/(ho2_rate + ro2_rate + no_rate + no3_rate + small)
        fragbound = min(1. - no_br,0.75)
      
      END SUBROUTINE
      
      REAL FUNCTION KVBS(A, B, C, temp, y_lonox, y_hinox)
        real, intent(in) :: A, B, C
        real, intent(in) :: temp
        real, intent(in) :: y_lonox, y_hinox
        
        KVBS = A*exp(-B/temp)*(temp/300.)**C &
          *(y_lonox*(1-no_br) + &
            y_hinox*no_br)
      
      END FUNCTION KVBS
      
      REAL FUNCTION KVBS3(A, B, C, temp, y)
        real, intent(in) :: A, B, C
        real, intent(in) :: temp
        real, intent(in) :: y
        
        KVBS3 = y*(A * EXP(- B /temp ) &
               *(temp/300.)**C)
       END FUNCTION KVBS3
       
       REAL FUNCTION KVBS4(A)
       
         real, intent(in) :: A

         KVBS4 = A*max(1. - fragbound - 0.1333*fragbound, 0.)*1.15
       
       END FUNCTION KVBS4
       
       REAL FUNCTION KVBS5(A)
         real, intent(in) :: A
         KVBS5 = A*fragbound
       
       END FUNCTION KVBS5
      
      END MODULE

