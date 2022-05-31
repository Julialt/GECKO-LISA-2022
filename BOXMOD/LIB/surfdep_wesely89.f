! This subroutine provides wesely (1989) resistance values for
! the 4 different types of surfaces, for 2 different seasons

      subroutine surfdep_wesely89(lout, iseas,
     &           rik,rlu0k,rack,rgsSk,rgsOk,rclSk,rclOk,z0k)

      USE akparameter_module
      implicit none

!inputs
      integer, intent(in)   :: lout  ! file number for error reportin g
      integer, intent(in)   :: iseas ! season index(1 = summer, 2 =winter)

!outputs
      real, intent(out) :: rik(msur)   !resistance stomatale minimum pour la vapeur d'eau
                                       !pour la surface (i)
      real, intent(out) :: rlu0k(msur) !resistance des feuilles de surface dans la canope
                                       !superieure ou resistance cuticulaire des feuilles
                                       !pour la surface (i)
      real, intent(out) :: rack(msur)  !resistance de transfert (depend de la hauteur et
                                       !la densite de la canope) pour la surface (i)
      real, intent(out) :: rgsSk(msur) !resistance de piegage au sol pour SO2
                                       !pour la surface (i)
      real, intent(out) :: rgsOk(msur) !resistance de piegage au sol pour O3
                                       !pour la surface (i)
      real, intent(out) :: rclSk(msur) !resistance des composants de la partie inferieure
                                       !de la canopee (ecorce, brindilles...) pour SO2
                                       !pour la surface (i)
      real, intent(out) :: rclOk(msur) !resistance des composants de la partie inferieure
                                       !de la canopee (ecorce, brindilles...) pour O3
                                       !pour la surface (i)
      real, intent(out) :: z0k(msur)   !hauteur de rugosite pour la surface (i)

! local
      integer, parameter :: i_cultivated = 2
      integer, parameter :: i_decidforest= 3
      integer, parameter :: i_conforest  = 4
      integer, parameter :: i_urban      = 1


      if (iseas == 1) then
        rik(i_cultivated)     = 9999.
        rik(i_decidforest)    = 9999.
        rik(i_conforest)      = 250.
        rik(i_urban)          = 9999.

        rlu0k(i_cultivated)   = 9999.
        rlu0k(i_decidforest)  = 9000.
        rlu0k(i_conforest)    = 4000.
        rlu0k(i_urban)        = 9999.

        rack(i_cultivated)    = 10.
        rack(i_decidforest)   = 1000.
        rack(i_conforest)     = 2000.
        rack(i_urban)         = 100.

        rgsSk(i_cultivated)   = 150.
        rgsSk(i_decidforest)  = 500.
        rgsSk(i_conforest)    = 500.
        rgsSk(i_urban)        = 400.

        rgsOk(i_cultivated)   = 150.
        rgsOk(i_decidforest)  = 200.
        rgsOk(i_conforest)    = 200.
        rgsOk(i_urban)        = 300.

        rclSk(i_cultivated)   = 9999.
        rclSk(i_decidforest)  = 9000.
        rclSk(i_conforest)    = 3000.
        rclSk(i_urban)        = 9999.

        rclOk(i_cultivated)   = 1000.
        rclOk(i_decidforest)  = 400.
        rclOk(i_conforest)    = 1000.
        rclOk(i_urban)        = 9999.

        z0k(i_cultivated)     = 0.05
        z0k(i_decidforest)    = 1.
        z0k(i_conforest)      = 1.
        z0k(i_urban)          = 1.

      else if (iseas == 2) then
        rik(i_cultivated)     = 60.
        rik(i_decidforest)    = 70.
        rik(i_conforest)      = 130.
        rik(i_urban)          = 9999.

        rlu0k(i_cultivated)   = 2000.
        rlu0k(i_decidforest)  = 2000.
        rlu0k(i_conforest)    = 2000.
        rlu0k(i_urban)        = 9999.

        rack(i_cultivated)    = 200.
        rack(i_decidforest)   = 2000.
        rack(i_conforest)     = 2000.
        rack(i_urban)         = 100.

        rgsSk(i_cultivated)   = 150.
        rgsSk(i_decidforest)  = 500.
        rgsSk(i_conforest)    = 500.
        rgsSk(i_urban)        = 400.

        rgsOk(i_cultivated)   = 150.
        rgsOk(i_decidforest)  = 200.
        rgsOk(i_conforest)    = 200.
        rgsOk(i_urban)        = 300.

        rclSk(i_cultivated)   = 2000.
        rclSk(i_decidforest)  = 2000.
        rclSk(i_conforest)    = 2000.
        rclSk(i_urban)        = 9999.

        rclOk(i_cultivated)   = 1000.
        rclOk(i_decidforest)  = 1000.
        rclOk(i_conforest)    = 1000.
        rclOk(i_urban)        = 9999.

        z0k(i_cultivated)     = 0.05
        z0k(i_decidforest)    = 1.
        z0k(i_conforest)      = 1.
        z0k(i_urban)          = 1.
      endif

      end subroutine
