      MODULE module_chamber_tools
      use akparameter_module, only: maxlsp
      use module_data_gecko_main, only: chrsp, numsp, rem
      use flags_module, only: inject_fg, lightsonoff_fg
      use time_mgmt_module, only: time
      use io_units_module, only: lout
      implicit none
      
      ! a dervied type that contains information about injected species (duration of injection and final amount)
      type injection
        ! user input
        character(maxlsp) :: code
        real :: final_conc
        real :: tstart, tend

        ! computed
        integer :: ind
        real    :: duration, em_rate
      end type injection


      type light
        ! user input
        real   :: ton, toff
      end type light

      type(light), allocatable        :: lights(:)
      integer                         :: mlights, nlights
      type(injection), allocatable    :: injections(:)
      integer                         :: minject,ninject

      real  :: A_V_ratio, vol
      contains

      subroutine init_injections(n)
        integer, intent(in)  :: n
        inject_fg = 1
        minject = n
        allocate(injections(minject))
        ninject = 0
      end subroutine

      subroutine set_injection( code, tstart,tend ,final_conc)
        real, intent(in) :: tstart, tend, final_conc
        character(maxlsp), intent(in) :: code

        write(lout,*) 'injecting ',final_conc, 'of ' ,code
        write(lout,* ) '   from ', tstart, ' to ', tend
        if (tend <= tstart) then
          write(lout,*) "error trying to inject with reversed time,", tend, "<", tstart
          stop
        endif

        ninject = ninject + 1
        if (ninject > minject) then
            write(lout,*) " too many injected species, ", ninject, ">",minject
            stop
        endif
        injections(ninject)%tstart = tstart
        injections(ninject)%tend = tend
        injections(ninject)%final_conc = final_conc
        injections(ninject)%code = code

        call akspnum(code, chrsp, numsp, injections(ninject)%ind)

        if (injections(ninject)%ind == 0) then
          write(lout,*) " could not find injected species: ", code
          stop
        endif

        injections(ninject)%duration = tend - tstart
        injections(ninject)%em_rate  = final_conc / (tend - tstart)

      end subroutine
      
      subroutine run_chamber_injections()
        integer :: i
        print*, ninject, time
        do i=1,ninject
           if (time < injections(i)%tstart .or. time > injections(i)%tend) cycle
           print *, 'injecting ', injections(i)%code
           print *, injections(i)%ind, injections(i)%em_rate
           rem(injections(i)%ind) = rem(injections(i)%ind) + injections(i)%em_rate
        enddo

      end subroutine 

      subroutine init_lights(n)
        integer, intent(in) :: n
        lightsonoff_fg = 1
        mlights = n
        allocate(lights(mlights))
        nlights = 0

      end subroutine init_lights

      subroutine set_light(ton, toff)
        real, intent(in) :: ton, toff
        if (toff <= ton) then
          write(lout,*) "error trying to setup light with reversed time,", toff, "<", ton
          stop
        endif
        nlights = nlights + 1
        if (nlights > mlights) then
          write(lout, *) " too many lights, ", nlights, ">", mlights
          stop
        endif
        lights(nlights)%ton = ton
        lights(nlights)%toff = toff

        write(lout, *) 'light is on from', ton, ' to ', toff

      end subroutine

      subroutine run_chamber_lights(xd)
        real, intent(out) :: xd

        integer :: i
        logical :: lighton
        lighton = .false.
        do i = 1,nlights
          if (time > lights(i)%ton .and. time < lights(i)%toff) then
            lighton = .true. !we are in a timeperiod when lights should be on
          endif
        enddo

        if (.not.lighton) then
          xd = 100. ! zen. angle > 90. = night
        endif
 

      end subroutine

      END MODULE  
