      subroutine pesticide_init

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calls subroutines which read input data for the 
!!    databases and the HRUs

        use hru_module, only : hru, sol_plt_ini
        use soil_module
        use organic_mineral_mass_module
        use constituent_mass_module
        use output_ls_pesticide_module
        use hydrograph_module, only : sp_ob
        use plant_module
        use pesticide_data_module
      
        implicit none 
        
        integer :: ihru = 0        !none          !counter       
        integer :: npmx = 0        !none          |total number of pesticides
        integer :: ly = 0          !none          |counter
        integer :: ipest = 0       !none          |counter    
        integer :: nly = 0         !none          |max soil layers
        integer :: npl = 0         !none          |max plants
        integer :: ipest_db = 0    !              | 
        integer :: isp_ini = 0     !              |
        integer :: ipl = 0         !none          |plant number
        real :: wt1 = 0.           !              |
        real :: solpst = 0.        !              |
        real :: pl_frac = 0.       !0-1           |fraction of pesticide applied to each plant
        
      !! allocate hru pesticides
      do ihru = 1, sp_ob%hru
        npmx = cs_db%num_pests
        if (npmx > 0) then
          nly = soil(ihru)%nly
          npl = pcom(ihru)%npl
          allocate (cs_soil(ihru)%ly(nly))
          allocate (cs_pl(ihru)%pl_in(npl))
          allocate (cs_pl(ihru)%pl_on(npl))
          allocate (cs_pl(ihru)%pl_up(npl))
          do ly = 1, nly
            allocate (cs_soil(ihru)%ly(ly)%pest(npmx), source = 0.)
            allocate (cs_pl(ihru)%pl_in(ipl)%pest(npmx), source = 0.)
            allocate (cs_pl(ihru)%pl_on(ipl)%pest(npmx), source = 0.)
            allocate (cs_pl(ihru)%pl_up(ipl)%pest(npmx), source = 0.)
            cs_soil(ihru)%ly(ly)%pest = 0.
          end do
          do ipl = 1, npl
            cs_pl(ihru)%pl_in(ipl)%pest = 0.
            cs_pl(ihru)%pl_on(ipl)%pest = 0.
            cs_pl(ihru)%pl_up(ipl)%pest = 0.
          end do
          allocate (cs_irr(ihru)%pest(npmx))
          cs_irr(ihru)%pest = 0.
        end if

        isp_ini = hru(ihru)%dbs%soil_plant_init
        ipest_db = sol_plt_ini(isp_ini)%pest
        do ipest = 1, npmx
          hpestb_d(ihru)%pest(ipest)%plant = pest_soil_ini(ipest_db)%plt(ipest)
          do ipl = 1, pcom(ihru)%npl
            if (pcom(ihru)%lai_sum > 1.e-6) then
              pl_frac = pcom(ihru)%plg(ipl)%lai / pcom(ihru)%lai_sum
            else
              pl_frac = 0.
            end if
            pl_frac = Min (pl_frac, 1.)
            cs_pl(ihru)%pl_on(ipl)%pest(ipest) = cs_pl(ihru)%pl_on(ipl)%pest(ipest) + pl_frac * pest_soil_ini(ipest_db)%plt(ipest)
          end do
          
          solpst = pest_soil_ini(ipest_db)%soil(ipest)
          do ly = 1, soil(ihru)%nly
            wt1 = soil(ihru)%phys(ly)%bd * soil(ihru)%phys(ly)%thick / 100.      !! mg/kg => kg/ha
            cs_soil(ihru)%ly(ly)%pest(ipest) = solpst * wt1
          end do
        end do

      end do    ! hru loop
                                   
      return
      end subroutine pesticide_init