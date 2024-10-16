      subroutine wet_salt_output(j) !rtb salt
    
      use output_ls_pesticide_module
      use res_pesticide_module
      use res_salt_module
      use plant_module
      use plant_data_module
      use time_module
      use basin_module
      use output_landscape_module
      use constituent_mass_module
      use hydrograph_module, only : sp_ob1, ob
      
      implicit none
      
      integer :: isalt = 0
      integer :: j
      integer :: iob = 0
      real :: const = 0.
                         
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine outputs salt ion mass in wetlands (by HRU)
      
      !objective number of HRU
      iob = sp_ob1%hru + j - 1
       
      !add daily values to monthly values
      do isalt = 1, cs_db%num_salts
        wetsalt_m(j)%salt(isalt)%inflow = wetsalt_m(j)%salt(isalt)%inflow + wetsalt_d(j)%salt(isalt)%inflow
        wetsalt_m(j)%salt(isalt)%outflow = wetsalt_m(j)%salt(isalt)%outflow + wetsalt_d(j)%salt(isalt)%outflow
        wetsalt_m(j)%salt(isalt)%seep = wetsalt_m(j)%salt(isalt)%seep + wetsalt_d(j)%salt(isalt)%seep
        wetsalt_m(j)%salt(isalt)%fert = wetsalt_m(j)%salt(isalt)%fert + wetsalt_d(j)%salt(isalt)%fert
        wetsalt_m(j)%salt(isalt)%irrig = wetsalt_m(j)%salt(isalt)%irrig + wetsalt_d(j)%salt(isalt)%irrig
        wetsalt_m(j)%salt(isalt)%div = wetsalt_m(j)%salt(isalt)%div + wetsalt_d(j)%salt(isalt)%div
        wetsalt_m(j)%salt(isalt)%mass = wetsalt_m(j)%salt(isalt)%mass + wetsalt_d(j)%salt(isalt)%mass
        wetsalt_m(j)%salt(isalt)%conc = wetsalt_m(j)%salt(isalt)%conc + wetsalt_d(j)%salt(isalt)%conc
      enddo
      wetsalt_m(j)%salt(1)%volm = wetsalt_m(j)%salt(1)%volm + wetsalt_d(j)%salt(1)%volm
      
      !daily print
      if (pco%salt_res%d == "y") then
        write (5090,100) time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                         (wetsalt_d(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                         (wetsalt_d(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                          wetsalt_d(j)%salt(1)%volm
        if (pco%csvout == "y") then
          write (5091,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                                        (wetsalt_d(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                                        (wetsalt_d(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                                         wetsalt_d(j)%salt(1)%volm
        endif
      endif

      !monthly print
      if (time%end_mo == 1) then
        !add monthly values to yearly values
        do isalt = 1, cs_db%num_salts
          wetsalt_y(j)%salt(isalt)%inflow = wetsalt_y(j)%salt(isalt)%inflow + wetsalt_m(j)%salt(isalt)%inflow
          wetsalt_y(j)%salt(isalt)%outflow = wetsalt_y(j)%salt(isalt)%outflow + wetsalt_m(j)%salt(isalt)%outflow
          wetsalt_y(j)%salt(isalt)%seep = wetsalt_y(j)%salt(isalt)%seep + wetsalt_m(j)%salt(isalt)%seep
          wetsalt_y(j)%salt(isalt)%fert = wetsalt_y(j)%salt(isalt)%fert + wetsalt_m(j)%salt(isalt)%fert
          wetsalt_y(j)%salt(isalt)%irrig = wetsalt_y(j)%salt(isalt)%irrig + wetsalt_m(j)%salt(isalt)%irrig
          wetsalt_y(j)%salt(isalt)%div = wetsalt_y(j)%salt(isalt)%div + wetsalt_m(j)%salt(isalt)%div
          wetsalt_y(j)%salt(isalt)%mass = wetsalt_y(j)%salt(isalt)%mass + wetsalt_m(j)%salt(isalt)%mass
          wetsalt_y(j)%salt(isalt)%conc = wetsalt_y(j)%salt(isalt)%conc + wetsalt_m(j)%salt(isalt)%conc
        enddo
        wetsalt_y(j)%salt(1)%volm = wetsalt_y(j)%salt(1)%volm + wetsalt_m(j)%salt(1)%volm
        const = float (ndays(time%mo + 1) - ndays(time%mo))
        do isalt=1,cs_db%num_salts !average mass and concentration
          wetsalt_m(j)%salt(isalt)%mass = wetsalt_m(j)%salt(isalt)%mass / const
          wetsalt_m(j)%salt(isalt)%conc = wetsalt_m(j)%salt(isalt)%conc / const
        enddo
        wetsalt_m(j)%salt(1)%volm = wetsalt_m(j)%salt(1)%volm / const
        if (pco%salt_res%m == "y") then
          write (5092,100) time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                           (wetsalt_m(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                           (wetsalt_m(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                            wetsalt_m(j)%salt(1)%volm
          if (pco%csvout == "y") then
            write (5093,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                                          (wetsalt_m(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                                          (wetsalt_m(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                                           wetsalt_m(j)%salt(1)%volm
          endif
        endif
        !zero out
        do isalt = 1, cs_db%num_salts
          wetsalt_m(j)%salt(isalt)%inflow = 0.
          wetsalt_m(j)%salt(isalt)%outflow = 0.
          wetsalt_m(j)%salt(isalt)%seep = 0.
          wetsalt_m(j)%salt(isalt)%fert = 0.
          wetsalt_m(j)%salt(isalt)%irrig = 0.
          wetsalt_m(j)%salt(isalt)%div = 0.
          wetsalt_m(j)%salt(isalt)%mass = 0.
          wetsalt_m(j)%salt(isalt)%conc = 0.
        enddo
        wetsalt_m(j)%salt(1)%volm = 0.
      endif
      
      !yearly print
      if (time%end_yr == 1) then
        !add yearly values to total values
        do isalt = 1, cs_db%num_salts
          wetsalt_a(j)%salt(isalt)%inflow = wetsalt_a(j)%salt(isalt)%inflow + wetsalt_y(j)%salt(isalt)%inflow
          wetsalt_a(j)%salt(isalt)%outflow = wetsalt_a(j)%salt(isalt)%outflow + wetsalt_y(j)%salt(isalt)%outflow
          wetsalt_a(j)%salt(isalt)%seep = wetsalt_a(j)%salt(isalt)%seep + wetsalt_y(j)%salt(isalt)%seep
          wetsalt_a(j)%salt(isalt)%fert = wetsalt_a(j)%salt(isalt)%fert + wetsalt_y(j)%salt(isalt)%fert
          wetsalt_a(j)%salt(isalt)%irrig = wetsalt_a(j)%salt(isalt)%irrig + wetsalt_y(j)%salt(isalt)%irrig
          wetsalt_a(j)%salt(isalt)%div = wetsalt_a(j)%salt(isalt)%div + wetsalt_y(j)%salt(isalt)%div
          wetsalt_a(j)%salt(isalt)%mass = wetsalt_a(j)%salt(isalt)%mass + wetsalt_y(j)%salt(isalt)%mass
          wetsalt_a(j)%salt(isalt)%conc = wetsalt_a(j)%salt(isalt)%conc + wetsalt_y(j)%salt(isalt)%conc
        enddo
        wetsalt_a(j)%salt(1)%volm = wetsalt_a(j)%salt(1)%volm + wetsalt_y(j)%salt(1)%volm
        const = time%day_end_yr
        do isalt=1,cs_db%num_salts !average mass and concentration
          wetsalt_y(j)%salt(isalt)%mass = wetsalt_y(j)%salt(isalt)%mass / const
          wetsalt_y(j)%salt(isalt)%conc = wetsalt_y(j)%salt(isalt)%conc / const
        enddo
        wetsalt_y(j)%salt(1)%volm = wetsalt_y(j)%salt(1)%volm / const
        if (pco%salt_res%y == "y") then
          write (5094,100) time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                           (wetsalt_y(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                           (wetsalt_y(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                            wetsalt_y(j)%salt(1)%volm
          if (pco%csvout == "y") then
            write (5095,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                                          (wetsalt_y(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                                          (wetsalt_y(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                                           wetsalt_y(j)%salt(1)%volm
          endif
        endif
        !zero out
        do isalt = 1, cs_db%num_salts
          wetsalt_y(j)%salt(isalt)%inflow = 0.
          wetsalt_y(j)%salt(isalt)%outflow = 0.
          wetsalt_y(j)%salt(isalt)%seep = 0.
          wetsalt_y(j)%salt(isalt)%fert = 0.
          wetsalt_y(j)%salt(isalt)%irrig = 0.
          wetsalt_y(j)%salt(isalt)%div = 0.
          wetsalt_y(j)%salt(isalt)%mass = 0.
          wetsalt_y(j)%salt(isalt)%conc = 0.
        enddo
        wetsalt_y(j)%salt(1)%volm = 0.
      endif
      
      !average annual print
      if (time%end_sim == 1 .and. pco%salt_res%a == "y") then
        !calculate average annual values
        do isalt = 1, cs_db%num_salts
          wetsalt_a(j)%salt(isalt)%inflow = wetsalt_a(j)%salt(isalt)%inflow / time%nbyr
          wetsalt_a(j)%salt(isalt)%outflow = wetsalt_a(j)%salt(isalt)%outflow / time%nbyr 
          wetsalt_a(j)%salt(isalt)%seep = wetsalt_a(j)%salt(isalt)%seep / time%nbyr 
          wetsalt_a(j)%salt(isalt)%fert = wetsalt_a(j)%salt(isalt)%fert / time%nbyr 
          wetsalt_a(j)%salt(isalt)%irrig = wetsalt_a(j)%salt(isalt)%irrig / time%nbyr
          wetsalt_a(j)%salt(isalt)%div = wetsalt_a(j)%salt(isalt)%div / time%nbyr
          wetsalt_a(j)%salt(isalt)%mass = wetsalt_a(j)%salt(isalt)%mass / time%nbyr 
          wetsalt_a(j)%salt(isalt)%conc = wetsalt_a(j)%salt(isalt)%conc / time%nbyr
        enddo
        wetsalt_a(j)%salt(1)%volm = wetsalt_a(j)%salt(1)%volm / time%nbyr
        write (5096,100) time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                         (wetsalt_a(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                         (wetsalt_a(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                          wetsalt_a(j)%salt(1)%volm
        if (pco%csvout == "y") then
          write (5097,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, j, ob(iob)%gis_id, & 
                                        (wetsalt_a(j)%salt(isalt)%inflow,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%outflow,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%seep,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%fert,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%irrig,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%div,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%mass,isalt=1,cs_db%num_salts), &
                                        (wetsalt_a(j)%salt(isalt)%conc,isalt=1,cs_db%num_salts), &
                                         wetsalt_a(j)%salt(1)%volm
        endif
      endif

      return
     
100   format (4i6,2i8,500e15.4)      

      end subroutine wet_salt_output