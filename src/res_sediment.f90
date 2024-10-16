      subroutine res_sediment

      use reservoir_data_module
      use reservoir_module
      use conditional_module
      use climate_module
      use time_module
      use hydrograph_module
      use water_body_module
      
      implicit none

      real :: trapres                       !              |
      real :: velofl                        !              |  
      real :: sed_ppm, sil_ppm, cla_ppm 

      if (wbody%flo < 1.e-6) then
        ! reservoir is empty
        wbody = hz
      else

        !! compute concentrations
	      if (wbody%flo > 0.) then
            sed_ppm = 1000000. * wbody%sed / wbody%flo
            sed_ppm = Max(1.e-6, sed_ppm)
            sil_ppm = 1000000. * wbody%sil / wbody%flo
            sil_ppm = Max(1.e-6, sil_ppm)
            cla_ppm = 1000000. * wbody%cla / wbody%flo
            cla_ppm = Max(1.e-6, cla_ppm)
	      else
            sed_ppm = 1.e-6
            sil_ppm = 1.e-6
            cla_ppm = 1.e-6
	      endif

        !! compute change in sediment concentration due to settling 
        if (sed_ppm > wbody_prm%sed%nsed) then
          sed_ppm = (sed_ppm - wbody_prm%sed%nsed) * wbody_prm%sed_stlr_co + wbody_prm%sed%nsed
          sed_ppm = Max (sed_ppm, wbody_prm%sed%nsed)
          wbody%sed = sed_ppm * wbody%flo / 1000000.      ! ppm -> t Jaehak 2024
          ht2%sed = sed_ppm * ht2%flo / 1000000.
          wbody%sed = Max(0.,wbody%sed - ht2%sed)
          
          sil_ppm = (sil_ppm - wbody_prm%sed%nsed) * wbody_prm%sed%sed_stlr + wbody_prm%sed%nsed
          wbody%sil = sil_ppm * wbody%flo / 1000000.      ! ppm -> t
          
          cla_ppm = (cla_ppm - wbody_prm%sed%nsed) * wbody_prm%sed%sed_stlr + wbody_prm%sed%nsed
          wbody%cla = cla_ppm * wbody%flo / 1000000.      ! ppm -> t

          !! assume all sand aggregates and gravel settles
          wbody%san = 0.
          wbody%sag = 0.
          wbody%lag = 0.
          wbody%grv = 0.
        end if

        !! compute sediment leaving reservoir - ppm -> t
        ht2%sed = sed_ppm * ht2%flo / 1000000.
        wbody%sed = wbody%sed - ht2%sed
        
        wbody%sed = max(0.,wbody%sed - ht2%sed)
        ht2%sil = sil_ppm * ht2%flo / 1000000.
        wbody%sil = max(0.,wbody%sil - ht2%sil)
        ht2%cla = cla_ppm * ht2%flo / 1000000.
        wbody%cla = max(0.,wbody%cla - ht2%cla)

      end if

      return
      end subroutine res_sediment