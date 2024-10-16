      subroutine recall_read

      use hydrograph_module
      use input_file_module
      use organic_mineral_mass_module
      use constituent_mass_module
      use maximum_data_module
      use time_module
      use exco_module
      
      implicit none      
 
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      character(len=16) :: ob_name = ""
      character(len=8) :: ob_typ = ""
      integer :: imax = 0             !none       |end of loop
      integer :: iyr = 0              !           |
      integer :: jday = 0             !           |
      integer :: mo = 0               !           |
      integer :: day_mo = 0           !           |
      integer :: eof = 0              !           |end of file
      logical :: i_exist              !none       |check to determine if file exists
      integer :: nbyr = 0             !none       !number of years the land use occurred 
      integer :: k = 0                !           |
      integer :: iyrs = 0             !           | 
      integer :: iyr_prev = 0         !none       |previous year
      integer :: istep = 0            !           | 
      integer :: ipestcom_db = 0      !none       !pointer to pestcom_db - fix*** ?? 
      integer :: ipc = 0              !none       |counter
      integer :: ii = 0               !none       |counter
      integer :: i = 0                !           |
      integer :: iexco_om = 0
      integer :: iexo_allo = 0
      integer :: idaystep = 0
      integer :: jday1 = 0
      integer :: mo1 = 0
      integer :: iyr1 = 0
      
      eof = 0
      imax = 0

      !read all recall files
      inquire (file=in_rec%recall_rec, exist=i_exist)
      if (i_exist .or. in_rec%recall_rec /= "null") then
      do
        open (107,file=in_rec%recall_rec)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
        imax = 0
          do while (eof == 0)
            read (107,*,iostat=eof) i
            if (eof < 0) exit
            imax = Max(imax,i) 
          end do
          db_mx%recall_max = imax
          
      allocate (recall(0:imax))
      allocate (rec_d(imax))
      allocate (rec_m(imax))
      allocate (rec_y(imax))
      allocate (rec_a(imax))
      
      rewind (107)
      read (107,*,iostat=eof) titldum
      if (eof < 0) exit
      read (107,*,iostat=eof) header
      if (eof < 0) exit
      
      do ii = 1, imax
        read (107,*,iostat=eof) i
        if (eof < 0) exit
        backspace (107)
        read (107,*,iostat = eof) k, recall(i)%name, recall(i)%typ, recall(i)%filename
        if (eof < 0) exit
        
        if (recall(i)%typ /= 4) then
          do 
            open (108,file = recall(i)%filename)
            read (108,*,iostat=eof) titldum
            if (eof < 0) exit
            read (108,*,iostat=eof) nbyr
            if (eof < 0) exit
            read (108,*,iostat=eof) header
            exit 
          end do
        
        select case (recall(i)%typ)
            
        case (0) !! subdaily
            allocate (recall(i)%hyd_flo(time%step*366,time%nbyr), source = 0.)
            allocate (recall(i)%hd(366,time%nbyr))
            
          case (1) !! daily
            allocate (recall(i)%hd(366,time%nbyr))
            
          case (2) !! monthly
            allocate (recall(i)%hd(12,time%nbyr))
            
          case (3) !! annual
            allocate (recall(i)%hd(1,time%nbyr))

          end select 
        
        !! save starting year of recall data
        read (108,*,iostat=eof) jday, mo, day_mo, iyr
        recall(i)%start_yr = iyr
        backspace (108)
        
        !! set start year if recall starts before start of simulation
        if (recall(i)%start_yr <= time%yrc) then
          iyrs = 1
          do
            read (108,*,iostat=eof) jday, mo, day_mo, iyr
            if (iyr == time%yrc)  then
              exit
            end if
          end do
          backspace (108)
        else
          !! seet star year if recall starts after start of  simulation
          iyrs = recall(i)%start_yr - time%yrc + 1
        end if
        
        !! read and store data
        do 
          iyr1 = iyr
          read (108,*,iostat=eof) jday1, mo1, day_mo, iyr
          if (eof < 0) exit
          if (iyr > time%yrc_end) exit
          backspace (108)
          
            !! increment iyrs (sequential year of recall data) if next year
            if (iyr1 /= iyr) then
              iyrs = iyrs + 1
            end if
            iyr1 = iyr
          
          !! read data for each time step
          select case (recall(i)%typ)
            case (0) !! subdaily
                 
              !! convert m3/s -> m3
              recall(i)%hyd_flo(istep,iyrs) = ht1%flo * 86400. / time%step
              
              !! reset daily step and sum the daily hyd
              if (istep > idaystep * time%step) then
                !! convert daily flow m3/s -> m3 -- other subdaily inputs are in t and kg
                !recall(i)%hd(idaystep,iyrs)%flo = recall(i)%hd(idaystep,iyrs)%flo * 86400. / time%step
                idaystep = idaystep + 1
                recall(i)%hd(idaystep,iyrs) = recall(i)%hd(idaystep,iyrs) + ht1
              else
                recall(i)%hd(idaystep,iyrs) = recall(i)%hd(idaystep,iyrs) + ht1
              end if
           
            case (1) !! daily
              read (108,*,iostat=eof) jday, mo, day_mo, iyr, ob_typ, ob_name,    &
                                                      recall(i)%hd(jday1,iyrs)
            case (2) !! monthly
              read (108,*,iostat=eof) jday, mo, day_mo, iyr, ob_typ, ob_name,    &
                                                      recall(i)%hd(mo1,iyrs)
            case (3) !! annual
              read (108,*,iostat=eof) jday, mo, day_mo, iyr, ob_typ, ob_name, ht1
              recall(i)%hd(1,iyrs) = ht1
            end select
            
        end do
        
        !! save end year of recall data
        recall(i)%end_yr = iyr
        close (108)
          
      else
          
     if (recall(i)%typ == 4) then
        iexo_allo = 1
        allocate (recall(i)%hd(1,1))
        !! xwalk with exco file to get sequential number
        do iexco_om = 1, db_mx%exco_om
          if (exco_db(iexco_om)%name == recall(i)%filename) then
            recall(i)%hd(1,1) = exco(iexco_om)
            exit
          end if
        end do
     end if
     
    end if
      
    end do
      close (107)
      exit
      enddo
      endif
      
      !read all rec_pest files
      inquire (file="pest.com", exist=i_exist)
      if (i_exist ) then
      do
        open (107,file="pest.com")
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
        imax = 0
          do while (eof == 0)
            read (107,*,iostat=eof) i
            if (eof < 0) exit
            imax = Max(imax,i) 
          end do
          
      allocate (rec_pest(0:imax))
      rewind (107)
      read (107,*,iostat=eof) titldum
      if (eof < 0) exit
      read (107,*,iostat=eof) header
      if (eof < 0) exit
      
      do ipc = 1, db_mx%pestcom
        read (107,*,iostat=eof) ipestcom_db   !pointer to pestcom_db - fix***
        if (eof < 0) exit

         do ii = 1, imax
           read (107,*,iostat=eof) i
           if (eof < 0) exit
           backspace (107)
           read (107,*,iostat = eof) k, rec_pest(i)%name, rec_pest(i)%typ, rec_pest(i)%filename
           if (eof < 0) exit
           open (108,file = rec_pest(i)%filename)
           read (108,*,iostat=eof) titldum
           if (eof < 0) exit
           read (108,*,iostat=eof) nbyr
           if (eof < 0) exit
           read (108,*,iostat=eof) header
           if (eof < 0) exit
        
        select case (rec_pest(i)%typ)
           case (1) !! daily
            allocate (rec_pest(i)%hd_pest(366,nbyr))
            
           case (2) !! monthly
            allocate (rec_pest(i)%hd_pest(12,nbyr))
            
           case (3) !! annual
            allocate (rec_pest(i)%hd_pest(1,nbyr))
        end select
           
        ! read and store entire year
       do 
         read (108,*,iostat=eof) iyr, istep
         if (eof < 0) exit
         if (iyr == time%yrc) exit
       end do
       
       backspace (108)
       iyr_prev = iyr
       iyrs = 1
       
       do 
         read (108,*,iostat=eof) iyr, istep, recall(i)%hd(istep,iyrs)
         if (eof < 0) exit
         !call hyd_convert_mass (recall(i)%hd(istep,iyrs))
         if (iyr /= iyr_prev) then
           iyr_prev = iyr
           iyrs = iyrs + 1
         endif
       end do
       close (108)
         end do 
        end do
        close (107)
      end do
      end if     
      
      return
      end subroutine recall_read