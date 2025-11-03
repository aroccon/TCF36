subroutine read_input
 use param 
 implicit none
 integer :: ios


nx=512
ny=256
nz=384
begin=10000
finish=200000
step=2000
rootpath='../../multi/output/'
  
 ! Example for 'in_post.txt'
 !512             ! Line 1: nx
 !0               ! Line 2: Begin (first saving to be processed)
 !140000             ! Line 3: End (last saving to be processed)
 !2000               ! Line 4: Step (step for processing: 1= every saving, 2=every other saving, ...)
 !"/leonardo_work/IscrB_SONORA/RUN1/multi/output"              !Line 4 : rootpath for results

 close(1)

 return
end subroutine
!*********************************************
subroutine print_start()
 use param
 implicit none
 
 write(*,'(a)') repeat('-', 70)
 write(*,'(a)') '           Count Droplets and Compute Their Size'
 write(*,'(a)') repeat('-', 70)
 
 write(*,'(a,i6)') 'Grid size      : nx = ', nx
 write(*,'(a,f10.4)') 'Domain size  :  lx = ', lx
 write(*,'(a,i6)') 'Grid size      : ny = ', ny
 write(*,'(a,f10.4)') 'Domain size  :  ly = ', ly
 write(*,'(a,i6)') 'Grid size      : nz = ', nz
 write(*,'(a,f10.4)') 'Domain size  :  lz = ', lz

 write(*,'(a,i8,a,i8,a,i6)') 'Time steps     : from ', begin, ' to ', finish, ' with step ', step
 
 write(*,'(a)') repeat('-', 70)
 write(*,*)
 
 return
end subroutine print_start
!*********************************************
