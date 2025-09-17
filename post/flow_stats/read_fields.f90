subroutine read_and_stats(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir,namefile
character(len=8) :: numfile
character(len=3) :: setnum
logical :: check


namedir='../../multi/output/'
write(numfile,'(i8.8)') nstep


allocate(u(nx,ny,nz))
allocate(v(nx,ny,nz))
allocate(w(nx,ny,nz))
allocate(phi(nx,ny,nz))

write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow'
if (uflag .eq. 1) then
  !reading u
  namefile=trim(namedir)//'u_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) u
  close(666,status='keep')
  !reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  !reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')
endif
if (phiflag .eq. 1) then
  !reading phi
  namefile=trim(namedir)//'phi_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) phi
  close(668,status='keep')
endif

!--- compute statistics ---
allocate(mean(nz,3),rms(nz,3),skw(nz,3),flt(nz,3))

mean=0.0d0
rms=0.0d0
skw=0.0d0
flt=0.0d0

! mean
do k=1,nz
  do i=1,nx
    do j=1,ny
      mean(k,1) = mean(k,1) + u(i,j,k)
      mean(k,2) = mean(k,2) + v(i,j,k)
      mean(k,3) = mean(k,3) + w(i,j,k)
    enddo
  enddo
enddo
mean=mean/(real(nx*ny))


! rms
do k=1,nz
  do i=1,nx
    do j=1,ny
      rms(k,1) = (u(i,j,k)-mean(k,1))**2
      rms(k,2) = (v(i,j,k)-mean(k,2))**2
      rms(k,3) = (w(i,j,k)-mean(k,3))**2
    enddo
  enddo
enddo
rms=rms/(real(nx*ny))
rms=sqrt(rms)

! to be implemented
! flt


! skw


! write output in a file (WIP)
! this part is copied from FLOW36, must be adapted

!open(66,status='replace',file='./output/statistics.dat',form='formatted')

!write(66,'(a,i8,a,3(i5),a)') 'Statistics gathered on ',counter,' flow fields, on a ',nx,ny,nz,' grid (nx,ny,nz)'
!write(66,'(13(a12,2x))') 'z','u mean','v mean','w mean','u rms','v rms','w rms','u skw','v skw','w skw','u flt','v flt','w flt'
!write(66,*)

!do i=1,nz
! write(66,'(f12.5,2x,12(es12.5,2x))') stats(i,1:13)
!enddo

deallocate(u,v,w,phi)
deallocate(mean,rms,skw,flt)

return
end
