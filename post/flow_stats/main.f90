program flow_stats
use commondata
implicit none

integer :: i,j,k
logical :: check
character(len=40) :: namefile

!! read input
open(10,file='input_par.inp',form='formatted')
read(10,*) nx
read(10,*) ny
read(10,*) nz
read(10,*) nstart
read(10,*) nend
read(10,*) dump
read(10,*) uflag
read(10,*) phiflag
read(10,*) lx
read(10,*) ly
read(10,*) lz

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

dx=lx/real(nx)
dy=ly/real(ny)
dz=lz/real(nz)

x=0.0d0
y=0.0d0
z=0.0d0
z(1)=dz/2

do i=1,nx-1
  x(i+1)=x(i)+dx
enddo
do j=1,ny-1
  y(j+1)=y(j)+dy
enddo
do k=1,nz-1
  z(k+1)=z(k)+dz
enddo

! read fluid data
do i=nstart,nend,dump
 call read_and_stats(i)
enddo

deallocate(x,y,z)

end program flow_stats
