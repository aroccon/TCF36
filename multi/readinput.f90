
!##########################################################################
!###########################################################################
subroutine readinput
use velocity
use phase
use param
use mpivar
implicit none

open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) restart
read(55,*) tstart
read(55,*) tfin
read(55,*) dump
!Flow parameters
read(55,*) inflow
read(55,*) inphi
read(55,*) dt
read(55,*) mu
read(55,*) rho
! forcing parameters
read(55,*) gradpx
read(55,*) gradpy
! phase-field parameters
read(55,*) radius
read(55,*) sigma
read(55,*) epsr   


! compute pre-defined constant 
twopi=8.0_8*atan(1.0_8)
lx=2.d0!*twopi
ly=2.d0!twopi
lz=2.d0
dx = lx/nx
dy = ly/ny
dz = lz/nz
dxi = 1.d0/dx
dyi = 1.d0/dy
dzi = 1.d0/dz
ddxi = 1.d0/dx/dx
ddyi = 1.d0/dy/dy
ddzi = 1.d0/dz/dz
rhoi=1.d0/rho
eps=epsr*dx
epsi=1.d0/eps
enum=1.e-16

if (rank .eq. 0) then
    !enable/disable for debug check parameters
    
    write(*,*) "------------------------------------------"
    write(*,*) "████████  ██████ ███████ ██████   ██████  "  
    write(*,*) "   ██    ██      ██           ██ ██       "       
    write(*,*) "   ██    ██      █████    █████  ███████  "  
    write(*,*) "   ██    ██      ██           ██ ██    ██ " 
    write(*,*) "   ██     ██████ ██      ██████   ██████  "
    write(*,*) "------------------------------------------"
    write(*,*) 'Grid:', nx, 'x', nx, 'x', nx
    write(*,*) "Restart ", restart
    write(*,*) "Tstart  ", tstart
    write(*,*) "Tfin    ", tfin
    write(*,*) "Dump    ", dump
    write(*,*) "Inflow  ", inflow
    write(*,*) "Deltat  ", dt
    write(*,*) "Mu      ", mu
    write(*,*) "Rho     ", rho
    write(*,*) "Gradpx  ", gradpx
    write(*,*) "Gradpy  ", gradpy
    write(*,*) "Radius  ", radius
    write(*,*) "Sigma   ", sigma
    write(*,*) "Eps     ", eps
    write(*,*) "Epsi    ", epsi
    write(*,*) "Lx      ", lx
    write(*,*) "dx", dx
    write(*,*) "dxi", dxi
    write(*,*) "ddxi", ddxi
    write(*,*) "dy", dx
    write(*,*) "dyi", dyi
    write(*,*) "ddyi", ddyi
    write(*,*) "dz", dz
    write(*,*) "dzi", dzi
    write(*,*) "ddzi", ddzi
    write(*,*) "rhoi", rhoi
endif

end subroutine


