module param
    integer, parameter :: nx=128
    integer, parameter :: ny=128
    integer, parameter :: nz=128
    double precision :: pi,lx,dx,dxi,ddxi,rhoi,twopi,dy,dyi,ddyi,dz,dzi,ddzi
    integer :: restart,tstart,tfin,dump
    double precision :: gamma, normod
    double precision :: dt,mu,rho !flow parameters
    integer :: inflow, inphi
    double precision :: f1,f2,f3,k0 ! forcing parameters
    double precision :: radius, sigma, epsr, eps, pos, val, epsi, enum ! phase-field parameters
    double precision :: times,timef
end module param


module mpivar
   ! MPI variables
   integer :: rank, ranks, ierr
   integer :: localRank, localComm
end module mpivar


module cudecompvar
   use cudecomp
   integer :: npx, npy, npz
   type(cudecompHandle) :: handle
   type(cudecompGridDesc) :: grid_desc
   type(cudecompGridDescConfig) :: config
   type(cudecompGridDescAutotuneOptions) :: options
   integer :: pdims(2) ! pr x pc pencils
   integer :: gdims(3) ! global grid dimensions
   integer :: halo(3) ! halo extensions
   integer :: halo_ext ! 0 no halo, 1 means 1 halo
   type(cudecompPencilInfo) :: piX, piY, piZ  ! size of the pencils in x- y- and z-configuration
   integer(8) :: nElemX, nElemY, nElemZ, nElemWork, nElemWork_halo
   logical :: halo_periods(3)
end module cudecompvar


module velocity
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
   double precision, allocatable :: rhsu_o(:,:,:), rhsv_o(:,:,:), rhsw_o(:,:,:)
   complex(8), allocatable :: rhsp_complex(:,:,:)
   double precision, allocatable :: rhsp(:,:,:), p(:,:,:)
   double precision, allocatable :: div(:,:,:)
   double precision :: uc, vc, wc, umax, gumax=1.0d0, cou, alpha, beta
   double precision :: h11, h12, h13, h21, h22, h23, h31, h32, h33
   double precision :: umean, vmean, wmean, gumean, gvmean, gwmean
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:), psidi(:,:,:), rhsphi_o(:,:,:)
   double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
   double precision, allocatable :: chempot(:,:,:)
   double precision, allocatable :: fxst(:,:,:), fyst(:,:,:), fzst(:,:,:)
end module phase




