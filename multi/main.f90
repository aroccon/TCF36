#define CHECK_CUDECOMP_EXIT(f) if (f /= CUDECOMP_RESULT_SUCCESS) call exit(1)

program main
use cudafor
use cudecomp
use cufft
use mpi
use velocity 
use phase
use param
use mpivar
use cudecompvar
use nvtx


implicit none
! timer for scaling test
real :: t_start, t_end, elapsed
! grid dimensions
integer :: comm_backend
integer :: pr, pc
! cudecomp
! cuFFT
integer :: planXf, planXb
integer :: planY, planZ
integer :: batchsize
integer :: status
! other variables (wavenumber, grid location)
real(8), allocatable :: x(:), y(:), z(:), kx(:), ky(:)
integer :: i,j,k,il,jl,kl,ig,jg,kg,t
integer :: im,ip,jm,jp,km,kp,last,idx
double complex :: a(nz), b(nz), c(nz), d(nz), sol(nz)
real(8), device, allocatable :: kx_d(:), ky_d(:)
! working arrays
complex(8), allocatable :: psi(:)
real(8), allocatable :: ua(:,:,:)
real(8), allocatable :: uaa(:,:,:)
real(8), allocatable :: psi_real(:)
! real(8), device, allocatable :: psi_real_d(:)
complex(8), device, allocatable :: psi_d(:)
complex(8), pointer, device, contiguous :: work_d(:), work_halo_d(:), work_d_d2z(:), work_halo_d_d2z(:)
character(len=40) :: namefile
character(len=4) :: itcount
! Code variables
real(8)::err,maxErr
complex(8), device, pointer :: psi3d(:,:,:)
real(8) :: k2
!integer :: il, jl, ig, jg
integer :: offsets(3), xoff, yoff
integer :: np(3)

! Enable or disable phase field (acceleration eneabled by default)
#define phiflag 1

!########################################################################################################################################
! 1. INITIALIZATION OF MPI AND cuDECOMP AUTOTUNING : START
!########################################################################################################################################
! MPI initialization, put in rank the local MPI rank number and ranks total number
! Same procedura defined in the cuDecomp documentation
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, ranks, ierr)

call mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localComm, ierr)
call mpi_comm_rank(localComm, localRank, ierr)
ierr = cudaSetDevice(localRank) !assign GPU to MPI rank

! Define grid and decomposition
call readinput

! hard coded
pr = 0
pc = 0
halo_ext=1
! comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
comm_backend = 0 ! Enable full autotuning

CHECK_CUDECOMP_EXIT(cudecompInit(handle, MPI_COMM_WORLD))
! config is a struct and pr and pc are the number of pencils along the two directions
! gdims is the global grid
! create an uninitialized configuration struct and initialize it to defaults using cudecompGridDescConfigSetDefaults. 
! Initializing to default values is required to ensure no entries are left uninitialized.
CHECK_CUDECOMP_EXIT(cudecompGridDescConfigSetDefaults(config))
pdims = [pr, pc] !pr and pc are the number of pencil along the different directions
config%pdims = pdims
! gdims = [nx, ny, nz]
! config%gdims = gdims
halo = [0, halo_ext, halo_ext] ! no halo along x neeed because is periodic and in physical space i have x-pencil
! for transpositions
config%transpose_comm_backend = comm_backend
config%transpose_axis_contiguous = .true.
! for halo exchanges
config%halo_comm_backend = CUDECOMP_HALO_COMM_MPI
! Setting for periodic halos in all directions (non required to be in config)
halo_periods = [.true., .true., .true.]

! create spectral grid descriptor first to select pdims for optimal transposes
gdims = [nx/2+1, ny, nz]
config%gdims = gdims

! Set up autotuning options for spectral grid (transpose related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_transpose_backend = .true.
   options%autotune_halo_backend = .false.
endif
options%transpose_use_inplace_buffers = .true.
options%transpose_input_halo_extents(:, 1) = halo
options%transpose_output_halo_extents(:, 4) = halo

CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_descD2Z, config, options))

! create physical grid descriptor
! take previous config and modify the global grid (nx instead of nx/2+1)
! reset transpose_comm_backend to default value to avoid picking up possible nvshmem
! transpose backend selection (this impacts how workspaces are allocated)
gdims = [nx, ny, nz]
config%gdims = gdims
config%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P

! Set up autotuning options for physical grid (halo related settings)
CHECK_CUDECOMP_EXIT(cudecompGridDescAutotuneOptionsSetDefaults(options))
options%dtype = CUDECOMP_DOUBLE_COMPLEX
if (comm_backend == 0) then
   options%autotune_halo_backend = .true.
endif
options%halo_extents(:) = halo
options%halo_periods(:) = halo_periods
options%halo_axis = 1
CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_desc, config, options))

! Print information on configuration
!if (rank == 0) then
!   write(*,"(' Running on ', i0, ' x ', i0, ' process grid ...')") config%pdims(1), config%pdims(2)
!   write(*,"(' Using ', a, ' transpose backend ...')") &
!            cudecompTransposeCommBackendToString(config%transpose_comm_backend)
!   write(*,"(' Using ', a, ' halo backend ...')") &
!            cudecompHaloCommBackendToString(config%halo_comm_backend)
!endif


! Get pencil info for the grid descriptor in the physical space
! This function returns a pencil struct (piX, piY or piZ) that contains the shape, global lower and upper index bounds (lo and hi), 
! size of the pencil, and an order array to indicate the memory layout that will be used (to handle permuted, axis-contiguous layouts).
! Additionally, there is a halo_extents data member that indicates the depth of halos for the pencil, by axis.
! Side note:  ! cudecompGetPencilInfo(handle, grid_desc, pinfo_x, 1, [1, 1, 1]) <- in this way the x-pencil also have halo elements
! If no halo regions are necessary, a NULL pointer can be provided in place of this array (or omitted)
! Pencil info in x-configuration present in PiX (shape,lo,hi,halo_extents,size)
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piX, 1, halo))
nElemX = piX%size !<- number of total elments in x-configuratiion (including halo)
! Pencil info in Y-configuration present in PiY
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piY, 2))
nElemY = piY%size
! Pencil info in Z-configuration present in PiZ
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_desc, piZ, 3))
nElemZ = piZ%size

! Get workspace sizes for transpose (1st row, not used) and halo (2nd row, used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_desc, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_desc, 1, halo, nElemWork_halo))






! Get pencil info for the grid descriptor in the complex space 
!gdims = [nx/2+1, ny, nz]
!config%gdims = gdims
!CHECK_CUDECOMP_EXIT(cudecompGridDescCreate(handle, grid_descD2Z, config, options))

CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piX_d2z, 1,halo))
nElemX_d2z = piX_d2z%size !<- number of total elments in x-configuratiion (include halo)
! Pencil info in Y-configuration present in PiY
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piY_d2z, 2))
nElemY_d2z = piY_d2z%size
! Pencil info in Z-configuration present in PiZ
CHECK_CUDECOMP_EXIT(cudecompGetPencilInfo(handle, grid_descD2Z, piZ_d2z, 3))
nElemZ_d2z = piZ_d2z%size
! Get workspace sizes for transpose (1st row,used) and halo (2nd row, not used)
CHECK_CUDECOMP_EXIT(cudecompGetTransposeWorkspaceSize(handle, grid_descD2Z, nElemWork_d2z))
CHECK_CUDECOMP_EXIT(cudecompGetHaloWorkspaceSize(handle, grid_descD2Z, 1, halo, nElemWork_halo_d2z))


!write(*,*) "piZ_d2z_lo", piZ_d2z%lo




! CUFFT initialization -- Create Plans
! Forward 1D FFT in X: D2Z
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXf, nx, CUFFT_D2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Forward'

! Backward 1D FFT in X: Z2D
batchSize = piX_d2z%shape(2)*piX_d2z%shape(3) !<- number of FFT (from x-pencil dimension)
status = cufftPlan1D(planXb, nx, CUFFT_Z2D, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating X plan Backward'

! it's always 2 and 3 because y-pencil have coordinates y,z,x
batchSize = piY_d2z%shape(2)*piY_d2z%shape(3)
status = cufftPlan1D(planY, ny, CUFFT_Z2Z, batchSize)
if (status /= CUFFT_SUCCESS) write(*,*) rank, ': Error in creating Y plan Forward & Backward'

! Z-plan removed (not neeeded)

! define grid
allocate(x(nx),y(ny),z(nz),kx(nx),ky(ny))
! location of the pressure nodes (cell centers)
x(1)=dx/2
do i = 2, nx
   x(i) = x(i-1) + dx
enddo
y(1)=dy/2
do i = 2, ny
   y(i) = y(i-1) + dy
enddo
z(1)=dz/2
do i = 2, nz
   z(i) = z(i-1) + dz
enddo
do i = 1, nx/2
   kx(i) = (i-1)*(twopi/lx)
enddo
do i = nx/2+1, nx
   kx(i) = (i-1-nx)*(twopi/lx)
enddo
do j = 1, ny/2
   ky(j) = (j-1)*(twopi/ly)
enddo
do j = ny/2+1, ny
   ky(j) = (j-1-ny)*(twopi/ly)
enddo
! allocate kx_d and ky_d on the device 
allocate(kx_d, source=kx)
allocate(ky_d, source=ky)

!########################################################################################################################################
! 1. INITIALIZATION AND cuDECOMP AUTOTUNING : END
!########################################################################################################################################





!########################################################################################################################################
! START STEP 2: ALLOCATE ARRAYS
!########################################################################################################################################
! allocate arrays
allocate(psi(max(nElemX, nElemY, nElemZ))) !largest among the pencil
allocate(psi_real(max(nElemX, nElemY, nElemZ))) !largest among the pencil
allocate(psi_d(max(nElemX_d2z, nElemY_d2z, nElemZ_d2z))) ! phi on device
allocate(ua(nx, piX%shape(2), piX%shape(3)))

! Pressure variable
allocate(rhsp(piX%shape(1), piX%shape(2), piX%shape(3))) 
allocate(p(piX%shape(1), piX%shape(2), piX%shape(3))) 

!allocate variables
!NS variables
allocate(u(piX%shape(1),piX%shape(2),piX%shape(3)),v(piX%shape(1),piX%shape(2),piX%shape(3)),w(piX%shape(1),piX%shape(2),piX%shape(3))) !velocity vector
! allocate(ustar(piX%shape(1),piX%shape(2),piX%shape(3)),vstar(piX%shape(1),piX%shape(2),piX%shape(3)),wstar(piX%shape(1),piX%shape(2),piX%shape(3))) ! provisional velocity field
allocate(rhsu(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(rhsu_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsv_o(piX%shape(1),piX%shape(2),piX%shape(3)),rhsw_o(piX%shape(1),piX%shape(2),piX%shape(3))) ! right hand side u,v,w
allocate(div(piX%shape(1),piX%shape(2),piX%shape(3)))
!PFM variables
#if phiflag == 1
allocate(phi(piX%shape(1),piX%shape(2),piX%shape(3)),rhsphi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(psidi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(tanh_psi(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(normx(piX%shape(1),piX%shape(2),piX%shape(3)),normy(piX%shape(1),piX%shape(2),piX%shape(3)),normz(piX%shape(1),piX%shape(2),piX%shape(3)))
allocate(fxst(piX%shape(1),piX%shape(2),piX%shape(3)),fyst(piX%shape(1),piX%shape(2),piX%shape(3)),fzst(piX%shape(1),piX%shape(2),piX%shape(3))) ! surface tension forces
#endif

! allocate arrays for transpositions and halo exchanges 
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_d, nElemWork))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_desc, work_halo_d, nElemWork_halo))
! allocate arrays for transpositions
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_d_d2z, nElemWork_d2z))
CHECK_CUDECOMP_EXIT(cudecompMalloc(handle, grid_descD2Z, work_halo_d_d2z, nElemWork_halo_d2z))
!########################################################################################################################################
! END STEP2: ALLOCATE ARRAYS
!########################################################################################################################################








!########################################################################################################################################
! START STEP 3: FLOW AND PHASE FIELD INIT
!########################################################################################################################################
! 3.1 Read/initialize from data without halo grid points (avoid out-of-bound if reading usin MPI I/O)
! 3.2 Call halo exchnages along Y and Z for u,v,w and phi
if (restart .eq. 0) then !fresh start Taylor Green or read from file in init folder
if (rank.eq.0) write(*,*) "Initialize velocity field (fresh start)"
   if (inflow .eq. 0) then
      if (rank.eq.0) write(*,*) "Initialize Taylor-green"
      do k = 1+halo_ext, piX%shape(3)-halo_ext
         kg = piX%lo(3) + k - 1 
         do j = 1+halo_ext, piX%shape(2)-halo_ext
            jg = piX%lo(2) + j - 1 
            do i = 1, piX%shape(1)
               call random_number(noise)
               u(i,j,k) =  10.d0 + 2.d0*sin(twopi/lx*x(i))*cos(twopi/ly*y(jg))*(1-z(kg)*z(kg))
               v(i,j,k) =  0.d0  - ly/lx*2.d0*cos(twopi/lx*x(i))*sin(twopi/ly*y(jg))*(1-z(kg)*z(kg))
               w(i,j,k) =  0.d0  
            enddo
         enddo
      enddo
   endif
   if (inflow .eq. 1) then
   if (rank.eq.0)  write(*,*) "Initialize from data"
         call readfield(1)
         call readfield(2)
         call readfield(3)
      endif
   endif
if (restart .eq. 1) then !restart, ignore inflow and read the tstart field 
   if (rank.eq.0)  write(*,*) "Initialize velocity field (from output folder), iteration:", tstart
   call readfield_restart(tstart,1)
   call readfield_restart(tstart,2)
   call readfield_restart(tstart,3)
endif


! update halo cells along y and z directions (enough only if pr and pc are non-unitary)
!$acc host_data use_device(u)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(v)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
!$acc host_data use_device(w)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 

! initialize phase-field
#if phiflag == 1
if (restart .eq. 0) then
if (rank.eq.0) write(*,*) 'Initialize phase field (fresh start)'
   if (inphi .eq. 0) then
   if (rank.eq.0) write(*,*) 'Spherical drop'
      do k = 1+halo_ext, piX%shape(3)-halo_ext
      kg = piX%lo(3) + k - 1 
         do j = 1+halo_ext, piX%shape(2)-halo_ext
         jg = piX%lo(2) + j - 1 
            do i = 1, piX%shape(1)
                pos=(x(i)-lx/2)**2d0 +  (y(jg)-ly/2)**2d0 + (z(kg)-lz/2)**2d0
                phi(i,j,k) = 0.5d0*(1.d0-tanh((sqrt(pos)-radius)/2/eps))
            enddo
        enddo
    enddo
   endif
   if (inphi .eq. 1) then
      if (rank.eq.0)  write(*,*) "Initialize phase-field from data"
      call readfield(5)
   endif
endif
if (restart .eq. 1) then
    write(*,*) "Initialize phase-field (restart, from output folder), iteration:", tstart
    call readfield_restart(tstart,5)
endif
! update halo
!$acc host_data use_device(phi)
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
!$acc end host_data 
#endif

!Save initial fields (only if a fresh start)
if (restart .eq. 0) then
   if (rank.eq.0) write(*,*) "Save initial fields"
   call writefield(tstart,1)
   call writefield(tstart,2)
   call writefield(tstart,3)
   call writefield(tstart,4)
   #if phiflag == 1
   call writefield(tstart,5)
   #endif
endif
!########################################################################################################################################
! END STEP 3: FLOW AND PHASE FIELD INIT
!########################################################################################################################################






! ########################################################################################################################################
! START TEMPORAL LOOP: STEP 4 to 9 REPEATED AT EVERY TIME STEP
! ########################################################################################################################################
! First step use Euler
alpha=1.0d0
beta=0.0d0
gumax=1.d0
tstart=tstart+1
gamma=1.d0*gumax
!$acc data copyin(piX)
!$acc data create(rhsu_o, rhsv_o, rhsw_o)
!$acc data copyin(mysin, mycos)
call cpu_time(t_start)
! Start temporal loop
do t=tstart,tfin
    ! Create custom label for each marker
    write(itcount,'(i4)') t
    ! Range with custom  color
    call nvtxStartRange("Iteration "//itcount,t)

    if (rank.eq.0) write(*,*) "Time step",t,"of",tfin
    call cpu_time(times)

   call nvtxStartRange("Phase-field")
   !########################################################################################################################################
   ! START STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   #if phiflag == 1
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            ! compute distance function psi (used to compute normals)
            val = min(phi(i,j,k),1.0d0) ! avoid machine precision overshoots in phi that leads to problem with log
            psidi(i,j,k) = eps*log((val+enum)/(1.d0-val+enum))
            ! compute here the tanh of distance function psi (used in the sharpening term) to avoid multiple computations of tanh
            tanh_psi(i,j,k) = tanh(0.5d0*psidi(i,j,k)*epsi)
         enddo
      enddo
   enddo
   !$acc end kernels

   gamma=1.d0*gumax
   !$acc parallel loop tile(16,4,2)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ! 4.1 RHS computation
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            ! convective (first three lines) and diffusive (last three lines)
            rhsphi(i,j,k) =   &
                  - (u(ip,j,k)*0.5d0*(phi(ip,j,k)+phi(i,j,k)) - u(i,j,k)*0.5d0*(phi(i,j,k)+phi(im,j,k)))*dxi  &  
                  - (v(i,jp,k)*0.5d0*(phi(i,jp,k)+phi(i,j,k)) - v(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,jm,k)))*dxi  &  
                  - (w(i,j,kp)*0.5d0*(phi(i,j,kp)+phi(i,j,k)) - w(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,j,km)))*dxi  &  
                        + gamma*(eps*(phi(ip,j,k)-2.d0*phi(i,j,k)+phi(im,j,k))*ddxi + &                   
                                 eps*(phi(i,jp,k)-2.d0*phi(i,j,k)+phi(i,jm,k))*ddxi + &                   
                                 eps*(phi(i,j,kp)-2.d0*phi(i,j,k)+phi(i,j,km))*ddxi)                      
            ! 4.1.3. Compute normals for sharpening term (gradient)
            normx(i,j,k) = (psidi(ip,j,k) - psidi(im,j,k))
            normy(i,j,k) = (psidi(i,jp,k) - psidi(i,jm,k))
            normz(i,j,k) = (psidi(i,j,kp) - psidi(i,j,km))
         enddo
      enddo
   enddo

   ! Update normx,normy and normz halos, required to then compute normal derivative
   !$acc host_data use_device(normx)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normx, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(normy)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normy, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(normz)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, normz, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   ! 4.1.3. Compute Sharpening term (gradient)
   ! Substep 2: Compute normals (1.e-16 is a numerical tollerance to avodi 0/0)
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            normod = 1.d0/(sqrt(normx(i,j,k)*normx(i,j,k) + normy(i,j,k)*normy(i,j,k) + normz(i,j,k)*normz(i,j,k)) + 1.0E-16)
            ! normod = 1.d0/(sqrt(normx(i,j,k)**2d0 + normy(i,j,k)**2d0 + normz(i,j,k)**2d0) + 1.0E-16)
            normx(i,j,k) = normx(i,j,k)*normod
            normy(i,j,k) = normy(i,j,k)*normod
            normz(i,j,k) = normz(i,j,k)*normod
         enddo
      enddo
   enddo
   !$acc end kernels

   ! Compute sharpening term
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
               ip=i+1
               jp=j+1
               kp=k+1
               im=i-1
               jm=j-1
               km=k-1
               if (ip .gt. nx) ip=1
               if (im .lt. 1) im=nx
               rhsphi(i,j,k)=rhsphi(i,j,k)-gamma*((0.25d0*(1.d0-tanh_psi(ip,j,k)*tanh_psi(ip,j,k))*normx(ip,j,k) - &
                                                      0.25d0*(1.d0-tanh_psi(im,j,k)*tanh_psi(im,j,k))*normx(im,j,k))*0.5*dxi + &
                                                     (0.25d0*(1.d0-tanh_psi(i,jp,k)*tanh_psi(i,jp,k))*normy(i,jp,k) - &
                                                      0.25d0*(1.d0-tanh_psi(i,jm,k)*tanh_psi(i,jm,k))*normy(i,jm,k))*0.5*dxi + &
                                                     (0.25d0*(1.d0-tanh_psi(i,j,kp)*tanh_psi(i,j,kp))*normz(i,j,kp) - &
                                                      0.25d0*(1.d0-tanh_psi(i,j,km)*tanh_psi(i,j,km))*normz(i,j,km))*0.5*dxi)
            enddo
        enddo
    enddo
    !$acc end kernels

   ! 4.2 Get phi at n+1 
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
            do i=1,nx
                phi(i,j,k) = phi(i,j,k) + dt*rhsphi(i,j,k)
            enddo
        enddo
    enddo
   !$acc end kernels

   ! 4.3 Call halo exchnages along Y and Z for phi 
   !$acc host_data use_device(phi)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, phi, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   #endif
   !########################################################################################################################################
   ! END STEP 4: PHASE-FIELD SOLVER (EXPLICIT)
   !########################################################################################################################################
   call nvtxEndRange






   call nvtxStartRange("Projection")
   !########################################################################################################################################
   ! START STEP 5: USTAR COMPUTATION (PROJECTION STEP)
   !########################################################################################################################################
   ! 5.1 compute rhs 
   ! 5.2 obtain ustar and store old rhs in rhs_o
   ! 5.3 Call halo exchnages along Y and Z for u,v,w

   ! Projection step, convective terms
   ! 5.1a Convective terms NS
   ! Loop on inner nodes
   !$acc parallel loop tile(16,4,2) 
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            ! Manual periodicity ony along x (x-pencil), along y and z directions use halos
            if (ip .gt. nx) ip=1  
            if (im .lt. 1) im=nx
            ! compute the products (conservative form)
            h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k))
            h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h13 = (u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k))
            h21 = (u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
            h22 = (v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k))
            h23 = (w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km))
            h31 = (w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km))
            h32 = (v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k))
            h33 = (w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km))
            ! compute the derivative
            h11=h11*0.25d0*dxi
            h12=h12*0.25d0*dyi
            h13=h13*0.25d0*dzi
            h21=h21*0.25d0*dxi
            h22=h22*0.25d0*dyi
            h23=h23*0.25d0*dzi
            h31=h31*0.25d0*dxi
            h32=h32*0.25d0*dyi
            h33=h33*0.25d0*dzi
            ! add to the rhs
            rhsu(i,j,k)=-(h11+h12+h13)
            rhsv(i,j,k)=-(h21+h22+h23)
            rhsw(i,j,k)=-(h31+h32+h33)
            ! viscous term
            h11 = mu*(u(ip,j,k)-2.d0*u(i,j,k)+u(im,j,k))*ddxi
            h12 = mu*(u(i,jp,k)-2.d0*u(i,j,k)+u(i,jm,k))*ddyi
            h13 = mu*(u(i,j,kp)-2.d0*u(i,j,k)+u(i,j,km))*ddzi
            h21 = mu*(v(ip,j,k)-2.d0*v(i,j,k)+v(im,j,k))*ddxi
            h22 = mu*(v(i,jp,k)-2.d0*v(i,j,k)+v(i,jm,k))*ddyi
            h23 = mu*(v(i,j,kp)-2.d0*v(i,j,k)+v(i,j,km))*ddzi
            h31 = mu*(w(ip,j,k)-2.d0*w(i,j,k)+w(im,j,k))*ddxi
            h32 = mu*(w(i,jp,k)-2.d0*w(i,j,k)+w(i,jm,k))*ddyi
            h33 = mu*(w(i,j,kp)-2.d0*w(i,j,k)+w(i,j,km))*ddzi
            rhsu(i,j,k)=rhsu(i,j,k)+(h11+h12+h13)*rhoi
            rhsv(i,j,k)=rhsv(i,j,k)+(h21+h22+h23)*rhoi
            rhsw(i,j,k)=rhsw(i,j,k)+(h31+h32+h33)*rhoi
            ! Pressure driven
            rhsu(i,j,k)=rhsu(i,j,k) - gradpx
            rhsv(i,j,k)=rhsv(i,j,k) - gradpy
         enddo
      enddo
   enddo

   ! Re-add here ST forces

   ! Surface tension forces
   #if phiflag == 1
   !$acc kernels
   !Obtain surface tension forces evaluated at the center of the cell (same as where phi is located)
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            im=i-1
            jm=j-1
            km=k-1
            if (ip .gt. nx) ip=1
            if (im .lt. 1) im=nx
            chempot=phi(i,j,k)*(1.d0-phi(i,j,k))*(1.d0-2.d0*phi(i,j,k))*epsi-eps*(phi(ip,j,k)+phi(im,j,k)+phi(i,jp,k)+phi(i,jm,k)+phi(i,j,kp)+phi(i,j,km)- 6.d0*phi(i,j,k))*ddxi
            ! chempot*gradphi
            fxst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(ip,j,k)-phi(im,j,k))*dxi
            fyst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dxi
            fzst(i,j,k)=6.d0*sigma*chempot*0.5d0*(phi(i,j,kp)-phi(i,j,km))*dxi
         enddo
      enddo
   enddo
   !$acc end kernels

   ! Update halo of fxst, fyst and fzst (required then to interpolate at velocity points)
   !$acc host_data use_device(fxst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fxst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(fyst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fyst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(fzst)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, fzst, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   
   ! Interpolate force at velocity points
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            im=i-1
            jm=j-1
            km=k-1
            if (im .lt. 1) im=nx
            rhsu(i,j,k)=rhsu(i,j,k) + 0.5d0*(fxst(im,j,k)+fxst(i,j,k))*rhoi
            rhsv(i,j,k)=rhsv(i,j,k) + 0.5d0*(fyst(i,jm,k)+fyst(i,j,k))*rhoi
            rhsw(i,j,k)=rhsw(i,j,k) + 0.5d0*(fzst(i,j,km)+fzst(i,j,k))*rhoi
            u(i,j,k) = u(i,j,k) + dt*(alpha*rhsu(i,j,k)-beta*rhsu_o(i,j,k))
            v(i,j,k) = v(i,j,k) + dt*(alpha*rhsv(i,j,k)-beta*rhsv_o(i,j,k))
            w(i,j,k) = w(i,j,k) + dt*(alpha*rhsw(i,j,k)-beta*rhsw_o(i,j,k))
            rhsu_o(i,j,k)=rhsu(i,j,k)
            rhsv_o(i,j,k)=rhsv(i,j,k)
            rhsw_o(i,j,k)=rhsw(i,j,k)
          enddo
      enddo
   enddo
   !$acc end kernels

   #else
   ! 5.2 find u, v and w star (explicit AB2), only in the inner nodes 
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
          do i=1,nx
              u(i,j,k) = u(i,j,k) + dt*(alpha*rhsu(i,j,k)-beta*rhsu_o(i,j,k))
              v(i,j,k) = v(i,j,k) + dt*(alpha*rhsv(i,j,k)-beta*rhsv_o(i,j,k))
              w(i,j,k) = w(i,j,k) + dt*(alpha*rhsw(i,j,k)-beta*rhsw_o(i,j,k))
              rhsu_o(i,j,k)=rhsu(i,j,k)
              rhsv_o(i,j,k)=rhsv(i,j,k)
              rhsw_o(i,j,k)=rhsw(i,j,k)
          enddo
      enddo
   enddo
   !$acc end kernels
   #endif

   ! store rhs* in rhs*_o 
   ! After first step move to AB2 
   alpha=1.5d0
   beta= 0.5d0
   ! !$acc kernels


   ! 5.3 update halos (y and z directions), required to then compute the RHS of Poisson equation because of staggered grid
   !$acc host_data use_device(u)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(v)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(w)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !########################################################################################################################################
   ! END STEP 5: USTAR COMPUTATION 
   !########################################################################################################################################
   call nvtxEndRange







   call nvtxStartRange("Poisson")
   !########################################################################################################################################
   ! START STEP 6: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################
   ! initialize rhs and analytical solution
   ! 6.1 Compute rhs of Poisson equation div*ustar: divergence at the cell center 
   ! I've done the halo updates so to compute the divergence at the pencil border i have the *star from the halo
   call nvtxStartRange("compute RHS")

   ! impose BC on u,v and w at k=1 and kg=nz (for u and v) and kg=nz+1
   ! can be improved as the bottom one
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            ! bottom wall
            kg = piX%lo(3)  + k -2
            if (kg .eq. 1) then
            u(i,j,k)=0.d0
            v(i,j,k)=0.d0
            w(i,j,k)=0.d0
            endif
            if (kg .eq. nz) then
            u(i,j,k)=0.d0
            v(i,j,k)=0.d0
            endif
            if (kg .eq. nz+1) then
            w(i,j,k)=0.d0
            endif
         enddo
      enddo
   enddo
   !$acc end kernels

   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            ip=i+1
            jp=j+1
            kp=k+1
            if (ip > nx) ip=1
            rhsp(i,j,k) =               (rho*dxi/dt)*(u(ip,j,k)-u(i,j,k))
            rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(v(i,jp,k)-v(i,j,k))
            rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(w(i,j,kp)-w(i,j,k))
            !rhsp(i,j,k) = 0.d0
         enddo
      enddo
   enddo
   !$acc end kernels
   call nvtxEndRange


   call nvtxStartRange("FFT forward w/ transpositions")

   !$acc host_data use_device(rhsp)
   status = cufftExecD2Z(planXf, rhsp, psi_d)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X forward error: ', status
   !$acc end host_data
   ! psi(kx,y,z) -> psi(y,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeXToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX,piX_d2z%halo_extents, [0,0,0]))
   ! psi(y,z,kx) -> psi(ky,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_FORWARD)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y forward error: ', status
   ! psi(ky,z,kx) -> psi(z,kx,ky)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToZ(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX)) 

   call nvtxEndRange


   np(piZ_d2z%order(1)) = piZ_d2z%shape(1)
   np(piZ_d2z%order(2)) = piZ_d2z%shape(2)
   np(piZ_d2z%order(3)) = piZ_d2z%shape(3)
   
   call c_f_pointer(c_devloc(psi_d), psi3d, piZ_d2z%shape)

   offsets(piZ_d2z%order(1)) = piZ_d2z%lo(1) - 1
   offsets(piZ_d2z%order(2)) = piZ_d2z%lo(2) - 1
   offsets(piZ_d2z%order(3)) = piZ_d2z%lo(3) - 1


   xoff = offsets(1)
   yoff = offsets(2)
   npx = np(1)
   npy = np(2)
   !write(*,*) "xoff, yoff", xoff, yoff

   call nvtxStartRange("Solution")

   !$acc parallel loop gang private(a, b, c, d, sol, factor)   
   !!$acc kernels ! kernels is safer but serialized the TDMA, rememeber a-d private if using parallel loop
   do jl = 1, npy
      jg = yoff + jl
      do il = 1, npx
         ig = xoff + il
         ! Set up tridiagonal system for each kx
         ! The system is: (A_j) * pc(kx,j-1) + (B_j) * pc(kx,j) + (C_j) * pc(kx,j+1) = rhs(kx,j)
         ! FD2 in z: -pc(k-1) + 2*pc(k) - pc(k+1)  --> Laplacian in z
         ! Neumann BC: d/dz pc = 0 at k=1 and k=nz
         ! Fill diagonals and rhs for each k
         do k = 1, nz
            a(k) =  1.0d0*dzi*dzi
            b(k) = -2.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
            c(k) =  1.0d0*dzi*dzi
            d(k) =  psi3d(k,il,jl)
         enddo

         ! Neumann BC at k=1 (bottom)
         b(1) = -2.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
         c(1) =  2.0d0*dzi*dzi
         a(1) =  0.0d0

         ! Neumann BC at j=ny (top)
         a(nz) =  2.0d0*dzi*dzi
         b(nz) = -2.0d0*dzi*dzi - kx_d(ig)*kx_d(ig) - ky_d(jg)*ky_d(jg)
         c(nz) =  0.0d0

         ! Special handling for kx=0 and ky=0 (mean mode): ix pressure on one point
         if ((ig .eq. 1) .and. (jg .eq. 1)) then
            b(1) = 1.0d0
            c(1) = 0.0d0
            d(1) = 0.0d0
         endif

         ! Thomas algorithm (TDMA) for tridiagonal system 
         ! Forward sweep
         do k=2,nz
            factor = a(k)/b(k-1)
            b(k) = b(k) - factor * c(k-1)
            d(k) = d(k) - factor * d(k-1)
         enddo

         ! Back substitution
         sol(nz) = d(nz)/b(nz)
         do k = nz-1, 1, -1
            sol(k) = (d(k) - c(k)*sol(k+1))/b(k)
         end do

         ! Store solution in array that do the back FFT
         do k=1,nz
            !sol(k)=cmplx(0.d0,0.d0)
            psi3d(k,il,jl) = sol(k)
         enddo      
      enddo
   enddo
   !!$acc end kernels


   call nvtxStartRange("FFT backwards along x and y w/ transpositions")

   ! psi(z,kx,ky) -> psi(ky,z,kx)
   CHECK_CUDECOMP_EXIT(cudecompTransposeZToY(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX))
   ! psi(ky,z,kx) -> psi(y,z,kx)
   status = cufftExecZ2Z(planY, psi_d, psi_d, CUFFT_INVERSE)
   if (status /= CUFFT_SUCCESS) write(*,*) 'Y inverse error: ', status
   ! psi(y,z,kx) -> psi(kx,y,z)
   CHECK_CUDECOMP_EXIT(cudecompTransposeYToX(handle, grid_descD2Z, psi_d, psi_d, work_d_d2z, CUDECOMP_DOUBLE_COMPLEX,[0,0,0], piX_d2z%halo_extents))
   !$acc host_data use_device(p)
   ! psi(kx,y,z) -> psi(x,y,z)
   status = cufftExecZ2D(planXb, psi_d, p)
   if (status /= CUFFT_SUCCESS) write(*,*) 'X inverse error: ', status
   !$acc end host_data

   ! normalize pressure
   !$acc kernels
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i=1,nx
            p(i,j,k) = p(i,j,k)/nx/ny
         enddo
      enddo
   enddo
   !$acc end kernels
      
   call nvtxEndRange



   ! update halo nodes with pressure 
   ! Update X-pencil halos 
   !$acc host_data use_device(p)
    CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, p, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
    CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, p, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
    !$acc end host_data 

   


   !########################################################################################################################################
   ! END STEP 7: POISSON SOLVER FOR PRESSURE
   !########################################################################################################################################
   call nvtxEndRange







   call nvtxStartRange("Correction")
   !########################################################################################################################################
   ! START STEP 8: VELOCITY CORRECTION
   ! ########################################################################################################################################
   ! 8.1 Correct velocity 
   ! 8.2 Call halo exchnages along Y 
   ! Correct velocity, pressure has also the halo
   !$acc kernels 
   do k=1+halo_ext, piX%shape(3)-halo_ext
      do j=1+halo_ext, piX%shape(2)-halo_ext
         do i = 1, piX%shape(1) ! equal to nx (no halo on x)
              im=i-1
              jm=j-1
              km=k-1
              if (im < 1) im=nx
              u(i,j,k)=u(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi
              v(i,j,k)=v(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dxi
              w(i,j,k)=w(i,j,k) - dt/rho*(p(i,j,k)-p(i,j,km))*dxi
          enddo
      enddo
   enddo
   !$acc end kernels 

   ! 8.3 update halos (y direction), required to then compute the RHS of Poisson equation because of staggered grid
   !$acc host_data use_device(u)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, u, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(v)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, v, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 
   !$acc host_data use_device(w)
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 2))
   CHECK_CUDECOMP_EXIT(cudecompUpdateHalosX(handle, grid_desc, w, work_halo_d, CUDECOMP_DOUBLE, piX%halo_extents, halo_periods, 3))
   !$acc end host_data 

   ! Re-impose BC on u,v and w at k=1 and kg=nz (for u and v) and kg=nz+1
   ! can be improved accessing directly kg?
   !$acc kernels
   do k=1, piX%shape(3)
      do j=1, piX%shape(2)
         do i=1,nx
            ! bottom wall
            kg = piX%lo(3)  + k - 2
            if (kg .eq. 1) then
            u(i,j,k)=0.d0
            v(i,j,k)=0.d0
            w(i,j,k)=0.d0
            endif
            if (kg .eq. nz) then
            u(i,j,k)=0.d0
            v(i,j,k)=0.d0
            endif
            if (kg .eq. nz+1) then
            w(i,j,k)=0.d0
            endif
         enddo
      enddo
   enddo
   !$acc end kernels

   ! find local maximum velocity
   uc=maxval(u)
   vc=maxval(v)
   wc=maxval(w)
   umax=max(wc,max(uc,vc))
   call MPI_Allreduce(umax,gumax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
   cou=gumax*dt*dxi
   if (rank.eq.0) then
      write(*,*) "CFL (max among tasks)", cou
      if (cou .gt. 7) stop
   endif


   call cpu_time(timef)
   if (rank.eq.0) print '(" Time elapsed = ",f6.1," ms")',1000*(timef-times)
   !########################################################################################################################################
   ! END STEP 8: VELOCITY CORRECTION  
   !########################################################################################################################################
   call nvtxEndRange


   !########################################################################################################################################
   ! START STEP 9: OUTPUT FIELDS 
   ! ########################################################################################################################################
   if (mod(t,dump) .eq. 0) then
      if (rank .eq. 0) write(*,*) "Saving output files"
         ! write velocity and pressure fiels (1-4)
         call writefield(t,1)
         call writefield(t,2)
         call writefield(t,3)
         call writefield(t,4)
         #if phiflag == 1
         ! write phase-field (5)
         call writefield(t,5)
         #endif
   endif
   !########################################################################################################################################
   ! END STEP 9: OUTPUT FIELDS N  
   !########################################################################################################################################

call nvtxEndRange
!call nvtxEndRange
enddo
call cpu_time(t_end)
elapsed = t_end-t_start
if (rank .eq. 0) write(*,*)  'Elapsed time (seconds):', elapsed
!$acc end data
!$acc end data
!$acc end data

! Remove allocated variables (add new)
deallocate(u,v,w)
deallocate(tanh_psi, mysin, mycos)
deallocate(rhsu,rhsv,rhsw)
deallocate(rhsu_o,rhsv_o,rhsw_o)
deallocate(phi,rhsphi,normx,normy,normz)

call mpi_finalize(ierr)

end program main