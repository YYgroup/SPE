!=============================================================================================
! Version 1.0
! Author: Zhaoyuan Meng
! Solve the Schrodinger-Pauli formulation of 2D incompressible viscous flows with constant density rho_0=1 using the pseudo-spectral method.
!     Solution domain: a periodic square of side 2pi;
!     Time integrating: first-order Euler scheme.
! 
! Update logs:
!   [v1][240115]: Implement the overall simulation process;
!=============================================================================================

program DNS
    implicit none
    include 'mpif.h'
    include "fftw3.f"
    integer status(MPI_STATUS_SIZE)
    integer, parameter :: nx=1, ny=1024, nz=1024
    real(8), parameter :: nu = 1d0
    real(8), parameter :: CFL = 0.5d0
    real(8), parameter :: time_stop = 10d0
    real(8), parameter :: xstart=-3.14159265359d0, ystart=-3.14159265359d0, zstart=-3.14159265359d0
    ! real(8), parameter :: xstart=0d0, ystart=0d0, zstart=0d0
    real(8), parameter :: lx=6.283185307179586d0, ly=6.283185307179586d0, lz=6.283185307179586d0
    real(8), parameter :: outputdt = 0.01d0
    real(8), parameter :: outputdt_field = 0.01d0
    integer, parameter :: initial_type = 2  ! =1: 2D-TG  =2: 2D vortex using rational map
    real(8) :: dx, dy, dz, time_t, dt, tempt, tempt_field
    integer i, j, k, ierr, outputstep
    integer(4) :: time_step
    integer(8) :: planxf, planxb, planyf, planyb, planzf, planzb   !!! fft plans
    integer id, nproc, nzp, nyp
    real(8) :: starttime, endtime
    complex(8), allocatable, dimension(:,:,:) :: psi1, psi2
    real(8), allocatable, dimension(:,:,:) :: velx, vely, velz, vorx, vory, vorz
    real(8), allocatable, dimension(:,:,:) :: meshx, meshy, meshz
    integer, allocatable, dimension(:) :: kx, ky, kz
    integer, allocatable, dimension(:,:,:) :: k2, kd, ik2
    integer, allocatable, dimension(:) :: counts, displs, countsb, displsb
    integer :: indexa, indexb, indexc, indexd

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, id, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    starttime = MPI_WTIME()
    nzp = nz/nproc
    nyp = ny/nproc
    call output_parameters(nu, CFL, initial_type, nproc, nzp, id)

    allocate(counts(nproc), displs(nproc), countsb(nproc), displsb(nproc))
    do i=1,nproc
        counts(i) = nx*nyp*nzp
        displs(i) = (i-1)*counts(i)
    end do
    do i=1,nproc
        countsb(i) = nx*ny*nzp
        displsb(i) = (i-1)*countsb(i)
    end do

    dx = lx / real(nx)
    dy = ly / real(ny)
    dz = lz / real(nz)
    allocate(psi1(nx,ny,nzp), psi2(nx,ny,nzp))
    allocate(velx(nx,ny,nzp), vely(nx,ny,nzp), velz(nx,ny,nzp), vorx(nx,ny,nzp), vory(nx,ny,nzp), vorz(nx,ny,nzp))
    allocate(meshx(nx,ny,nzp), meshy(nx,ny,nzp), meshz(nx,ny,nzp))
    allocate(k2(nx,ny,nzp), kd(nx,ny,nzp), ik2(nx,ny,nzp))
    allocate(kx(nx), ky(ny), kz(nzp))

    call set_fft_plans(nx, ny, nz, planxf, planxb, planyf, planyb, planzf, planzb)
    call wavenumber(nx, ny, nz, nzp, kx, ky, kz, k2, kd, id)
    ik2 = int(sqrt(k2 + 1d-8) + 0.5d0)

    call initialize_mesh(nx, ny, nzp, xstart, ystart, zstart, dx, dy, dz, meshx, meshy, meshz, id)
    call initialize_psi(nx, ny, nz, nzp, meshx, meshy, meshz, velx, vely, velz, vorx, vory, vorz, psi1, psi2, kx, ky, kz, k2, kd, ik2, initial_type, id, nproc, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs)

    call openfiles(id)

    time_t = 0d0
    time_step = 0
    tempt = time_t
    tempt_field = time_t
    do while (time_t <= time_stop + 0.01d0)
        if (time_t >= tempt) then
            tempt = tempt + outputdt
            call get_statistics(nx, ny, nzp, velx, vely, velz, vorx, vory, vorz, kx, ky, kz, k2, ik2, kd, planxf, planyf, planzf, counts, displs, id, nproc, time_t)
        end if
      
        if (time_t >= tempt_field) then
            tempt_field = tempt_field + outputdt_field
            outputstep = floor(1000*time_t)
            call output_field(outputstep, nx, ny, nzp, meshx, meshy, meshz, vely, velz, vorx, psi1, psi2, id, nproc, countsb, displsb)
        end if

        ! call get_dt(CFL, velx, vely, velz, nx, ny, nzp, id, dt, dx, dy, dz, nproc)
        dt = CFL*dz
        call time_advance_psi(psi1, psi2, nu, dt, nx, ny, nz, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
        call compute_fluid_quantities(vorx, vory, vorz, velx, vely, velz, psi1, psi2, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
        ! call test_div_vel(velx, vely, velz, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)

        if (id == 0) then
            print *, 'time step = ', time_step
            print *, 'dt = ', dt
            print *, 'time = ', time_t
            print *, '***************************************************************'
            write(45,*) time_t, dt
        end if

        time_t = time_t + dt
        time_step = time_step + 1
    end do

    deallocate(counts, displs, countsb, displsb)
    deallocate(psi1, psi2)
    deallocate(velx, vely, velz, vorx, vory, vorz, meshx, meshy, meshz)
    deallocate(kx, ky, kz, k2, kd)

    call destroy_fft_plans(planxf, planxb, planyf, planyb, planzf, planzb)
    call closefiles(id)
    
    endtime = MPI_WTIME()
    if (id==0) print *, 'computational time is', endtime - starttime
    call MPI_Finalize(ierr)
end program DNS