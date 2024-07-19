subroutine set_fft_plans(nx, ny, nz, planxf, planxb, planyf, planyb, planzf, planzb)
	implicit none
    include "fftw3.f"
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nz
    complex(8), allocatable, dimension (:) :: tempx, tempy, tempz

    allocate(tempx(nx), tempy(ny), tempz(nz))
    call dfftw_plan_dft_1d(planxf,nx,tempx,tempx,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planxb,nx,tempx,tempx,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planyf,ny,tempy,tempy,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planyb,ny,tempy,tempy,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planzf,nz,tempz,tempz,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(planzb,nz,tempz,tempz,FFTW_BACKWARD,FFTW_ESTIMATE)
    deallocate(tempx, tempy, tempz)
end subroutine set_fft_plans

subroutine destroy_fft_plans(planxf,planxb,planyf,planyb,planzf,planzb)
	implicit none
	include "fftw3.f"
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb   
    call dfftw_destroy_plan(planxf)
    call dfftw_destroy_plan(planxb)
    call dfftw_destroy_plan(planyf)
    call dfftw_destroy_plan(planyb)
    call dfftw_destroy_plan(planzf)
    call dfftw_destroy_plan(planzb)
end subroutine destroy_fft_plans

subroutine output_parameters(nu, CFL, initial_type, nproc, nzp, id)
    implicit none
    integer :: initial_type, nproc, nzp, id
    real(8) :: nu, CFL
    if (id == 0) then
        print *, '======================================================================='
        print *, 'DNS of the Schrodinger-Pauli formulation of 2D incompressible viscous flows with constant density rho_0=1.'

        print *, 'nu = ', nu
        print *, 'CFL = ', CFL

        select case (initial_type)
        case (1)
            print *, 'Initial condition: 2D Taylor-Green vortices'
        case (2)
            print *, 'Initial condition: 2D Taylor-Green vortices with random noise'
        end select

        print * , 'nproc=', nproc
        print * , 'nzp=', nzp
        print *, '======================================================================='
    end if
end subroutine output_parameters

subroutine openfiles(id)
    implicit none
    integer :: id
    if (id == 0) then
        open(31, file='./dns_statistics/enstrophy.dat', status='unknown')
        open(32, file='./dns_statistics/total_energy.dat', status='unknown')
        open(34, file='./dns_statistics/mean_velocity.dat', status='unknown')
        open(40, file='./dns_statistics/integral_length.dat', status='unknown')
        open(42, file='./dns_statistics/enstrophy_spectrum.dat', status='unknown')
        open(43, file='./dns_statistics/energy_spectrum.dat', status='unknown')
        open(45, file='./dns_statistics/time.dat', status='unknown')
    end if
end subroutine openfiles

subroutine closefiles(id)
    implicit none
    integer :: id
    if (id == 0) then
        close(31)
        close(32)
        close(34)
        close(40)
        close(42)
        close(43)
        close(45)
    end if
end subroutine closefiles

SUBROUTINE wavenumber(nx,ny,nz,nzp,kx,ky,kz,k2,kd,id)
    implicit none
    integer nx, ny, nz, nzp, i, j, k, id
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    do i=1,nx
        kx(i) = mod(i-1+nx/2,nx) - nx/2
    end do
    do j=1,ny
        ky(j) = mod(j-1+ny/2,ny) - ny/2
    end do
    do k=1,nzp
        kz(k) = mod(k-1+nz/2+id*nzp,nz) - nz/2
    end do
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                k2(i,j,k) = kx(i)**2 + ky(j)**2 + kz(k)**2
            end do
        end do
    end do
    kd(:,:,:) = 1
    do k=1,nzp
        do j=1,ny
            do i=1,nx
            	if (k2(i,j,k) > (nx**2 + ny**2 + nz**2)/18) then
            		kd(i,j,k) = 0
            	end if
            end do
        end do
    end do
END SUBROUTINE wavenumber

subroutine initialize_mesh(nx,ny,nzp,xstart,ystart,zstart,dx,dy,dz,meshx,meshy,meshz,id)
	implicit none
    integer nx, ny, nzp, i, j, k, id
    real(8) :: xstart, ystart ,zstart
    real(8) :: dx, dy, dz
    real(8), dimension (nx,ny,nzp) :: meshx, meshy, meshz 
    do i=1,nx
        meshx(i,:,:) = xstart + real(i-1)*dx
    end do
    do j=1,ny
        meshy(:,j,:) = ystart + real(j-1)*dy
    end do
    do k=1,nzp
        meshz(:,:,k) = zstart + real(k-1+nzp*id)*dz
    end do
end subroutine initialize_mesh 

subroutine initialize_psi(nx, ny, nz, nzp, meshx, meshy, meshz, velx, vely, velz, vorx, vory, vorz, psi1, psi2, kx, ky, kz, k2, kd, ik2, initial_type, id, nproc, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs)
	implicit none
    include 'mpif.h'
    integer nx, ny, nz, nzp, id, nproc, ierr, i, j, k
    integer initial_type
    real(8), dimension (nx,ny,nzp) :: meshx, meshy, meshz, velx, vely, velz, vorx, vory, vorz
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd, ik2
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb
    integer, dimension(nproc) :: counts, displs
 
    real(8), parameter :: pi = 3.14159265359d0
    real(8), allocatable, dimension (:,:,:) :: Hy, r, fr
    complex(8), allocatable, dimension (:,:,:) :: u, v, P, Q
    real(8), parameter :: sigma = 3

    select case (initial_type)
    case (1) ! 2D-TG
        allocate(Hy(nx, ny, nzp))
        do k = 1, nzp
            do j = 1, ny
                do i = 1, nx
                    if (meshy(i, j, k) <= pi) then
                        Hy(i, j, k) = meshy(i, j, k)/2d0
                    else
                        Hy(i, j, k) = pi - meshy(i, j, k)/2d0
                    end if
                end do
            end do
        end do
        psi1 = cos(Hy) * exp(cmplx(0d0, 1d0)*cos(meshz)*(2d0 - cos(meshy)))
        psi2 = sin(Hy) * exp(-cmplx(0d0, 1d0)*cos(meshz)*(2d0 + cos(meshy)))
        deallocate(Hy)
    case (2) ! 2D vortex using rational map
        allocate(r(nx, ny, nzp), fr(nx, ny, nzp), u(nx, ny, nzp), v(nx, ny, nzp), P(nx, ny, nzp), Q(nx, ny, nzp))
        r = sqrt(meshy**2 + meshz**2)
        fr = exp(-(r/sigma)**4)
        u = 2d0*cmplx(meshy, meshz)*fr / (1d0 + r**2)
        v = cmplx(0d0, 1d0) * (r**2 + 1d0 - 2d0*fr) / (1d0 + r**2)
        P = u
        Q = v**2
        psi1 = P / sqrt(abs(P)**2 + abs(Q)**2)
        psi2 = Q / sqrt(abs(P)**2 + abs(Q)**2)
        deallocate(r, fr, u, v, P, Q)
        call non_divergence_projection(psi1, psi2, nx, ny, nz, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    end select

    call compute_fluid_quantities(vorx, vory, vorz, velx, vely, velz, psi1, psi2, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
end subroutine initialize_psi

subroutine transpose_yz(mx,my,mz,nproc,spec_tran,spectemp_tran,id,counts,displs)
    implicit none
    include 'mpif.h'
    integer mx, my, mz, nproc, j, k, imp, mzp, myp, id, ierr
    integer, dimension (nproc) :: counts, displs
    complex(8), dimension (mx,my,mz/nproc) :: spec_tran
    complex(8), dimension (mx,mz,my/nproc) :: spectemp_tran
    complex(8), allocatable, dimension (:,:,:) :: specb, spectemp1

    allocate(spectemp1(mx,my/nproc,mz), specb(mx,my/nproc,mz/nproc))
    mzp = mz/nproc
    myp = my/nproc
    do imp=1,nproc
        do k=1,mzp
            do j=1,myp
                specb(:,j,k) = spec_tran(:, j+(imp-1)*myp, k)
            enddo
        enddo
        call mpi_gatherv(specb,counts(id+1),MPI_DOUBLE_COMPLEX,spectemp1,counts,displs,MPI_DOUBLE_COMPLEX,imp-1,mpi_comm_world,ierr)
    enddo
    do j=1,myp 
        do k=1,mz
            spectemp_tran(:,k,j) = spectemp1(:,j,k)
        enddo
    enddo
    deallocate(specb, spectemp1)
end subroutine transpose_yz

!! 3D forward fast fourier transformation
!! mpi for z direction, 
!!    number of blocks: nproc=nz/nzp 
!!    id = 0, 1, ... , nproc-1
!! input:
!!    phy: real, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!	  fft plans: planxf, planyf, planzf
!!    counts, displs: integer, dimension (nproc) 
!!    counts(i)=nx*nyp*nzp 
!!    displs(i)=(i-1)*counts(i)
!! output:
!!	  spec: complex,dimension (nx,ny,nzp)	
subroutine fourier_forward(phy,spec_ff,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
	implicit none
	include 'mpif.h'
	include "fftw3.f"
	integer(8) planxf, planyf, planzf
    !!!MPI
    integer :: id, nproc, ierr
	integer nx, ny, nz, nzp, nyp, i, j, k
    real(8), dimension (nx,ny,nz/nproc) :: phy
    complex(8), dimension (nx,ny,nz/nproc) :: spec_ff
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:) :: tempx, tempy, tempz
    complex(8), allocatable, dimension (:,:,:) :: spectemp

    allocate(tempx(nx), tempy(ny), tempz(nz))
    allocate(spectemp(nx,nz,ny/nproc))
    nyp = ny/nproc
    nzp = nz/nproc
    do k=1,nzp
        do j=1,ny
        	do i=1,nx
            	tempx(i) = cmplx(phy(i,j,k)/real(nx), 0d0)
            enddo
            call dfftw_execute_dft(planxf, tempx, tempx)
            spec_ff(:,j,k) = tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy = spec_ff(i,:,k)/real(ny)
            call dfftw_execute_dft(planyf, tempy, tempy)
            spec_ff(i,:,k) = tempy
        enddo
    enddo
    call transpose_yz(nx, ny, nz, nproc, spec_ff, spectemp, id, counts, displs)
    do j=1,nyp
        do i=1,nx
            tempz = spectemp(i,:,j)/real(nz)
            call dfftw_execute_dft(planzf, tempz, tempz)
            spectemp(i,:,j) = tempz
        enddo
    enddo
    call transpose_yz(nx, nz, ny, nproc, spectemp, spec_ff, id, counts, displs)
    deallocate(tempx, tempy, tempz)
    deallocate(spectemp)
end subroutine fourier_forward


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D backward fast fourier transformation
!! mpi for z direction, 
!!    number of blocks: nproc=nz/nzp 
!!    id = 0, 1, ... , nproc-1
!! input:
!!    spec: complex, dimension (nx,ny,nzp)
!!    mesh size: nx, ny, nz
!!	  fft plans: planxf, planyf, planzf
!!    counts, displs: integer, dimension (nproc) 
!!    counts(i)=nx*nyp*nzp 
!!    displs(i)=(i-1)*counts(i)
!! output:
!!	  phy: real, dimension (nx,ny,nzp)	
subroutine fourier_backward(phy,spec_fb,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
	implicit none
	include 'mpif.h'
	include "fftw3.f"
	integer(8) planxb, planyb, planzb
    integer :: id, nproc, ierr
	integer nx, ny, nz, nzp, nyp, i, j, k
    real(8), dimension (nx,ny,nz/nproc) :: phy
    integer, dimension (nx,ny,nz/nproc) :: kd
    complex(8), dimension (nx,ny,nz/nproc) :: spec_fb
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:) :: tempx, tempy, tempz
    complex(8), allocatable, dimension (:,:,:) :: spectemp

    allocate(tempx(nx), tempy(ny), tempz(nz))
    allocate(spectemp(nx,nz,ny/nproc))
    nyp = ny/nproc
    nzp = nz/nproc
    spec_fb = spec_fb*kd
    do k=1,nzp
        do j=1,ny
            tempx = spec_fb(:,j,k)
            call dfftw_execute_dft(planxb,tempx,tempx)
            spec_fb(:,j,k) = tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy = spec_fb(i,:,k)
            call dfftw_execute_dft(planyb,tempy,tempy)
            spec_fb(i,:,k) = tempy
        enddo
    enddo
    call transpose_yz(nx,ny,nz,nproc,spec_fb,spectemp,id,counts,displs)
    do j=1,nyp
        do i=1,nx
            tempz = spectemp(i,:,j)
            call dfftw_execute_dft(planzb,tempz,tempz)
            spectemp(i,:,j) = tempz
        enddo
    enddo
    call transpose_yz(nx,nz,ny,nproc,spectemp,spec_fb,id,counts,displs)
    phy = real(spec_fb)
    deallocate(tempx, tempy, tempz)
    deallocate(spectemp)
end subroutine fourier_backward

subroutine dx_dy_dz_dp_dm(phy,dphy,switch_d,nx,ny,nzp,kx,ky,kz,k2,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    implicit none
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    integer :: id, nproc, switch_d
    real(8), dimension (nx,ny,nzp) :: phy, dphy
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:,:,:) :: spec

    allocate(spec(nx,ny,nzp))
    nz = nzp*nproc
    call fourier_forward(phy,spec,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    if(switch_d==1) then
        do i=1,nx
            spec(i,:,:) = spec(i,:,:)*cmplx(0d0, real(kx(i)))
        end do
    else if(switch_d==2) then
        do j=1,ny
            spec(:,j,:) = spec(:,j,:)*cmplx(0d0, real(ky(j)))
        end do
    else if(switch_d==3) then
        do k=1,nzp
            spec(:,:,k) = spec(:,:,k)*cmplx(0d0, real(kz(k)))
        end do
    else if(switch_d==6) then
        spec = -real(k2)*spec
    else if(switch_d==-6) then
        spec = -spec/real(k2)
        if(id==0) then
            spec(1,1,1) = cmplx(0d0, 0d0)
        end if
    end if
    call fourier_backward(dphy,spec,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    deallocate(spec)
end subroutine dx_dy_dz_dp_dm

! if (switch_d==-1) (dphy1,dphy2,dphy3) is divergence free
subroutine cross_vector(phy1,phy2,phy3,dphy1,dphy2,dphy3,switch_d,nx,ny,nzp,kx,ky,kz,k2,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    implicit none
    integer*8 planxf,planyf,planzf,planxb,planyb,planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz, kd
    integer, dimension (nx,ny,nzp) :: k2
    integer :: id, nproc, switch_d
    real(8), dimension (nx,ny,nzp) :: phy1, phy2, phy3, dphy1, dphy2, dphy3
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:,:,:) :: spec1, spec2, spec3, spec4, spec5, spec6

    allocate(spec1(nx,ny,nzp), spec2(nx,ny,nzp), spec3(nx,ny,nzp), spec4(nx,ny,nzp), spec5(nx,ny,nzp), spec6(nx,ny,nzp))
    nz = nzp*nproc
    call fourier_forward(phy1,spec1,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(phy2,spec2,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(phy3,spec3,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    if(switch_d==1) then
    	do k=1,nzp
    		do j=1,ny
    			do i=1,nx
    				spec4(i,j,k) = cmplx(0d0,1d0)*(real(ky(j))*spec3(i,j,k) - real(kz(k))*spec2(i,j,k))
    				spec5(i,j,k) = cmplx(0d0,1d0)*(real(kz(k))*spec1(i,j,k) - real(kx(i))*spec3(i,j,k))
    				spec6(i,j,k) = cmplx(0d0,1d0)*(real(kx(i))*spec2(i,j,k) - real(ky(j))*spec1(i,j,k))
    			enddo
    		enddo
    	enddo
    elseif(switch_d==-1) then
    	do k=1,nzp
    		do j=1,ny
    			do i=1,nx
    				spec4(i,j,k) = cmplx(0d0,1d0)*(real(ky(j))*spec3(i,j,k) - real(kz(k))*spec2(i,j,k)) / real(k2(i,j,k))
    				spec5(i,j,k) = cmplx(0d0,1d0)*(real(kz(k))*spec1(i,j,k) - real(kx(i))*spec3(i,j,k)) / real(k2(i,j,k))
    				spec6(i,j,k) = cmplx(0d0,1d0)*(real(kx(i))*spec2(i,j,k) - real(ky(j))*spec1(i,j,k)) / real(k2(i,j,k))
    			enddo
    		enddo
    	enddo
    	if(id==0) then
            spec4(1,1,1) = cmplx(0d0, 0d0)
            spec5(1,1,1) = cmplx(0d0, 0d0)
            spec6(1,1,1) = cmplx(0d0, 0d0)
        endif
    endif
    deallocate(spec1, spec2, spec3)
    call fourier_backward(dphy1,spec4,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call fourier_backward(dphy2,spec5,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call fourier_backward(dphy3,spec6,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    deallocate(spec4, spec5, spec6)
end subroutine cross_vector

subroutine divergence(phy1,phy2,phy3,dphy,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    implicit none
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, dimension (nx,ny,nzp) :: kd
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer :: id, nproc
    real(8), dimension (nx,ny,nzp) :: phy1, phy2, phy3, dphy
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:,:,:) :: spec1, spec2, spec3

    allocate(spec1(nx,ny,nzp), spec2(nx,ny,nzp), spec3(nx,ny,nzp))
    nz = nzp*nproc
    call fourier_forward(phy1,spec1,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(phy2,spec2,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(phy3,spec3,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    do k=1,nzp
    	do j=1,ny
    		do i=1,nx
    			spec1(i,j,k) = cmplx(0d0,1d0)*(real(kx(i))*spec1(i,j,k) + real(ky(j))*spec2(i,j,k) + real(kz(k))*spec3(i,j,k))
    		enddo
    	enddo
    enddo
    call fourier_backward(dphy,spec1,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    deallocate(spec1, spec2, spec3)
end subroutine divergence

subroutine gradient(phy,dphy1,dphy2,dphy3,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    implicit none
    integer nx, ny, nz, nzp, i, j, k
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: kd
    integer :: id, nproc
    real(8), dimension (nx,ny,nzp) :: dphy1, dphy2, dphy3, phy
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:,:,:) :: spec1, spec2, spec3, spec4

    nz = nzp*nproc
    allocate(spec1(nx,ny,nzp), spec2(nx,ny,nzp), spec3(nx,ny,nzp), spec4(nx,ny,nzp))
    call fourier_forward(phy,spec4,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    do k=1,nzp
    	do j=1,ny
    		do i=1,nx
    			spec1(i,j,k) = cmplx(0d0,1d0)*real(kx(i))*spec4(i,j,k)
    			spec2(i,j,k) = cmplx(0d0,1d0)*real(ky(j))*spec4(i,j,k)
    			spec3(i,j,k) = cmplx(0d0,1d0)*real(kz(k))*spec4(i,j,k)
    		enddo
    	enddo
    enddo
    deallocate(spec4)
    call fourier_backward(dphy1,spec1,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call fourier_backward(dphy2,spec2,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call fourier_backward(dphy3,spec3,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    deallocate(spec1, spec2, spec3)
end subroutine gradient

subroutine cross_product(nx,ny,nzp,phy1,phy2,phy3,phy11,phy12,phy13,phy21,phy22,phy23)
	implicit none
    integer nx, ny, nzp
    real(8), dimension (nx,ny,nzp) :: phy1, phy2, phy3, phy11, phy12, phy13, phy21, phy22, phy23
    phy21 = phy2*phy13 - phy3*phy12
    phy22 = phy3*phy11 - phy1*phy13
    phy23 = phy1*phy12 - phy2*phy11
end subroutine cross_product

subroutine output_field(outputstep, nx, ny, nzp, meshx, meshy, meshz, vely, velz, vorx, psi1, psi2, id, nproc, countsb, displsb)
    implicit none
    include 'mpif.h'
    integer outputstep, nx, ny, nzp, nz, id, nproc, ierr, i
    integer, dimension (nproc) :: countsb, displsb
    real(8), dimension (nx,ny,nzp) :: meshx, meshy, meshz, vely, velz, vorx
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer indexa, indexb, indexc, indexd, indexe, indexf
    character*200 name
    integer, parameter :: nbox = 10
    character*40, dimension (nbox) :: varname
    real(4), allocatable, dimension (:,:,:,:) :: data_box, data_box_core

    allocate(data_box_core(nx,ny,nzp,nbox))
    data_box_core(:,:,:,1) = meshx
    data_box_core(:,:,:,2) = meshy
    data_box_core(:,:,:,3) = meshz
    data_box_core(:,:,:,4) = vely
    data_box_core(:,:,:,5) = velz
    data_box_core(:,:,:,6) = vorx
    data_box_core(:,:,:,7) = real(psi1)
    data_box_core(:,:,:,8) = aimag(psi1)
    data_box_core(:,:,:,9) = real(psi2)
    data_box_core(:,:,:,10) = aimag(psi2)
    if (id == 0) then
        nz = nzp*nproc
        allocate(data_box(nx,ny,nz,nbox))
        indexa = mod(outputstep,1000000)/100000
        indexb = mod(outputstep,100000)/10000
        indexc = mod(outputstep,10000)/1000
        indexd = mod(outputstep,1000)/100
        indexe = mod(outputstep,100)/10
        indexf = mod(outputstep,10)/1
        name = './field_output/field'//char(indexa+48)//char(indexb+48)//char(indexc+48)//char(indexd+48)//char(indexe+48)//char(indexf+48)//'.plt'
        varname(1) = 'x'
        varname(2) = 'y'
        varname(3) = 'z'
        varname(4) = 'vely'
        varname(5) = 'velz'
        varname(6) = 'vorx'
        varname(7) = 'Re(psi1)'
        varname(8) = 'Im(psi1)'
        varname(9) = 'Re(psi2)'
        varname(10) = 'Im(psi2)'
    end if

    do i = 1, nbox
        call mpi_gatherv(data_box_core(:,:,:,i), countsb(id+1), MPI_real, data_box(:,:,:,i), countsb, displsb, MPI_real, 0, mpi_comm_world, ierr)
    end do

    deallocate(data_box_core)
    if (id==0) then
        call output_3d_tecplot_bin(data_box, nx, ny, nz, nbox, name, varname)
        deallocate(data_box)
    end if
    call MPI_Barrier(mpi_comm_world, ierr)
end subroutine output_field


! ================= output binary data for tecplot =====================
! input:
!    v_box: real, dimension (nx,ny,nz,nbox)
!    nbox: number of output varible
!    mesh size: nx, ny, nz
!	  name: character*200 name of the output data
!    varname: character*40, dimension (nbox) names of varible
subroutine output_3d_tecplot_bin(v_box,nx,ny,nz,nbox,name,varname)
	implicit none
    integer nx, ny, nz, nbox
    real(4), dimension (nx,ny,nz,nbox) :: v_box
    real(8), dimension(nbox) :: min_value, max_value
    real(4) ZONEMARKER, EOHMARKER
    integer len, i
    character*40 Title, var
    character*200 name
    character*40, dimension(nbox) :: varname
    character*40 Zonename
    character(40) instring
    ZONEMARKER= 299.0
    EOHMARKER = 357.0
    do i=1,nbox
    	min_value(i) = minval(v_box(:,:,:,i))
        max_value(i) = maxval(v_box(:,:,:,i))
    enddo
    open(unit=99, file=name, form="BINARY")
    !I. The header section.
    !1.1 Magic number, Version number
    write(99) "#!TDV112"
    !1.2. Integer value of 1.
    write(99) 1
    !1.3. Title and variable names.
    !Filetype
    write(99) 0
    !1.3.1. The TITLE.
    Title = ""
    call dumpstring(Title)
    !1.3.2 Number of variables (NumVar) in the datafile.
    write(99) nbox
    !1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]
    do i=1,nbox
    	call dumpstring(varname(i))
    enddo
    !1.4. Zones
    !Zone marker. Value = 299.0
    write(99) ZONEMARKER
    !Zone name.
    Zonename = "ZONE 001"
    call dumpstring(Zonename)
    !ParentZone
    write(99) -1
    !StrandID
    write(99) -1
    !solution time
    write(99) 0
    write(99) 0
    !not used
    write(99) -1
    !ZoneType  
    write(99) 0
    !DataPacking 0=Block, 1=Point
    write(99) 0
    !Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
    write(99) 0
    !Number of user defined face neighbor connections (value >= 0)
    write(99) 0
    !IMax,JMax,KMax
    write(99) nx
    write(99) ny
    write(99) nz
    !1=Auxiliary name/value pair to follow   0=No more Auxiliar name/value pairs.
    write(99) 0
    !I HEADER OVER
    !EOHMARKER, value=357.0
    write(99) EOHMARKER
    !II. Data section
    !2.1 zone
    write(99) Zonemarker
    !variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit                
    do i=1,nbox
        write(99) 1                                                               
    enddo
    !Has variable sharing 0 = no, 1 = yes.       
    write(99) 0
    !Has passive variables 0 = no, 1 = yes.       
    write(99) 0                                  
    !Zone number to share connectivity list with (-1 = no sharing).
    write(99) -1
    !min value
    !max value
    do i = 1, nbox
        write(99) min_value(i)
        write(99) max_value(i)                                                               
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Zone Data. Each variable is in data format asspecified above.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(99) v_box
    close(99)
end subroutine output_3d_tecplot_bin

subroutine dumpstring(instring)
    !!!for binary output
    character(40) instring
    integer len
    len = LEN_TRIM(instring)
    do i = 1, len
        ii = ICHAR(instring(i:i))
        write(99) ii
    enddo
    write(99) 0
end

subroutine time_advance_psi(psi1, psi2, nu, dt, nx, ny, nz, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
	include 'mpif.h'
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer :: nx, ny, nz, nzp, id, nproc
    real(8) :: nu, dt
    integer, dimension (nproc) :: counts, displs
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer, dimension (nx) :: kx
    integer, dimension (ny) :: ky
    integer, dimension (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd

    complex(8), allocatable, dimension (:,:,:) :: psi1_spec, psi2_spec
    real(8), allocatable, dimension (:,:,:) :: Vr, Vi, P
    complex(8), allocatable, dimension (:,:,:) :: Q

    allocate(psi1_spec(nx,ny,nzp), psi2_spec(nx,ny,nzp))
    call fourier_forward_cmplx(psi1, psi1_spec, nx, ny, nz, planxf, planyf, planzf, counts, displs, id, nproc)
    call fourier_forward_cmplx(psi2, psi2_spec, nx, ny, nz, planxf, planyf, planzf, counts, displs, id, nproc)
    psi1_spec = psi1_spec * exp(-cmplx(nu, 0.5d0)*k2*dt)
    psi2_spec = psi2_spec * exp(-cmplx(nu, 0.5d0)*k2*dt)
    call fourier_backward_cmplx(psi1, psi1_spec, nx, ny, nz, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call fourier_backward_cmplx(psi2, psi2_spec, nx, ny, nz, planxb, planyb, planzb, counts, displs, id, nproc, kd)

    allocate(Vr(nx,ny,nzp), Vi(nx,ny,nzp), P(nx,ny,nzp), Q(nx,ny,nzp))
    call compute_potential(Vr, Vi, P, Q, psi1, psi2, nu, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    psi1_spec = exp(cmplx(Vi, -Vr-P)*dt)*psi1 + (exp(-cmplx(0d0, 1d0)*Q*dt) - 1d0)*psi2
    psi2_spec = exp(cmplx(Vi, -Vr+P)*dt)*psi2 + (exp(-cmplx(0d0, 1d0)*conjg(Q)*dt) - 1d0)*psi1

    psi1 = psi1_spec / sqrt(psi1_spec*conjg(psi1_spec) + psi2_spec*conjg(psi2_spec))
    psi2 = psi2_spec / sqrt(psi1_spec*conjg(psi1_spec) + psi2_spec*conjg(psi2_spec))
    call non_divergence_projection(psi1, psi2, nx, ny, nz, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    deallocate(psi1_spec, psi2_spec)
end subroutine time_advance_psi

subroutine compute_potential(Vr, Vi, P, Q, psi1, psi2, nu, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
	include 'mpif.h'
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer :: nx, ny, nz, nzp, id, nproc
    real(8) :: nu
    integer, dimension (nproc) :: counts, displs
    real(8), dimension (nx,ny,nzp) :: Vr, Vi, P
    complex(8), dimension (nx,ny,nzp) :: Q
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer, dimension (nx) :: kx
    integer, dimension (ny) :: ky
    integer, dimension (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd

    real(8), allocatable, dimension(:,:,:) :: a, b, c, d, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz, ddx, ddy, ddz
    real(8), allocatable, dimension(:,:,:) :: s1, s2, s3, f1, f2, f3
    real(8), allocatable, dimension(:,:,:) :: GradPsi2, Laplacian_s1, Laplacian_s2, Laplacian_s3

    allocate(a(nx,ny,nzp), b(nx,ny,nzp), c(nx,ny,nzp), d(nx,ny,nzp), dax(nx,ny,nzp), day(nx,ny,nzp), daz(nx,ny,nzp), dbx(nx,ny,nzp), dby(nx,ny,nzp), dbz(nx,ny,nzp), dcx(nx,ny,nzp), dcy(nx,ny,nzp), dcz(nx,ny,nzp), ddx(nx,ny,nzp), ddy(nx,ny,nzp), ddz(nx,ny,nzp))
    a = real(psi1)
    b = aimag(psi1)
    c = real(psi2)
    d = aimag(psi2)
    call gradient(a, dax, day, daz, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call gradient(b, dbx, dby, dbz, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call gradient(c, dcx, dcy, dcz, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call gradient(d, ddx, ddy, ddz, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    deallocate(a, b, c, d)
    allocate(GradPsi2(nx,ny,nzp))
    GradPsi2 = dax**2 + day**2 + daz**2 + dbx**2 + dby**2 + dbz**2 + dcx**2 + dcy**2 + dcz**2 + ddx**2 + ddy**2 + ddz**2
    Vr = GradPsi2/2d0
    Vi = nu*GradPsi2
    deallocate(dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz, ddx, ddy, ddz)

    allocate(s1(nx,ny,nzp), s2(nx,ny,nzp), s3(nx,ny,nzp))
    s1 = 1d0 - 2d0*psi2*conjg(psi2)
    s2 = 2d0*(aimag(psi1)*real(psi2) - real(psi1)*aimag(psi2))
    s3 = 2d0*(real(psi1)*real(psi2) + aimag(psi1)*aimag(psi2))
    allocate(Laplacian_s1(nx,ny,nzp), Laplacian_s2(nx,ny,nzp), Laplacian_s3(nx,ny,nzp))
    call dx_dy_dz_dp_dm(s1, Laplacian_s1, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call dx_dy_dz_dp_dm(s2, Laplacian_s2, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call dx_dy_dz_dp_dm(s3, Laplacian_s3, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    allocate(f1(nx,ny,nzp), f2(nx,ny,nzp), f3(nx,ny,nzp))
    call compute_f(f1, f2, f3, GradPsi2, s1, s2, s3, psi1, psi2, nu, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    deallocate(s1, s2, s3, GradPsi2)

    ! call test_data_max(f1, nx, ny, nzp, id, nproc)
    ! call test_data_max(f2, nx, ny, nzp, id, nproc)
    ! call test_data_max(f3, nx, ny, nzp, id, nproc)

    P = Laplacian_s1/4d0 - f1
    Q = cmplx(Laplacian_s3, Laplacian_s2)/4d0 - cmplx(f3, f2)
    deallocate(Laplacian_s1, Laplacian_s2, Laplacian_s3, f1, f2, f3)
end subroutine compute_potential

subroutine compute_f(f1, f2, f3, GradPsi2, s1, s2, s3, psi1, psi2, nu, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nz, nzp, id, nproc, i, j, k
    real(8) :: nu
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    real(8), dimension (nx,ny,nzp) :: GradPsi2, f1, f2, f3, s1, s2, s3
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb
    integer, dimension(nproc) :: counts, displs

    real(8), allocatable, dimension(:,:,:) :: vorx, vory, vorz, velx, vely, velz
    real(8), allocatable, dimension(:,:,:) :: ds1x, ds1y, ds1z, ds2x, ds2y, ds2z, ds3x, ds3y, ds3z
    real(8) :: A(2, 3), AAT(2, 2), IAAT(2, 2), Ap(3, 2), b(2), det

    allocate(vorx(nx,ny,nzp), vory(nx,ny,nzp), vorz(nx,ny,nzp), velx(nx,ny,nzp), vely(nx,ny,nzp), velz(nx,ny,nzp))
    allocate(ds1x(nx,ny,nzp), ds1y(nx,ny,nzp), ds1z(nx,ny,nzp), ds2x(nx,ny,nzp), ds2y(nx,ny,nzp), ds2z(nx,ny,nzp), ds3x(nx,ny,nzp), ds3y(nx,ny,nzp), ds3z(nx,ny,nzp))
    call compute_fluid_quantities(vorx, vory, vorz, velx, vely, velz, psi1, psi2, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    call gradient(s1, ds1x, ds1y, ds1z, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call gradient(s2, ds2x, ds2y, ds2z, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call gradient(s3, ds3x, ds3y, ds3z, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)

    do k = 1, nzp
        do j = 1, ny
            do i = 1, nx
                A(1, 1) = ds1y(i,j,k)
                A(1, 2) = ds2y(i,j,k)
                A(1, 3) = ds3y(i,j,k)
                A(2, 1) = ds1z(i,j,k)
                A(2, 2) = ds2z(i,j,k)
                A(2, 3) = ds3z(i,j,k)
                AAT(1, 1) = A(1, 1)*A(1, 1) + A(1, 2)*A(1, 2) + A(1, 3)*A(1, 3)
                AAT(1, 2) = A(1, 1)*A(2, 1) + A(1, 2)*A(2, 2) + A(1, 3)*A(2, 3)
                AAT(2, 1) = A(2, 1)*A(1, 1) + A(2, 2)*A(1, 2) + A(2, 3)*A(1, 3)
                AAT(2, 2) = A(2, 1)*A(2, 1) + A(2, 2)*A(2, 2) + A(2, 3)*A(2, 3)
                det = AAT(1, 1)*AAT(2, 2) - AAT(1, 2)*AAT(2, 1)
                IAAT(1, 1) = AAT(2, 2) / det
                IAAT(1, 2) = -AAT(1, 2) / det
                IAAT(2, 1) = -AAT(2, 1) / det
                IAAT(2, 2) = AAT(1, 1) / det
                Ap(1, 1) = A(1, 1)*IAAT(1, 1) + A(2, 1)*IAAT(2, 1)
                Ap(1, 2) = A(1, 1)*IAAT(1, 2) + A(2, 1)*IAAT(2, 2)
                Ap(2, 1) = A(1, 2)*IAAT(1, 1) + A(2, 2)*IAAT(2, 1)
                Ap(2, 2) = A(1, 2)*IAAT(1, 2) + A(2, 2)*IAAT(2, 2)
                Ap(3, 1) = A(1, 3)*IAAT(1, 1) + A(2, 3)*IAAT(2, 1)
                Ap(3, 2) = A(1, 3)*IAAT(1, 2) + A(2, 3)*IAAT(2, 2)
                b(1) = 2*nu*GradPsi2(i,j,k)*vely(i,j,k)
                b(2) = 2*nu*GradPsi2(i,j,k)*velz(i,j,k)
                f1(i,j,k) = Ap(1, 1)*b(1) + Ap(1, 2)*b(2)
                f2(i,j,k) = Ap(2, 1)*b(1) + Ap(2, 2)*b(2)
                f3(i,j,k) = Ap(3, 1)*b(1) + Ap(3, 2)*b(2)
            end do
        end do
    end do

    call test_eq_rhs(f1, f2, f3, ds1y, ds2y, ds3y, ds1z, ds2z, ds3z, vely, velz, GradPsi2, nu, nx, ny, nzp, id, nproc)
    deallocate(vorx, vory, vorz, velx, vely, velz, ds1x, ds1y, ds1z, ds2x, ds2y, ds2z, ds3x, ds3y, ds3z)
end subroutine compute_f

subroutine compute_fluid_quantities(vorx, vory, vorz, velx, vely, velz, psi1, psi2, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nzp, id, nproc
    real(8), dimension (nx,ny,nzp) :: velx, vely, velz, vorx, vory, vorz
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb
    integer, dimension(nproc) :: counts, displs, countsb, displsb

    real(8), allocatable, dimension(:,:,:) :: s1, s2, s3
    real(8), allocatable, dimension(:,:,:) :: ds1x, ds1y, ds1z, ds2x, ds2y, ds2z, ds3x, ds3y, ds3z

    allocate(s1(nx,ny,nzp), s2(nx,ny,nzp), s3(nx,ny,nzp))
    s1 = 1d0 - 2d0*psi2*conjg(psi2)
    s2 = 2d0*(aimag(psi1)*real(psi2) - real(psi1)*aimag(psi2))
    s3 = 2d0*(real(psi1)*real(psi2) + aimag(psi1)*aimag(psi2))

    allocate(ds1x(nx,ny,nzp), ds1y(nx,ny,nzp), ds1z(nx,ny,nzp), ds2x(nx,ny,nzp), ds2y(nx,ny,nzp), ds2z(nx,ny,nzp), ds3x(nx,ny,nzp), ds3y(nx,ny,nzp), ds3z(nx,ny,nzp))
    call gradient(s1,ds1x,ds1y,ds1z,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call gradient(s2,ds2x,ds2y,ds2z,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    call gradient(s3,ds3x,ds3y,ds3z,nx,ny,nzp,kx,ky,kz,planxf,planyf,planzf,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    vorx = 0.5d0 * (s1*(ds2y*ds3z - ds2z*ds3y) + s2*(ds3y*ds1z - ds3z*ds1y) + s3*(ds1y*ds2z - ds1z*ds2y))
    vory = 0.5d0 * (s1*(ds2z*ds3x - ds2x*ds3z) + s2*(ds3z*ds1x - ds3x*ds1z) + s3*(ds1z*ds2x - ds1x*ds2z))
    vorz = 0.5d0 * (s1*(ds2x*ds3y - ds2y*ds3x) + s2*(ds3x*ds1y - ds3y*ds1x) + s3*(ds1x*ds2y - ds1y*ds2x))
    deallocate(s1, s2, s3)
    deallocate(ds1x, ds1y, ds1z, ds2x, ds2y, ds2z, ds3x, ds3y, ds3z)
    call cross_vector(vorx, vory, vorz, velx, vely, velz, -1, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
end subroutine compute_fluid_quantities

subroutine non_divergence_projection(psi1, psi2, nx, ny, nz, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nz, nzp, id, nproc
    complex(8), dimension (nx,ny,nzp) :: psi1, psi2
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    integer(8) :: planxf, planyf, planzf, planxb, planyb, planzb
    integer, dimension(nproc) :: counts, displs

    real(8), allocatable, dimension(:,:,:) :: a, b, c, d, Laplacian_a, Laplacian_b, Laplacian_c, Laplacian_d
    real(8), allocatable, dimension(:,:,:) :: fs, q
    complex(8), allocatable, dimension(:,:,:) :: spec

    allocate(a(nx,ny,nzp), b(nx,ny,nzp), c(nx,ny,nzp), d(nx,ny,nzp))
    a = real(psi1)
    b = aimag(psi1)
    c = real(psi2)
    d = aimag(psi2)
    allocate(Laplacian_a(nx,ny,nzp), Laplacian_b(nx,ny,nzp), Laplacian_c(nx,ny,nzp), Laplacian_d(nx,ny,nzp))
    call dx_dy_dz_dp_dm(a, Laplacian_a, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call dx_dy_dz_dp_dm(b, Laplacian_b, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call dx_dy_dz_dp_dm(c, Laplacian_c, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    call dx_dy_dz_dp_dm(d, Laplacian_d, 6, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    deallocate(a, b, c, d)
    allocate(fs(nx,ny,nzp))
    fs = real(-cmplx(0d0,1d0)*conjg(psi1)*cmplx(Laplacian_a, Laplacian_b) - cmplx(0d0,1d0)*conjg(psi2)*cmplx(Laplacian_c, Laplacian_d))
    deallocate(Laplacian_a, Laplacian_b, Laplacian_c, Laplacian_d)
    allocate(spec(nx,ny,nzp))
    call fourier_forward(fs, spec, nx, ny, nz, planxf, planyf, planzf, counts, displs, id, nproc)
    deallocate(fs)
    spec = -spec / real(k2)
    if (id == 0) then
        spec(1,1,1) = cmplx(0d0, 0d0)
    end if
    allocate(q(nx,ny,nzp))
    call fourier_backward(q, spec, nx, ny, nz, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    psi1 = psi1 * exp(-cmplx(0d0, 1d0) * q)
    psi2 = psi2 * exp(-cmplx(0d0, 1d0) * q)
    deallocate(q, spec)
end subroutine non_divergence_projection

subroutine get_dt(CFL, velx, vely, velz, nx, ny, nzp, id, dt, dx, dy, dz, nproc)
    implicit none
    include 'mpif.h'
    real(8) :: dx, dy, dz, dt
    real(8) :: CFL
    real(8), dimension (nx,ny,nzp) :: velx, vely, velz
    integer :: nx, ny, nzp, id, ierr, nproc
    real(8) :: lam11, lam21, lam31, lam1, lam2, lam3
    real(8) :: dt1, dt2, dt3

    lam11 = maxval(abs(velx))
    lam21 = maxval(abs(vely))
    lam31 = maxval(abs(velz))
    call MPI_REDUCE(lam11,lam1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(lam21,lam2,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(lam31,lam3,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    if (id == 0) then
    	dt1 = CFL*dx/lam1
    	dt2 = CFL*dy/lam2
    	dt3 = CFL*dz/lam3
    	dt = min(dt1, dt2, dt3)
    end if
    call mpi_bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
end subroutine get_dt

subroutine get_statistics(nx, ny, nzp, velx, vely, velz, vorx, vory, vorz, kx, ky, kz, k2, ik2, kd, planxf, planyf, planzf, counts, displs, id, nproc, time_t)
    implicit none
    include 'mpif.h'
    integer(8) planxf, planyf, planzf
    integer, dimension (nx,ny,nzp) :: k2, kd, ik2
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer :: id, nproc, ierr
    integer nx, ny, nzp, nz, i, j, k, nek
    real(8), parameter :: pi=3.141592653589793d0
    real(8) time_t
    real(8), dimension (nx,ny,nzp) :: velx, vely, velz, vorx, vory, vorz
    integer, dimension (nproc) :: counts, displs
    real(8) :: enstrophy, total_energy, mean_velocity, integral_length
    real(8) :: enstrophy1, total_energy1, integral_length1
    complex(8), allocatable, dimension (:,:,:) :: spec1, spec2, spec3, spec4, spec5, spec6
    real(8), allocatable, dimension (:,:,:) :: spectral_k

    nz = nzp*nproc
    nek = int(sqrt((nx**2+ny**2+nz**2)/18d0 + 0d0))
	
    ! get total energy and energy spectrum:
    allocate(spec1(nx,ny,nzp), spec2(nx,ny,nzp), spec3(nx,ny,nzp), spec4(nx,ny,nzp), spec5(nx,ny,nzp), spec6(nx,ny,nzp), spectral_k(nx,ny,nzp))
    call fourier_forward(velx,spec1,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(vely,spec2,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    call fourier_forward(velz,spec3,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    spectral_k = 0.5d0 * real(spec1*conjg(spec1) + spec2*conjg(spec2) + spec3*conjg(spec3))
    call get_spectrum(time_t,nx,ny,nzp,nek,spectral_k,43,k2,ik2,id)
    total_energy1 = sum(spectral_k)
    call MPI_REDUCE(total_energy1, total_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    do i = 1, nek
        integral_length1 = integral_length1 + sum(spectral_k, mask=(ik2 == i))/real(i)
    end do
    call MPI_REDUCE(integral_length1, integral_length, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! get total enstrophy and enstrophy spectrum:
    do k=1,nzp
        do j=1,ny
            do i=1,nx
                spec4(i,j,k) = cmplx(0d0,1d0)*(real(ky(j))*spec3(i,j,k) - real(kz(k))*spec2(i,j,k))
                spec5(i,j,k) = cmplx(0d0,1d0)*(real(kz(k))*spec1(i,j,k) - real(kx(i))*spec3(i,j,k))
                spec6(i,j,k) = cmplx(0d0,1d0)*(real(kx(i))*spec2(i,j,k) - real(ky(j))*spec1(i,j,k))
            end do
        end do
    end do
    deallocate(spec1, spec2, spec3)
    spectral_k = 0.5d0 * real(spec4*conjg(spec4) + spec5*conjg(spec5) + spec6*conjg(spec6))
    deallocate(spec4, spec5, spec6)
    call get_spectrum(time_t, nx, ny, nzp, nek, spectral_k, 42, k2, ik2, id)
    enstrophy1 = sum(spectral_k)
    call MPI_REDUCE(enstrophy1, enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    deallocate(spectral_k)

    if (id == 0) then
        integral_length = 3d0*pi*integral_length/4d0/total_energy
        mean_velocity = sqrt(2d0*total_energy/3d0)

        total_energy = total_energy * 4d0 * pi**2
        enstrophy = enstrophy * 4d0 * pi**2

        print *, 'Energy = ', total_energy
        print *, 'Enstropy = ', enstrophy

        write(31,*) time_t, enstrophy
        write(32,*) time_t, total_energy
        write(34,*) time_t, mean_velocity
        write(40,*) time_t, integral_length
    end if
end subroutine get_statistics

subroutine get_spectrum(time_t, nx, ny, nzp, nek, spec_k, num_data, k2, ik2, id)
    implicit none
    include 'mpif.h'
    integer nx, ny, nzp, num_data, nek, i, id, ierr
    real(8) :: time_t
    real(8), dimension (nx,ny,nzp) :: spec_k
    integer, dimension (nx,ny,nzp) :: k2, ik2
    real(8) spectrum_k1, spectrum_k

    if (id == 0) write(num_data,*) 'zone t=','"',time_t,'"'
    do i = 1, nek
        spectrum_k1 = sum(spec_k, mask=(ik2 == i))
        call MPI_REDUCE(spectrum_k1, spectrum_k, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (id == 0) write(num_data,*) i, spectrum_k
    end do
end subroutine get_spectrum

subroutine fourier_forward_cmplx(phy,spec_ff,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
	implicit none
	include 'mpif.h'
	include "fftw3.f"
	integer(8) planxf, planyf, planzf
    integer :: id, nproc, ierr
	integer nx, ny, nz, nzp, nyp, i, j, k
    complex(8), dimension (nx,ny,nz/nproc) :: phy
    complex(8), dimension (nx,ny,nz/nproc) :: spec_ff
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:) :: tempx, tempy, tempz
    complex(8), allocatable, dimension (:,:,:) :: spectemp

    allocate(tempx(nx), tempy(ny), tempz(nz))
    allocate(spectemp(nx,nz,ny/nproc))
    nyp = ny/nproc
    nzp = nz/nproc
    do k=1,nzp
        do j=1,ny
        	do i=1,nx
            	tempx(i) = phy(i,j,k)/real(nx)
            enddo
            call dfftw_execute_dft(planxf,tempx,tempx)
            spec_ff(:,j,k) = tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy = spec_ff(i,:,k)/real(ny)
            call dfftw_execute_dft(planyf,tempy,tempy)
            spec_ff(i,:,k) = tempy
        enddo
    enddo
    call transpose_yz(nx,ny,nz,nproc,spec_ff,spectemp,id,counts,displs)
    do j=1,nyp
        do i=1,nx
            tempz = spectemp(i,:,j)/real(nz)
            call dfftw_execute_dft(planzf,tempz,tempz)
            spectemp(i,:,j) = tempz
        enddo
    enddo
    call transpose_yz(nx,nz,ny,nproc,spectemp,spec_ff,id,counts,displs)
    deallocate(tempx, tempy, tempz)
    deallocate(spectemp)
end subroutine fourier_forward_cmplx

subroutine fourier_backward_cmplx(phy,spec_fb,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
	implicit none
	include 'mpif.h'
	include "fftw3.f"
	integer(8) planxb, planyb, planzb
    integer :: id, nproc, ierr
	integer nx, ny, nz, nzp, nyp, i, j, k
    complex(8), dimension (nx,ny,nz/nproc) :: phy
    integer, dimension (nx,ny,nz/nproc) :: kd
    complex(8), dimension (nx,ny,nz/nproc) :: spec_fb
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:) :: tempx, tempy, tempz
    complex(8), allocatable, dimension (:,:,:) :: spectemp

    allocate(tempx(nx), tempy(ny), tempz(nz))
    allocate(spectemp(nx,nz,ny/nproc))
    nyp = ny/nproc
    nzp = nz/nproc
    spec_fb = spec_fb*kd
    do k=1,nzp
        do j=1,ny
            tempx = spec_fb(:,j,k)
            call dfftw_execute_dft(planxb,tempx,tempx)
            spec_fb(:,j,k) = tempx
        enddo
    enddo
    do k=1,nzp
        do i=1,nx
            tempy = spec_fb(i,:,k)
            call dfftw_execute_dft(planyb,tempy,tempy)
            spec_fb(i,:,k) = tempy
        enddo
    enddo
    call transpose_yz(nx,ny,nz,nproc,spec_fb,spectemp,id,counts,displs)
    do j=1,nyp
        do i=1,nx
            tempz = spectemp(i,:,j)
            call dfftw_execute_dft(planzb,tempz,tempz)
            spectemp(i,:,j) = tempz
        enddo
    enddo
    call transpose_yz(nx,nz,ny,nproc,spectemp,phy,id,counts,displs)
    deallocate(tempx, tempy, tempz)
    deallocate(spectemp)
end subroutine fourier_backward_cmplx

subroutine Laplacian_cmplx(phy, dphy, switch_d, nx, ny, nzp, kx, ky, kz, k2, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    implicit none
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nz, nzp, i, j, k
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    integer :: id, nproc, switch_d
    complex(8), dimension (nx,ny,nzp) :: phy, dphy
    integer, dimension (nproc) :: counts, displs
    complex(8), allocatable, dimension (:,:,:) :: spec

    allocate(spec(nx,ny,nzp))
    nz = nzp*nproc
    call fourier_forward_cmplx(phy,spec,nx,ny,nz,planxf,planyf,planzf,counts,displs,id,nproc)
    if(switch_d==1) then
        do i=1,nx
            spec(i,:,:) = spec(i,:,:)*cmplx(0d0, real(kx(i)))
        end do
    else if(switch_d==2) then
        do j=1,ny
            spec(:,j,:) = spec(:,j,:)*cmplx(0d0, real(ky(j)))
        end do
    else if(switch_d==3) then
        do k=1,nzp
            spec(:,:,k) = spec(:,:,k)*cmplx(0d0, real(kz(k)))
        end do
    else if(switch_d==6) then
        spec = -real(k2)*spec
    else if(switch_d==-6) then
        spec = -spec/real(k2)
        if(id==0) then
            spec(1,1,1) = cmplx(0d0, 0d0)
        end if
    end if
    call fourier_backward_cmplx(dphy,spec,nx,ny,nz,planxb,planyb,planzb,counts,displs,id,nproc,kd)
    deallocate(spec)
end subroutine Laplacian_cmplx

subroutine gaussian(u, nx, ny, nzp)
    implicit none
    complex(8), dimension(nx,ny,nzp) :: u
    integer :: nx, ny, nzp
    integer :: i, j, k
    real(8) :: t1, t2
    real(8), parameter :: pi = 3.1415926535898d0
    
    u = (0d0, 0d0)
    do k = 1, nzp
        do j = 1, ny
            do i = 1, nx
                call random_number(t1) 
                call random_number(t2) 
                if (t1 <= 1.e-10) t1 = 1.e-10
                if (t2 <= 1.e-10) t2 = 1.e-10
                t2 = 2d0*pi*t2
                u(i,j,k) = sqrt(-2.0*log(t1))*cmplx(cos(t2),sin(t2))
            end do
        end do
    end do
end subroutine gaussian

subroutine test_div_vel(velx, vely, velz, nx, ny, nzp, kx, ky, kz, k2, kd, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc)
    implicit none
    include 'mpif.h'
    integer(8) planxf, planyf, planzf, planxb, planyb, planzb
    integer nx, ny, nzp, id, nproc, ierr
    integer, DIMENSION (nx) :: kx
    integer, DIMENSION (ny) :: ky
    integer, DIMENSION (nzp) :: kz
    integer, dimension (nx,ny,nzp) :: k2, kd
    real(8), dimension (nx,ny,nzp) :: velx, vely, velz
    integer, dimension (nproc) :: counts, displs

    real(8), allocatable, dimension(:,:,:) :: theta
    real(8) :: theta_max_id, theta_max

    allocate(theta(nx,ny,nzp))
    call divergence(velx, vely, velz, theta, nx, ny, nzp, kx, ky, kz, planxf, planyf, planzf, planxb, planyb, planzb, counts, displs, id, nproc, kd)
    theta_max_id = maxval(abs(theta))
    call MPI_REDUCE(theta_max_id, theta_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (id == 0) print *, theta_max
    deallocate(theta)
end subroutine test_div_vel

subroutine test_data_max(f, nx, ny, nzp, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nzp, id, nproc, ierr
    real(8), dimension (nx,ny,nzp) :: f

    real(8) :: f_max_id, f_max

    f_max_id = maxval(abs(f))
    call MPI_REDUCE(f_max_id, f_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (id == 0) print *, f_max
end subroutine test_data_max

subroutine test_data_min(f, nx, ny, nzp, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nzp, id, nproc, ierr
    real(8), dimension (nx,ny,nzp) :: f

    real(8) :: f_min_id, f_min

    f_min_id = minval(abs(f))
    call MPI_REDUCE(f_min_id, f_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    if (id == 0) print *, f_min
end subroutine test_data_min

subroutine test_eq_rhs(f1, f2, f3, ds1y, ds2y, ds3y, ds1z, ds2z, ds3z, vely, velz, GradPsi2, nu, nx, ny, nzp, id, nproc)
    implicit none
    include 'mpif.h'
    integer nx, ny, nzp, id, nproc, ierr
    real(8) :: nu
    real(8), dimension (nx,ny,nzp) :: f1, f2, f3, ds1y, ds2y, ds3y, ds1z, ds2z, ds3z, vely, velz, GradPsi2

    real(8), allocatable, dimension (:,:,:) :: rhs_y, rhs_z, rhs
    real(8) :: rhs_max_id, rhs_max, rhs_mean_id, rhs_mean

    allocate(rhs_y(nx,ny,nzp), rhs_z(nx,ny,nzp), rhs(nx,ny,nzp))
    rhs_y = ds1y*f1 + ds2y*f2 + ds3y*f3 - 2*nu*GradPsi2*vely
    rhs_z = ds1z*f1 + ds2z*f2 + ds3z*f3 - 2*nu*GradPsi2*velz
    rhs = sqrt(rhs_y**2 + rhs_z**2)
    rhs_max_id = maxval(rhs)
    rhs_mean_id = sum(rhs)
    call MPI_REDUCE(rhs_max_id, rhs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(rhs_mean_id, rhs_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (id == 0) then
        print *, 'The maximum rhs in Eq:', rhs_max
        print *, 'The average rhs in Eq:', rhs_mean / dble(nx*ny*nzp*nproc)
    end if
    deallocate(rhs_y, rhs_z, rhs)
end subroutine test_eq_rhs

subroutine solve_linear_system(A, b, x, n)
    implicit none
    integer :: n, i, j, k, p
    real(8) :: A(n, n), b(n), x(n)
    real(8) :: temp, max, factor

    do k = 1, n-1
        ! Find the pivot row with the largest element in column k
        max = abs(A(k,k))
        p = k
        do i = k+1, n
            if (abs(A(i,k)) > max) then
                max = abs(A(i,k))
                p = i
            end if
        end do

        ! Swap the pivot row with the current row
        if (p /= k) then
            do j = k, n
                temp = A(k,j)
                A(k,j) = A(p,j)
                A(p,j) = temp
            end do
            temp = b(k)
            b(k) = b(p)
            b(p) = temp
        end if

        ! Eliminate the elements below the pivot
        do i = k+1, n
            factor = A(i,k) / A(k,k)
            do j = k+1, n
                A(i,j) = A(i,j) - factor * A(k,j)
            end do
            b(i) = b(i) - factor * b(k)
        end do
    end do

    ! Perform back substitution to find the solution vector x
    x(n) = b(n) / A(n,n)
    do i = n-1, 1, -1
        temp = b(i)
        do j = i+1, n
            temp = temp - A(i,j) * x(j)
        end do
        x(i) = temp / A(i,i)
    end do
end subroutine
