program onda_2d

!##################################################
implicit none

!Velocity model array
real, allocatable, dimension(:,:) :: mod_vel, mod_vs

!Density model array
real, allocatable, dimension(:,:) :: mod_den

!Velocity field arrays
real, allocatable, dimension(:,:) :: U,V

!Displacement field arrays
real, allocatable, dimension(:,:) :: Xx,Zz,Tt

!Inverse density array
real, allocatable, dimension(:,:) :: b

!Elastic Properties
real, allocatable, dimension(:,:) :: L, M

!Cerjan factors array
real, allocatable, dimension(:) :: vetor_cerjan

!Seismogram array
real, allocatable, dimension(:,:) :: seismogram

!#################################################

!MODELLING VARIABLES

!Velocity model file
character(len=30) :: file_vp_model, file_vs_model

!Density model file
character(len=30) :: file_density_model

!Seismic source variables
integer :: nfx, nfz
real :: A, fcorte, fc, t0, tm, fonte, passos_fonte, fx, fz 

!Numerical solution stability variables
real :: alfa, beta

!Grid, time & spatial operators variables
integer :: Nx, Nz, h, passos, tmax
real :: x, z, dt, delta, aux 

!Cerjan layers variables
integer :: pontos_cerjan
real :: fator_cerjan


!Snapshot variables
character(len=3) :: snap_char, string
integer :: nsnap,c_snap, dt_snap
real :: passos_snap 

!Seismogram variables
integer :: receptor_z
real :: passos_seismo

!Counter variables
integer :: i,j,k,n,t 

!Modelling constants
integer :: comp_byte=4
real :: pi=3.14159
real :: start_time,final_time

!#################################################

call CPU_TIME(start_time)
call read_verify_input()
call allocate_variables()
call read_verify_velocity_model()
call read_verify_density_model
call calculate_source_modelling_parameters()

write(*,*)"Number of snapshots:",nsnap

open(unit=50,file="fonte.txt")
do i=1,passos 
    write(50,*) i*dt,staggered_grid_source(A,pi,fc,k,dt,t0)
enddo

close(50)
!LOOP DO TEMPO
do k=1, passos
    if (k<=passos_fonte) then
    Xx(nfz,nfx)=Xx(nfz,nfx)+staggered_grid_source(A,pi,fc,k,dt,t0)
    Zz(nfz,nfx)=Zz(nfz,nfx)+staggered_grid_source(A,pi,fc,k,dt,t0)
    end if
    
    call elastic_second_order_staggered_grid_operator_loop()
    call anti_reflection_cerjan_conditions()

    if(mod(k,dt_snap)==0) then
       c_snap=c_snap+1
       write(*,*) "Gerando snapshot:", c_snap
       write(snap_char,"(i3)") c_snap
       call take_snapshot()
    end if

    do t=1,Nx
        seismogram(k,t)=Xx(receptor_z,t)
    end do

end do

open(25,file="seismogram.bin",status="unknown",access="direct",form="unformatted",recl=comp_byte*Nx*k)
write(25,rec=1) seismogram 
close(25)

call CPU_TIME(final_time)
write(*,*) "Tempo de processamento decorrido:",final_time-start_time,"segundos."

contains
!#################################################

subroutine read_verify_input()
implicit none

open (30, file="input.txt")
read(30,*) x, z, h, tmax
read(30,*) file_vp_model
read(30,*) file_density_model
read(30,*) file_vs_model
read(30,*) fx, fz, A, fcorte
read(30,*) alfa, beta
read(30,*) fator_cerjan, pontos_cerjan
read(30,*) nsnap
read(30,*) receptor_z
Nx=nint(x/h)+1
Nz=nint(z/h)+1
nfx=nint(fx/h)+1
nfz=nint(fz/h)+1
write(*,*)"Nx=", Nx, "Nz=", Nz, "h=", h, "tmax=", tmax
write(*,*) "Modelo de velocidades:", file_vp_model
write(*,*) " Modelo de densidades:", file_density_model
write(*,*)"fx=", nfx, "fz=", nfz, "A=", A, "fcorte=", fcorte
write(*,*)"alfa=", alfa, "beta=", beta
write(*,*)"Fator cerjan=", fator_cerjan, "Tamanho da camada cerjan=",pontos_cerjan
close(30)

end subroutine read_verify_input

!#################################################

subroutine read_verify_velocity_model()
implicit none

open(100,file=file_vp_model,status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
read(100,rec=1) ((mod_vel(i,j),i=1,Nz),j=1,Nx)
close(100)

open(95,file=file_vs_model,status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
read(95,rec=1) ((mod_vs(i,j),i=1,Nz),j=1,Nx)
close(95)

open(60,file="model_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(60,rec=1) ((mod_vel(i,j),i=1,Nz),j=1,Nx)
close(60)

end subroutine read_verify_velocity_model

!#################################################

subroutine read_verify_density_model()
implicit none

open(99,file=file_density_model,status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
read(99,rec=1) ((mod_den(i,j),i=1,Nz),j=1,Nx)
close(99)

open(61,file="density_model_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(61,rec=1) ((mod_den(i,j),i=1,Nz),j=1,Nx)
close(61)

end subroutine read_verify_density_model

!#################################################

subroutine allocate_variables()
implicit none

allocate (U(Nz,Nx),V(Nz,Nx))
allocate (b(Nz,Nx),L(Nz,Nx),M(Nz,Nx))
allocate (mod_vel(Nz, Nx),mod_den(Nz,Nx),mod_vs(Nz,Nx))
allocate (Xx(Nz,Nx),Zz(Nz,Nx),Tt(Nz,Nx))
allocate (vetor_cerjan(pontos_cerjan-1))


end subroutine allocate_variables

!#################################################

subroutine calculate_source_modelling_parameters()
implicit none
print*, " "
print*, " "
write(*,*)"Parametros da fonte e de modelagem:"
fc=fcorte/(3*sqrt(pi))
write(*,*)"fc=", fc
t0=2*sqrt(pi)/fcorte
write(*,*)"t0=", t0
tm=2*t0
write(*,*)"tm=", tm
dt=h/(beta*2500) !######################################################ATENCAO VP MAX####################################
write(*,*) "dt=", dt
delta=dt/h 
write(*,*) "dt/h=", delta
U=0
V=0
P=0
L=0
M=0
Xx=0
Zz=0
Tt=0
b=0
passos=nint(tmax/dt)
passos_snap=passos
dt_snap=nint(passos_snap/nsnap)
write(*,*) "Numero de passos=", passos 
passos_fonte=nint(tm/dt)
write(*,*) "Fonte aplicada ate", passos_fonte, "passos"
allocate (seismogram(passos,Nx))
seismogram=0
c_snap=0
vetor_cerjan=0

do j=1, Nx
    do i=1, Nz
        b(i,j)=1/(mod_den(i,j))
    end do
end do

M=(mod_vs**2)*mod_den
L=(mod_vel**2)*mod_den-2*M

open(42,file="b_model_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(42,rec=1) ((b(i,j),i=1,Nz),j=1,Nx)
close(42)

open(37,file="vs_model_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(37,rec=1) ((mod_vs(i,j),i=1,Nz),j=1,Nx)
close(37)

open(39,file="M_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(39,rec=1) ((M(i,j),i=1,Nz),j=1,Nx)
close(39)

open(38,file="L_verification.bin",status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
write(38,rec=1) ((L(i,j),i=1,Nz),j=1,Nx)
close(38)

end subroutine calculate_source_modelling_parameters

!#######################################################

real function staggered_grid_source(A,pi,fc,k,dt,t0) 
implicit none
integer :: k
real :: A, pi, fc, dt, t0

staggered_grid_source=-A*(k*dt-t0)*exp(-pi*(pi*fc*(k*dt-t0))**2)

end function staggered_grid_source

!#######################################################

subroutine elastic_second_order_staggered_grid_operator_loop()
implicit none

do j=2, Nx
    do i=2,Nz
        U(i,j)=U(i,j)+delta*b(i,j)*(Xx(i,j+1)-Xx(i,j)+Tt(i,j)-Tt(i-1,j))
        V(i,j)=V(i,j)+delta*b(i,j)*(Tt(i,j)-Tt(i,j-1)+Zz(i+1,j)-Zz(i,j))
    end do
end do

do j=2, Nx-1
    do i=2, Nz-1
        Xx(i,j)=Xx(i,j)+delta*(L(i,j)+2*M(i,j))*(U(i,j)-U(i,j-1))+delta*L(i,j)*(V(i,j)-V(i-1,j))
        Zz(i,j)=Zz(i,j)+delta*(L(i,j)+2*M(i,j))*(V(i,j)-V(i-1,j))+delta*L(i,j)*(U(i,j)-U(i,j-1))
        Tt(i,j)=Tt(i,j)+delta*M(i,j)*(V(i,j+1)-V(i,j)+U(i+1,j)-U(i,j))
    end do
end do

end subroutine elastic_second_order_staggered_grid_operator_loop

!#######################################################

subroutine anti_reflection_cerjan_conditions()
implicit none

do i=1,pontos_cerjan-1
    vetor_cerjan(i)=exp(-((fator_cerjan*(pontos_cerjan-1-i))**2))
end do


!Applying the cerjan conditions on the left region
do j=1,pontos_cerjan-1
    do i=1,Nz
        Xx(i,j)=vetor_cerjan(j)*Xx(i,j)
        Zz(i,j)=vetor_cerjan(j)*Zz(i,j)
        Tt(i,j)=vetor_cerjan(j)*Tt(i,j)
    end do
end do
close(23)

!Applying the cerjan conditions on the right region
n=pontos_cerjan
do j=Nx-pontos_cerjan+2, Nx 
    n=n-1
    do i=1,Nz
        Xx(i,j)=vetor_cerjan(n)*Xx(i,j)
        Zz(i,j)=vetor_cerjan(n)*Zz(i,j)
        Tt(i,j)=vetor_cerjan(n)*Tt(i,j)
    end do
end do

!Applying the cerjan conditions to the bottom region
n=pontos_cerjan
do i=Nz-pontos_cerjan+2,Nz
    n=n-1
    do j=1,Nx
        Xx(i,j)=vetor_cerjan(n)*Xx(i,j)
        Zz(i,j)=vetor_cerjan(n)*Zz(i,j)
        Tt(i,j)=vetor_cerjan(n)*Tt(i,j)
    end do
end do

end subroutine anti_reflection_cerjan_conditions

!#######################################################
subroutine take_snapshot()
implicit none
open (20,file="snap"//trim(adjustl(snap_char))//".bin",status="replace",access="direct",form="unformatted", &
& recl=comp_byte*Nx*Nz)
open (27,file="snap"//trim(adjustl(snap_char))//"_U.bin",status="replace",access="direct",form="unformatted", &
& recl=comp_byte*Nx*Nz)
open (28,file="snap"//trim(adjustl(snap_char))//"_V.bin",status="replace",access="direct",form="unformatted", &
& recl=comp_byte*Nx*Nz)
open (33,file="snap"//trim(adjustl(snap_char))//"_T.bin",status="replace",access="direct",form="unformatted", &
& recl=comp_byte*Nx*Nz)
write(20,rec=1) ((Xx(i,j),i=1,Nz),j=1,Nx)
write(27,rec=1) ((U(i,j),i=1,Nz),j=1,Nx)
write(28,rec=1) ((V(i,j),i=1,Nz),j=1,Nx)
write(33,rec=1) ((Tt(i,j),i=1,Nz),j=1,Nx)
close(20)
close(27)
close(28)

end subroutine take_snapshot

!#######################################################
end program onda_2d
