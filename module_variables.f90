module module_variables

!#################################################

!ALLOCATABLE ARRAYS

implicit none

!Velocity model array
real, allocatable, dimension(:,:) :: mod_vel, mod_vs

!Density model array
real, allocatable, dimension(:,:) :: mod_den

!Pressure field arrays
real, allocatable, dimension(:,:) :: P1, P2, P3, P

!Velocity field arrays
real, allocatable, dimension(:,:) :: U,V

!Displacement field arrays
real, allocatable, dimension(:,:) :: Xx,Zz,Tt

!Inverse density array
real, allocatable, dimension(:,:) :: b

!Elastic Properties
real, allocatable, dimension(:,:) :: K_modulus, L, M

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

!#################################################

end module module_variables



