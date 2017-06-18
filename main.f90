!  Main_Programm
!
!

!****************************************************************************
!
!****************************************************************************

program Main

implicit none

integer jmax,j, n_picsel, ind_calc, pol_ind, shade_type ! flags & params
real,allocatable :: rut(:, :, :)! geom
complex, allocatable :: j_e(:, :)! surface currents
real, allocatable :: s(:)! area
real, allocatable :: norm(:, :)! vet norm
real, allocatable :: rkt(:, :)! tochki collokacii
real, allocatable :: tau(:, :, :) ! local basis
real h_max ! partition diametr

real k  !volnovoe chislo
real k_v(3) ! volnovoy vektor
real k_v_isx(3) ! volnovoy vektor in read file
real k_v_0(3) ! ort_of_volnovoy vektor
real ssum ! summarnaya ploshad
real e_0_enter(3) ! vector of polarization

integer m,n !parameters or cels partition
real rb ! vortex radius
real dl, g_abs

real r_s ! radius-vector to source
real h0 ! polarization vector
integer, allocatable :: sha(:) ! shaded part
complex, allocatable :: j0(:,:), je0(:,:) ! array for even currents j0
!complex, allocatable :: sigma(:)

integer, allocatable :: soot(:,:,:), ind(:,:,:) ! arrays for element neighbours
real, allocatable :: shad_cont(:,:,:)
real :: taut(3),w,mu,mu0,c,c0,e0(3),eps,eps0,r(3),zr,sigma,pi
integer :: i,jj,num_shad,alpha

integer i_beg, i_end, INFO, i_dummy

print*, "==========================="
!======= input parameters =====================================================!
print *,"Reading init parameters..."
call read_Init_Param_single_solve(ind_calc,k,k_v_0,e_0_enter,pol_ind)
k_v = k_v_0*k

if(ind_calc.eq.1) then
	print*,"*** Direct diagramm ***" 
else
	print*,"*** Inverse diagramm ***"
endif

if (pol_ind.eq.1) then
	print*, "- vertical polarization -"
elseif (pol_ind.eq.2) then
	print*, "- horizontal polarization -"
else
	print*, "tutu - 1"
endif
print*, "-------------------------"

!======= constants ============================================================!
! physical
pi=3.141592653589793
 c0 = 299792458	! light speed
mu0 = 4*pi*1.d-7
eps0 = 1/(mu0*c0**2)
mu = 1
eps = 1
 c = c0/sqrt(eps*mu)
w = k*c
print*,"w, Hz ->", w/(2*pi)
print*, "Wavelength, m ->", 2*pi/k

! solution & tolerances
n_picsel = 500 ! number of directions to calculate RCS
zr = 1.d-4  ! zero value dependending on the problem's size

!======= input geometry =======================================================!
print *, "Reading geometry..."
call readGeom_0(jmax) ! bulk geometry file creation
allocate(rut(3,4,jmax))
call readGeom(rut,jmax) ! full geometry creation
print*, "Number of cells ->", jmax
!======= 
call diam_partition(rut,jmax,h_max)
print*, "h_max ->", h_max
print*, "Wavelength to element ratio ->", 2*pi/k/h_max
print*, "-------------------------"

!======== calculate norms, tau, etc ===========================================!
print *, "Initializing geometry..."
allocate(s(jmax),norm(3,jmax),rkt(3,jmax))
allocate(tau(3,2,jmax))
call Init_Geom(rut,jmax,rkt,s,norm,tau)

!======== calculate shaded part ===============================================!
print *, "Find shaded part..."
call SYSTEM_CLOCK(i_beg, i_dummy, i_dummy)
allocate(sha(jmax))
shade_type = 2	! define the type of shading procedure 
if (shade_type.eq.1) then  
	! convex body
	call shade_geom(norm,k_v,sha,jmax)
else	
	! ray-tracing method
	call shade_geom_rt(k_v,sha,jmax,rkt,rut)
endif
! Write to file
call sha2file(jmax,sha)

!======== calculate shaded contour ============================================!
print *, "The Shaded Contour find we should (c)"
! find all neighbours
allocate(soot(4,5,jmax),ind(4,5,jmax))
call sootv(rut,jmax,soot,ind,zr) 

allocate(shad_cont(3,2,jmax))
if(ind_calc.eq.1) then
	! Direct RCS -> need to calculate shaded part
	call shade_contour(sha,jmax,rut,soot,shad_cont,num_shad)
	! write to geodat-file
	call geom2geodat(rut,shad_cont,jmax,1,1,50,30,num_shad)
	call SYSTEM_CLOCK(i_end, i_dummy, i_dummy)
	print *, "shade elements&contour time = ", i_end-i_beg
elseif (ind_calc.eq.2) then
	! Inverse RCS -> no need in it
	print *,"shaded contour will be calculated later"
endif

!======== form currents j0 ====================================================!
print *, "Calculate currents j0"
allocate(j0(3,jmax),je0(3,jmax))

if(ind_calc.eq.1) then
	! Direct RCS
	call SYSTEM_CLOCK(i_beg, i_dummy, i_dummy)
	call calc_j0(k_v,rkt,w,mu,mu0,e_0_enter,jmax,pol_ind,c,norm,sha,j0,je0)
	! Write to file
	call j02file(jmax,je0)
	call calc_im_re_j0(je0,jmax)
	call SYSTEM_CLOCK(i_end, i_dummy, i_dummy)
	print *, "current find time  time = ", i_end-i_beg
elseif (ind_calc.eq.2) then		
	! Inverse RCS
	print*, 'currents saved for later'	
else
	! Reserve for future
	print*,"tutu"
endif

!======== calc EPR (RCS) ======================================================!
print *, "Calculate RCS"
call SYSTEM_CLOCK(i_beg, i_dummy, i_dummy)

if(ind_calc.eq.1) then
	! Direct RCS
	! eps1 = eps0*eps, eps=1 -> eps1 = eps0
	print*, 'direct RCS computation...'
	call calc_rcs_direct(k,w,jmax,eps0,je0,s,rkt,n_picsel)
elseif (ind_calc.eq.2) then
	! Inverse RCS
	print*, 'inverse RCS computation...'
	call calc_rcs_inverse(k,w,jmax,eps0,j0,je0,s,rkt,rut,n_picsel,mu,mu0,&
 pol_ind,c,norm,sha)	
else
	! Reserved for future
	print*,"tutu"
endif
call SYSTEM_CLOCK(i_end, i_dummy, i_dummy)
print *, "RCS find  time = ", i_end-i_beg

print*, "Solution COMPLETED"
stop
end program Main
 
