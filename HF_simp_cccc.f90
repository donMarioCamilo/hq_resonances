!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM FOR COMPUTING THE GROUND STATE OF A MANY BODY SYSTEM !!
!!  VIA HARTREE FOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!In this demo the same-flavoured tetraquark ground state will be computed
PROGRAM HF_simp_cccc
    IMPLICIT none
    
    DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
    !mass: quark mass (charm)
    DOUBLE PRECISION, PARAMETER         ::  mass = 1.575d3
    !alpha: strong coupling constant at the mass scale
    DOUBLE PRECISION, PARAMETER         ::  alpha = 0.34
    !Bohr radius
    DOUBLE PRECISION, PARAMETER         ::  bohr_r = 1.d0/(mass*alpha)

    !radial coordinates
    ! 1fm = 5.068 GeV-1
    DOUBLE PRECISION, ALLOCATABLE       ::  r(:)
    ! The array r will start at r0 and increase every spacing until rmax is reached (arange)
    DOUBLE PRECISION, PARAMETER         ::  r0 = 1.d-4*bohr_r, rmax = 1.0d1*bohr_r! * (5.068d0)!1.d-5, rmax = 6.d-3
    !Orbital spacing    
    DOUBLE PRECISION, PARAMETER         ::  gridsize = 1000
    DOUBLE PRECISION, PARAMETER         ::  spacing = (rmax-r0)/gridsize

    !Initial guess function
    DOUBLE PRECISION, ALLOCATABLE       ::  phi_init(:)
    !Non-iterable matrix: one particle hamiltonian (laplacian)
    DOUBLE PRECISION, ALLOCATABLE       ::  h(:,:)

    !Iterables
    !Eigenvector resulting function, Fock matrix
    DOUBLE PRECISION, ALLOCATABLE       ::  phi(:), F_matrix(:,:)
    !Eigenvalues
    DOUBLE PRECISION                    ::  energy_new, energy_old
    !Final energy                       
    DOUBLE PRECISION                    ::  binding_energy
    !Orbital array length
    integer                             ::  n

    !Max iterations in scf procedure
    INTEGER, PARAMETER                  ::  max_step = 50
    !Energy threshold for convergence
    DOUBLE PRECISION, PARAMETER         ::  energy_threshold = 1.d-3

    !Iterators
    INTEGER                             ::  i,step

    !Debugging                          
    DOUBLE PRECISION                    ::   prenormcheck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !The independent variable array is obtained as the integration grid
    r = arange(r0,rmax,spacing)
    n = size(r)
    
    !Initial guess for the orbital (hydrogen-like)
    phi_init = (2.d0/sqrt((bohr_r)**3)) * r * EXP(-r/(bohr_r))
    
    !the grid (r) and phi_init are written into phi_init.dat for later debuging
    OPEN(unit=1,file='results/phi_init.dat',status='replace',action='write')

    DO i=1,n
        WRITE(1,'(2f15.8)') r(i),phi_init(i)
    END DO
    CLOSE(1)

    !Integration is checked by phi_init normalization
    !Should be close to 1.0
    write(*,*) 'Norma de phi_init**2'
    write(*,*) sqrt(integrate(r,phi_init*phi_init))
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Set one-particle hamiltonian dimensions
    ALLOCATE(h(n,n))

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!         BEGINNING OF THE SCF PROCEDURE           !!!!!!!!!!!!!
    WRITE(*,*)  'Starting SCF procedure'

    !Here lies the definition of the one-particle hamiltonian
    !Note that the factor behind /mass should take into account centre-of-mass corrections
    h = -0.375d0/mass * laplacian(r)! + alpha_e * V_coulomb(2.d0,r)

    !Value for the energy initialized
    energy_old = 1000.d0
    
    !Energy results for every step are saved into log.dat file
    OPEN(unit=13,file='results/log.dat',status='replace',action='write')
    !Normalized eigenvectors for each iteration are saved into phi.dat file
    OPEN(unit=15,file='results/phi.dat',status='replace',action='write')
    
    DO step = 1,max_step
        ALLOCATE(F_matrix(n,n))
        !Definition of the Fock matrix
        !The color factor lies behind the Hartree term
        F_matrix = h - 4.d0*alpha/3.d0*Hartree(r,phi_init)
        !diagonalize subroutine is called
        !Lowest eigenvalue is saved in energy_new
        !Corresponding eigenvector saved as phi
        CALL diagonalize(F_matrix,energy_new,phi)

        WRITE(13,*) 'Iteration, Eigenvalue:'
        WRITE(13,*) step, energy_new

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!        ASK FOR CONVERGENCE             !!!!!!!!!!!!
        IF (ABS((energy_new-energy_old)/energy_old) < energy_threshold) THEN
            WRITE(*,*) 'Convergence met in',step,'steps'
            exit
        END IF
        IF (step == max_step) THEN
            WRITE(*,*) 'WARNING: CONVERGENCE NOT MET. Last two eigenvalues are:'
            WRITE(*,'(2f8.4)')  energy_old, energy_new
        END IF
        DEALLOCATE(F_matrix)

        !Normalization of the obtained eigenvector
        prenormcheck = sqrt(integrate(r,phi*phi))
        phi_init = phi/prenormcheck

        !The values of the normalized eigenvector are saved in phi.dat
        WRITE(15,*) phi_init(:)
       
        !The new energy will be the old energy in the next SCF iteration
        energy_old = energy_new

        !The energy is printed via console
        !Note that this formula varies with the number of bodies
        !This is not needed for any calculations. Might be commented out
        binding_energy = 2.d0*(energy_new + integrate(r,phi_init*MATMUL(h,phi_init)))
        WRITE(*,*) 'Step:',step,'completed with energy:', binding_energy
    END DO
    WRITE(13,*)'Ground State Energy:', binding_energy
    WRITE(*,*)
    WRITE(*,*)'Ground State Energy:', binding_energy
    WRITE(*,*)'System mass:', binding_energy+4.d0*mass

    close(13)
    close(15)

    CONTAINS

    FUNCTION laplacian(x)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(in)        ::  x(:)
        double precision                    ::  space
        DOUBLE PRECISION, ALLOCATABLE       ::  laplacian(:,:)

        INTEGER     ::  ilap   !Iterator
        

        ALLOCATE(laplacian(size(x),size(x)))
        !The matrix is initialized
        laplacian = 0.0d0
        space = x(2)-x(1)
        !For the first 2 elements, the advanced derivative is calculated
        ilap=1
        laplacian(ilap,ilap) = -2.0d0/space**2
        laplacian(ilap,ilap+1) = 1.d0/space**2
        ! laplacian(i,i) = 1.d0/space**2
        ! laplacian(i,i+1) = -2.d0/space**2
        ! laplacian(i,i+2) = 1.d0/space**2
        
        !For the central n-4 elements, the centered derivative is calculated
        DO ilap = 2,n-1
            laplacian(ilap,ilap) = -2.0d0/space**2
            laplacian(ilap,ilap+1) = 1.d0/space**2
            laplacian(ilap,ilap-1) = 1.d0/space**2 
        END DO
        !For the last 2 elements, the retarded derivative is calculated
        ilap=n
        laplacian(ilap,ilap) = -2.0d0/space**2
        laplacian(ilap,ilap-1) = 1.d0/space**2
        ! laplacian(i,i) = 1.d0/space**2
        ! laplacian(i,i-1) = -2.d0/space**2
        ! laplacian(i,i-2) = 1.d0/space**2
        RETURN
    END FUNCTION
        
    FUNCTION integrate(rgrid,rfun)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)        ::  rgrid(:),rfun(:)

        ! DOUBLE PRECISION, ALLOCATABLE       ::  cosgrid1(:),cosgrid2(:)
        ! DOUBLE PRECISION, ALLOCATABLE       ::  int_r(:)!,int_cos1(:),int_cos2(:)
        ! DOUBLE PRECISION                    ::  auxr!,auxcos1,auxcos2,auxcos
        DOUBLE PRECISION                    ::  integrate
        INTEGER                             ::  Npartr
        double precision                    ::  space
        ! DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
        INTEGER                             ::  ir!,  ntrash

        
        space  = rgrid(2)-rgrid(1)
        Npartr = size(rgrid)
        
        integrate = rfun(1) + rfun(Npartr)+4*rfun(Npartr-1)
        DO ir = 1,Npartr/2
            integrate = integrate + 4*rfun(2*ir)+2*rfun(1+2*ir)
        END DO
        integrate = integrate * space/3.d0

        ! DEALLOCATE(cosgrid1,cosgrid2)
        !DEALLOCATE(int_r)!,int_cos1,int_cos2)
        RETURN
    END FUNCTION
    



    SUBROUTINE diagonalize(matrix,eigenvalue,eigenvector)

        IMPLICIT NONE

        !Matrix
        DOUBLE PRECISION, INTENT(IN)    :: matrix(:,:)
        DOUBLE PRECISION, INTENT(OUT)   :: eigenvalue
        DOUBLE PRECISION, INTENT(OUT), ALLOCATABLE  ::  eigenvector(:)
        INTEGER                         :: n_size

        !Eigenvalue Problem
        DOUBLE PRECISION, ALLOCATABLE   :: rvalues(:), ivalues(:)
        DOUBLE PRECISION, ALLOCATABLE   :: vectors(:,:) ! Matriz para almacenar los autovectores

    
        !Iteratives
        DOUBLE PRECISION, ALLOCATABLE   ::  iv1(:),fv1(:)
        INTEGER                         ::  ierr,idiag
    
        n_size = SIZE(matrix, DIM = 1)
        ierr = 0
        ALLOCATE (rvalues(n_size))
        ALLOCATE (ivalues(n_size))
        ALLOCATE (vectors(n_size,n_size))
        ALLOCATE (iv1(n_size),fv1(n_size))
        
            CALL rg(n_size,n_size,matrix,rvalues,ivalues,1,vectors,iv1,fv1,ierr)!nm,n,a,wr,wi,matz,z,iv1,fv1,ierr

            !Imposition: eigenvalues should not be greater than -2000
            !This might vary for each problem
            !It is necessary, since the laplacian introduces divergences near the origin
            DO idiag = 1,size(rvalues)
                if (rvalues(idiag)<-2000.d0) then
                    rvalues(idiag) = 0.d0
                end if
            end do
            eigenvalue = MINVAL(rvalues)
    
            eigenvector = vectors(:,MINLOC(rvalues,DIM=1))

            DEALLOCATE (rvalues)
            DEALLOCATE (ivalues)
            DEALLOCATE (vectors)
            DEALLOCATE (iv1,fv1)
        END SUBROUTINE

    SUBROUTINE orbital_save(rgrid,phi_init_save)
        DOUBLE PRECISION, INTENT(IN)        ::  phi_init_save(:),rgrid(:)!,cos_int(:,:)
        ! INTEGER, INTENT(IN)                 ::  min_index, max_index
        DOUBLE PRECISION                    ::  coul_phi,lap_phi
        DOUBLE PRECISION, ALLOCATABLE       ::  lap_m(:,:),coul_m(:,:)!lap_phi(:)!,hart_m(:,:),fock_m(:,:)
        INTEGER                             ::  n_size

        n_size = size(phi_init_save)
        ! ALLOCATE(lap_phi(n),coul_phi(n),hart_phi(n),fock_phi(n))
        ALLOCATE(lap_m(n_size,n_size),coul_m(n_size,n_size))!,lap_phi(n))!,hart_m(n,n),fock_m(n,n))
        
        lap_m = laplacian(rgrid)
        lap_phi = integrate(rgrid,phi_init_save*MATMUL(lap_m,phi_init_save))
        ! coul_m = V_coulomb(Z,rgrid)
        ! coul_phi =  integrate(rgrid,phi_init_save*MATMUL(coul_m,phi_init_save))
        ! hart_m = Hartree(r,phi_init,cos_int,1)
        ! hart_phi =  DOT_PRODUCT(phi_init,MATMUL(hart_m,phi_init))
        ! fock_m = Fock(r,phi_init,cos_int,1)
        ! fock_phi =  DOT_PRODUCT(phi_init,MATMUL(fock_m,phi_init))
        OPEN(unit=20,file='results/laplacian_orb.dat',status='replace',action='write')
        OPEN(unit=21,file='results/laplacian_matrix.dat',status='replace',action='write')
        OPEN(unit=30,file='results/coulomb_orb.dat',status='replace',action='write')
        OPEN(unit=31,file='results/coulomb_matrix.dat',status='replace',action='write')
        ! OPEN(unit=40,file='results/hartree_orb.dat',status='replace',action='write')
        ! OPEN(unit=41,file='results/hartree_matrix.dat',status='replace',action='write')
        ! OPEN(unit=50,file='results/fock_orb.dat',status='replace',action='write')
        ! OPEN(unit=51,file='results/fock_matrix.dat',status='replace',action='write')

        WRITE(20,'(f12.8)') lap_phi
        WRITE(30,'(f12.8)') coul_phi
        ! WRITE(40,'(f12.8)') hart_phi
        ! WRITE(50,'(f12.8)') fock_phi

        DO i=1,n_size
            WRITE(21,*) lap_m(i,:)
            WRITE(31,*) coul_m(i,:)
            ! WRITE(41,*) hart_m(i,:)
            ! WRITE(51,*) fock_m(i,:)
        END DO
        CLOSE(20)
        CLOSE(21)
        CLOSE(30)
        CLOSE(31)
        
        DEALLOCATE(lap_m,coul_m)!,lap_phi)!,hart_m,fock_m)
        WRITE(*,*) 'Matrices and expected values saved'
        END SUBROUTINE

        FUNCTION arange(a,b,astep)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN)    ::  a,b,astep
            
            DOUBLE PRECISION, ALLOCATABLE   ::  arange(:)
    
            INTEGER                         ::  iarange,len
            
            ! write(*,*) a,b,step
            len = INT((b-a)/astep)
            
            ALLOCATE(arange(len))
            arange(1) = a
            DO iarange = 2,len
                arange(iarange) = arange(iarange-1) + astep
                
            END DO
            
            RETURN
        END FUNCTION


    FUNCTION Hartree(rgrid,phi1)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)        ::  rgrid(:),phi1(:)
        DOUBLE PRECISION                    ::  r1,rfun
        ! DOUBLE PRECISION, ALLOCATABLE       ::  cosgrid1(:),cosgrid2(:)
        DOUBLE PRECISION, ALLOCATABLE       ::  int_r(:)!,int_cos1(:),int_cos2(:)
        DOUBLE PRECISION, ALLOCATABLE       ::  Hartree(:,:)
        INTEGER                             ::  Npartr!, Npartcos

        ! DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
        INTEGER                             ::  icol, ir
        Npartr = size(rgrid)
        ALLOCATE(int_r(Npartr))
        

        ALLOCATE(Hartree(Npartr,Npartr))
        
        Hartree = 0.d0
        !icol is the index for the Hartree matrix
        !ir is the index for the integration grid
        do icol = 1,Npartr
                do ir = 1,Npartr
                    r1 = rgrid(ir)
                    rfun = phi1(ir)**2 / MAX(r1,rgrid(icol))!*r**2 revisar en funcion del calculo

                    int_r(ir) = rfun
                end do

                Hartree(icol,icol) = integrate(rgrid,int_r) + Hartree(icol,icol)
        end do

        DEALLOCATE(int_r)
        RETURN
    END FUNCTION

end program HF_simp_cccc
