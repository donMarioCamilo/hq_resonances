!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM FOR COMPUTING THE GROUND STATE OF A MANY BODY SYSTEM !!
!!  VIA HARTREE FOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!In this demo muonium-like system ground state will be computed
!For same-flavoured mesons, set all masses equal
PROGRAM HF_muonium
    IMPLICIT none
    
    DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
    !Electron mass (eV)
    DOUBLE PRECISION, PARAMETER         ::  mass_charm = 1575.d0
    DOUBLE PRECISION, PARAMETER         ::  mass_bottom = 4768.d0
    !Fine structure constant
    DOUBLE PRECISION, PARAMETER         ::  alpha_charm = 0.34d0
    DOUBLE PRECISION, PARAMETER         ::  alpha_bottom = 0.21d0
    DOUBLE PRECISION, PARAMETER         ::  alpha = sqrt(alpha_charm*alpha_bottom)
    !Bohr radius
    DOUBLE PRECISION, PARAMETER         ::  bohr_r = 1.d0/(sqrt(mass_charm*mass_bottom)*alpha)!Cambiar la masa y la alpha

    !radial coordinates
    DOUBLE PRECISION, ALLOCATABLE       ::  r(:)
    DOUBLE PRECISION, PARAMETER         ::  r0 = 0.001d0*bohr_r, rmax = 1.0d1*bohr_r
    
    !Initial guess function
    DOUBLE PRECISION, ALLOCATABLE       ::  phi_init_charm(:),phi_init_bottom(:)

    !Non-iterable matrix
    DOUBLE PRECISION, ALLOCATABLE       ::  h_charm(:,:),h_bottom(:,:)
    !Iterables
    !Eigenvalue resulting function
    DOUBLE PRECISION, ALLOCATABLE       ::  phi_charm(:),phi_bottom(:), F_matrix_charm(:,:),F_matrix_bottom(:,:)
    DOUBLE PRECISION                    ::  energy_new, energy_old, eps_charm, eps_bottom
    !Parameters
    !Orbital spacing    
    double precision,parameter          ::  gridsize = 1000
    DOUBLE PRECISION, PARAMETER         ::  spacing = (rmax-r0)/gridsize
    !Orbital array length
    integer                             ::  n

    !Max iterations in scf procedure
    INTEGER, PARAMETER                  ::  max_step = 50
    !Energy threshold for convergence
    DOUBLE PRECISION, PARAMETER         ::  energy_threshold = 1.d-3

    !Iterators
    INTEGER                             ::  i,step

    !Checkers
    ! DOUBLE PRECISION                    ::  intcheck

    !Debugging         
    DOUBLE PRECISION, allocatable       ::  fockcheck(:)                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    r = arange(r0,rmax,spacing)
    n = size(r)
    ALLOCATE(fockcheck(n))
    ! allocate(phi_init(n))
    !Initial guess for the orbital
    
    !Initial guess for the orbital
    phi_init_charm = (2.d0*sqrt((mass_charm*alpha_charm)**3)) * r * EXP(-r*(mass_charm*alpha_charm))
    phi_init_bottom = (2.d0*sqrt((mass_bottom*alpha_bottom)**3)) * r * EXP(-r*(mass_bottom*alpha_bottom))

    ! phi_init_mu = 0.d0
    ! phi_init_mu(1) = 1.d0
    ! phi_init_mu = phi_init_mu/sqrt(integrate(r,phi_init_mu**2))
    !Both are written into phi_init.dat for later debuging
    OPEN(unit=1,file='results/phi_init_e.dat',status='replace',action='write')
    OPEN(unit=2,file='results/phi_init_mu.dat',status='replace',action='write')
    DO i=1,n
        WRITE(1,'(2f15.8)') r(i),phi_init_charm(i)
        WRITE(2,'(2f15.8)') r(i),phi_init_bottom(i)
    END DO
    CLOSE(1)
    CLOSE(2)

    !Integration is checked by phi_init normalization
    write(*,*) 'Norma de phi_init_e**2'
    write(*,*) sqrt(integrate(r,phi_init_charm**2))
    write(*,*) 'Norma de phi_init_mu**2'
    write(*,*) sqrt(integrate(r,phi_init_bottom**2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(h_charm(n,n),h_bottom(n,n))
    
    ! CALL orbital_save(r,phi_init)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!         BEGINNING OF THE SCF PROCEDURE           !!!!!!!!!!!!!
    WRITE(*,*)  'Starting SCF procedure'
    h_charm = -0.25d0/mass_charm * laplacian(r)
    h_bottom = -0.25d0/mass_bottom * laplacian(r)
    energy_old = 1000.d0
    !OPEN(unit=12,file='results/scf_matrix.dat',status='replace',action='write')
    OPEN(unit=13,file='results/log.dat',status='replace',action='write')
    OPEN(unit=14,file='results/F_matmul_phi.dat',status='replace',action='write')
    OPEN(unit=15,file='results/phi_e.dat',status='replace',action='write')
    OPEN(unit=16,file='results/phi_mu.dat',status='replace',action='write')
    !OPEN(unit=23,file='results/operator.dat',status='replace',action='write')
    WRITE(13,*) 'N = ', n
    WRITE(13,*) 'r0, rmax = ', r0,rmax
    
    DO step = 1,max_step
        ALLOCATE(F_matrix_charm(n,n),F_matrix_bottom(n,n))
        F_matrix_bottom = h_bottom - 4.d0/3.d0*alpha*Hartree(r,phi_init_charm)
        CALL diagonalize(F_matrix_bottom,eps_bottom,phi_bottom)
        phi_init_bottom = phi_bottom/sqrt(integrate(r,phi_bottom**2))


        F_matrix_charm = h_charm - 4.d0/3.d0*alpha*Hartree(r,phi_init_bottom)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!         ORBITAL CHECK           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fockcheck = MATMUL(F_matrix_charm,phi_init_charm)

        fockcheck = fockcheck/sqrt(integrate(r,fockcheck**2))
        WRITE(14,*) fockcheck(:)
        ! ! WRITE(23,*) integrate(r,phi_init*fockcheck)+integrate(r,phi_init*MATMUL(h,phi_init))
        ! IF(step==1) then
        !     DO i = 1,n
        !         WRITE(12,*) F_matrix_e(i,:)
        !     END DO
        ! end if  
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL diagonalize(F_matrix_charm,eps_charm,phi_charm)
        write(*,*) eps_bottom,eps_charm
        
        phi_init_charm = phi_charm/sqrt(integrate(r,phi_charm**2))

        energy_new = 0.5d0*(eps_charm+eps_bottom + &
                integrate(r,phi_init_charm*MATMUL(h_charm,phi_init_charm)) + &
                integrate(r,phi_init_bottom*MATMUL(h_bottom,phi_init_bottom)))

        
        WRITE(13,*) 'Iteration, Energy:'
        WRITE(13,*) step, energy_new
        WRITE(13,*) 'Charm Energy:',0.5d0*(eps_charm + integrate(r,phi_init_charm*MATMUL(h_charm,phi_init_charm)))
        WRITE(13,*) 'Bottom Energy:',0.5d0*(eps_bottom + integrate(r,phi_init_bottom*MATMUL(h_bottom,phi_init_bottom)))

        IF (ABS(energy_new-energy_old)/ABS(energy_old) < energy_threshold) THEN
            WRITE(*,*) 'Convergence met in',step,'steps'
            exit
        END IF
        IF (step == max_step) THEN
            WRITE(*,*) 'WARNING: CONVERGENCE NOT MET. Last two energy values are:'
            WRITE(*,'(2f8.4)')  energy_old, energy_new
        END IF
        
        ! Solution orbital is saved in phi.dat
        WRITE(15,*) phi_init_charm(:)
        WRITE(16,*) phi_init_bottom(:)

        energy_old = energy_new
        WRITE(*,*) 'Step:',step,'completed with energy (MeV):', energy_new
        
        deallocate(F_matrix_charm,F_matrix_bottom)
    END DO
    ! write(*,*) integrate(r,phi_init,MATMUL(F_matrix,phi_init))
    ! stop
    !energy_new = 0.5d0*(energy_new + integrate(r,phi_init_e*MATMUL(h,phi_init_e)))
    
    WRITE(13,*)'Ground State Energy (MeV):', energy_new
    WRITE(*,*)
    WRITE(*,*)'Ground State Energy (MeV):', energy_new

    WRITE(13,*)'Ground State Mass (MeV):', energy_new+1.d0*mass_charm+1.d0*mass_bottom
    WRITE(*,*)
    WRITE(*,*)'Ground State Mass (MeV):', energy_new+1.d0*mass_charm+1.d0*mass_bottom
    ! WRITE(*,*)'Ground State Energy:', 2.d0*integrate(r,phi_init*MATMUL(h,phi_init))+2.d0* &
    !     integrate(r,phi_init*MATMUL(Hartree(r,phi_init,1),phi_init))-integrate(r,phi_init*&
    !     MATMUL(Fock(r,phi_init,weights,1),phi_init))
    close(12)
    close(13)
    close(14)
    close(15)


    CONTAINS

        
    SUBROUTINE get_grid(a,b,n_size,r_out,weight_out)
        IMPLICIT NONE

        INTEGER, INTENT(IN)             ::  n_size
        DOUBLE PRECISION, INTENT(IN)    ::  a,b

        INTEGER                         ::  ntrash
        DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT)   ::  r_out(:),weight_out(:)

        ALLOCATE(r_out(20*n_size),weight_out(20*n_size))
        CALL DSG20Rww(a,b,n_size,r_out,weight_out,ntrash)
    END SUBROUTINE

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

    FUNCTION V_coulomb(Z_coulomb,x)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN)        ::  Z_coulomb
        DOUBLE PRECISION, INTENT(IN)        ::  x(:)
        DOUBLE PRECISION, ALLOCATABLE       ::  V_coulomb(:,:)

        INTEGER     ::  i_cou   !Iterator
        

        ALLOCATE(V_coulomb(size(x),size(x)))
        !The matrix is initialized
        V_coulomb = 0.0d0
        
        DO i_cou = 1,size(x)
            V_coulomb(i_cou,i_cou) = -Z_coulomb*1.d0/(x(i_cou))
        END DO

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
        OPEN(unit=21,file='results/integrate_debug.dat',status='replace',action='write')
        write(21,*) 'Norm Val:',integrate
                close(21)
        RETURN
    END FUNCTION


  
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
    ! INTEGER, INTENT(IN)                 ::  l
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
    
    do icol = 1,Npartr
            do ir = 1,Npartr
                r1 = rgrid(ir)
                rfun = phi1(ir)**2 / MAX(r1,rgrid(icol))!*r**2 revisar en funcion del calculo

                int_r(ir) = rfun
                
            end do

            Hartree(icol,icol) = integrate(rgrid,int_r)
            
    end do

    DEALLOCATE(int_r)
    RETURN
END FUNCTION

    ! FUNCTION Fock(rgrid,phi1,w,l)
    !     IMPLICIT NONE
    !     !Grid,orbital and weights
    !     DOUBLE PRECISION, INTENT(IN)        ::  rgrid(:),phi1(:),w(:)
    !     INTEGER, INTENT(IN)                 ::  l

    !     DOUBLE PRECISION, ALLOCATABLE       ::  Fock(:,:)
    !     INTEGER                             ::  Npartr

    !     ! DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
    !     INTEGER                             ::  icol,irow, il!, ntrash

    !     Npartr = size(rgrid)/20

    !     ALLOCATE(Fock(size(rgrid),size(rgrid)))
    !     Fock = 0.d0
    !     do irow = 1,size(rgrid)
    !         do icol = 1,size(rgrid)
    !             do il = 1,l
    !                 Fock(irow,icol) = Fock(irow,icol) + w(icol)*phi1(icol)*phi1(irow)/MAX(rgrid(irow),rgrid(icol))
    !             end do
    !         end do
    !     end do

    !     RETURN
    ! END FUNCTION
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
            ! do i=1,n
            !     write(*,*) rvalues(i)
            ! end do
            DO idiag = 1,size(rvalues)
                if (abs(rvalues(idiag))>1000.0d0) then
                    rvalues(idiag) = 100.d0
                end if
            end do
            eigenvalue = MINVAL(rvalues)!(min_index:max_index))
    
            eigenvector = vectors(:,MINLOC(rvalues,DIM=1))

            !write(*,*) 'VALUE:',eigenvalue
            OPEN(unit=22,file='results/eigenvectors.dat',status='replace',action='write')
            do idiag=1,n_size
                write(22,*) vectors(:,idiag)
                ! write(*,*) rvalues(i)
            end do
            close(22)
            OPEN(unit=24,file='results/eigenvalues.dat',status='replace',action='write')
            do idiag=1,n_size
                write(24,*) rvalues(idiag)
            end do
            close(24)
            ! OPEN(unit=23,file='results/foundeigen.dat',status='replace',action='write')
            ! DO i = 1,n_size
            !     foundeigen = vectors(:,i)
            !     foundeigen = foundeigen(min_index:max_index)
            !     write(23,*) foundeigen/sqrt(integrate(r0,rmax,r,foundeigen,foundeigen))
            ! END DO
            ! close(23)
            ! write(*,*) 'Eigenvalue for the 91st eigenvector:'
            ! write(*,*) rvalues(91)
            
            ! write(*,*) size(rvalues)
            ! write(*,*) MINLOC(rvalues,DIM=1)
            DEALLOCATE (rvalues)
            DEALLOCATE (ivalues)
            DEALLOCATE (vectors)
            DEALLOCATE (iv1,fv1)
        END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!            DEBUGGING       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FUNCTION int_check(rgrid,phi,check)
    !     IMPLICIT NONE
    !     INTEGER,intent(in)                  ::  check
    !     DOUBLE PRECISION, INTENT(IN)        ::  rgrid(:),phi(:)
    !     DOUBLE PRECISION                    ::  r,costheta1,costheta2,rfun
    !     DOUBLE PRECISION, ALLOCATABLE       ::  cosgrid1(:),cosgrid2(:)
    !     DOUBLE PRECISION, ALLOCATABLE       ::  int_r(:),int_cos(:)
    !     DOUBLE PRECISION                    ::  int_check
    !     DOUBLE PRECISION                    ::  auxr,auxcos,auxcos1,auxcos2

    !     INTEGER                             ::  Npartr, Npartcos

    !     DOUBLE PRECISION, PARAMETER         ::  pi = 3.14159265d0
    !     INTEGER                             ::  ir, icos, ntrash

    !     IF (size(rgrid) > 2000) THEN
    !         WRITE(*,*)  'WARNING: dobint.f is not prepared for', &
    !          'computing more than 2000 gauss points'
    !         WRITE(*,*)
    !     END IF

    !     Npartcos = 10
    !     ALLOCATE(cosgrid1(20*Npartcos),cosgrid2(20*Npartcos))
    !     CALL DSG20R(-1.d0,-1d-1,Npartcos,cosgrid1,Ntrash)
    !     CALL DSG20R(1.d-1,1.d0,Npartcos,cosgrid2,Ntrash)
    !     ALLOCATE(int_r(size(rgrid)), int_cos(20*Npartcos))
        
    !     Npartr = size(rgrid)/20

    !     do ir = 1,20*Npartr
    !         r = rgrid(ir)
    !         rfun = r**2 * phi(ir)**2
    !         do icos = 1,20*Npartcos
    !             costheta1 = cosgrid1(icos)
    !             costheta2 = cosgrid2(icos)
    !             int_cos(icos) = 1.d0
    !         end do

    !         CALL DRG20R(cosgrid1(1),cosgrid1(size(cosgrid1)),Npartcos,int_cos,auxcos1)
    !         CALL DRG20R(cosgrid2(1),cosgrid2(size(cosgrid2)),Npartcos,int_cos,auxcos2)
    !         auxcos = auxcos1 + auxcos2
    !         int_r(ir) = auxcos * rfun
    !     end do

    !     CALL DRG20R(rgrid(1),rgrid(size(rgrid)),Npartr,int_r,auxr)
    !     int_check = 2*pi*auxr
    !     IF (check==1) THEN
    !         WRITE(*,*) 'Numerical integration checked. Precision(%):', 100-100*abs(int_check-1)
    !     END IF
    !     IF ((check==1).or.(check==2)) THEN
    !         IF (ABS(int_check-1.d0) > 5.d-2) THEN
    !             WRITE(*,*) 'Problem while checking phi_init normalization'
    !             WRITE(*,*) 'Please check the integration method.'
    !         END IF
    !         WRITE(*,*)
    !     END IF
    !     DEALLOCATE(cosgrid1,cosgrid2)
    !     DEALLOCATE(int_r,int_cos)
    !     RETURN
    ! END FUNCTION

    ! SUBROUTINE orbital_save(rgrid,phi_init_save)
    !     DOUBLE PRECISION, INTENT(IN)        ::  phi_init_save(:),rgrid(:)!,cos_int(:,:)
    !     ! INTEGER, INTENT(IN)                 ::  min_index, max_index
    !     DOUBLE PRECISION                    ::  coul_phi,lap_phi
    !     DOUBLE PRECISION, ALLOCATABLE       ::  lap_m(:,:),coul_m(:,:)!lap_phi(:)!,hart_m(:,:),fock_m(:,:)
    !     INTEGER                             ::  n_size

    !     n_size = size(phi_init_save)
    !     ! ALLOCATE(lap_phi(n),coul_phi(n),hart_phi(n),fock_phi(n))
    !     ALLOCATE(lap_m(n_size,n_size),coul_m(n_size,n_size))!,lap_phi(n))!,hart_m(n,n),fock_m(n,n))
        
    !     lap_m = laplacian(rgrid)
    !     lap_phi = integrate(rgrid,phi_init_save*MATMUL(lap_m,phi_init_save))
    !     coul_m = V_coulomb(Z,rgrid)
    !     coul_phi =  integrate(rgrid,phi_init_save*MATMUL(coul_m,phi_init_save))
    !     ! hart_m = Hartree(r,phi_init,cos_int,1)
    !     ! hart_phi =  DOT_PRODUCT(phi_init,MATMUL(hart_m,phi_init))
    !     ! fock_m = Fock(r,phi_init,cos_int,1)
    !     ! fock_phi =  DOT_PRODUCT(phi_init,MATMUL(fock_m,phi_init))
    !     OPEN(unit=20,file='results/laplacian_orb.dat',status='replace',action='write')
    !     OPEN(unit=21,file='results/laplacian_matrix.dat',status='replace',action='write')
    !     OPEN(unit=30,file='results/coulomb_orb.dat',status='replace',action='write')
    !     OPEN(unit=31,file='results/coulomb_matrix.dat',status='replace',action='write')
    !     ! OPEN(unit=40,file='results/hartree_orb.dat',status='replace',action='write')
    !     ! OPEN(unit=41,file='results/hartree_matrix.dat',status='replace',action='write')
    !     ! OPEN(unit=50,file='results/fock_orb.dat',status='replace',action='write')
    !     ! OPEN(unit=51,file='results/fock_matrix.dat',status='replace',action='write')

    !     WRITE(20,'(f12.8)') lap_phi
    !     WRITE(30,'(f12.8)') coul_phi
    !     ! WRITE(40,'(f12.8)') hart_phi
    !     ! WRITE(50,'(f12.8)') fock_phi

    !     DO i=1,n_size
    !         WRITE(21,*) lap_m(i,:)
    !         WRITE(31,*) coul_m(i,:)
    !         ! WRITE(41,*) hart_m(i,:)
    !         ! WRITE(51,*) fock_m(i,:)
    !     END DO
    !     CLOSE(20)
    !     CLOSE(21)
    !     CLOSE(30)
    !     CLOSE(31)
    !     ! CLOSE(40)
    !     ! CLOSE(41)
    !     ! CLOSE(50)
    !     ! CLOSE(51)
    !     ! DEALLOCATE(lap_phi,coul_phi,hart_phi,fock_phi)
    !     DEALLOCATE(lap_m,coul_m)!,lap_phi)!,hart_m,fock_m)
    !     WRITE(*,*) 'Matrices and expected values saved'
    !     END SUBROUTINE

end program HF_muonium
