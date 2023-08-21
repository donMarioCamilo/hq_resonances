PROGRAM main
    IMPLICIT none
    
    ! Helium Hartree-Fock
    
    ! Parameters
    INTEGER, PARAMETER  :: n_orbitals = 4
    INTEGER, PARAMETER  :: n_atoms = 1
    INTEGER, PARAMETER  :: Z(n_atoms) = [2]
    INTEGER, PARAMETER  :: n_occupation = 1
    INTEGER, PARAMETER  :: max_it = 30
    DOUBLE PRECISION, PARAMETER :: scf_threshold = 1E-5
    
    DOUBLE PRECISION    :: alpha(n_orbitals,1)
    DOUBLE PRECISION    :: coeff(n_orbitals,1)
    DOUBLE PRECISION    :: coord(n_orbitals,3)
    INTEGER             :: l(n_orbitals,3)
    DOUBLE PRECISION    :: nuc_nuc_repulsion_E
    
    DOUBLE PRECISION    :: S_matrix(n_orbitals,n_orbitals),&
                           K_matrix(n_orbitals,n_orbitals),&
                           nuc_atraction(n_orbitals,n_orbitals),&
                           Two_e(n_orbitals,n_orbitals,n_orbitals,n_orbitals)
    DOUBLE PRECISION    :: electronic_energy, energy
    
    DOUBLE PRECISION    :: distance_iterator, distance_i
    
    INTEGER             :: iread
    
    nuc_nuc_repulsion_E = 0.0
    
    CALL gaussian(n_orbitals, alpha, coeff, l)
    
    distance_iterator = 0.0d0
    
    
    coord(1,:) = [0.0d0, 0.0d0,0.0d0]
    coord(2,:) = [0.0d0, 0.0d0,0.0d0]
    coord(3,:) = [0.0d0, 0.0d0,0.0d0]
    coord(4,:) = [0.0d0, 0.0d0,0.0d0]
    !Uncomment next lines for calculating the relevant matrices. Check comments.txt
    CALL matrix_calculator(n_orbitals,n_atoms,Z,alpha,coeff,coord)
    !CALL V_ne_calculator(n_orbitals,n_atoms,Z,alpha,coeff,coord) !not necessary if matrix_calculator is included
    CALL Two_e_calculator(n_orbitals,alpha,coeff,coord)    
    
    !Get matrices from the .dat files
    S_matrix = give_matrix(n_orbitals,'overlap.dat')
    K_matrix = give_matrix(n_orbitals,'kinetic.dat')
    nuc_atraction = give_matrix(n_orbitals,'nuclear_atraction.dat')
    Two_e = give_4matrix(n_orbitals,"two_electron.dat")
    IF (n_atoms > 1) THEN
        nuc_nuc_repulsion_E = nuclear_nuclear_repulsion(n_atoms,Z,coord)
    END IF
    
    electronic_energy = scf_loop(n_orbitals, n_occupation, S_matrix,&
        K_matrix, nuc_atraction, Two_e, scf_threshold, max_it)
    energy = electronic_energy + nuc_nuc_repulsion_E
    WRITE(*,*) 'Energy:'
    WRITE(*,'(f8.4)') energy

    
    
    
    CONTAINS
    
    FUNCTION give_matrix(n_orbitals,name)
            IMPLICIT NONE
            
            INTEGER, INTENT(IN)     :: n_orbitals
            CHARACTER(len=*), INTENT(IN)   :: name
            DOUBLE PRECISION        :: give_matrix(n_orbitals,n_orbitals)
            
            INTEGER                 :: i, j, iread, jread
            DOUBLE PRECISION        :: element
            
            OPEN(unit = 10, file = name, status = 'OLD', action = 'READ')
            
            DO iread = 1,n_orbitals
                DO jread = 1,n_orbitals
                    READ(10, '(i3, i3, f8.4)') i,j, element
                    give_matrix(i,j) = element
                END DO
            END DO
            CLOSE(10)
            RETURN
    END FUNCTION
    
    FUNCTION give_4matrix(n_orbitals,name)
            IMPLICIT NONE
            
            INTEGER, INTENT(IN)     :: n_orbitals
            CHARACTER(len=*), INTENT(IN)   :: name
            DOUBLE PRECISION        :: give_4matrix(n_orbitals,n_orbitals,&
                    n_orbitals,n_orbitals)
            
            INTEGER                 :: i, j, k, l, iread, jread, kread, lread
            DOUBLE PRECISION        :: element
            
            OPEN(unit = 15, file = name, status = 'OLD', action = 'READ')
            
            DO iread = 1,n_orbitals
                DO jread = 1,n_orbitals
                    DO kread = 1,n_orbitals
                        DO lread = 1,n_orbitals
                        READ(15, '(i3, i3, i3, i3, f8.4)') i,j,k,l, element
                            give_4matrix(i,j,k,l) = element
                        END DO
                    END DO
                END DO
            END DO
            CLOSE(15)
            RETURN
    END FUNCTION
    
    FUNCTION nuclear_nuclear_repulsion(n_atoms,Z,coord)
            IMPLICIT NONE
            
            INTEGER, INTENT(IN)             :: n_atoms
            INTEGER, INTENT(IN)             :: Z(n_atoms)
            DOUBLE PRECISION, INTENT(IN)    :: coord(n_atoms,3)
            
            DOUBLE PRECISION                :: Rx2,Ry2,Rz2
            DOUBLE PRECISION                :: nuc_distance
            DOUBLE PRECISION                :: nuclear_nuclear_repulsion
            
            !ITERATIVES
            INTEGER                         :: p,q
            
            nuclear_nuclear_repulsion = 0.0
            
            DO p = 1,n_atoms
                DO q = 1,n_atoms
                    IF (q > p) THEN
                        Rx2 = (coord(p,1) - coord(q,1))**2
                        Ry2 = (coord(p,2) - coord(q,2))**2
                        Rz2 = (coord(p,3) - coord(q,3))**2
                        nuc_distance = SQRT(Rx2+Ry2+Rz2)
                        nuclear_nuclear_repulsion = nuclear_nuclear_repulsion + (Z(p)*Z(q))/nuc_distance
                    END IF
                END DO
            END DO
            RETURN
    END FUNCTION
            
    
    FUNCTION scf_loop(n_orbitals, n_occupation, S_matrix, K_matrix, nuc_atraction, Two_e, threshold, max_it) RESULT(energy_new)
            IMPLICIT NONE
            
            INTEGER, INTENT(IN)                 ::  n_orbitals
            INTEGER, INTENT(IN)                 ::  n_occupation
            INTEGER, INTENT(IN)                 ::  max_it
            DOUBLE PRECISION, INTENT(IN)        ::  threshold
            
            DOUBLE PRECISION, INTENT(IN)        ::  S_matrix(n_orbitals,n_orbitals), K_matrix(n_orbitals,n_orbitals), &
                                                    nuc_atraction(n_orbitals,n_orbitals)
            DOUBLE PRECISION, INTENT(IN)        ::  Two_e(n_orbitals,n_orbitals,n_orbitals,n_orbitals)
            
            DOUBLE PRECISION                    ::  energy_new, energy_old
            DOUBLE PRECISION                    ::  density_matrix(n_orbitals,n_orbitals)
            DOUBLE PRECISION                    ::  G_matrix(n_orbitals,n_orbitals) !Coulomb and exchange contributions to Fock Matrix
            DOUBLE PRECISION                    ::  Fock_matrix(n_orbitals,n_orbitals)
            DOUBLE PRECISION                    ::  X_trans(n_orbitals,n_orbitals)
            DOUBLE PRECISION                    ::  Fock_unitary(n_orbitals,n_orbitals)
            DOUBLE PRECISION                    ::  molecular_orbitals(n_orbitals,n_orbitals)
            
            DOUBLE PRECISION                    ::  eigen_values(n_orbitals), eigen_vectors(n_orbitals,n_orbitals)
            
            INTEGER                             ::  it, orbital_it, iter, i, j
            energy_new = 0.0d0
            density_matrix = 0.0d0
            
            DO it = 1,max_it
                energy_old = energy_new
                G_matrix = G_calculator(n_orbitals,density_matrix,Two_e)
                !DO i=1,n_orbitals
                !    WRITE(*,'(2f8.4)') G_matrix(:,i)
                !END DO
                eigen_vectors = 0.0d0
                Fock_matrix = K_matrix + nuc_atraction + G_matrix
                !WRITE(*,*) 'Fock matrix:'
                !DO iter = 1,n_orbitals
                !    WRITE(*,*) Fock_matrix(:,iter)
                !END DO
                !Diagonalise Overlap matrix
                X_trans = diagonalize('T',n_orbitals,S_matrix)
                !WRITE(*,*) 'X_trans_'
                !DO iter = 1,n_orbitals
                !    WRITE(*,*) X_trans(:,iter)
                !END DO
                
                Fock_unitary = MATMUL(TRANSPOSE(X_trans),(MATMUL(Fock_matrix,X_trans)))
                
                !WRITE(*,*) '0Fock_unitary:'
                !DO orbital_it = 1, n_orbitals
                !    WRITE(*,*) Fock_unitary(:,orbital_it)
                !END DO
                
                eigen_vectors = diagonalize('F',n_orbitals,Fock_unitary)
                
                !WRITE(*,*) '1X_trans:'
                !DO orbital_it = 1, n_orbitals
                !    WRITE(*,*) X_trans(orbital_it,:)
                !END DO
                !
                !WRITE(*,*) '2eigen_vectors:'
                !DO orbital_it = 1, n_orbitals
                !    WRITE(*,'(i3,3f8.4)') orbital_it, eigen_vectors(orbital_it,:)
                !END DO
                !WRITE (*,*) 'matmul'
                !WRITE(*,*) matmul(X_trans,eigen_vectors(:,1))
                !WRITE(*,*) 'Molecular:'
                molecular_orbitals = MATMUL(X_trans,eigen_vectors)
                !DO orbital_it = 1, n_orbitals
                !    molecular_orbitals(:,orbital_it) = MATMUL(X_trans,eigen_vectors(orbital_it,:))
                !END DO
                
                density_matrix = density_matrix_calculator(n_orbitals, n_occupation, molecular_orbitals)
                !WRITE(*,*) '3mos matrix:'
                !DO i=1,n_orbitals
                !    WRITE(*,'(2f8.4)') molecular_orbitals(i,:)
                !END DO
                !WRITE(*,*) '4Density matrix:'
                !DO i=1,n_orbitals
                !    WRITE(*,'(2f8.4)') density_matrix(i,:)
                !END DO
                energy_new = 0.0d0
                DO i = 1, n_orbitals
                    DO j = 1, n_orbitals
                        energy_new = energy_new + density_matrix(i,j) * (K_matrix(i,j) + nuc_atraction(i,j) + 0.5*G_matrix(i,j))
                    END DO
                END DO
                
                        
                IF (ABS(energy_new - energy_old) < threshold) THEN
                  WRITE(*,*) 'Convergence met'
                  !WRITE(*,'(f8.4)') energy_new
                  RETURN
                END IF
    
            END DO
            WRITE(*,*) 'CONVERGENCE NOT MET. Last two energy values are:'
            WRITE(*,'(2f8.4)')  energy_old, energy_new
            RETURN
    END FUNCTION
    
    FUNCTION G_calculator(n, P, ee_matrix)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)                     ::  n
        DOUBLE PRECISION, INTENT(IN)            ::  P(n,n)
        DOUBLE PRECISION, INTENT(IN)            ::  ee_matrix(n,n,n,n)
    
        DOUBLE PRECISION                        ::  G_calculator(n,n)
        DOUBLE PRECISION                        ::  J_term, K_term
        
        INTEGER                                 ::  i,j,k,l
        
        G_calculator = 0.0d0
        
        DO i = 1,n
            DO j = 1,n
                DO k = 1,n
                    DO l = 1,n
                        J_term = ee_matrix(i,j,k,l)
                        K_term = ee_matrix(i,l,k,j)
                        G_calculator(i,j) = G_calculator(i,j) + P(k,l)*(J_term - 0.5*K_term)
                    END DO
                END DO
            END DO
        END DO
        
        RETURN
    END FUNCTION
    
    FUNCTION density_matrix_calculator(n,n_occupation, molecular_orbitals) RESULT(P)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)                 ::  n !This should be the n of molecular orbitals, not atomic!!!!
        INTEGER, INTENT(IN)                 ::  n_occupation
        DOUBLE PRECISION, INTENT(IN)        ::  molecular_orbitals(n,n)
        
        
        DOUBLE PRECISION                    ::  P(n,n)
        DOUBLE PRECISION                    ::  C, C_dagger
        INTEGER                             ::  i,j, orbital_i, occupation
        
        occupation = 2
        P = 0.0d0
        
        DO i = 1,n
            DO j = 1,n
                DO orbital_i = 1,n_occupation
                    C = molecular_orbitals(i,orbital_i)
                    C_dagger = molecular_orbitals(j,orbital_i)
                    P(i,j) = P(i,j) + occupation * C * C_dagger
                END DO
                !WRITE(*,*) 'C'
                !WRITE(*,*) C
                !WRITE(*,*) C_dagger
                !WRITE(*,*) P(i,j)
            END DO
        END DO
        
        RETURN
    
    END FUNCTION
    
    FUNCTION diagonalize(transform_out,N,matrix)
    
        IMPLICIT NONE

        INTEGER, INTENT(IN)             :: N
        CHARACTER(len=*), INTENT(IN)    :: transform_out
        !Matrix
        DOUBLE PRECISION, INTENT(IN)    :: matrix(N,N)
        DOUBLE PRECISION, ALLOCATABLE   :: matrixdiag(:,:)
        DOUBLE PRECISION, ALLOCATABLE   :: s_halfinv(:,:)
   
        !Eigenvalue Problem
        DOUBLE PRECISION, ALLOCATABLE   :: rvalues(:), ivalues(:), dummyvectors(:,:)
        DOUBLE PRECISION                :: vectors(N,N) ! Matriz para almacenar los autovectores
        DOUBLE PRECISION                :: X_trans(N,N)
        DOUBLE PRECISION                :: diagonalize(N,N)
        DOUBLE PRECISION, ALLOCATABLE   :: work(:)
    
        !Iteratives
        INTEGER ::  it, i, j, info, lwork
    

        lwork=8*N

        ALLOCATE (s_halfinv(N,N))
        ALLOCATE (rvalues(N))
        ALLOCATE (ivalues(N))
        ALLOCATE (dummyvectors(N,N))
        ALLOCATE (work(lwork))
        ALLOCATE (matrixdiag(N,N))
    
        matrixdiag=matrix
    
        IF (transform_out == 'F') THEN
            CALL DSYEV('V', 'U', N, matrixdiag, N, rvalues, work, lwork, info)!, ivalues, dummyvectors, N, vectors, N, work, lwork, info)

            DO i = 1,N
                vectors(i,:) = matrixdiag(i,:)
            END DO
    
            diagonalize = vectors
            DEALLOCATE (s_halfinv)
            DEALLOCATE (rvalues)
            DEALLOCATE (ivalues)
            DEALLOCATE (dummyvectors)
            DEALLOCATE (work)
            DEALLOCATE (matrixdiag)
            RETURN
        
        ELSEIF (transform_out == 'T') THEN
            CALL DGEEV('N', 'V', N, matrixdiag, N, rvalues, ivalues, dummyvectors, N, vectors, N, work, lwork, info)

            s_halfinv=0.0
            DO i=1,N
                s_halfinv(i,i)=(matrixdiag(i,i))**(-0.5)
            END DO
    
            !transformation matrix
            X_trans=MATMUL(vectors,MATMUL(s_halfinv,TRANSPOSE(vectors)))
            
            diagonalize = X_trans
            DEALLOCATE (s_halfinv)
            DEALLOCATE (rvalues)
            DEALLOCATE (ivalues)
            DEALLOCATE (dummyvectors)
            DEALLOCATE (work)
            DEALLOCATE (matrixdiag)
            RETURN
        END IF
        
    END FUNCTION
END PROGRAM
