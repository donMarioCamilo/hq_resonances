SUBROUTINE matrix_calculator(n,n_atoms,Z,alpha,coeff,coord)
    IMPLICIT NONE

    !Constants
    DOUBLE PRECISION, PARAMETER ::  pi = 4.D0*DATAN(1.D0)
    
    INTEGER, INTENT(IN)             :: n
    INTEGER, INTENT(IN)             :: n_atoms
    INTEGER, INTENT(IN)             :: Z(n_atoms)
    DOUBLE PRECISION, INTENT(IN)    :: alpha(n,3)
    DOUBLE PRECISION, INTENT(IN)    :: coeff(n,3)
    DOUBLE PRECISION, INTENT(IN)    :: coord(n,3)
    
    !Matrix
    !DOUBLE PRECISION, ALLOCATABLE   :: S_matrix(:,:)!, K_matrix(:,:), V_matrix(:,:), Two_e(:,:,:,:)
    !preguntar si tiene sentido no definir matrices
    
    DOUBLE PRECISION   :: S_val, S_element, K_element, V_ne ! S_matrix, K_matrix, V_nuclear_electron_atraction_matrix, Two_e
    
    !Iteratives
    INTEGER ::  p,q,r,s, iread, atom_i
    
    !For calculations
    DOUBLE PRECISION        :: c1c2, alphasum, alphaprod, expfactor, normfactor, distance(3)  !For calculations
    DOUBLE PRECISION        :: alphacoord(3), PG(3), PG2, PGx2, PGy2, PGz2, PG_ne(3), PG2_ne   !For kinetic energy and nuclear electron atraction
    
    
    OPEN(unit = 1, file = 'kinetic.dat', status = 'replace', action = 'write')  !Kinetic Matrix
    OPEN(unit = 2, file = 'nuclear_atraction.dat', status = 'replace', action = 'write')  !Potential Matrix
    OPEN(unit = 3, file = 'overlap.dat', status = 'replace', action = 'write')  !Overlap Matrix
    !OPEN(unit = 4, file = 'two_electron.dat', status = 'replace', action = 'write')  !Two electron Matrix
    

    !Overlap, Hamilton, Q and F matrix elements
    DO p = 1,n
        DO q = 1,n
            S_element = 0.0
            K_element = 0.0
            V_ne = 0.0
            DO r = 1, SIZE(alpha(p,:))
                DO s = 1, SIZE(alpha(q,:))
                    !Factors for calculations
                    c1c2 = coeff(p,r)*coeff(q,s)
                    alphasum = alpha(p,r)+alpha(q,s)
                    alphaprod = alpha(p,r)*alpha(q,s)
                    expfactor = alphaprod/alphasum
                    distance = coord(p,:)-coord(q,:)
                    normfactor = (2*alpha(p,r)/pi)**0.75*(2*alpha(q,s)/pi)**0.75
                    
                    !Factors for kinetic energy
                    alphacoord = alpha(p,r)*coord(p,:) + alpha(q,s)*coord(q,:)
                    PG = alphacoord/alphasum - coord(q,:)
                    PG2 = SUM(PG*PG)
                    PGx2 = PG(1)*PG(1)
                    PGy2 = PG(2)*PG(2)
                    PGz2 = PG(3)*PG(3)
                    
                    S_val = normfactor * c1c2 * (pi/alphasum)**(3.0/2.0) * EXP(-expfactor*SUM(distance*distance))
                
                    S_element = S_element + S_val
                    
                    K_element = K_element + 3.0*alpha(q,s)*S_val
                    K_element = K_element - 2.0*alpha(q,s)*alpha(q,s)*S_val * (PGx2 + 0.5/alphasum)
                    K_element = K_element - 2.0*alpha(q,s)*alpha(q,s)*S_val * (PGy2 + 0.5/alphasum)
                    K_element = K_element - 2.0*alpha(q,s)*alpha(q,s)*S_val * (PGz2 + 0.5/alphasum)
                    
                    !Nuclear electron atraction
                    DO atom_i = 1,n_atoms
                        PG_ne = alphacoord/alphasum - coord(atom_i,:)
                        PG2_ne = SUM(PG_ne*PG_ne)
                        V_ne = V_ne - Z(atom_i) * normfactor * c1c2 * EXP(-expfactor*SUM(distance*distance)) *&
                             (2.0*pi/alphasum) * F0(alphasum*PG2_ne)
                    END DO
                END DO
            END DO
            WRITE (1,'(i3, i3, f8.4)') p, q, K_element
            WRITE (2,'(i3, i3, f8.4)') p, q, V_ne
            WRITE (3,'(i3, i3, f8.4)') p, q, S_element
        END DO
    END DO
    !WRITE (*,*) 'Overlap matrix S:'
    !DO iread=1,n
    !    WRITE (*,'(4f8.4)') S_matrix(iread,:)
    !END DO
    
    
    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    !CLOSE(4)
    
    !DEALLOCATE(S_matrix)
    !DEALLOCATE(K_matrix)
    !DEALLOCATE(V_matrix)
    !DEALLOCATE(Two_e)
    
    
    !WRITE(*,'(f10.8)') lower_incomplete_gamma(1.5,1.0d0)
    CONTAINS 
        FUNCTION F0(x)
                DOUBLE PRECISION, INTENT(IN)    ::  x
                DOUBLE PRECISION                ::  F0
                DOUBLE PRECISION, PARAMETER     ::  pi = 4.D0*DATAN(1.D0)
            
                IF (ABS(x) < 1E-5) THEN
                    F0 = 1.0D0
                    RETURN
                END IF
                F0 = x**(-0.5) * SQRT(pi)/2 * ERF(x**(0.5))
                RETURN
            END FUNCTION

    END SUBROUTINE