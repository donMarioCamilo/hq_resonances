SUBROUTINE gaussian(n_orbitals,alpha,coeff,l)
    IMPLICIT NONE
    
    !File for loading the GTOs
    !To be improved via modules
    
    INTEGER, INTENT(IN)   :: n_orbitals
    
    
    !STO-3G
    DOUBLE PRECISION, INTENT(OUT)    :: alpha(n_orbitals,1)
    DOUBLE PRECISION, INTENT(OUT)    :: coeff(n_orbitals,1)
    !DOUBLE PRECISION, INTENT(OUT)    :: coord(n_orbitals,3)
    INTEGER, INTENT(OUT)             :: l(n_orbitals,3)
    
    alpha(:,1) = [0.297104, 1.236745, 5.749982, 38.216677]

    coeff(:,1) = [1,1,1,1]
    
    !coord(1,:) = [0., 0.,0.]
    !coord(2,:) = [0., 0., 1.4]
    
    l(1,:) = [0, 0, 0]
    l(2,:) = [0, 0, 0]
    l(3,:) = [0, 0, 0]
    l(4,:) = [0, 0, 0]
    
END SUBROUTINE