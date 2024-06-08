SUBROUTINE gaussian(n_orbitals,alpha,coeff,l)
    IMPLICIT NONE
    
    !File for loading the GTOs
    !To be improved via modules
    
    INTEGER, INTENT(IN)   :: n_orbitals
    
    
    !STO-3G
    DOUBLE PRECISION, INTENT(OUT)    :: alpha(n_orbitals,3)
    DOUBLE PRECISION, INTENT(OUT)    :: coeff(n_orbitals,3)
    !DOUBLE PRECISION, INTENT(OUT)    :: coord(n_orbitals,3)
    INTEGER, INTENT(OUT)             :: l(n_orbitals,3)
    
    alpha(1,:) = [0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00]
    alpha(2,:) = [0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00]

    coeff(1,:) = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00]
    coeff(2,:) = [0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00]
    
    !coord(1,:) = [0., 0.,0.]
    !coord(2,:) = [0., 0., 1.4]
    
    l(1,:) = [0, 0, 0]
    l(2,:) = [0, 0, 0]
    
END SUBROUTINE