SUBROUTINE Two_e_calculator(n,alpha,coeff,coord)
    IMPLICIT NONE

    !Constants
    DOUBLE PRECISION, PARAMETER ::  pi = 4.D0*DATAN(1.D0)
    
    INTEGER, INTENT(IN)             :: n
    DOUBLE PRECISION, INTENT(IN)    :: alpha(n,1)
    DOUBLE PRECISION, INTENT(IN)    :: coeff(n,1)
    DOUBLE PRECISION, INTENT(IN)    :: coord(n,3)
    
    !Matrix
    DOUBLE PRECISION, ALLOCATABLE   :: Two_e(:,:,:,:)
    
    !Iteratives
    INTEGER ::  p,q,r,s, pp,qq,rr,ss
    
    !For calculations
    DOUBLE PRECISION        :: c1c2c3c4, alphasum_pq, alphasum_rs, alphaprod_pq, alphaprod_rs, normfactor
    DOUBLE PRECISION        :: position_pq(3), position_rs(3), distance2_pq, distance2_rs
    DOUBLE PRECISION        :: alphacoord_rs(3), alphacoord_pq(3), PG_pqrs(3), PG2_pqrs, denom
    DOUBLE PRECISION        :: term1, term2, term3, term4
    
    OPEN(unit = 4, file = 'two_electron.dat', status = 'replace', action = 'write')  !Two electron Matrix
    ALLOCATE (Two_e(n,n,n,n))
    Two_e = 0.0

    !Two electron integral
    DO p = 1,n
    DO q = 1,n
    DO r = 1, n
    DO s = 1, n
                DO pp = 1, SIZE(alpha(p,:))
                DO qq = 1, SIZE(alpha(q,:))
                DO rr = 1, SIZE(alpha(r,:))
                DO ss = 1, SIZE(alpha(s,:))
                        !Factors for calculations
                        c1c2c3c4 = coeff(p,pp)*coeff(q,qq)*coeff(r,rr)*coeff(s,ss)
                        alphasum_pq = alpha(p,pp)+alpha(q,qq)
                        alphasum_rs = alpha(r,rr)+alpha(s,ss)
                        denom = 1.0/alphasum_pq + 1.0/alphasum_rs
                        alphacoord_pq = alpha(p,pp)*coord(p,:) + alpha(q,qq)*coord(q,:)
                        alphacoord_rs = alpha(r,rr)*coord(r,:) + alpha(s,ss)*coord(s,:)
                        normfactor = (2*alpha(p,pp)/pi)**0.75*(2*alpha(q,qq)/pi)**0.75*&
                                    (2*alpha(r,rr)/pi)**0.75*(2*alpha(s,ss)/pi)**0.75
                    
                        !alphaprod = alpha(p,r)*alpha(q,s)
                        !expfactor = alphaprod/alphasum
                        !distance = coord(p,:)-coord(q,:)
                    
                    
                        !Factors for kinetic energy
                    
                        PG_pqrs = alphacoord_pq/alphasum_pq - alphacoord_rs/alphasum_rs
                        PG2_pqrs = SUM(PG_pqrs*PG_pqrs)
                        alphaprod_pq = alpha(p,pp)*alpha(q,qq) / alphasum_pq
                        alphaprod_rs = alpha(r,rr)*alpha(s,ss) / alphasum_rs
                        position_pq = coord(p,:) - coord(q,:)
                        position_rs = coord(r,:) - coord(s,:)
                        distance2_pq = SUM(position_pq*position_pq)
                        distance2_rs = SUM(position_rs*position_rs)
                    
                        term1 = 2.0*PI**2/(alphasum_pq*alphasum_rs)
                        term2 = SQRT(PI/(alphasum_pq+alphasum_rs))
                        term3 = EXP(-alphaprod_pq*distance2_pq)
                        term4 = EXP(-alphaprod_rs*distance2_rs)
                        Two_e(p,q,r,s) = Two_e(p,q,r,s) + normfactor * c1c2c3c4 *&
                                         term1 * term2 * term3 * term4 * F0(PG2_pqrs/denom) 
                END DO
                END DO
                END DO
                END DO
                WRITE (4,'(i3, i3, i3, i3, f8.4)') p, q, r, s, Two_e(p,q,r,s)
    END DO
    END DO
    END DO
    END DO
 
    
    CLOSE(4)
    DEALLOCATE(Two_e)
    
    CONTAINS 
        FUNCTION F0(x)
                DOUBLE PRECISION, INTENT(IN)    ::  x
                DOUBLE PRECISION                ::  F0
                DOUBLE PRECISION, PARAMETER     ::  pi = 4.D0*DATAN(1.D0)
            
                IF (ABS(x) < 1E-4) THEN
                    F0 = 1.0D0
                    RETURN
                END IF
                F0 = x**(-0.5) * SQRT(pi)/2 * ERF(x**(0.5))
                RETURN
        END FUNCTION
    !    function boys(x,n)
    !            implicit none
    !            DOUBLE PRECISION, intent(in) :: x
    !            integer, intent(in) :: n
    !            real :: boys
    !            
    !            IF (x==0.) THEN
    !                boys = 1.0/(2.0*n+1)
    !                RETURN
    !            ELSE
    !                boys =  lower_incomplete_gamma(n+0.5,x)*GAMMA(n+0.5) * (1.0/(2*x**(n+0.5)))
    !            END IF
    !            
    !            
    !        
    !     end function boys
    !            
    !FUNCTION lower_incomplete_gamma(s, x) 
    !    IMPLICIT NONE
    !    
    !    REAL, INTENT(IN) :: s
    !    DOUBLE PRECISION, INTENT(IN) :: x
    !    DOUBLE PRECISION             :: lower_incomplete_gamma
    !    INTEGER, PARAMETER :: n = 1e6  ! Number of intervals
    !    REAL :: dx
    !    INTEGER :: i
    !    REAL :: t
    !
    !    ! Compute step size
    !    dx = x / REAL(n)
    !
    !    ! Initialize integral
    !    lower_incomplete_gamma = 0.0d0
    !
    !    ! Perform summation
    !    DO i = 1, n-1
    !        t = i * dx
    !        lower_incomplete_gamma = lower_incomplete_gamma + t**(s-1) * EXP(-t)
    !    END DO
    !
    !    ! Apply trapezoidal rule
    !    lower_incomplete_gamma = dx * (0.5 * (1.0 + x**(s-1) * EXP(-x)) + lower_incomplete_gamma)
    !    lower_incomplete_gamma = lower_incomplete_gamma/GAMMA(s)
    !    RETURN
    !END FUNCTION
    
        
    END SUBROUTINE
