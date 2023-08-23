!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!this function prepares the two-body potential in NRQCD
!in coordinate space. The second version below is meant to be used with BLM scale fixing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function VpNRQCD(order,r,mred)
!r is the radial variable to evaluate the potential
!order can be 0,1,2 for LO,NLO,NNLO
!See Antonio Vairo' hep-ph/0010191, for the 1/m pieces Kniehl et al prd65, 091503 2002
      implicit none
      integer,intent(in) :: order
      double precision, parameter ::  pi = 4.D0*DATAN(1.D0)
      double precision, parameter ::  gammaE = 0.57721d0
	  double precision, intent(in)::  r, mred
      double precision :: a1,a2,beta0,Riemann3,alphas,alpha, VpNRQCD
      double precision :: VLO,VNLO,VNNLO,musq,beta1 &
                        ,freezingscale,b2
      integer :: switch1overM,Coulombswitch,Nf
      ! common/coupling/alphas,musq,freezingscale
      ! common/betafunction/beta0,beta1
      ! common/lightflavors/Nf
      ! common/pinumber/pi
      ! common/orderpert/order
      ! common/recoil/switch1overM
      ! common/electm/Coulombswitch
      external alpha

      OPEN(unit = 10, file = "instructions.dat", status = 'OLD', action = 'READ')
      READ(10, '(i1)') Nf
      !WRITE(*,*) 'Nf'
      !WRITE(*,*) Nf
      
      READ(10, '(i1)') switch1overM
      !WRITE(*,*) 'switch1m'
      !WRITE(*,*) switch1overM

      READ(10, '(i1)') Coulombswitch
      !WRITE(*,*) 'Coulomb'
      !WRITE(*,*) Coulombswitch
      CLOSE(10)
      if(order.lt.0.or.order.gt.2) then
       WRITE(*,*) 'The potential is only given to LO,NLO,NNLO',order
       
       stop
      endif
! c The Leading Order potential. Since r is in GeV**-1 can pass to the
! c functions that compute the running coupling constant
! c For consistency alpha_s is computed at the correct order for each piece
! c except the LO potential where we use alpha ran at NLO at the corresponding
! c scale to avoid nonsense answers
      if(order.eq.0) then
       VLO=-(4d0/3d0)*alpha(r**(-2),0)/r
       VpNRQCD=VLO
       return
      endif

      
      a1=31d0/3d0-10d0/9d0*Nf
! c the NLO potential, first the LO piece with alpha_s evaluated at NLO
      beta0=11 - Nf*2.d0/3
      if(order.eq.1) then
       VNLO=-(4d0/3d0)*alpha(r**(-2),1)/r
!c then the NLO piece with alpha_s evaluated at LO
       VLO=-(4d0/3d0)*alpha(r**(-2),0)/r
       VNLO=VNLO + VLO*&
           alpha(r**(-2),0)/(4d0*pi)*(a1+2d0*gammaE*beta0)
       VpNRQCD=VNLO
!c if active, add the 1/m correction

      beta1=102 - Nf*38.d0/3

      
       if(switch1overM.eq.1) then
        VpNRQCD=VpNRQCD-alpha(r**(-2),1)**2/(mred*r**2)
!c *(7d0/9d0) This factor for Smirnov's convention, without it, Nora's
       endif
       
       return
      endif  
      Riemann3=1.202d0
! c      a2=(4343d0/162d0+4d0*pi**2-pi**4/4d0+22d0/3d0*Riemann3)*9d0
! c     &  -(1798d0/81d0+56d0/3d0*Riemann3)*3d0/2d0*Nf
! c     & +(100d0/81d0*Nf-(55d0/3d0-16d0*Riemann3)*(2d0/3d0))*Nf
!c a simpler expression for a2:
      a2=4343d0/18d0+36d0*pi**2-9d0*pi**4/4d0+66d0*Riemann3&
        -(1229d0/27d0+52d0/3d0*Riemann3)*Nf+100d0/81d0*Nf**2
! c the NNLO potential. No if_then needed, order=2 automatically,
! c the former if_then  have excluded the other possible values
! c alpha_s is evaluated at NNLO only for the VLO part,
      VNNLO=-(4d0/3d0)*alpha(r**(-2),2)/r
!c then at NLO for the NLO part,
      VLO=-(4d0/3d0)*alpha(r**(-2),1)/r
      VNNLO=VNNLO + VLO*&
           alpha(r**(-2),1)/(4d0*pi)*(a1+2d0*gammaE*beta0)
!c and at LO for the NNLO part of the potential
      VLO=-(4d0/3d0)*alpha(r**(-2),0)/r
      VNNLO=VNNLO+VLO*(alpha(r**(-2),0)/(4d0*pi))**2 * ( a2 +&
       beta0**2*(pi**2/3d0+4d0*gammaE**2)+&
       gammaE*(4d0*beta0*a1+2d0*beta1)               )
!c at this order add the Coulomb potential. Charges are 2/3 for c, -1/3 for b
       
      if(Coulombswitch.eq.1) then
       if(Nf.eq.4) then
        VNNLO=VNNLO-1./(9.*137.*r)
       elseif(Nf.eq.3) then
        VNNLO=VNNLO-4./(9.*137.*r)
       endif
      endif
      if(switch1overM.eq.1) then
       if(Nf.eq.3) then
        b2=-20.836d0
       elseif(Nf.eq.4) then
        b2=-18.943d0
       elseif(Nf.eq.5) then
        b2=-17.049d0
       else
        write(*,*) 'Potential rigged for Nf 3 to 5', Nf
        
        stop
       endif
       !READ(10,'(i1)') musq
       musq=r**2
       VNNLO=VNNLO-alpha(r**(-2),2)**2/(mred*r**2)*(7d0/9d0)&
             -alpha(r**(-2),2)**3/(3d0*pi*mred*r**2)*( b2&
            +(7d0/6d0*beta0 + 68d0/3d0)*dlog(musq*r**2))
      endif
      VpNRQCD=VNNLO
      
      return


      ! contains 

      ! function alpha(r,order)
      !       double precision :: alpha
      !       double precision :: r
      !       integer          :: order
      !       alpha = 0.4d0
      !       return
      ! end function

      ! function beta0(Nf)
      !       double precision :: beta0
      !       integer          :: Nf
      !       beta0=11-Nf*2.d0/3
      !       return
      ! end function

      ! function beta1(Nf)
      !       double precision :: beta1
      !       integer          :: Nf
      !       beta1=102 - Nf*38.d0/3
      !       return
      ! end function

      end
