program plot_VpNRQCD
	implicit none
	double precision				::	r,Vp1,rlog,Vp2,Vp3
	double precision				::	mred
	double precision				::	rmax,rmin,rscale
	integer							::	iterator,rlength
	interface
          double precision function VpNRQCD(order,r,mred)
            double precision, intent(in) :: r,mred
			integer, intent(in)			 :: order
          end function VpNRQCD
        end interface
	
	rlength = 1e4

	rmin = -0.7d0
	r = rmin
	rmax = 3.d0
	rscale = (rmax-rmin)/rlength
	OPEN(unit = 3, file = 'Vp.dat',status = 'replace', action = 'write')

	DO iterator = 1,rlength
		r = r+rscale
		!rlog = 10**r
		Vp1 = VpNRQCD(0,10**r,1.d0)
		Vp2 = VpNRQCD(1,10**r,1.d0)
		Vp3 = VpNRQCD(2,10**r,1.d0)
		!WRITE(*,*) Vp
		WRITE(3,'(4f12.4)') 10**r, Vp1, Vp2, Vp3
	END DO

	CLOSE(3)
end program
