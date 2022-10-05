
	  subroutine ellipticity_check( defoGrad_F, dS_dC, stress_PK2,
     & 								idele )
c		
c @todo Extend or change this, so we can call it from utan
c @todo Option to choose dS_dC or E and PK2 or Cauchy_stress
c
      use Tensor
c
	  implicit none
	  include 'iounits.inc'
c
c INPUT
		 type(Tensor2), intent(in) :: defoGrad_F, stress_PK2
		 type(Tensor4), intent(in) :: dS_dC
		 integer, intent(in) :: idele
c LOCAL
		 type(Tensor4) :: dP_dF
		 type(Tensor2) :: delta, q
		 type(Tensor1) :: normal_vector
		 integer dim, n_angle_steps
		 integer i, I_, J_, k, K_, L_, i_alpha, i_beta
		 real*8, parameter :: PI = 4. * atan(1.)
		 real*8 :: det_q, det_q_min
c
		delta = identity2(delta)
c
c Transform dS_dC to dP_dF
	  do i=1,3
		do J_=1,3
			do k=1,3
				do L_=1,3
					do I_=1,3
						do K_=1,3
						  dP_dF%abcd(i,J_,k,L_) = dP_dF%abcd(i,J_,k,L_)
     &					   + 2.*dS_dC%abcd(I_,J_,K_,L_)
     &					     * defoGrad_F%ab(i,I_)
     &						 * defoGrad_F%ab(k,K_)
						enddo
					enddo
				enddo
			enddo
		enddo
	  enddo
	  do i=1,3
		do J_=1,3
			do k=1,3
				do L_=1,3
					dP_dF%abcd(i,J_,k,L_) = dP_dF%abcd(i,J_,k,L_)
     &					+ stress_PK2%ab(L_,J_) * delta%ab(i,k)
				enddo
			enddo
		enddo
	  enddo
c
c Choose the number of steps to discretise the normal vectors
	  dim = 3;
		select case( dim )
		case ( 2 )
		n_angle_steps = 1000
		case ( 3 )
		n_angle_steps = 20!100
		case default
	   write( *, * ) "ellipticity_check<<
     & Dimension (",dim,") is not available,"
		stop
	  end select 
c
c Loop over all possible normal vectors and compute the det(q) for each
	  det_q_min = 9.9999e20
	  do i_alpha=0,n_angle_steps
		do i_beta=0,n_angle_steps
		 select case( dim )
		  case ( 2 )
		   normal_vector%a(1) = cos(i_alpha*PI/20.)
		   normal_vector%a(2) = sin(i_alpha*PI/20.)
		   normal_vector%a(3) = 0
		 case ( 3 )
		   normal_vector%a(1) = sin(i_alpha*PI/20.) * cos(i_beta*PI/20.)
		   normal_vector%a(2) = sin(i_alpha*PI/20.) * sin(i_beta*PI/20.)
		   normal_vector%a(3) = cos(i_alpha*PI/20.)
	     end select 
c
		do i=1,3
			do J_=1,3
				do k=1,3
					do L_=1,3
						q%ab(i,k) = q%ab(i,k)
     &								+ normal_vector%a(J_)
     &					 			  * dP_dF%abcd(i,J_,k,L_)
     &					  			  * normal_vector%a(L_)
					enddo
				enddo
			enddo
		enddo
c
		 det_q = det( q )
		 det_q_min = min(det_q_min,det_q)
c
		if ( dim==2 ) exit ! inner loop over i_beta
c		 
		enddo ! i_beta
	enddo ! i_alpha
c
c
	    if ( det_q_min < -1e-14 ) then
			write(iomsg,*) "### loss of ellipticity_loE detected",
     &		det_q_min," for idele ",idele," ###"
		endif	
c
	  end subroutine ellipticity_check
