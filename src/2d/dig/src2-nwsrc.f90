

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module, only: grav, dry_tolerance,deg2rad,friction_depth
      use geoclaw_module, only: manning_coefficient,friction_forcing

      use digclaw_module, only: alpha,alpha_seg,bed_normal,curvature
      use digclaw_module, only: entrainment,entrainment_rate
      use digclaw_module, only: i_ent,i_fsphi,i_phi,i_theta
      use digclaw_module, only: mu,rho_f,sigma_0
      use digclaw_module, only: admissibleq,auxeval
      use digclaw_module, only: calc_pmtanh

      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
<<<<<<< HEAD
      real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
      
=======
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

>>>>>>> main
      ! Output
      real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
<<<<<<< HEAD
      real(kind=8) :: gz,gx,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      real(kind=8) :: rhoh,manning_n
      real(kind=8) :: D,tau,sigbed,kperm,compress,chi,tol
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,hvnorm0
      real(kind=8) :: shear,sigebar,chitanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,gacc,beta
=======
      real(kind=8) :: gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      real(kind=8) :: D,tau,sigbed,kperm,compress,pm,coeff
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,taucf,fsphi,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,chi,beta
>>>>>>> main
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining,dtcoeff
      real(kind=8) :: gamma,zeta,krate,p_eq,dgamma

<<<<<<< HEAD
      integer :: i,j,ii,jj,jjend,icount
=======
      integer :: i,j,ii,jj,icount
      logical :: ent
>>>>>>> main


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

<<<<<<< HEAD
      manning_n = manning_coefficient(1) ! Current implementation of friction has manning as an array 
      ! take the first element for now. If only one value is provided to geo_data.manning_coefficient 
=======
      gmod=grav

      if (friction_forcing) then
         coeff = manning_coefficient(1)
      else
         coeff = 0.d0
      endif

      ! Current implementation of friction has manning as an array
      ! take the first element for now. If only one value is
      ! provided to geo_data.manning_coefficient
>>>>>>> main
      ! it will be manning_coefficient(1)
      ! DIG: Decide if this should be handled in some other way.

<<<<<<< HEAD
      tol = dry_tolerance !# to prevent divide by zero in gamma
      
      gz = grav  !needed later for bed-normal direction gravity
      gx = 0.d0
      theta=0.d0 

      do i=1-mbc+1,mx+mbc-1
         do j=1-mbc+1,my+mbc-1
            
            if (q(1,i,j)<=dry_tolerance) cycle
            
=======
      if (entrainment>0) then
         ent = .true.
      else
         ent = .false.
      endif

      do j=1-mbc+1,my+mbc-1
         do i=1-mbc+1,mx+mbc-1
         ! DIG: 1/12/24: KRB and MJB notice that here we are looping over ghost cells.
         ! These ghost cells have not been updated by the riemann solver? Are they used
         ! meaningfully (e.g., theta is diffed below.)
         ! 1/30/2024 - Leaving this as is for the moment, this is something to evaluate later.

            ! adjust gravity if bed_normal = 1
            theta = 0.d0
            dtheta = 0.d0
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gmod = grav*cos(theta)
               dtheta = -(aux(i_theta,i+1,j) - theta)/dx
            endif

            ! Get state variable
>>>>>>> main
            h = q(1,i,j)
            hu = q(2,i,j)
            hv = q(3,i,j)
            hm = q(4,i,j)
            p =  q(5,i,j)
<<<<<<< HEAD
            phi = aux(ia_phi,i,j)
            chi = q(iq_seg,i,j)/h
            chi = max(0.0,chi)
            chi = min(1.0,chi)
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)

            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(ia_theta,i,j)
               gz = grav*cos(theta)
               gx = grav*sin(theta)
            endif
            if (curvature==1) then
               b_xx=(aux(1,i+1,j)-2.d0*aux(1,i,j)+aux(1,i-1,j))/(dx**2)
               b_yy=(aux(1,i,j+1)-2.d0*aux(1,i,j)+aux(1,i,j-1))/(dy**2)
               b_xy=(aux(1,i+1,j+1)-aux(1,i-1,j+1) -aux(1,i+1,j-1)+aux(1,i-1,j-1))/(4.0*dx*dy)
               dtheta = -(aux(ia_theta,i+1,j) - theta)/dx
               gacc = max((u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy + u**2*dtheta,0.d0)!max:currently only consider enhancement not reduction of gz (ie. basin not a hump)
               gz = gz + gacc
=======
            phi = aux(i_phi,i,j)
            pm = q(6,i,j)/h
            pm = max(0.0d0,pm)
            pm = min(1.0d0,pm)
            fsphi = aux(i_fsphi,i,j)
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)

            !integrate momentum source term
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            ! calculate total velocity
            vnorm = sqrt(u**2 + v**2)
            hvnorm = sqrt(hu**2 + hv**2)
            hvnorm0 = hvnorm

            !integrate friction
            hvnorm = dmax1(0.d0,hvnorm - dt*tau/rho)
            hvnorm = hvnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(rho*h**2))
            if (hvnorm<1.d-16) hvnorm = 0.d0

            ! adjust based on curvature
            if (hvnorm>0.d0.and.curvature==1) then
               b_xx=(aux(1,i+1,j)-2.d0*aux(1,i,j)+aux(1,i-1,j))/(dx**2)
               b_yy=(aux(1,i,j+1)-2.d0*aux(1,i,j)+aux(1,i,j-1))/(dy**2)
               b_xy=(aux(1,i+1,j+1)-aux(1,i-1,j+1) -aux(1,i+1,j-1)+aux(1,i-1,j-1))/(4.d0*dx*dy)
               chi = (u**2*b_xx + v**2*b_yy + 2.0d0*u*v*b_xy)/gmod
               chi = max(chi,-1.d0)
               taucf = chi*tau
               hvnorm = dmax1(0.d0,hvnorm - dt*taucf/rho)
               taucf = u**2*dtheta*tau/gmod
               hvnorm = dmax1(0.d0,hvnorm - dt*taucf/rho)
>>>>>>> main
            endif

            call vareval(h,u,v,m,p,rho,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation)

            rhoh = h*rho !this is invariant in src and always >0 below
            hvnorm0 = sqrt(hu**2 + hv**2)
            vnorm = hvnorm0/h

            if (hvnorm0>0.d0) then
               !integrate dynamic friction !DIG: TO DO - move dynamic friction to Riemann solver
               vnorm = dmax1(0.d0,vnorm - dt*tau/rhoh) !exact solution for Coulomb friction
               vnorm = vnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(h*rhoh)) !exact solution (prior to h change) for effective viscous friction
               ! velocity determined, calculate directions etc. from vnorm
               hvnorm = h*vnorm
               hu = hvnorm*hu/hvnorm0 + gx*h*dt !gx=0 unless bed-normal !DIG: last term should ultimately be in Riemann solver
               hv = hvnorm*hv/hvnorm0
               u = hu/h
               v = hv/h
               ! velocity now constant for remainder of src2
            endif

<<<<<<< HEAD
            if (p_initialized==0) cycle !DIG: deprecate?

            !determine m and p evolution from multiple integration steps maintaining physically admissable m,p. h = rhoh/rho
            ! Note: if m = 0, it cannot increase as dm/dt - m D. In that case, rho=rho_f ==> D = 0 for all t, and dp/dt =0.
            if (m>0.d0) then
               dt_remaining = dt
               if (D==0.d0 & (vnorm==0.d0)) then !at a critical point already, rhs = 0
                     dt_remaining = 0.d0
               endif
               do while (dt_remaining>0.d0)
                  call integrate_mp(h,u,v,m,p,rho,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation,dt_remaining,dt_taken)
                  dt_remaining = max(0.d0,dt_remaining-dt_taken)
               enddo
            endif




=======
            ! call admissible q and auxeval before moving on to shear induced
            ! dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            ! calculate velocity
            vnorm = sqrt(u**2 + v**2)
>>>>>>> main

            !integrate shear-induced dilatancy
            shear = 2.d0*vnorm/h
            krate = 1.5d0*shear*m*tanpsi/alpha
<<<<<<< HEAD
            !sigebar = sigebar*exp(krate*dt)
            !p = rho*gmod*h + sigma_0 - sigebar
=======
>>>>>>> main
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dt*3.d0*vnorm*tanpsi/(h*compress)
            endif

<<<<<<< HEAD
            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               chitanh01=0.d0

=======
            ! update pmtanh01 and rho_fp for segregation
            ! DIG: if segregation is compartmentalized
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               pmtanh01=0.d0
>>>>>>> main
      		else
         		seg = 1.d0
               call calc_chitanh(chi,seg,chitanh01)
               rho_fp = max(0.d0,(1.d0-chitanh01))*rho_f
      		endif
<<<<<<< HEAD
            !chitanh01 = seg*(0.5*(tanh(20.0*(chi-0.80))+1.0))
            !chitanh01 = seg*(0.5*(tanh(40.0*(chi-0.90))+1.0))
            
            !integrate pressure relaxation
            !if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
            !   zeta = 3.d0/(compress*h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !else
            !   zeta = (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !endif
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.d0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
            p_hydro = h*rho_fp*gmod
            p_litho = (rho_s*m + (1.d0-m)*rho_fp)*gmod*h

            !if (abs(compress*krate)>0.0) then
            !   p_eq = p_hydro + 3.0*vnorm*tanpsi/(compress*h*krate)
            !else
            !   p_eq = p_hydro
            !endif
            !if (abs(chi-.5)>.49) then
            !chitanh01 = 0.5*(tanh(20.0*(chi-0.80))+1.0)
            p_eq = p_hydro !*(1.0-chitanh01)
            !p_eq = max(p_eq,0.0)
            !p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dt)

=======

            ! integrate pressure relaxation
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.d0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
            p_eq = h*rho_fp*gmod
            p = p_eq + (p-p_eq)*exp(krate*dt)
>>>>>>> main

            ! call admissible q and auxeval before moving on to dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
<<<<<<< HEAD
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
            
=======
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
>>>>>>> main

            ! calculate rate of change
            krate = D*(rho-rho_fp)/rho
<<<<<<< HEAD
            hu = hu*exp(dt*krate/h)
            hv = hv*exp(dt*krate/h)
            hm = hm*exp(-dt*D*rho_fp/(h*rho))
            h = h + krate*dt
=======
>>>>>>> main

            ! integrate hu, hv, hm, and h.
            hu = hu*exp(dt*krate/h)
            hv = hv*exp(dt*krate/h)
            hm = hm*exp(-dt*D*rho_fp/(h*rho))
            h = h + krate*dt


            !======================mass entrainment===========================

            ! call admissible q and auxeval before moving on to mass entrainment.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

            vnorm = sqrt(u**2 + v**2)
<<<<<<< HEAD
            vlow = 0.1d0

            if (ent.and.vnorm.gt.vlow.and.(aux(ia_theta,i,j)>0.d0)) then
               b_x = (aux(1,i+1,j)+q(7,i+1,j)-aux(1,i-1,j)-q(7,i-1,j))/(2.d0*dx)
               b_y = (aux(1,i,j+1)+q(7,i,j+1)-aux(1,i,j-1)-q(7,i,j-1))/(2.d0*dy)
               dbdv = (u*b_x+v*b_y)/vnorm
               slopebound = 1.d10
               b_eroded = q(7,i,j)
               if (dbdv<slopebound.and.b_eroded<aux(ia_theta,i,j)) then
                  b_remaining = aux(ia_theta,i,j)-b_eroded
                  m2 = 0.6d0
                  rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0
                  beta2 = 0.66d0
                  t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                  !write(*,*) '------------'
                  !write(*,*) 'vu',t1bot
                  beta = 1.d0-m!tanh(10.d0*m) !tan(1.5*p/(rho*gmod*h))/14.0
                  gamma= rho*beta2*(vnorm**2)*(beta*gmod*manning_n**2)/(tanh(h+1.d-2)**(1.0/3.0))
                  !write(*,*) 'gamma', gamma
                  t1bot = t1bot + gamma
                  t1bot = t1bot + tau!+p*tan(phi)
                  !write(*,*) 'tau',tau
                  t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))
                  !write(*,*) 't2top',t2top
                  prat = p/(rho*h)
                  !dh = dt*(t1bot-t2top)/(beta2*tanh(vnorm+1.d-2)*rho2)
                  vreg = ((vnorm-vlow)**2/((vnorm-vlow)**2+1.d0))
                  dtcoeff = entrainment_rate*dt*vreg/(beta2*(vnorm+vlow)*rho2)
                  !dh = dtcoeff*t1bot/(1.d0 + dtcoeff*tan(phi))
                  dh = dtcoeff*(t1bot-t2top)
                  dh = entrainment_rate*dt*(t1bot-t2top)/(rho2*beta2*vnorm)
                  !write(*,*) 'dh',dh
                  !write(*,*) 'dh/dt', dh/dt
                  dh = min(dh,b_remaining)
                  h = h + dh
                  hm = hm + dh*m2
                  q(7,i,j) = q(7,i,j) + dh

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
                  p = prat*rho*h
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
=======
            vlow = 0.1d0 ! minimum velocity for entrainment to occur. ! DIG: should this be a user
            ! specified variable.

            if (ent.and.vnorm.gt.vlow) then
               if (aux(i_ent,i,j)>0.d0) then
                  b_x = (aux(1,i+1,j)+q(7,i+1,j)-aux(1,i-1,j)-q(7,i-1,j))/(2.d0*dx)
                  b_y = (aux(1,i,j+1)+q(7,i,j+1)-aux(1,i,j-1)-q(7,i,j-1))/(2.d0*dy)
                  dbdv = (u*b_x+v*b_y)/vnorm
                  slopebound = 1.d10
                  b_eroded = q(7,i,j)

                  if (dbdv<slopebound.and.b_eroded<aux(i_ent,i,j)) then
                     b_remaining = aux(i_ent,i,j)-b_eroded

                     ! value for m for entrained material
                     m2 = 0.6d0
                     ! eventually make this a user defined variable or an aux value.
                     ! should there also be a substrate value for chi?

                     ! calculate entrained material density, using default values
                     ! for rho_f and rho_s
                     rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0

                     beta2 = 0.66d0

                     ! calculate top and bottom shear stress.
                     t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d0-2.d0))
                     beta = 1.d0-m

                     t1bot = t1bot + tau

                     t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))

                     ! calculate pressure ratio
                     prat = p/(rho*h)

                     ! regularize v
                     vreg = ((vnorm-vlow)**2/((vnorm-vlow)**2+1.d0))

                     dtcoeff = entrainment_rate*dt*vreg/(beta2*(vnorm+vlow)*rho2)

                     ! calculate dh
                     dh = entrainment_rate*dt*(t1bot-t2top)/(rho2*beta2*vnorm)
                     dh = min(dh,b_remaining)

                     ! increment h based on dh
                     h = h + dh
                     hm = hm + dh*m2

                     ! store amount eroded in q7
                     q(7,i,j) = q(7,i,j) + dh

                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                     call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                     ! update pressure based on prior pressure ratio.
                     p = prat*rho*h

                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  endif
>>>>>>> main
               endif
            endif

            !===================================================================
            ! end of entrainment, put state variables back in q.

            q(1,i,j) = h
            q(2,i,j) = hu
            q(3,i,j) = hv
            q(4,i,j) = hm
            q(5,i,j) = p
            q(6,i,j) = chi*h

         enddo
      enddo



      ! Manning friction------------------------------------------------
      if (friction_forcing) then
<<<<<<< HEAD
      if (manning_n>0.d0.and.friction_depth>0.d0) then
         do i=1,mx
            do j=1,my
=======
      if (coeff>0.d0.and.friction_depth>0.d0) then

         do j=1,my
            do i=1,mx

               if (bed_normal==1) gmod = grav*cos(aux(i_theta,i,j))
>>>>>>> main
                  h=q(1, i,j)
               if (h<=friction_depth) then
                 !# apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)
                  hm = q(4,i,j)
                  p =  q(5,i,j)
                  phi = aux(ia_phi,i,j)
                  theta = aux(ia_theta,i,j)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<dry_tolerance) cycle
<<<<<<< HEAD
                  chi = q(6,i,j)/h
                  chi = max(0.0,chi)
                  chi = min(1.0,chi)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

                  if (h.lt.dry_tolerance) then
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     beta = 1.d0-m !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*manning_n**2)/(h**(7.0/3.0))
=======
                  pm = q(6,i,j)/h
                  pm = max(0.0d0,pm)
                  pm = min(1.0d0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                  if (h.lt.dry_tolerance) then
                     q(1,i,j)=0.d0
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     beta = 1.d0-m
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*coeff**2)/(h**(7.0d0/3.0d0))
>>>>>>> main
                     dgamma=1.d0 + dt*gamma
                     q(2,i,j)= q(2,i,j)/dgamma
                     q(3,i,j)= q(3,i,j)/dgamma
                  endif
               endif
            enddo
         enddo
         endif
      endif
     ! ----------------------------------------------------------------

      return
      end subroutine src2

   !====================================================================
   ! subroutine integrate_mp: integrate portion of rhs for m and p
   !     dm/dt ~ -Dm
   !     dp/dt ~ D + (m-m_eqn)
   !
   !     integrated for p_tilde = p - p_eq
   !     dm/dt ~ -p_tilde m
   !     dp_tilde/dt ~ p_tilde + (m-m_eqn)
   !     note: p_eq is h(t) dependent. dp_tilde/dt = dp/dt - rho_f gz *dh/dt ~ dp/dt + p_tilde
   !====================================================================

   subroutine integrate_mp(h,u,v,m,p,rhoh,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation,dt_remaining,dt)

      use geoclaw_module, only: grav, dry_tolerance
      use digclaw_module, only: rho_f, rho_s, m_crit
      implicit none

      use geoclaw_module, only: grav, dry_tolerance

      !i/o
      real(kind=8), intent(in) :: vnorm,rhoh,phi_bed,gz
      real(kind=8), intent(inout)  :: h,m,p,kperm,dt_remaining,dt_taken
      real(kind=8), intent(out) ::

      !local
      real(kind=8) :: p_eq, p_tilde
      real(kind=8) :: 

      call vareval(h,u,v,m,p,rhoh,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation)
      shear = vnorm/h !

      p_eq = rho_f*gz*h*M_saturation !Note: M_saturation =1, unless experimenting with segregation models
      p_tilde = p - p_eq
      m_eq = !DIG: WIP - maybe this and above calculated in vareval

      if (D==0.d0 & (sigma_e==0.d0|shear==0.d0|tanpsi==0.d0)) then
         dt = dt_remaining
         return
      endif

      if (D==0.d0) then !integrate only shear induced dilatancy. no change in h or m.
         rhs_p = -(3.d0/alpha_c)*m*sigma_e*shear*tanpsi
         if (rhs_p>0.d0) then !m<m_eq => increasing p
            dt = min(abs((rhoh*gz-p)/rhs_p),dt_remaining) !enforce dt st p1<= rho g h
            !forward euler. rhs_p is nonlinear function of p through tanpsi(m_eq)
            !DIG: - update to RK. 
            p = p + rhs_p*dt 
            return
         elseif (rhs_p<0.d0) then
            dt = min(abs(p/rhs_p),dt_remaining) !enforce dt st p1>=0
            !forward euler. rhs_p is nonlinear function of p through tanpsi(m_eq)
            !DIG: - update to RK. 
            p = p + rhs_p*dt 
            return
         else
            dt = dt_remaining
            return
         endif 
      elseif (D>0.d0) then
         !note: by integrating m (=> rho and hence h), dh/dt terms on rhs of p are avoided by using new p_eq(h).

         rhs_ptilde =  -(3.d0*alphainv/h)*(vnorm*tanpsi +0.5d0*

         
      else

      endif

      if (dabs(alpha_seg-1.0d0)<1.d-6) then
         seg = 0.0d0
         rho_fp = rho_f
         pmtanh01=0.0d0
      else
         seg = 1.0d0
         call calc_pmtanh(pm,seg,pmtanh01)
         rho_fp = (1.0d0-pmtanh01)*rho_f
      endif
      !pmtanh01 = seg*(0.5*(tanh(20.0*(pm-0.80))+1.0))
      !pmtanh01 = seg*(0.5*(tanh(40.0*(pm-0.90))+1.0))

      if (bed_normal.eq.1) gmod=grav*dcos(theta)
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_fp*(1.d0-m)
      shear = 2.0d0*vnorm/hbounded
      sigbed = dmax1(0.d0,rho*gmod*h - p)
      sigbedc = rho_s*(shear*delta)**2 + sigbed
      if (sigbedc.gt.0.0d0) then
         S = (mu*shear/(sigbedc))
      else
         S = 0.d0
      endif
      !Note: m_eqn = m_crit/(1+sqrt(S))
      !From Boyer et. al
      !S = 0.0
      !m_eqn = m_crit/(1.d0 + sqrt(S))
      !if (m.gt.m_eqn) write(*,*) 'm,m_eqn,S:',m,m_eqn,S,sigbed,shear
      !tanpsi = c1*(m-m_eqn)*tanh(shear/0.1)
      !pmlin = seg*2.0*(pm-0.5)
      !pmtan = seg*0.06*(tan(3.*(pm-0.5)))
      !pmtanh = seg*tanh(3.*pmlin)
      !pmtanh01 = seg*0.5*(tanh(8.0*(pm-0.75))+1.0)
      !pmtanh01 = seg*0.5*(tanh(20.0*(pm-0.80))+1.0)
      !pmtanh01s = seg*4.0*(tanh(8.0*(pm-0.95))+1.0)

   
      kperm = kappita*exp(-(m-m0)/(0.04d0))!*(10**(pmtanh01))
      !m_crit_pm - max(pm-0.5,0.0)*(0.15/0.5) - max(0.5-pm,0.0)*(0.15/0.5)
      !m_crit_pm =  max(pm-0.7,0.0)*((m_crit- 0.55)/0.5) + max(0.3-pm,0.0)*((m_crit-0.55)/0.5)
      m_crit_pm =  0.d0! max(pm-0.6,0.0)*((m_crit- 0.55)/0.4) + max(0.3-pm,0.0)*((m_crit-0.55)/0.3)
      !m_crit_pm = max(pm-0.9,0.0)*((m_crit- 0.55)/0.1) + max(0.1-pm,0.0)*((m_crit-0.55)/0.1);

      m_crit_pm = pmtanh01*0.09d0
      m_crit_m = m_crit - m_crit_pm
      m_eqn = m_crit_m/(1.d0 + sqrt(S))
      tanpsi = c1*(m-m_eqn)*tanh(shear/0.1)

      !kperm = kperm + 1.0*pm*kappita
      !compress = alpha/(sigbed + 1.d5)

      !if (m.le.0.04) then
          !eliminate coulomb friction for hyperconcentrated/fluid problems
          !klugey, but won't effect debris-flow problems
         !sigbed = sigbed*0.5d0*(tanh(400.d0*(m-0.02d0)) + 1.0d0)
      !endif

      if (m.le.1.d-16) then
         compress = 1.d16
         kperm = 0.0d0
         tanpsi = 0.0d0
         sigbed=0.0d0
      else
         compress = alpha/(m*(sigbed +  sigma_0))
      endif

      !if (m.le.0.4) then
      !   kperm = tanh(10.d0*m)*kperm
      !   tanpsi = tanh(10.d0*m)*tanpsi
      !   sigbed= tanh(10.d0*m)*sigbed
      !endif

      if (p_initialized.eq.0.and.vnorm.le.0.d0) then
      !if (vnorm.le.0.d0) then
         tanpsi = 0.d0
         D = 0.d0
      elseif (h*mu.gt.0.d0) then
         D = 2.0d0*(kperm/(mu*h))*(rho_fp*gmod*h - p)
      else
         D = 0.d0
      endif

      tanphi = dtan(phi_bed + datan(tanpsi))! + phi_seg_coeff*pmtanh01*dtan(phi_bed)
      !if (S.gt.0.0) then
      !   tanphi = tanphi + 0.38*mu*shear/(shear + 0.005*sigbedc)
      !endif

      tau = dmax1(0.d0,sigbed*tanphi)

      !tau = (grav/gmod)*dmax1(0.d0,sigbed*tanphi)
      !kappa: earth pressure coefficient
      !if (phi_int.eq.phi_bed) then
      !   sqrtarg = 0.d0
      !else
      !   sqrtarg = 1.d0-(dcos(phi_int)**2)*(1.d0 + dtan(phi_bed)**2)
      !endif

      !kappa = (2.d0 - pm*2.d0*dsqrt(sqrtarg))/(dcos(phi_int)**2)
      !kappa = kappa - 1.d0
      kappa = 1.d0
      !kappa = 0.4d0

   end subroutine integrate_mp