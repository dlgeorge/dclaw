

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
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

      ! Output
      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      real(kind=8) :: D,tau,sigbed,kperm,compress,pm,coeff
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,taucf,fsphi,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,chi,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining,dtcoeff
      real(kind=8) :: gamma,zeta,krate,p_eq,dgamma

      integer :: i,j,ii,jj,icount
      logical :: ent


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

      gmod=grav

      if (friction_forcing) then
         coeff = manning_coefficient(1)
      else
         coeff = 0.d0
      endif

      ! Current implementation of friction has manning as an array
      ! take the first element for now. If only one value is
      ! provided to geo_data.manning_coefficient
      ! it will be manning_coefficient(1)
      ! DIG: FIX.

      if (entrainment>0) then
         ent = .true.
      else
         ent = .false.
      endif

      do j=1-mbc+1,my+mbc-1
         do i=1-mbc+1,mx+mbc-1
         ! DIG: 1/12/24: KRB and MJB notice that here we are looping over ghost cells.
         ! These ghost cells have not been updated by the reimann solver? Are they used
         ! meaningfully (e.g., theta is diffed below. )

            ! adjust gravity if bed_normal = 1
            theta = 0.d0
            dtheta = 0.d0
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gmod = grav*cos(theta)
               dtheta = -(aux(i_theta,i+1,j) - theta)/dx
            endif

            ! Get state variable
            h = q(1,i,j)
            if (h<=dry_tolerance) cycle
            hu = q(2,i,j)
            hv = q(3,i,j)
            hm = q(4,i,j)
            p =  q(5,i,j)
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
            endif

            if (hvnorm0>0.d0) then
               hu = hvnorm*hu/hvnorm0
               hv = hvnorm*hv/hvnorm0
            endif

            ! call admissible q and auxeval before moving on to shear induced
            ! dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            ! calculate velocity
            vnorm = sqrt(u**2 + v**2)

            !integrate shear-induced dilatancy
            shear = 2.d0*vnorm/h
            krate = 1.5d0*shear*m*tanpsi/alpha
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dt*3.d0*vnorm*tanpsi/(h*compress)
            endif

            !DIG move segregation to segeval
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               pmtanh01=0.d0
      		else
         		seg = 1.d0
               call calc_pmtanh(pm,seg,pmtanh01)
               rho_fp = max(0.d0,(1.d0-pmtanh01))*rho_f
      		endif

            !integrate pressure relaxation
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.d0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
            p_eq = h*rho_fp*gmod
            p = p_eq + (p-p_eq)*exp(krate*dt)

            ! call admissible q and auxeval before moving on to dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            ! calculate rate of change
            krate = D*(rho-rho_fp)/rho

            ! integrate hu, hv, hm, and h.
            hu = hu*exp(dt*krate/h)
            hv = hv*exp(dt*krate/h)
            hm = hm*exp(-dt*D*rho_fp/(h*rho))
            h = h + krate*dt


            !======================mass entrainment===========================

            ! call admissible q and auxeval before moving on to mass entrainment.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            vnorm = sqrt(u**2 + v**2)
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
                     m2 = 0.6d0
                     rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0
                     beta2 = 0.66d0
                     t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d0-2.d0))
                     !write(*,*) '------------'
                     !write(*,*) 'vu',t1bot
                     beta = 1.d0-m!tanh(10.d0*m) !tan(1.5d0*p/(rho*gmod*h))/14.0d0
                     !! DIG not used. segfault. gamma= rho*beta2*(vnorm**2)*(beta*gmod*coeff**2)/(tanh(h+1.d0-2.d0)**(1.0d0/3.0d0))
                     !write(*,*) 'gamma', gamma
                     !! DIG not used. segfault. t1bot = t1bot + gamma
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
                     !write(*,*) 'dh/dt', dh/dt
                     dh = min(dh,b_remaining)
                     h = h + dh
                     hm = hm + dh*m2
                     q(7,i,j) = q(7,i,j) + dh

                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                     call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
                     p = prat*rho*h
                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  endif
               endif
            endif

            !===================================================================
            ! end of entrainment, put state variables back in q.

            q(1,i,j) = h
            q(2,i,j) = hu
            q(3,i,j) = hv
            q(4,i,j) = hm
            q(5,i,j) = p
            q(6,i,j) = pm*h

         enddo
      enddo



      ! Manning friction------------------------------------------------
      if (friction_forcing) then
      if (coeff>0.d0.and.friction_depth>0.d0) then

         do j=1,my
            do i=1,mx

               if (bed_normal==1) gmod = grav*cos(aux(i_theta,i,j))
                  h=q(1, i,j)
               if (h<=friction_depth) then
                 !# apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)
                  hm = q(4,i,j)
                  p =  q(5,i,j)
                  phi = aux(i_phi,i,j)
                  theta = aux(i_theta,i,j)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<dry_tolerance) cycle
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
      end
