c======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                 ql,qr,auxl,auxr,fwave,s,amdq,apdq)
c======================================================================
c
c Solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      # This Riemann solver is for debris flow equations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module, only: g => grav, drytol => dry_tolerance
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      use digclaw_module

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local
      integer m,i,mw,maxiter,mhu,nhv,mcapa,icom,jcom,waves
      double precision dtcom,dxcom,dycom,tcom
      double precision wall(3),fw(6,3),sw(3),wave(6,3)
      double precision lamL(3),lamR(3),beta(3)
      !logical entropy(5)
      logical rare1,rare2,wallprob,drystate
      !logical entropycorr1,entropycorr2

      double precision drytol,gmod,veltol
      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL
      double precision pR,pL,hmL,hmR,mL,mR,phi_bedL,phi_bedR
      double precision hstar,hstartest,s1m,s2m,bL,bR
      double precision dxdc,dx,taudirL,taudirR
      double precision theta,thetaL,thetaR
      double precision h1M,h2M,hu1M,hu2M,u1M,u2M,heR,heL
      double precision sE1,sE2
      double precision chiHL,chiHR,chiL,chiR,fsL,fsR

      gmod=grav
      veltol = 1.d-3
      waves = 3

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

      !-----------------------Initializing-----------------------------------
         !inform of a bad Riemann problem from the start
c         if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
c            write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
c         endif

         !Initialize Riemann problem for grid interface

         do mw=1,mwaves
            s(i,mw)=0.d0
            !entropy(mw)=.false.
            do m=1,meqn
               fwave(i,m,mw)=0.d0
            enddo
         enddo
         do mw=1,waves
            sw(mw) = 0.d0
            do m=1,6
               wave(m,mw) = 0.d0
               fw(m,mw) = 0.d0
            enddo
         enddo

         !skip problem if in a completely dry area
         if (qr(i-1,1).le.drytol.and.ql(i,1).le.drytol) then
            go to 30
         endif

c        !set normal direction
         if (ixy.eq.1) then
            mhu=2
            nhv=3
            dx = dxcom
            taudirR = auxl(i,i_taudir_x)
            taudirL = auxr(i-1,i_taudir_x)
         else
            mhu=3
            nhv=2
            dx = dycom
            taudirR = auxl(i,i_taudir_y)
            taudirL = auxr(i-1,i_taudir_y)
         endif

         fsL = auxr(i-1,i_fsphi)
         fsR = auxl(i,i_fsphi)

         if (bed_normal.eq.1) then
            thetaL = auxr(i-1,i_theta)
            thetaR = auxl(i,i_theta)
            theta = 0.5d0*(thetaL+thetaR)
            gmod = grav*dcos(0.5d0*(thetaL+thetaR))
         else
            thetaL = 0.d0
            thetaR = 0.d0
            theta = 0.d0
         endif

         !zero (small) negative values if they exist and set velocities
         call admissibleq(ql(i,1),ql(i,mhu),ql(i,nhv),
     &            ql(i,4),ql(i,5),uR,vR,mR,thetaR)

         call admissibleq(qr(i-1,1),qr(i-1,mhu),qr(i-1,nhv),
     &            qr(i-1,4),qr(i-1,5),uL,vL,mL,thetaL)


         !Riemann problem variables
         hL = qr(i-1,1)
         hR = ql(i,1)
         huL = qr(i-1,mhu)
         huR = ql(i,mhu)
         hvL=qr(i-1,nhv)
         hvR=ql(i,nhv)
         hmL = qr(i-1,4)
         hmR = ql(i,4)
         pL = qr(i-1,5)
         pR = ql(i,5)
         bL = auxr(i-1,1)  - qr(i-1,7)
         bR = auxl(i,1) - ql(i,7)
         phi_bedL = auxr(i-1,i_phi)
         phi_bedR = auxl(i,i_phi)
         chiHL = qr(i-1,6)
         chiHR = ql(i,6)

         if (hL.ge.drytol) then
            chiL = chiHL/hL
         endif
         if (hR.ge.drytol) then
            chiR = chiHR/hR
         endif

         !test for wall problem vs. inundation problem
         do mw=1,waves
            wall(mw) = 1.d0
         enddo
         drystate=.false.
         wallprob = .false.
         if (hR.le.drytol) then
            hR = 0.d0
            pR = 0.d0
            chiR = chiL
            drystate=.true.
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                 rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               wallprob=.true.
               hR=hL
               huR=-huL
               hvR=hvL
               hmR=hmL
               bR=bL
               uR=-uL
               vR=vL
               mR=mL
               pR=pL
               chiHR=chiHL
               chiR = chiL
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hL+bL.lt.bR) then
               !bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            hL = 0.d0
            pL = 0.d0
            chiL= chiR
            drystate=.true.
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               wallprob=.true.
               hL=hR
               huL=-huR
               hvL=hvR
               hmL=hmR
               mL = mR
               bL=bR
               uL=-uR
               vL=vR
               pL=pR
               chiHL=chiHR
               chiL = chiR
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hR+bR.lt.bL) then
               !bL=hR+bR
            endif
         endif
         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

          ! current dclaw Riemann solver
          call riemann_dig2_aug_sswave_ez(ixy,6,3,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         thetaL,thetaR,phi_bedL,phi_bedR,dx,sw,fw,wave,wallprob,
     &         taudirL,taudirR,chiL,chiR,fsL,fsR)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,6
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

         s(i,1) = sw(1)
         s(i,2) = sw(2)
         s(i,3) = sw(2)
         s(i,4) = sw(2)
         s(i,5) = sw(3)

         fwave(1,i,1) =   fw(1,1)
         fwave(1,i,mhu) = fw(2,1)
         fwave(1,i,nhv) = fw(3,1)
         fwave(1,i,4)   = fw(4,1)
         fwave(1,i,5) =   fw(5,1)
         fwave(1,i,6) =   fw(6,1)

         fwave(5,i,1) =   fw(1,3)
         fwave(5,i,mhu) = fw(2,3)
         fwave(5,i,nhv) = fw(3,3)
         fwave(5,i,4)   = fw(4,3)
         fwave(5,i,5) =   fw(5,3)
         fwave(5,i,6) =   fw(6,3)

         fwave(2,i,1) =   fw(1,2)
         fwave(2,i,mhu) = fw(2,2)
         fwave(2,i,nhv) = fw(3,2)
         fwave(2,i,4)   = 0.0
         fwave(2,i,5) =  0.0
         fwave(2,i,6) = fw(6,2)

         fwave(3,i,1) =   0.0
         fwave(3,i,mhu) = 0.0
         fwave(3,i,nhv) = 0.0
         fwave(3,i,4)   = fw(4,2)
         fwave(3,i,5) =  0.0
         fwave(3,i,6) =  0.0

         fwave(4,i,1) =   0.0
         fwave(4,i,mhu) = 0.0
         fwave(4,i,nhv) = 0.0
         fwave(4,i,4)   = 0.0
         fwave(4,i,5) =  fw(5,2)
         fwave(4,i,6) =  0.0

 30      continue
      enddo


c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxl(3,i))*deg2rad
          endif

          do mw=1,mwaves
c             if (s(mw,i) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
c                endif
               s(mw,i)=dxdc*s(mw,i)
               do m=1,meqn
                  fwave(m,mw,i)=dxdc*fwave(m,mw,i)
               enddo
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         amdq(1:meqn,:) = 0.d0
         apdq(1:meqn,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
               if (s(mw,i) < 0.d0) then
                     amdq(1:meqn,i) = amdq(1:meqn,i) + fwave(1:meqn,mw,i)
               else if (s(mw,i) > 0.d0) then
                  apdq(1:meqn,i)  = apdq(1:meqn,i) + fwave(1:meqn,mw,i)
               else
                 amdq(1:meqn,i) = amdq(1:meqn,i) + 0.5d0 * fwave(1:meqn,mw,i)
                 apdq(1:meqn,i) = apdq(1:meqn,i) + 0.5d0 * fwave(1:meqn,mw,i)
               endif
            enddo
         enddo
!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo

      return
      end subroutine
