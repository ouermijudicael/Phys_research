
      program monte

c======================================================================
c
c     Version: 2015 September 30
c
c     current(i):   photon position
c     direction(i): photon direction cosines
c
c     x(j),z(k):    radial and vertical grid boundaries
c     xmfp(nq,j,k): photon mean free path in cell (nq,j,k)
c
c======================================================================

      parameter(jmax=100,jmax1=jmax+1,kmax=50,kmax1=kmax+1)
      parameter(nmu=30)
      real*8 current,direction,xmfp
      common/c2/x(jmax1),y(jmax1),z(kmax1)
      common/c1/xmfp(jmax1,jmax1,kmax1),xint(jmax1,jmax1,kmax1),
     &          absal(jmax1,jmax1,kmax1)
      dimension current(3),direction(3)
      
      ! tajo
      dimension S(4)
      
      common/c3/nref(nmu),ntrans(nmu)
 
      data nref/nmu*0/
      data ntrans/nmu*0/

c**********************************************************************

      pi=acos(-1.0d0)
      twopi=2.0d0*pi

c**********************************************************************
c      call traveltest
c      call gridtest
      call scattest
      write(6,*) 'Enter nphot(photons=2**nphot), output interval, seed'
      read(5,10)   nphot,nout,nseed
      print*, "nphot, nout, nseed =",nphot,", ",nout,", ",nseed
      nphot=2**nphot
 10   format(3i12)
      write(6,10) nphot
      rseed=rand(nseed)

      call grid

c======================================================================
c======================================================================

c   Loop over individual photon trajectories. Run nphot photons.

      do nset = 1,nphot
      
      print*, "photon =", nset
c  photon: initial position and direction cosines 

c    initial direction 
 
      w= 2.0*rand()-1.0
      angle=twopi*rand()
      direction(1)=sqrt(1.0-w*w)*cos(angle)
      direction(2)=sqrt(1.0-w*w)*sin(angle)
      direction(3)=w

      xcos=0.0d0
      do i=1,3
        xcos=xcos+direction(i)**2
      enddo
      do i=1,3
        direction(i)=direction(i)/sqrt(xcos)
      enddo

c    initial position 

      nx=jmax/2
      ny=jmax/2
      nz=kmax/2+1
      current(1)=x(nx)
      current(2)=y(ny)
      current(3)=z(nz)
      lunit=8
c     save info for debugging purposes
c  start random walk
      print*, "absal =", absal(nx, ny, nz)
      print*, "internal energy =", xint(nx, ny, nz)

      nabs=0
      n=0
 100  n=n+1

c      print*,"pos=(",current(1),", ",current(2),", ",current(3),")"
c      print*,"dir=(",direction(1),",",direction(2),",",direction(3),")"
c  step: distance to next interaction in units of optical depth
 
      rand_num=rand() 
      step=-alog(rand_num)

c    check: how many cells does the photon cross?
      call travel(nx,ny,nz,step,nesc,current,direction)
      print*,"nesc = ", nesc
      if(nesc.eq.1) go to 1000

c  check: was photon absorbed or scattered?
c      absal(nx, ny, nz) = .01
      rand_num=rand()
c      print*, "rand_num=", rand_num
c      print*, "absal =", absal(nx, ny, nz)
c8      write(8,*) n,rand_num,absal(nx,ny,nz)
      if(rand_num.lt.absal(nx,ny,nz)) then

c    absorbed:
        print*, "the photon was absorbed"
        xint(nx,ny,nz)=xint(nx,ny,nz)+1.0
c8        write(8,*) n,rand_num,absal(nx,ny,nz)
        nabs=1
        go to 1000

      else 

c    scattered:

        u=direction(1)
        v=direction(2)
        w=direction(3)
c        print*, "Before (u, v, w) = (", u, ", ",v, ", ", w, ")"
        !call scat(u,v,w)
        S = (/1.00, 0.00, 0.00, 0.00/)
c        call scat2(u,v,w, S)
        call im_scat(u,v,w)
        direction(1)=u
        direction(2)=v
        direction(3)=w
        print*, "After (u, v, w) = (", u, ", ",v, ", ", w, ")"

      endif

      go to 100 

 1000 continue

c8      write(8,*) nset,n,nabs
c8      write(8,*) nx,ny,nz

      end do

c======================================================================
c======================================================================

      nabs=0
      do i=1,jmax1
        do j=1,jmax1
          do k=1,kmax1
            nabs=nabs+xint(i,j,k)
          enddo
        enddo
      enddo

      nescref=0
      nesctrans=0
      xfrac=1.0/nmu
      do mu=1,nmu
        nescref=nescref+nref(mu)
        nesctrans=nesctrans+ntrans(mu)
        xmu=xfrac*mu
        cosdcos=0.5*(xmu**2-(xmu-xfrac)**2)
        xinv=1.0/cosdcos
        xmu=(180.0/pi)*acos(xmu)
        write(7,3000) xmu,xinv*nref(mu),xinv*ntrans(mu)
 3000   format(1p3e12.4,3i8)
      enddo
      write(lunit,3100) nphot,nescref,nesctrans,nescref+nesctrans,nabs
 3100 format(5i8)

      end 


c====================================================================
c  This subroutine is to test the sactering routine
c====================================================================
       subroutine scattest
       u = 1/sqrt(3.00)
       v = 1/sqrt(3.00)
       w = 1/sqrt(3.00)
c       S = (/1.00, 0.00, 0.00, 0.00/)
       do i=0, 10000
           call scat2(u, v, w)
           print*,"u, v, w = (", u, v, w,")"
       enddo
       end

c======================================================================
!
!      subroutine scat(u,v,w)
!
c  subroutine scat applies the scattering phase function, P(theta), 
c  for isotropic scattering and Rayleigh scattering to find 
c  the photon scattering angles, theta and phi

c  on input, the variables u,v,w are the initial direction cosines,
c  on exit, u,v,w are the direction cosines of the scattered photon.
!
!      pi=acos(-1.0d0)
!      twopi=2.0d0*pi
!
!      p12=1.0d0/2.0d0
!      p13=1.0d0/3.0d0
!      p14=1.0d0/4.0d0
!
c  set nscat = 0 ===> isotropic scattering
c      nscat = 1 ===> Rayleigh scattering
c      nscat = 2 ===> Henyey-Greenstein scattering
!
!      nscat=1
!
c  isotropic scattering: probabiliy of scattering into a given solid
c                        angle element, dOmega; given by cos(theta)
!
!      if(nscat.eq.0) then 
!        a=2.0d0*rand()-1.0d0
!!      endif
!
c  rayleigh scattering: probability of scattering into solid
c                       angle element, d(Omega), given by 
c                       (3/4)(cos[theta]+(1/3)cos^3[theta])
!          
!      if(nscat.eq.1) then                        
!        phase=2.0d0*(0.5d0-rand())
!        q=-4.0d0*phase 
!        disc=sqrt(1.0d0+p14*q**2)
!        a1=(-p12*q+disc)**p13
!        a2=-(p12*q+disc)**p13
!        a=a1+a2
!      endif
! 
!c  henyey-greenstein scattering: probability of scattering into 
!c                    solid angle element, d(Omega)
          
!      if(nscat.eq.2) then
!        g=0.1
!        g2=g*g
!        a=0.5*(1.0+g2-((1.0-g2)/(1.0-g+2.0*g*rand()))**2)/g
!      endif

c  axially symmetric scattering

!      delta=twopi*rand()

c  find the direction cosines for the scattered photon in the frame
c  where the initial direction of the photon's motion defines the z-axis. 
c  next, transform back to the lab frame to continue the walk.
c
c  perform two rotations to return to the laboratory frame. the first
c  transformation rotates the lab frame about the y-axis to return the
c  z-axis to its original direction. the second rotates the frame about 
c  the new (old) z-axis to return the x and y axes to their original 
c  directions.

!      b=sqrt(1.0d0-a*a)
!      c=cos(delta)
!      d=sin(delta)
!
!      bc=b*c
!      bd=b*d
!      aw=a*w
!      bcw=bc*w
!      sinlab=sqrt(1.0d0-w*w)

c  transform direction cosines to lab frame

!      uf=(bcw*u-bd*v)/sinlab+a*u
!      vf=(bcw*v+bd*u)/sinlab+a*v
!      wf=-bc*sinlab+a*w

c  make sure the direction cosines have length 1

!      coslength=sqrt(uf*uf+vf*vf+wf*wf)

!      u=uf/coslength
!      v=vf/coslength
!      w=wf/coslength

!      return
!      end


c======================================================================
c testing grid routine 
c=======================================================================
      subroutine gridtest
     
      parameter(jmax=100,jmax1=jmax+1,kmax=50,kmax1=kmax+1)
      common/c2/x(jmax1),y(jmax1),z(kmax1)
      common/c1/xmfp(jmax1,jmax1,kmax1),xint(jmax1,jmax1,kmax1),
     &          absal(jmax1,jmax1,kmax1)
      call grid
      nx=jmax/2 - 1
      ny=jmax/2 - 1
      nz=kmax/2+1 - 1
      print*, "(nx, ny, nz)= (",nx,",",ny,",",",",nz,")"
      print*, "(x,y,z) =",x(nx),",",y(ny),",",z(nz)
      print*, "xmfp =", xmfp(nx,ny,nz)
      print*, "xint =", xint(nx, ny, nz)
      print*, "absal =", absal(nx,ny,nz)
      end

c======================================================================

      subroutine grid
      parameter(jmax=100,jmax1=jmax+1,kmax=50,kmax1=kmax+1)
      common/c2/x(jmax1),y(jmax1),z(kmax1)
      common/c1/xmfp(jmax1,jmax1,kmax1),xint(jmax1,jmax1,kmax1),
     &          absal(jmax1,jmax1,kmax1)

c  set up the computational grid
      print*, "Entering grid routine"
      deltax=1.0
      deltay=1.0
      deltaz=1.0
      do j=1,jmax1
        x(j)=(j-1)*deltax
        y(j)=(j-1)*deltay
      enddo 
      do k=1,kmax1
        z(k)=(k-1-kmax/2)*deltaz
      enddo 

c  set mean free path, xmfp(q,j,k), in each compuational cell, (j,k)

      taumax=10.0
      epsilon=exp(2.0*alog(taumax)/kmax)
      write(6,*) taumax,epsilon
      do j=1,jmax1
!!
!        epsilon=1.0
!!
        do i=1,jmax1
          delta=1.0
          do k=kmax/2+1,kmax
            xmfp(i,j,k)=1.0*delta
            xmfp(i,j,kmax1-k)=xmfp(i,j,k)
            delta=delta*epsilon
            if(k.ge.3*kmax/4) delta=1.0e6
          enddo
          xmfp(i,j,kmax1)=xmfp(i,j,kmax)
        enddo
      enddo

c  carve a hole in the center of the slab

      ncut=5
      ncut=ncut*ncut
      nout=jmax/2
      nout=nout*nout
      do i=1,jmax1
        do j=1,jmax1
          nradius=(i-jmax/2)**2+(j-jmax/2)**2
          soften=1.0
          if(nradius.le.ncut) soften=1.0e06
          if(nradius.ge.nout) soften=1.0e04
          tauz=0.0
          xradius=nradius
          xout=nout
          xopen=(1.0+xradius/xout) 
          do k=1,kmax1
!!
           xmfp(i,j,k)=soften*xmfp(i,j,k)*xopen
!!
            if(i.eq.jmax/2+1)  then 
              write(10,*) j,k,deltax/xmfp(i,j,k)
              tauz=tauz+deltax/xmfp(i,j,k)
            endif
          enddo
          if(i.eq.jmax/2+1)  write(11,*) j,tauz
        enddo
      enddo
c
      write(99,100) (k,z(k),xmfp(jmax/2,jmax/2,k),k=1,kmax1)
 100  format(i8,1p2e12.4)

c  initialize the internal energy array

      do j=1,jmax1
        do i=1,jmax1
          do k=1,kmax1
            test1 = 20
            xint(j,i,k)=0.0d00
          enddo
        enddo
      enddo

c  set absorption-to-(absorption+scattering) fraction

      do j=1,jmax1
        do i=1,jmax1
          do k=1,kmax1
            test2 = 30
            absal(j,i,k)=0.01
          enddo
        enddo
      enddo
c
      print*, "Leaving grid routine"
      return
      end 

!==========================================================================!
! Testing the travel routine                                               !
!==========================================================================!

      subroutine traveltest
      parameter(jmax=100,jmax1=jmax+1,kmax=50,kmax1=kmax+1)
      parameter(nmu=30)
      real*8 current,direction,xmfp
      common/c2/x(jmax1),y(jmax1),z(kmax1)
      common/c1/xmfp(jmax1,jmax1,kmax1),xint(jmax1,jmax1,kmax1),
     &          absal(jmax1,jmax1,kmax1)
      dimension current(3),direction(3)
      common/c3/nref(nmu),ntrans(nmu)

      call grid

c     initial position 
      nx=jmax/2
      ny=jmax/2
      nz=kmax/2+1
      current(1)=x(nx)
      current(2)=y(ny)
      current(3)=z(nz)
      
      direction(1)= 1 / sqrt(3.00)
      direction(2)= 1 / sqrt(3.00)
      direction(3)= 1 / sqrt(3.00)
      step = -alog(rand())
      print*, " ------BEFORE TRAVEL------"
      print*, "position ="
      print*,current(1),", ",current(2),", ",current(3)
      print*, "direction ="
      print*,direction(1),",",direction(2),"," ,direction(3)
      print*,"(nx, ny, nz) ="
      print*,nx,", ", ny,", ", nz
      print*, "nesc =", nesc
      print*, "step =", step
      call travel(nx,ny,nz,step,nesc,current,direction)
      print*, " ------AFTER TRAVEL------"
      print*, "position ="
      print*,current(1),", ",current(2),", ",current(3)
      print*, "direction ="
      print*,direction(1),",",direction(2),"," ,direction(3)
      print*,"(nx, ny, nz) ="
      print*,nx,", ", ny,", ", nz
      print*, "nesc =", nesc
      print*, "step =", step
      
      end     



      subroutine travel(nx,ny,nz,step,nesc,current,direction)

c======================================================================
c
c     subroutine travel advances the photon in the grid, starting
c     from its current location (x,y,z), current, moving in direction, 
c     direction, (U,V,W).
c
c     x(nx),y(ny),z(nz):  radial and vertical grid boundaries
c     xmfp(nx,ny,nz):     photon mean free path in cell (nx,ny,nz)
c
c======================================================================

      parameter(jmax=100,jmax1=jmax+1,kmax=50,kmax1=kmax+1)
      parameter(nmu=30)
      real*8 current,direction,xmfp
      common/c2/x(jmax1),y(jmax1),z(kmax1)
      common/c1/xmfp(jmax1,jmax1,kmax1),xint(jmax1,jmax1,kmax1),
     &          absal(jmax1,jmax1,kmax1)
      dimension current(3),direction(3)
      common/c3/nref(nmu),ntrans(nmu)

c***********************************************************************
c      print*, "Entering travel "
      nesc=0
      stepsum=0.0

c   is the photon moving up or down?

      if(direction(3).gt.0.0) then

c     is the photon moving left or right?

        if(direction(1).gt.0.0) then

c       is the photon moving forward or backward?

          if(direction(2).gt.0.0) then 

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 110        continue                        ! up and right and forward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx+1)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny+1)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz+1)-current(3))/xmfp(nx,ny,nz)

            if(stepx.gt.taux .or. stepy.gt.tauy .or. stepz.gt.tauz) then
              sx=(x(nx+1)-current(1))/direction(1)
              sy=(y(ny+1)-current(2))/direction(2)
              sz=(z(nz+1)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz+1
                if(nz.gt.kmax1) then
                  mu=nmu*direction(3)+1
                  nref(mu)=nref(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx+1
                if(nx.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny+1
                if(ny.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 110

          else

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 111        continue                        ! up and right and backward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx+1)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz+1)-current(3))/xmfp(nx,ny,nz)

            if(stepx.gt.taux .or. stepy.lt.tauy .or. stepz.gt.tauz) then
              sx=(x(nx+1)-current(1))/direction(1)
              sy=(y(ny)-current(2))/direction(2)
              sz=(z(nz+1)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz+1
                if(nz.gt.kmax1) then
                  mu=nmu*direction(3)+1
                  nref(mu)=nref(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx+1
                if(nx.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny-1
                if(ny.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 111

          endif

        else

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c         is photon going forward or backward?

          if(direction(2).gt.0.0) then 

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 112        continue                        ! up and left and forward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny+1)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz+1)-current(3))/xmfp(nx,ny,nz)

            if(stepx.lt.taux .or. stepy.gt.tauy .or. stepz.gt.tauz) then
              sx=(x(nx)-current(1))/direction(1)
              sy=(y(ny+1)-current(2))/direction(2)
              sz=(z(nz+1)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz+1
                if(nz.gt.kmax1) then
                  mu=nmu*direction(3)+1
                  nref(mu)=nref(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx-1
                if(nx.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny+1
                if(ny.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 112

          else

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 113        continue                        ! up and left and backward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz+1)-current(3))/xmfp(nx,ny,nz)

            if(stepx.lt.taux .or. stepy.lt.tauy .or. stepz.gt.tauz) then
              sx=(x(nx)-current(1))/direction(1)
              sy=(y(ny)-current(2))/direction(2)
              sz=(z(nz+1)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz+1
                if(nz.gt.kmax1) then
                  mu=nmu*direction(3)+1
                  nref(mu)=nref(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx-1
                if(nx.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny-1
                if(ny.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 113

          endif

        endif

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      else

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c       is the photon moving left or right?

        if(direction(1).gt.0.0) then

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c         is the photon moving forward or backward?

          if(direction(2).gt.0.0) then 

 120        continue                        ! down and right and forward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx+1)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny+1)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz)-current(3))/xmfp(nx,ny,nz)

            if(stepx.gt.taux .or. stepy.gt.tauy .or. stepz.lt.tauz) then
              sx=(x(nx+1)-current(1))/direction(1)
              sy=(y(ny+1)-current(2))/direction(2)
              sz=(z(nz)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz-1
                if(nz.lt.1) then
                  mu=-nmu*direction(3)+1
                  ntrans(mu)=ntrans(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx+1
                if(nx.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny+1
                if(ny.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cisde                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 120

          else

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 121        continue                        ! down and right and backward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx+1)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz)-current(3))/xmfp(nx,ny,nz)

            if(stepx.gt.taux .or. stepy.lt.tauy .or. stepz.lt.tauz) then
              sx=(x(nx+1)-current(1))/direction(1)
              sy=(y(ny)-current(2))/direction(2)
              sz=(z(nz)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz-1
                if(nz.lt.1) then
                  mu=-nmu*direction(3)+1
                  ntrans(mu)=ntrans(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx+1
                if(nx.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny-1
                if(ny.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 121

          endif

        else

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c         is photon going forward or backward?

          if(direction(2).gt.0.0) then 

 122        continue                        ! down and left and forward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny+1)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz)-current(3))/xmfp(nx,ny,nz)

            if(stepx.lt.taux .or. stepy.gt.tauy .or. stepz.lt.tauz) then
              sx=(x(nx)-current(1))/direction(1)
              sy=(y(ny+1)-current(2))/direction(2)
              sz=(z(nz)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz-1
                if(nz.lt.1) then
                  mu=-nmu*direction(3)+1
                  ntrans(mu)=ntrans(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx-1
                if(nx.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cisde                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny+1
                if(ny.gt.jmax) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 122

          else

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 123        continue                        ! down and left and backward

            stepx=(step-stepsum)*direction(1)
            stepy=(step-stepsum)*direction(2)
            stepz=(step-stepsum)*direction(3)
            taux=(x(nx)-current(1))/xmfp(nx,ny,nz)
            tauy=(y(ny)-current(2))/xmfp(nx,ny,nz)
            tauz=(z(nz)-current(3))/xmfp(nx,ny,nz)

            if(stepx.lt.taux .or. stepy.lt.tauy .or. stepz.lt.tauz) then
              sx=(x(nx)-current(1))/direction(1)
              sy=(y(ny)-current(2))/direction(2)
              sz=(z(nz)-current(3))/direction(3)
c
              nsmin=3
              if(sx.lt.sy) then 
                if(sx.lt.sz) nsmin=1
              else
                if(sy.lt.sz) nsmin=2
              endif
c
              if(nsmin.eq.3) then
                stepsum=stepsum+tauz/direction(3)
                current(1)=current(1)+sz*direction(1)
                current(2)=current(2)+sz*direction(2)
                current(3)=current(3)+sz*direction(3)
                nz=nz-1
                if(nz.lt.1) then
                  mu=-nmu*direction(3)+1
                  ntrans(mu)=ntrans(mu)+1
                  goto 1000 
                endif
              endif

              if(nsmin.eq.1) then
                stepsum=stepsum+taux/direction(1)
                current(1)=current(1)+sx*direction(1)
                current(2)=current(2)+sx*direction(2)
                current(3)=current(3)+sx*direction(3)
                nx=nx-1
                if(nx.lt.1) then
                  mu=nmu*direction(3)
                  if(mu.gt.0) then 
cside                    nref(mu+1)=nref(mu+1)+1
                  else
cside                    ntrans(-mu+1)=ntrans(-mu+1)+1
                  endif
                  goto 1000 
                endif
              endif

              if(nsmin.eq.2) then
                stepsum=stepsum+tauy/direction(2)
                current(1)=current(1)+sy*direction(1)
                current(2)=current(2)+sy*direction(2)
                current(3)=current(3)+sy*direction(3)
                ny=ny-1
                if(ny.lt.1) then
                  mu=nmu*direction(3)+1
                  if(mu.gt.0) then
cside                    nref(mu)=nref(mu)+1
                  else
cside                    ntrans(-mu)=ntrans(-mu)+1
                  endif
                  goto 1000 
                endif
              endif

            else

              current(1)=current(1)+stepx*xmfp(nx,ny,nz)
              current(2)=current(2)+stepy*xmfp(nx,ny,nz)
              current(3)=current(3)+stepz*xmfp(nx,ny,nz)
              go to 200

            endif
            go to 123

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          endif

        endif

      endif

 200    return  
c      print*, "leaving travel"
c      print*,"pos=(",current(1),", ",current(2),", ",current(3),")"
c      print*,"dir=(",direction(1),",",direction(2),",",direction(3),")"
c      return

 1000 nesc=1
c      print*, "Leaving travel "
c      print*,"pos=(",current(1),", ",current(2),", ",current(3),")"
c      print*,"dir=(",direction(1),",",direction(2),",",direction(3),")"
      return

      end
