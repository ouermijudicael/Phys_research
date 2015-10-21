      subroutine im_scat(u, v, w)
c subroutine scat applies the scattering phase function, P(theta), 
c  for isotropic scattering and Rayleigh scattering to find 
c  the photon scattering angles, theta and phi

c  on input, the variables u,v,w are the initial direction cosines,
c  on exit, u,v,w are the direction cosines of the scattered photon.

      pi=acos(-1.0d0)
      twopi=2.0d0*pi

      p12=1.0d0/2.0d0
      p13=1.0d0/3.0d0
      p14=1.0d0/4.0d0

c  set nscat = 0 ===> isotropic scattering
c      nscat = 1 ===> Rayleigh scattering
c      nscat = 2 ===> Henyey-Greenstein scattering
      ! it was 1
      nscat=0

c  isotropic scattering: probabiliy of scattering into a given solid
c                        angle element, dOmega; given by cos(theta)

      if(nscat.eq.0) then 
        a=2.0d0*rand()-1.0d0
      endif

c  rayleigh scattering: probability of scattering into solid
c                       angle element, d(Omega), given by 
c                       (3/4)(cos[theta]+(1/3)cos^3[theta])
          

          
      if(nscat.eq.1) then                        
        phase=2.0d0*(0.5d0-rand())
        q=-4.0d0*phase 
        disc=sqrt(1.0d0+p14*q**2)
        a1=(-p12*q+disc)**p13
        a2=-(p12*q+disc)**p13
        a=a1+a2
      endif
 
c  henyey-greenstein scattering: probability of scattering into 
c                    solid angle element, d(Omega)
          
      if(nscat.eq.2) then
        g=0.1
        g2=g*g
        a=0.5*(1.0+g2-((1.0-g2)/(1.0-g+2.0*g*rand()))**2)/g
      endif

c  axially symmetric scattering

      delta=twopi*rand()

c  find the direction cosines for the scattered photon in the frame
c  where the initial direction of the photon's motion defines the z-axis. 
c  next, transform back to the lab frame to continue the walk.
c
c  perform two rotations to return to the laboratory frame. the first
c  transformation rotates the lab frame about the y-axis to return the
c  z-axis to its original direction. the second rotates the frame about 
c  the new (old) z-axis to return the x and y axes to their original 
c  directions.

      b=sqrt(1.0d0-a*a)
      c=cos(delta)
      d=sin(delta)
      bc=b*c
      bd=b*d
      aw=a*w
      bcw=bc*w
                 
      sinlab=sqrt(1.0d0-w*w)

c  transform direction cosines to lab frame

      uf=(bcw*u-bd*v)/sinlab+a*u
      vf=(bcw*v+bd*u)/sinlab+a*v
      wf=-bc*sinlab+a*w

c  make sure the direction cosines have length 1

      coslength=sqrt(uf*uf+vf*vf+wf*wf)

      u=uf/coslength
      v=vf/coslength
      w=wf/coslength

      return
      end

c======================================================================


