      subroutine partsim(tgam,tthet,dgam,dthet,istepz)
c     ==================================================================
c     define the system of ode of the canonic variables
c     t... are the values of the variables at current position
c     d... are the value of the differentioan equation
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'particle.cmn'
      include 'work.cmn'
c
      complex*16 ctmp
      integer ip,istepz,ih,nharmpart
      real*8 tgam,tthet,dgam,dthet,ztemp1,ztemp2,btper0
      dimension tgam(*),tthet(*)
      dimension dgam(*),dthet(*)
c
c     diff-eq. for longitudinal variables:
c        - energy (gamma)
c        - phase (theta) 
c     ------------------------------------------------------------------
c
c     
      ztemp1=-2.d0*xlamds/xlamd
      ztemp2=xlamd/xlamds
c
      nharmpart=1        !default - only fundamentalacts on electron
      if (iharmsc.ne.0) then   
       nharmpart=nharm   ! self-consistent equation of motion
      endif
c
c      nharmpart=1
c
      do ip=1,npart
c       
        ctmp=dcmplx(0.0)
        do ih =1,nharmpart ! loop over harmonics if self-consitent method is enabled
          ctmp=ctmp+cpart1(ip+(ih-1)*npart)
     +     *dcmplx(dcos(dble(ih)*tthet(ip)),-dsin(dble(ih)*tthet(ip)))
        enddo ! ih
c
        btper0=btper(ip)+ztemp1*dble(ctmp)  !perpendicular velocity
        btpar(ip)=dsqrt(1.d0-btper0/tgam(ip)**2)     !parallel velocity
c
        dthet(ip)=dthet(ip)+ztemp2*(1.-1./btpar(ip))+1.      !dtheta/dz
        dgam (ip)=dgam(ip) +dimag(ctmp)/btpar(ip)/tgam(ip)-ez(ip)  !dgamma/dz
      end do       ! ip
      return
      end     !partsim
c

      subroutine partsorc(istepz)
c     ==================================================================
c     calculates source term for gamma-theta integration 
c     when higher harmonic coupling is considered
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'field.cmn'
      include 'input.cmn'
      include 'magnet.cmn'
      include 'particle.cmn'
      include 'sim.cmn'
      include 'work.cmn'
c
      complex*16 clocal
      integer ip,idx,istepz,ih,ix,nharmpart
      integer idx1,idx2,idx3,idx4,ioffset
      real*8 wei1,wei2,wei3,wei4,rtmp
      real*8  zconst,wei,aw2,xi

c      
c     call space charge routine to calculate the local field Ez
c
c
      nharmpart=1        !default - only fundamentalacts on electron
      if (iharmsc.ne.0) then   
       nharmpart=nharm   ! self-consistent equation of motion
      endif

c
c     check for drifts
c
      if (awz(istepz).lt.tiny) then
         aw2=awdz(istepz)*awdz(istepz)
         do ip=1,npart
           btper(ip)=aw2+px(ip)*px(ip)+py(ip)*py(ip)+1.
         enddo
         do ip=1,npart*nharmpart
           cpart1(ip)=dcmplx(0.0)
         enddo
         return   
      endif
 
c
c      xpart(1)=0
c      ypart(1)=0

      call rpos(istepz,xpart,ypart)                          !position of particles  
c
c     tmp
c
      do ip=1,npart
        aw2=faw2(istepz,xporb(ip),yporb(ip))     !square of wiggler field amplitude
        aw2=aw2+awdz(istepz)*awdz(istepz)             !artificial delay in drifts 
        btper(ip)=aw2+px(ip)*px(ip)+py(ip)*py(ip)+1.       
        rtmp=sqrt(aw2)*xlamds/xlamd               ! effective K-parameter
c
c       interpolation to the grid (index and weight of the 4 surrounding grid points)
c  
        wei1=wx(ip)*wy(ip)
        idx1=ipos(1,ip)
        wei2=wx(ip)*(1.d0-wy(ip))    
        idx2=ipos(2,ip)
        wei3=(1.d0-wx(ip))*wy(ip)
        idx3=ipos(3,ip)
        wei4=(1.d0-wx(ip))*(1.d0-wy(ip))
        idx4=ipos(4,ip)
        do ih=1,nharmpart,2                               !sum over odd harmonics
           ioffset=(ih-1)*ncar*ncar        
           clocal=wei1*crfield(idx1+ioffset)              ! note Genesis field u is related to electric field by
           clocal=clocal+wei2*crfield(idx2+ioffset)       ! u = k*n/ k_u^2 * e/mc^2  E_0
           clocal=clocal+wei3*crfield(idx3+ioffset)
           clocal=clocal+wei4*crfield(idx4+ioffset)
           cpart1(ip+(ih-1)*npart)=conjg(clocal)*rtmp
     +           *besselcoupling(ih)/dble(ih)
        enddo
        do ih=2,nharmpart,2                               !sum over even harmonics
           xi=dble(ih)*xlamd/xlamds/gamma(ip)/gamma(ip)*dsqrt(aw2)      ! nk K/gamma/k_u * x' = coupling -> x' = px/gamma
           ioffset=(ih-1)*ncar*ncar        
           clocal=wei1*crfield(idx1+ioffset)              ! note Genesis field u is related to electric field by
           clocal=clocal+wei2*crfield(idx2+ioffset)       ! u = k*n/ k_u^2 * e/mc^2  E_0
           clocal=clocal+wei3*crfield(idx3+ioffset)
           clocal=clocal+wei4*crfield(idx4+ioffset)
           cpart1(ip+(ih-1)*npart)=conjg(clocal)*rtmp
     +           *besselcoupling(ih)/dble(ih)
     +           *xi*dcmplx(0.,sqrt(2.)*px(ip))          ! missing px here and add I*py for helical undulator 
           
        enddo
      enddo    
C      write(*,*) cpart1(1),rtmp,besselcoupling(1)
c
      return
      end     !partsim     
c

c 
      subroutine harmcoupling(awloc)
c     ============================================================
c     routine to calculate the coupling to higher modes
c     ------------------------------------------------------
c
      include 'genesis.def'
      include 'field.cmn' 
      include 'input.cmn'
c
      real*8 xi,awloc
      integer ih
c
      if (awloc.lt.tiny) then      !drift -> set coupling to zero
        do ih=1,nharm
          besselcoupling(ih)=0.
        enddo
        return
      endif
c
      xi = awloc**2/(1.d0+awloc**2)/2
c
      if (iwityp.ne.0) then         ! helical undulator 
        besselcoupling(1)=1.
        do ih=2,nharm               ! even harmonic is disabled due to the possible
          besselcoupling(ih)=0.     ! break in symmetry
        enddo
      else
         do ih=1,nharm
           if(mod(ih,2).eq.1) then  
            besselcoupling(ih) = (bessj((ih-1)/2,xi*dble(ih))
     +          -bessj((ih+1)/2,xi*dble(ih)))*((-1)**((ih-1)/2))
           else    
            besselcoupling(ih)=0.5*(bessj((ih-2)/2,xi*dble(ih))
     +          -bessj((ih+2)/2,xi*dble(ih)))*(-1)**((ih-2)/2) 
           endif                !note that for full coupling it has to 
                                !mulitplied with the transverse momentum 
         enddo
      endif
      return
      end        ! of harmcoupling
