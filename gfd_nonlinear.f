      program toy

c     thermal convection (nonlinear)

      parameter (nz=201,nn=2,nx=401)

      dimension psi(nz,0:nn),omg(nz,0:nn),tem(nz,0:nn),
     $domgdt1(nz,0:nn),dtemdt1(nz,0:nn),
     $domgdt2(nz,0:nn),dtemdt2(nz,0:nn),
     $z(nz),sub(nz),dia(nz),sup(nz),work1(nz),work2(nz)
      character ifin*1,ifmov*1

      dimension temmov(nx,nz),psimov(nx,nz),
     $cosa(0:nn,nx),sina(0:nn,nx),x(nx)
      real*4 temmov,psimov,ra,pr,aspect
      integer*4 nx4,nz4
      character movfile*10,movfileid*4,movstr*6

c    -set up

      write(6,*) 'aspect ratio'
      read(5,*) aspect
      write(6,*) 'rayleigh number'
      read(5,*) ra
      write(6,*) 'prandtl number'
      read(5,*) pr

      dz=1./float(nz-1)
      oodz2=1./(dz*dz)

      dtmin=min(dz*dz/(4.*pr),dz*dz/4.)
   10 write(6,*) 'timestep'
      read(5,*) dt
      if(dt .le. 0.) stop
      if(dt .ge. dtmin) then
         write(6,*) 'dt must be .lt. ',dtmin
         go to 10
      endif

      write(6,*) 'number of timesteps'
      read(5,*) nstep
      write(6,*) 'print max values every ? timesteps'
      read(5,*) iprnt

      pi=4.*atan(1.)
      c=pi/aspect
      c1=0.25*c/dz
      nx4=nx
      nz4=nz

      sub(1)=0.
      sup(1)=0.
      do k=2,nz-1
         sub(k)=-oodz2
         sup(k)=-oodz2
      enddo
      sub(nz)=0.
      sup(nz)=0.

      z(1)=0.
      do k=2,nz-1
         z(k)=float(k-1)*dz
      enddo
      z(nz)=1.

c    -initial conditions

      do n=0,nn
         do k=1,nz
            psi(k,n)=0.
            omg(k,n)=0.
            tem(k,n)=0.
            dtemdt1(k,n)=0.
            dtemdt2(k,n)=0.
            domgdt1(k,n)=0.
            domgdt2(k,n)=0.
         enddo
      enddo

      write(6,*) 'is there an input file "in"?'
      read(5,*) ifin
      if(ifin .eq. 'y') then
         open(8,file='in',status='old',form='unformatted')
         read(8) tem,omg,psi
         read(8) dtemdt2,domgdt2
      else
c       -initializing only the n = 0 and 1 modes
         amp=0.5
         do k=1,nz
            tem(k,0)=1.-z(k)
            tem(k,1)=amp*sin(pi*z(k))
         enddo
      endif

c    -movie stuff

      write(6,*) 'make movie?'
      read(5,*) ifmov
      if(ifmov .eq. 'y') then
         write(6,*) 'save movie data every ? steps'
         read(5,*) imovie
         movfileid='mov.'
         ooaspect=1./aspect
         dx=aspect/float(nx-1)
         x(1)=0.
         do i=2,nx-1
            x(i)=float(i-1)*dx
         enddo
         x(nx)=aspect
         do n=0,nn
            do i=1,nx
               sina(n,i)=sin(float(n)*pi*x(i)*ooaspect)
               cosa(n,i)=cos(float(n)*pi*x(i)*ooaspect)
            enddo
         enddo
      endif

c    -time integration

      do istep=1,nstep

c       -compute time derivatives
         do k=2,nz-1
            do n=0,nn

c             -linear terms
               dtemdt1(k,n)=oodz2*tem(k-1,n)-
     $            (2.*oodz2+float(n*n)*c*c)*tem(k,n)+
     $            oodz2*tem(k+1,n)
               domgdt1(k,n)=(oodz2*omg(k-1,n)-
     $            (2.*oodz2+float(n*n)*c*c)*omg(k,n)+
     $            oodz2*omg(k+1,n)+float(n)*c*ra*tem(k,n))*pr



c             -nonlinear terms
               do n1=0,nn
c                -for n2+n1=n
                  n2=n-n1
                  if((n2.ge.1) .and. (n2.le.nn)) then
                     dtemdt1(k,n)=dtemdt1(k,n)-
     $                  c1*(-n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                  tem(k,n1)+
     $                  n2*psi(k,n2)*(tem(k+1,n1)-tem(k-1,n1)))
                     domgdt1(k,n)=domgdt1(k,n)-
     $                  c1*(-n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                  omg(k,n1)+n2*psi(k,n2)*
     $                  (omg(k+1,n1)-omg(k-1,n1)))
                  endif
c                -for n2-n1=n
                  n2=n+n1
                  if((n2.ge.1) .and. (n2.le.nn)) then
                     dtemdt1(k,n)=dtemdt1(k,n)-
     $                  c1*(n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                  tem(k,n1)+
     $                  n2*psi(k,n2)*(tem(k+1,n1)-tem(k-1,n1)))
                     if(n.ne.0) then
                        domgdt1(k,n)=domgdt1(k,n)+
     $                     c1*(n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                     omg(k,n1)+n2*psi(k,n2)*
     $                     (omg(k+1,n1)-omg(k-1,n1)))
                     endif
                  endif
c                -for n1-n2=n
                  n2=n1-n 
                  if((n.ne.0) .and.
     $               (n2.ge.1) .and. (n2.le.nn)) then
                     dtemdt1(k,n)=dtemdt1(k,n)-
     $                  c1*(n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                  tem(k,n1)+
     $                  n2*psi(k,n2)*(tem(k+1,n1)-tem(k-1,n1)))
                     domgdt1(k,n)=domgdt1(k,n)-
     $                  c1*(n1*(psi(k+1,n2)-psi(k-1,n2))*
     $                  omg(k,n1)+n2*psi(k,n2)*
     $                  (omg(k+1,n1)-omg(k-1,n1)))
                  endif
               enddo

            enddo
         enddo

c       -update solution
         do k=2,nz-1
            tem(k,0)=tem(k,0)+
     $         0.5*dt*(3.*dtemdt1(k,0)-dtemdt2(k,0))
            do n=1,nn
               tem(k,n)=tem(k,n)+
     $           0.5*dt*(3.*dtemdt1(k,n)-dtemdt2(k,n))
               omg(k,n)=omg(k,n)+
     $           0.5*dt*(3.*domgdt1(k,n)-domgdt2(k,n))
            enddo
         enddo
         do n=1,nn
            c3=float(n*n)*c*c
            dia(1)=1.
            do k=2,nz-1
               dia(k)=2.*oodz2+c3
            enddo
            dia(nz)=1.
            call tridisl(nz,omg(1,n),psi(1,n),
     $         sub,dia,sup,work1,work2)
         enddo


c       -check values at mid-depth
         if(mod(istep,iprnt) .eq. 0) then
            k=(nz+1)/2
            write(6,*) 'step =',istep
            do n=1,nn
               write(6,2) n,tem(k,n),omg(k,n),psi(k,n)
    2          format(i3,3(2x,1pe10.3))
            enddo
         endif

c       -check max values
c         if(mod(istep,iprnt) .eq. 0) then
c            tmax=0.
c            omgmax=0.
c            psimax=0.
c            do n=1,nn
c               do k=1,nz
c                  tmax=max(tmax,abs(tem(k,n)))
c                  omgmax=max(omgmax,abs(omg(k,n)))
c                  psimax=max(psimax,abs(psi(k,n)))
c               enddo
c            enddo
c            write(6,2) istep,tmax,omgmax,psimax
c    2       format(i8,3(2x,1pe10.3))
c         endif

c       -movie data
         if((ifmov .eq. 'y') .and.
     $      (imovie.gt.0) .and. (mod(istep,imovie).eq.0)) then
            do i=1,nx
               do k=1,nz
                  temmov(i,k)=0.
                  psimov(i,k)=0.
                  do n=0,nn
                     temmov(i,k)=temmov(i,k)+tem(k,n)*cosa(n,i)
                     psimov(i,k)=psimov(i,k)+psi(k,n)*sina(n,i)
                  enddo
               enddo
            enddo
            write(movstr,"(i6)") istep/imovie
            do ii=1,6
               if(movstr(ii:ii) .eq. ' ') movstr(ii:ii)='0'
            enddo
            movfile=movfileid//movstr(1:6)
            open(2,file=movfile,form='formatted')
            write(2,22) nx4,nz4,aspect,ra,pr
   22       format(2i7,3(1x,1pe10.3))
            do k=1,nz
               write(2,23) (temmov(i,k),i=1,nx)
   23          format(7(1x,1pe10.3))
            enddo
            do k=1,nz
               write(2,23) (psimov(i,k),i=1,nx)
            enddo
            close(2)
         endif

c       -make previous time derivatives for next step
         do k=2,nz-1
            do n=1,nn
               dtemdt2(k,n)=dtemdt1(k,n)
               domgdt2(k,n)=domgdt1(k,n)
            enddo
         enddo 

      enddo

c    -store solution
      open(10,file='out',status='new',form='unformatted')
      write(10) tem,omg,psi
      write(10) dtemdt2,domgdt2

      stop
      end



      subroutine tridisl(n,rhs,sol,sub,dia,sup,a,g)

c     tridiagonal matrix solver

c     dia(1) sup(1)
c     sub(2) dia(2) sup(2)
c            sub(3) dia(3) sup(3)
c                   ...
c                          sub(n-2) dia(n-2) sup(n-2)
c                                   sub(n-1) dia(n-1) sup(n-1)
c                                            sub(n)   dia(n)

      dimension sub(n),dia(n),sup(n),a(n),g(n),rhs(n),sol(n)

c *** lu decomposition

      a(1)=1./dia(1)
      g(1)=sup(1)*a(1)
      do k=2,n-1
         a(k)=1./(dia(k)-sub(k)*g(k-1))
         g(k)=sup(k)*a(k)
      enddo
      a(n)=1./(dia(n)-sub(n)*g(n-1))

c *** solve

      sol(1)=rhs(1)*a(1)
      do k=2,n
         sol(k)=(rhs(k)-sub(k)*sol(k-1))*a(k)
      enddo
      do k=n-1,1,-1
         sol(k)=sol(k)-g(k)*sol(k+1)
      enddo

      return
      end
