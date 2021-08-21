      INTEGER N,m
      PARAMETER(N=4,m=20000)     ! N = Number od differential equation, m = Number of divisions between x1 and x2
      INTEGER:: i
      ! complex::al,be

      REAL h,x,y(N),dydx(N),yout(N)
      EXTERNAL derivs
	  open(10,file='e:\1\Rabi_Model_schrodinger_sigma(x).txt') ! Output program results --- for in sigma basis is easier
	  !open(20,file='e:\1\c.txt')

      x1=0.0           ! x1 = t1
	  x2=30            ! x2 = t2
      ! initial conditons
      y(1)=1/sqrt(2.)            ! real part of c1(t):= y(1) + iy(2)
      y(2)=0                      ! complex part of c1(t):= y(1) + iy(2)
      y(3)=1/sqrt(2.)             ! real part of c2(t):= y(3) + iy(4)
      y(4)=0                      ! complex part of c2(t):= y(3) + iy(4)

      ! calculate in sigma(x) base
	  x=x1
	  h=(x2-x1)/m      ! Number of divisions 
      do i=0,m-1
	       x=x1+i*h
           call derivs(x,y,dydx)
           call rk4(y,dydx,N,x,h,yout,derivs)
           x=x1+i*h
           if (mod(i,1)==0) then
	          A=y(1)*y(1)+y(2)*y(2)       ! |c1(t)|^2
              B=y(3)*y(3)+y(4)*y(4)       ! |c2(t)|^2
              rho_a_2_11 = A * A
              rho_a_2_22 = B * B
              tr_r_2  = rho_a_2_11 + rho_a_2_22 ! trace (rho_a) ^ 2 
              U = A * A
              O = B * B
              S= -((B*LOG(B)/LOG(2.))+(A*LOG(A)/LOG(2.)))   ! Von Neumann entropy
              G= 2*(1-U-O)                ! G = L (Linear entropy)
              C = SQRT(G) !Concurence
!	         et=y(1)*2*(y(3)+y(2)*y(4))
      	     write(10,'(22(f10.4,2x))')x,y(1),y(2),y(3),y(4),tr_r_2,S,G  ! Computational values to create in the file a.txt
           end if
	     do j=1,n
	        y(j)=yout(j) !/sqrt((r11+r22))
	     end do
	  end do

      END

      SUBROUTINE derivs(x,y,dydx)
      REAL x,y(*),dydx(*)
	  REAL la
      
      ! ------------------- State --------------------
      
      ! psi(t) = c1(t)|+x,n> + c2(t)|-x,n+1>

      ! --------------- initial values --------------- 
      
	  wa=1 ! atomic frequency  
	  wp=0.6 ! field frequency 
	  g=0.006 ! coupling strength; g =< (0.5)wp
      n=1  ! number of photons
      
      ! ----------------- Equations ------------------
      
      la = g / (wp + ( wa * exp((-2) * ( g / ((wp+wa) * (wp+wa)))))) ! la = lambda
      
      x = 4 * la * la
      
      Rr= 2 * (la * wp + g)    
      
      !Rr=la*wp+g-wa*la*exp(-2*la*la)/(n+1)*aalgu(n,1,x)    ! Rr= R[ar]
      
      G0n = exp(-2 * la * la) * algu(n,x) 
      
      alpha = wp * la * la + 2 * la * g
      
      beta =0.5 * wa * G0n
      
      gamma = Rr
      
      A = wp * n + alpha + beta
      
      B = gamma * sqrt(n+1.)
      
      C = (wp * (n + 1)) + alpha - beta !!!!!! I changed "n" to "n+1" in this file !!!!!!
      
      ! ----------- Differential equations -----------

      dydx(1)= A*y(2)+B*y(4)
      dydx(2)=-A*y(1)-B*y(3)
      dydx(3)= C*y(4)+B*y(2)
      dydx(4)=-C*y(3)-B*y(1)
      

      return
     
	  end

      ! ------------ Runge-Kutta Function ------------
      
      SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs) 
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i
      REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5
      h6=h/6.
	  x=0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))

14    continue
      return
      END


      ! ------- Laguerre polynomial Function -------

       function algu (n,x)
       al0=1
       al1=1-x
       if (n==0) then
          algu=al0 
	 return
       else if (n==1) then
	     algu=al1
	 return
       end if
       alk=al1
       alkm=al0
       do k=2,n
          alk1=((2*k-1.-x)*alk-(k-1)*alkm)/k
          alkm=alk
          alk=alk1
       end do
       algu=alk1
       return
       end
       
      ! -- Associate Laguerre polynomial Function --
      
       function aalgu (n,k,x)
       al0=1
       al1=1+k-x
       if (n==0) then
           aalgu=al0 
	 return
       else if (n==1) then
	      aalgu=al1
	 return
       end if
       alk=al1
       alkm=al0
       do i=2,n
           alk1=(k+2*i-1.-x)/i*alk-(k-1+i)/i*alkm
           alkm=alk
           alk=alk1
       end do
       aalgu=alk1
       return
       end

