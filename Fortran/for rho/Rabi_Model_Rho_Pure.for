      INTEGER N,m
      PARAMETER(N=32,m=20000)     ! N = Number od differential equation, m = Number of divisions between x1 and x2
      INTEGER:: i
      ! complex::al,be

      REAL h,x,y(N),dydx(N),yout(N)
      EXTERNAL derivs
	  open(10,file='e:\1\Rabi_Model_Rho_pure.txt') ! Output program results 
	  ! open(20,file='e:\1\c.txt')

      x1=0.0           ! x1 = t1
	  x2=30            ! x2 = t2
      ! initial conditons
      y(1)=  0
      y(2)=  0
      y(3)=  0
      y(4)=  0
      y(5)=  0
      y(6)=  0
      y(7)=  0
      y(8)=  0 
      y(9)=  0
      y(10)= 0
      y(11)= 0.5
      y(12)= 0
      y(13)= 0.5
      y(14)= 0
      y(15)= 0
      y(16)= 0 
      y(17)= 0
      y(18)= 0
      y(19)= 0.5
      y(20)= 0
      y(21)= 0.5
      y(22)= 0
      y(23)= 0
      y(24)= 0
      y(25)= 0
      y(26)= 0
      y(27)= 0
      y(28)= 0
      y(29)= 0
      y(30)= 0
      y(31)= 0
      y(32)= 0

      
	  x=x1
	  h=(x2-x1)/m      ! Number of divisions 
      do i=0,m-1
	       x=x1+i*h
           call derivs(x,y,dydx)
           call rk4(y,dydx,N,x,h,yout,derivs)
           x=x1+i*h
           if (mod(i,1)==0) then
              S=-((y(11)*LOG(y(11))/LOG(2.))+(y(21)*LOG(y(21))/LOG(2.)))     ! Von Neumann entropy
              G= 2*(1-y(11)*y(11)-y(21)*y(21))                ! G = L (Linear entropy)
              Concurence = SQRT(G)
              trace_rho = y(1)+y(2)+y(11)+y(12)+y(21)+y(22)+y(31)+y(32)
!	         et=y(1)*2*(y(3)+y(2)*y(4))
      	     write(10,'(45(f20.4,2x))')x,y(1)+y(2),y(3)+y(4),y(5)+y(6)   ! (#(f10.4,2x)) == # -> number of outputs ; f10 -> serial number for output file
     &         ,y(7)+y(8),y(9)+y(10),y(11)+y(12),y(13)+y(14),y(15)+y(16)
     &         ,y(17)+y(18),y(19)+y(20),y(21)+y(22),y(23)+y(24),y(25)
     &         +y(26),y(27)+y(28),y(29)+y(30),y(31)+y(32),S,G,Concurence
     &          ,trace_rho ! Computational values to create in the file a.txt
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
      n_p=1  ! number of photons
      
      ! ----------------- Equations ------------------
      
      la = g / (wp + ( wa * exp((-2.) * ( g / ((wp+wa) * (wp+wa)))))) ! la = lambda
      
      xn = 4 * la * la
      
      R= 2 * (la * wp + g)    
      
      !Rr=la*wp+g-wa*la*exp(-2*la*la)/(n+1)*aalgu(n,1,x)    ! Rr= R[ar]
      
      G0n = exp(-2. * la * la) * algu(n_p,xn) 
      
      beta = (wp * la * la) + (2 * la * g)
      
      A = 0.5 * wa * G0n ! = alpha
      
      Xm = n_p+1

      C = R * sqrt(Xm)
      ! print*, C
      
      
      ! ----------- Differential equations -----------

      dydx(1)=   0
      dydx(2)=   0
      dydx(3)=  -(wp*y(4))-(C*y(6))
      
      dydx(4)=  +(wp*y(3))+(C*y(5))
      
      dydx(5)=  -2*A*y(6)-C*y(4)
      
      dydx(6)=  +(2*A*y(5))+C*y(3)
      
      dydx(7)=  -((2*A)+wp)*y(8)
      
      dydx(8)=  +((2*A)+wp)*y(7)
      
      dydx(9)=  +wp*y(10)+C*y(18)
      
      dydx(10)= -wp*y(9)-C*y(17)
      
      dydx(11)= -C*(y(14)-y(20))
      
      dydx(12)= +C*(y(13)-y(19))
      
      dydx(13)= -(2*A-wp)*y(14)-C*(y(12)-y(22))
      
      dydx(14)= +(2*A-wp)*y(13)+C*(y(11)-y(21))
      
      dydx(15)= -2*A*y(16)+C*y(24)
      
      dydx(16)= +2*A*y(15)-C*y(23)
      
      dydx(17)= +2*A*y(18)+C*y(10)
      
      dydx(18)= -2*A*y(17)-C*y(9)
      
      dydx(19)= -(wp-2*A)*y(20)+C*(y(12)-y(22))
      
      dydx(20)= +(wp-2*A)*y(19)-C*(y(11)-y(21))
      
      dydx(21)= +C*(y(14)-y(20))
      
      dydx(22)= -C*(y(13)-y(19))
      
      dydx(23)= -wp*y(24)+C*y(16)
      
      dydx(24)= +wp*y(23)-C*y(15)
      
      dydx(25)= +(wp+2*A)*y(26)
      
      dydx(26)= -(wp+2*A)*y(25)
      
      dydx(27)= +2*A*y(28)-C*y(30)
      
      dydx(28)= -2*A*y(27)+C*y(29)
      
      dydx(29)= +wp*y(30)-C*y(28)
      
      dydx(30)= -wp*y(29)+C*y(27)
      
      dydx(31)=  0
      dydx(32)=  0
      

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

