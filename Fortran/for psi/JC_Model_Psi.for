      INTEGER N,m
      PARAMETER(N=4,m=100000) ! N = Number od differential equation, m = Number of divisions between x1 and x2
      INTEGER:: i
!	  complex::al,be

	  REAL h,x,y(N),dydx(N),yout(N)
      EXTERNAL derivs
      open(10,file='e:\1\JC_Model_Psi.txt')  ! Output program results
      ! -------------- initial conditons --------------
      x1=0.0                      ! x1 = t1
	  x2=30                       ! x2 = t2
      y(1)=1/sqrt(2.)             ! real part of c1(t):= y(1) + iy(2)
      y(2)=0                      ! complex part of c1(t):= y(1) + iy(2)
      y(3)=1/sqrt(2.)             ! real part of c2(t):= y(3) + iy(4)
      y(4)=0                      ! complex part of c2(t):= y(3) + iy(4)

	   x=x1
	   h=(x2-x1)/m                        ! Number of divisions
       do i=0,m-1
	      x=x1+i*h
          call derivs(x,y,dydx)
          call rk4(y,dydx,N,x,h,yout,derivs)
              x=x1+i*h
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
          if (mod(i,10)==0) then
              write(10,'(22(f10.4,2x))')x,y(1),y(2),y(3),y(4),tr_r_2,S,G   ! Computational values to create in the file a.txt
          end if
	      do j=1,n
	        y(j)=yout(j) !/sqrt((r11+r22))
	     end do
      end do

      END


      SUBROUTINE derivs(x,y,dydx)
      
      REAL x,y(N),dydx(N)
      
      ! ------------------- State --------------------
    
      ! psi(t) = c1(t)|e,n> + c2(t)|g,n+1>
      
      ! --------------- initial values --------------- 
      

	  wa=1      ! atomic frequency  
	  wp=0.6      ! field frequency
	  g=0.6    ! coupling strength 
      n = 1     ! number of photons
      
      ! ----------------- Equations ------------------    
   
      
      A = n*wp+(wa/2)
      B = (1+n)*wp-(wa/2)
      X = n+1
      C = SQRT(X)
      ! print*, A
      
      ! ----------- Differential equations -----------
      

      dydx(1)=  A*y(2)+g*C*y(4)
      dydx(2)= -A*y(1)-g*C*y(3)
      dydx(3)=  B*y(4)+g*C*y(2)
      dydx(4)= -B*y(3)-g*C*y(1)


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
