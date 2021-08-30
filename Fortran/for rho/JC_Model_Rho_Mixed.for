      INTEGER N,m
      PARAMETER(N=32,m=100000) ! N = Number od differential equation, m = Number of divisions between x1 and x2
      INTEGER:: i
!	  complex::al,be

	  REAL h,x,y(N),dydx(N),yout(N)
      EXTERNAL derivs
      open(10,file='e:\1\JC_Model_Rho_Mixed.txt')  ! Output program results
      ! -------------- initial conditons --------------
      x1=0.0                      ! x1 = t1
	  x2=10                       ! x2 = t2
      y(1)=  0.25
      y(2)=  0
      y(3)=  0
      y(4)=  0
      y(5)=  0
      y(6)=  0
      y(7)=  0
      y(8)=  0 
      y(9)=  0
      y(10)= 0
      y(11)= 0.25
      y(12)= 0
      y(13)= 0
      y(14)= 0
      y(15)= 0
      y(16)= 0 
      y(17)= 0
      y(18)= 0
      y(19)= 0
      y(20)= 0
      y(21)= 0.25
      y(22)= 0
      y(23)= 0
      y(24)= 0
      y(25)= 0
      y(26)= 0
      y(27)= 0
      y(28)= 0
      y(29)= 0
      y(30)= 0
      y(31)= 0.25
      y(32)= 0


	   x=x1
	   h=(x2-x1)/m                        ! Number of divisions
       do i=0,m-1
	      x=x1+i*h
          call derivs(x,y,dydx)
          call rk4(y,dydx,N,x,h,yout,derivs)
              x=x1+i*h
              trace_rho = y(1)+y(2)+y(11)+y(12)+y(21)+y(22)+y(31)+y(32)
          if (mod(i,10)==0) then
              write(10,'(48(f10.4,2x))')x,y(1),y(2),y(3),y(4),y(5),y(6)   ! (#(f10.4,2x)) == # -> number of outputs ; f10 -> serial number for output file
     &         ,y(7),y(8),y(9),y(10),y(11),y(12),y(13),y(14),y(15),y(16)
     &         ,y(17),y(18),y(19),y(20),y(21),y(22),y(23),y(24),y(25)
     &         ,y(26),y(27),y(28),y(29),y(30),y(31),y(32),trace_rho ! Computational values to create in the file a.txt
          end if
	      do j=1,n
	        y(j)=yout(j) !/sqrt((r11+r22))
	     end do
      end do

      END


      SUBROUTINE derivs(x,y,dydx)
      
      REAL x,y(*),dydx(*)
      

      
      ! --------------- initial values --------------- 
      

	  wa=1      ! atomic frequency  
	  wp=0.6      ! field frequency
	  g=0.8   ! coupling strength 
      n_p = 20     ! number of photons
      
      ! ----------------- Equations ------------------    
   
      
      Xm = n_p+1
      C = g*SQRT(Xm)

      
      
      ! ----------- Differential equations -----------
      

      dydx(1)=   0
      dydx(2)=   0
      dydx(3)=  -C*y(6)-wp*y(4)
      dydx(4)=  +C*y(5)+wp*y(3)
      dydx(5)=  -C*y(4)-wa*y(6)
      dydx(6)=  +C*y(3)+wa*y(5)
      dydx(7)=  -(wa+wp)*y(8)
      dydx(8)=  +(wa+wp)*y(7)
      dydx(9)=  +C*y(18)+wp*y(10)
      dydx(10)= -C*y(17)-wp*y(9)
      dydx(11)= +C*(-y(14)+y(20))
      dydx(12)= +C*(y(13)-y(19))
      dydx(13)= +C*(y(22)-y(12))-(wa-wp)*y(14)
      dydx(14)= -C*(y(21)-y(11))+(wa-wp)*y(13)
      dydx(15)= +C*y(24)-wa*y(16)
      dydx(16)= -C*y(23)+wa*y(15)
      dydx(17)= +C*y(10)+wa*y(18)
      dydx(18)= -C*y(9)-wa*y(17)
      dydx(19)= +C*(y(12)-y(22))-(-wa+wp)*y(20)
      dydx(20)= -C*(y(11)-y(21))+(-wa+wp)*y(19)
      dydx(21)= +C*(y(14)-y(20))
      dydx(22)= -C*(y(13)-y(19))
      dydx(23)= +C*y(16)-wp*y(24)
      dydx(24)= -C*y(15)+wp*y(23)
      dydx(25)= +(wa+wp)*y(26)
      dydx(26)= -(wa+wp)*y(25)
      dydx(27)= -C*y(30)+wa*y(28)
      dydx(28)= +C*y(29)-wa*y(27)
      dydx(29)= -C*y(28)+wp*y(30)
      dydx(30)= +C*y(27)-wp*y(29)
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
