module integration
implicit none
contains

subroutine intg(n,n1,h,kx,ky,Nxi,Nyi,Nxf,Nyf,Ntime,step,step1,r0,rf,t,Du,Dv,theta_dot,r_dot,w,u,v)

 real(kind=8), dimension(-600:600,-600:600),intent(inout) :: u, v
 real(kind=8),dimension(-600:600,-600:600),intent(inout) :: w
 real(kind=8), intent(in) :: r_dot
 real(kind=8), intent(in) :: theta_dot
 real(kind=8),intent(in) :: h, kx, ky, Du, Dv
 integer, intent(in) :: Nxi, Nyi, Nxf, Nyf
 integer, intent(in) :: t, Ntime
 integer, intent(in) :: step, step1
 real(kind=8), intent(in) :: r0, rf
 
 real(kind=8), dimension(-600:600,-600:600) :: u1, v1
 real(kind=8), dimension(10) :: xk1, xk2, xk3 ,xk4
 real(kind=8) :: u2,v2
 integer :: i,j
 integer, dimension(3) :: Mx, My
 integer, intent(inout) :: n, n1

   


do i = Nxi+1, Nxf-1
     
      do j = Nyi+1, Nyf-1
      

      Mx(1) = i + 1
      Mx(2) = i
      Mx(3) = i - 1
      
      My(1) = j + 1
      My(2) = j
      My(3) = j - 1
      
!     periodic boundary conditions      
!      if(Mx(1) .eq. Nxf) Mx(1) = Nxi+1
!      if(Mx(3) .eq. Nxi) Mx(3) = Nxf-1
!      if(My(1) .eq. Nyf) My(1) = Nyi+1
!      if(My(3) .eq. Nyi) My(3) = Nyf-1
!     neumann boundaty conditions 
      
      if(Mx(1) .eq. Nxf) Mx(1) = Nxf - 1
      if(Mx(3) .eq. Nxi) Mx(3) = Nxi + 1
      if(My(1) .eq. Nyf) My(1) = Nyf - 1
      if(My(3) .eq. Nyi) My(3) = Nyi + 1
       
  
       
      u1(i,j) = u(Mx(2),My(2)) + (Du * h / (kx*kx)) * (u(Mx(1),My(2)) - 2.0d0*u(Mx(2),My(2)) + u(Mx(3),My(2)))  &
     &                         + (Du * h / (ky*ky)) * (u(Mx(2),My(1)) - 2.0d0*u(Mx(2),My(2)) + u(Mx(2),My(3)))
      
      v1(i,j) = v(Mx(2),My(2)) + (Dv * h / (kx*kx)) * (v(Mx(1),My(2)) - 2.0d0*v(Mx(2),My(2)) + v(Mx(3),My(2)))  &
     &                         + (Dv * h / (ky*ky)) * (v(Mx(2),My(1)) - 2.0d0*v(Mx(2),My(2)) + v(Mx(2),My(3)))
      
      
 ! print the results    
    
      
      end do
    
    end do
    
    
 ! Runge-kutta
 
!------------------------------------------------------
! Counting used in the growth steps    


    if((t.eq.1) .and. (mod(t,step).eq.1)) n = 0
  
    if((t.ne.1) .and. (mod(t,step).eq.1)) n = n + 1
    
    if((t.eq.1) .and. (mod(t,step1).eq.1)) n1 = 0

    if((t.ne.1) .and. (mod(t,step1).eq.1)) n1 = n1 + 1
!------------------------------------------------------
  
  do i = Nxi+1, Nxf-1
     
      do j = Nyi+1, Nyf-1
    
  !----------------------------------------------------------   
    ! calculating w(i,j)
    if(t .le. int(rf/(r_dot*h)))then  
    
      if(mod(t,step1).eq.1)then
          call light(t,i,j,n,n1,Nxi,Nyi,Nxf,Nyf,r0,rf,step,step1,kx,ky,theta_dot,r_dot,w)
      endif
    
  endif
  !----------------------------------------------------------
    
           u2 = u1(i,j)
           v2 = v1(i,j)
   
           xk1(1)=h*f(i,j,w,u2,v2)
           xk1(2)=h*g(i,j,w,u2,v2)
           
           xk2(1)=h*f(i,j,w,u2+xk1(1)/2,v2+xk1(2)/2)
           xk2(2)=h*g(i,j,w,u2+xk1(1)/2,v2+xk1(2)/2)
           
           xk3(1)=h*f(i,j,w,u2+xk2(1)/2,v2+xk2(2)/2)
           xk3(2)=h*g(i,j,w,u2+xk2(1)/2,v2+xk2(2)/2)
           
           xk4(1)=h*f(i,j,w,u2+xk3(1),v2+xk3(2))
           xk4(2)=h*g(i,j,w,u2+xk3(1),v2+xk3(2))
           
	   u(i,j)=u2+xk1(1)/6+xk2(1)/3+xk3(1)/3+xk4(1)/6
	   v(i,j)=v2+xk1(2)/6+xk2(2)/3+xk3(2)/3+xk4(2)/6
    
    
     enddo
    enddo
      
end subroutine intg 

!-------------------------------------------------------------------

subroutine light(t,i,j,n,n1,Nxi,Nyi,Nxf,Nyf,r0,rf,step,step1,kx,ky,theta_dot,r_dot,w)
real(kind=8), dimension(-600:600,-600:600), intent(inout) :: w
real(kind=8), intent(in) :: r_dot
real(kind=8), intent(in) :: theta_dot
integer, intent(in) :: Nxi,Nyi,Nxf,Nyf
integer, intent(in) :: i, j
integer, intent(in) :: t
integer, intent(in) :: n, n1
integer, intent(in) :: step, step1
real(kind=8), intent(in) :: kx, ky
real(kind=8), intent(in) :: r0, rf



real, parameter :: Pi=3.14159265359
real(kind=8) :: vx0, vy0, vx, vy
real(kind=8) :: theta, theta0
real(kind=8) :: phi2, phi3_p, phi3_m, phi4, phi0_p, phi0_m
real(kind=8) :: nm1, nm2
real(kind=8) :: prod
real(kind=8) :: r
real(kind=8) :: z

! reference vector
vx0 = r0 
vy0 = 0

nm1 = sqrt(vx0*vx0 + vy0*vy0)

 if(i.eq.0 .and. j.eq.0)goto 10

  r = (r0 + r_dot*n) 
 
  if(r.gt.rf) r = rf

  z = (theta_dot)*1.0*(1.0d0 + real(n1))! It is multiplied by 10, because the growth occurs after a time interval equal to 10 u.t. 
  
  theta = mod(z,360.)
  
 
  
 if(i.eq.Nxi+1 .and. j.eq.Nyi+1)print*,theta,z
  
  if(theta.eq.0)theta=360
   
  theta0 = theta - (theta_dot)
  
  phi3_p = theta
  phi3_m = theta
  
  phi0_p = theta0
  phi0_m = theta0
  
  
  if(theta0 .lt. 0) then
  
  


  phi3_p = theta
  phi0_p = 0.0
  
  phi3_m = 360.0
  phi0_m = 360.0 + theta0
  
  
  endif

 if(i.eq.Nxi+1 .and. j.eq.Nyi+1)print*,t,phi3_p, phi3_m, phi0_p,phi0_m
 
  vx = real(i) * kx
  vy = real(j) * ky
  
  
  nm2 = sqrt(vx*vx + vy*vy)
  
  prod = vx0*vx+vy0*vy
  
  phi2 = acos(prod/(nm1*nm2))
  
  
  if(t.ne.1)then
  
   if(w(i,j).eq.0)goto 10
  
  endif

! First filter    
  if(vx*vx+vy*vy .le. r*r)then

  
! Second filter   
  if(vy .ge. 0)then
   
    phi4 = phi2*180/Pi   
    
    if((phi4 .lt. phi3_p) .and. (phi4 .ge. phi0_p))then
    
     w(i,j) = 0   
    else     
     w(i,j) = 1.5
    endif
  
  
  elseif(vy .lt. 0)then

   phi4 = 360 - phi2*180/Pi
      
   if(phi4 .le. phi3_m .and. phi4 .ge. phi0_m)then 
     w(i,j) = 0
   else
     w(i,j) = 1.5
   endif
  
 endif
 
  else
  
   
   w(i,j) = 1.5
 
 endif
 
10 continue 

if(i.eq.0 .and. j.eq.0)then
  w(i,j) = 0
endif





end subroutine light

!-------------------------------------------------------------------

      real function f(i,j,w,u2,v2)
      implicit none
      real(kind=8) :: u2,v2
      real(kind=8) :: a
      integer, intent(in) :: i,j
      real(kind=8), dimension(-600:600,-600:600), intent(in) :: w
               a=12.d0
                             
               f=a - u2 - (4*u2*v2)/(1+u2*u2) - w(i,j) 

               
      return
      END

      real function g(i,j,w,u2,v2)
      implicit none
      real(kind=8) :: u2,v2
      real(kind=8) :: b,s
      integer, intent(in) :: i,j
      real(kind=8), dimension(-600:600,-600:600), intent(in) :: w
               b=0.30
               s=50.0d0

               g=s*(b*(u2 - (u2*v2)/(1+u2*u2) + w(i,j)))

             
      return
      END
   end module integration
!##############################################################
!##############################################################
!##############################################################
program explicit
use integration
implicit none

 real(kind=8), dimension(-600:600,-600:600) :: u, v
 real(kind=8),dimension(-600:600,-600:600) :: w
 real(kind=8) :: h, kx, ky 
 real(kind=8) :: time, xmin, xmax, ymin, ymax
 real(kind=8) :: Du, Dv
 real(kind=8) :: r0, rf
 integer :: step, step1
 integer :: i, j, t
 integer :: Ntime, Nxi, Nyi, Nxf, Nyf 
 
 real(kind=8) :: theta_dot
 real(kind=8) :: r_dot
 
 integer :: cont, n, n1
 
 u = 0.0d0
 v = 0.0d0
 
 w = 1.5

! initial radius

 r0 = 5
 

 
! linear velocity (radius/time) 
 r_dot = 0.0268519 * 1.0d0 
 
! angular velocity (degree/time)
 theta_dot = (360.0 * r_dot)/(7.25d0 * 2.6d0) 
 
 ! Defining space
 
 xmin = -40.0d0 ! here is the continuous value of xmin  
 xmax =  40.0d0 ! here is the continuous value of xmax
 
 ymin = -40.0d0 ! here is the continuous value of ymin  
 ymax =  40.0d0 ! here is the continuous value of ymax
 
 ! Defining the size of the space step
 
 kx = 0.5d0
 
 ky = 0.5d0 
  
 ! Defining the discret number of steps
   
 Nxi = int(xmin/kx)
 Nyi = int(ymin/ky)
 
 Nxf = int(xmax/kx)
 Nyf = int(ymax/ky)
 
 print*,"Nxi = ", Nxi, "Nxf = ", Nxf
 
 ! final radius
 rf = xmax - 5.0
 
 ! time needed to achieve rf and to complete one aditional rotation

 ! Defining time
 
 time = rf/r_dot + 10.


 ! Defining the size of the time step
 
 h = 1.0d-3
   
 print*,"final radius = ",rf
   
 print*,"time requested to radius reach the final value = ", rf/r_dot 

 Ntime = time/h
 
 step = int(1./h)
 
 step1 = int(1./h)
 
 ! Defining the diffusion coefficient
 
 Du = 1.0d0
 Dv = 50.0d0*1.00d0
 
 ! initial condition for chemicals 
 
 do i = Nxi, Nxf-1
  do j = Nyi, Nyf-1
  
     u(i,j) = 2.0d0 + 0.15d0*rand()
     v(i,j) = 6.76d0 + 0.15d0*rand()
         
  enddo
 enddo 


 ! testing Neumann condition 
 print*,"h = ",h,"kx = ky =", kx
 print*,"total time = ", time, Ntime 
 print*,"D Dt / Dx Dx = ", (Dv*h)/(kx*kx)
 
 cont = 9
 
 do t = 1, Ntime
 
    call intg(n,n1,h,kx,ky,Nxi,Nyi,Nxf,Nyf,Ntime,step,step1,r0,rf,t,Du,Dv,theta_dot,r_dot,w,u,v)
    
    if(t.eq.Ntime)then
    
      cont = cont + 1
     
     open(unit = cont)
     
     do i = Nxi+1, Nxf - 1
      do j = Nyi+1, Nyf - 1         
        
       write(cont,*)i*kx,j*ky,u(i,j) 
        
      end do
     end do
     
     close(unit = cont)
    endif
 
   
 
  end do 

 end program explicit   
