module global
implicit none
contains


 subroutine in_cond(N_cells_x,N_cells_y,u,v)
 real(kind=8), allocatable, intent(inout) :: u(:,:), v(:,:)
 integer, intent(in) :: N_cells_x, N_cells_y
 
 integer :: i, j
 
  do i = 1, N_cells_x
    do j = 1, N_cells_y
    
     print*,rand()
     
    enddo
 enddo

 
 do i = 1, N_cells_x
    do j = 1, N_cells_y
    
     

     u(i,j) = 0.0622535d0 * (1.0d0 + 0.01d0 * rand()*(-1.0d0)**int(2*rand()))
     v(i,j) = 0.0622535d0 * (1.0d0 + 0.01d0 * rand()*(-1.0d0)**int(2*rand()))
     
    enddo
 enddo


 end subroutine in_cond
 
 subroutine coupling(N_cells_x,N_cells_y,kc)
 real(kind=8), allocatable, intent(inout) :: kc(:,:)
 integer, intent(in) :: N_cells_x, N_cells_y
  
 integer :: ix,jy
 integer :: i, j
 

       kc = 0.01d0        
    

 end subroutine coupling
 
 
 subroutine epsilonp(N_cells_x,N_cells_y,ff0,ff)
 real(kind = 8), allocatable, intent(inout) :: ff(:,:)
 real(kind=8), intent(in) :: ff0
 integer, intent(in) :: N_cells_x, N_cells_y
 
 integer :: ix,jy
 integer :: number_pert
 integer :: mx(2), my(2)
 integer, allocatable :: ci(:), cj(:)
 integer :: i, j, k
 
  
  do i = 1, N_cells_x
     do j = 1, N_cells_y
    
    if(i.le.1)then
    
       ff(i,j) = ff0 * (1.0d0 + 0.00d0 * rand()*(-1.0d0)**int(2*rand()))

    else
    
       ff(i,j) = 2.7d0 * (1.0d0 + 0.0d0 * rand()*(-1.0d0)**int(2*rand()))
    
    endif
    
     enddo
  enddo
     
 end subroutine epsilonp
 
 
 
 subroutine intg(h,length_cell_x,length_cell_y,N_cells_x,N_cells_y,kc,kc_cross,Du,Dv,ff,u,v)

 real(kind=8), allocatable,intent(inout) :: u(:,:), v(:,:)
 real(kind=8), allocatable,intent(in) :: kc(:,:)
 real(kind=8), allocatable, intent(in) :: ff(:,:)
 real(kind=8), intent(in) :: length_cell_x,length_cell_y
 real(kind=8), intent(in) ::  kc_cross
 real(kind=8), intent(in) :: h, Du, Dv
 integer, intent(in) :: N_cells_x, N_cells_y
 
 real(kind=8), allocatable :: u1(:,:), v1(:,:)
 real(kind=8) :: u2,v2
 real(kind=8), allocatable :: u3(:,:), v3(:,:)
 real(kind=8) :: kx, ky
 real(kind=8), dimension(10) :: xk1, xk2, xk3 ,xk4
 integer, dimension(3) :: Mx, My
 integer :: i,j
      
  allocate(u1(N_cells_x,N_cells_y))
  allocate(v1(N_cells_x,N_cells_y))
  
  allocate(u3(N_cells_x,N_cells_y))
  allocate(v3(N_cells_x,N_cells_y))
  
   kx = length_cell_x
   
   ky = length_cell_y
  
    
    do i = 1, N_cells_x
     do j = 1, N_cells_y
     
      Mx(1) = i + 1
      Mx(2) = i
      Mx(3) = i - 1
      
      My(1) = j + 1
      My(2) = j
      My(3) = j - 1
      
      
 !     neumann boundary conditions 
      if(Mx(1) .eq. N_cells_x+1) Mx(1) = N_cells_x
      if(Mx(3) .eq. 0) Mx(3) = 1
      if(My(1) .eq. N_cells_y+1) My(1) = N_cells_y
      if(My(3) .eq. 0) My(3) = 1
     
     u1(i,j) = u(Mx(2),My(2)) + (Du * h / (kx*kx)) * (u(Mx(1),My(2)) &
     &     - 2.0d0*u(Mx(2),My(2)) + u(Mx(3),My(2))) + (Du * h / (ky*ky)) * (u(Mx(2),My(1))&
     &     - 2.0d0*u(Mx(2),My(2)) + u(Mx(2),My(3)))
     
     v1(i,j) = v(Mx(2),My(2)) + (Dv * h / (kx*kx)) * (v(Mx(1),My(2)) &
     &     - 2.0d0*v(Mx(2),My(2)) + v(Mx(3),My(2))) + (Dv * h / (ky*ky)) * (v(Mx(2),My(1))&
     &     - 2.0d0*v(Mx(2),My(2)) + v(Mx(2),My(3)))
          
     
           u2 = u1(i,j)
           v2 = v1(i,j)
   
           xk1(1)=h*f(u2,v2,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y,ff)
           xk1(2)=h*g(u2,v2,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y)
           
           xk2(1)=h*f(u2+xk1(1)/2.0d0,v2+xk1(2)/2.0d0,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y,ff)
           xk2(2)=h*g(u2+xk1(1)/2.0d0,v2+xk1(2)/2.0d0,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y)
           
           xk3(1)=h*f(u2+xk2(1)/2.0d0,v2+xk2(2)/2.0d0,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y,ff)
           xk3(2)=h*g(u2+xk2(1)/2.0d0,v2+xk2(2)/2.0d0,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y)
           
           xk4(1)=h*f(u2+xk3(1),v2+xk3(2),u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y,ff)
           xk4(2)=h*g(u2+xk3(1),v2+xk3(2),u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y)
           
	   u3(i,j)=u2+xk1(1)/6.0d0+xk2(1)/3.0d0+xk3(1)/3.0d0+xk4(1)/6.0d0
	   v3(i,j)=v2+xk1(2)/6.0d0+xk2(2)/3.0d0+xk3(2)/3.0d0+xk4(2)/6.0d0
 
         if (isnan(u(i,j))) stop 'u is a NaN'
     enddo
    enddo
    
        u = u3
        
        v = v3

end subroutine intg 


      real function f(u2,v2,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y,ff)
      implicit none
      real(kind=8), allocatable ,intent(in) :: u(:,:), v(:,:)
      real(kind=8), allocatable ,intent(in) :: kc(:,:)
      real(kind=8), allocatable ,intent(in) :: ff(:,:)
      real(kind=8), intent(in) :: u2, v2
      real(kind=8), intent(in) :: kc_cross 
      integer, intent(in) :: N_cells_x,N_cells_y     
      integer, intent(in) :: i,j
      
      real(kind=8) :: q, e
      real(kind=8) :: R
     
            
            
               e = 0.01d0
               q = 0.002d0
               
               if((i.eq.1).and.(j.eq.1))then
                             
               R = kc(i,j) * u(i+1,j) + kc(i,j) * u(i,j+1) + kc_cross * u(i+1,j+1) - &
     &               u2 * (kc(i,j) + kc(i,j) + kc_cross)

               elseif((i.eq.1).and.(j.eq.N_cells_y))then
               
                R = kc(i,j) * u(i+1,j) + kc(i,j) * u(i,j-1) + kc_cross * u(i+1,j-1) - &
     &               u2 * (kc(i,j) + kc(i,j) + kc_cross )
                
                elseif((i.eq.N_cells_x).and.(j.eq.1))then
               
                R = kc(i,j) * u(i-1,j) + kc(i,j) * u(i,j+1) +kc_cross * u(i-1,j+1) - &
     &               u2 * (kc(i,j) + kc(i,j) + kc_cross)
               
                elseif((i.eq.N_cells_x).and.(j.eq.N_cells_y))then
               
                R = kc(i,j) * u(i-1,j) + kc(i,j) * u(i,j-1) + kc_cross * u(i-1,j-1) - &
     &               u2 * (kc(i,j) + kc(i,j) + kc_cross)
               
                elseif(((i.eq.1).and.(j.gt.1)).and.(j.lt.N_cells_y))then
               
                R = kc(i,j) * u(i+1,j) + kc(i,j) * (u(i,j-1)+u(i,j+1)) + & 
     &              kc_cross * (u(i+1,j-1)+u(i+1,j+1)) - &
     &              u2 * (kc(i,j) + kc(i,j) * 2.0d0 + kc_cross * 2.0d0)
                
                elseif(((i.eq.N_cells_x).and.(j.gt.1)).and.(j.lt.N_cells_y))then
               
                R = kc(i,j) * u(i-1,j) + kc(i,j) * (u(i,j-1)+u(i,j+1)) +  &
     &              kc_cross * (u(i-1,j-1)+u(i-1,j+1)) - &
     &              u2 * (kc(i,j) + kc(i,j) * 2.0d0 + kc_cross * 2.0d0 )
                
                elseif(((j.eq.1).and.(i.gt.1)).and.(i.lt.N_cells_x))then
               
                R = kc(i,j) * (u(i-1,j)+u(i+1,j)) + kc(i,j) * (u(i,j+1)) + &
     &              kc_cross * (u(i-1,j+1)+u(i+1,j+1)) - &
     &              u2 * (kc(i,j) * 2.0d0 + kc(i,j) + kc_cross * 2.0d0)
                
                elseif(((j.eq.N_cells_y).and.(i.gt.1)).and.(i.lt.N_cells_x))then
               
                R = kc(i,j) * (u(i-1,j)+u(i+1,j)) + kc(i,j) * (u(i,j-1)) + &
     &              kc_cross * (u(i-1,j-1)+u(i+1,j-1)) - &
     &              u2 * (kc(i,j) * 2.0d0 + kc(i,j) + kc_cross * 2.0d0 )
                
                else
                
                R = kc(i,j) * (u(i-1,j)+u(i+1,j)) + kc(i,j) * (u(i,j-1)+u(i,j+1)) + &
     &              kc_cross * (u(i-1,j-1)+u(i+1,j-1)+u(i-1,j+1)+u(i+1,j+1)) - &
     &              u2 * (kc(i,j) * 2.0d0 + kc(i,j) * 2.0d0 + kc_cross * 4.0d0)
     
                endif
                
                
                
                
               
               f = (u2 - u2*u2 - ff(i,j)*v2*(u2-q)/(u2+q)   + R )/e 
               
             
      return
      END

      real function g(u2,v2,u,v,i,j,kc,kc_cross,N_cells_x,N_cells_y)
      implicit none
      real(kind=8), allocatable ,intent(in) :: u(:,:), v(:,:)
      real(kind=8), allocatable ,intent(in) :: kc(:,:)
      real(kind=8), intent(in) :: u2, v2
      real(kind=8), intent(in) :: kc_cross 
      integer, intent(in) :: N_cells_x,N_cells_y     
      integer, intent(in) :: i,j
      
      real(kind=8) :: R
      
              if((i.eq.1).and.(j.eq.1))then
                             
               R = kc(i,j) * v(i+1,j) + kc(i,j) * v(i,j+1) + kc_cross * v(i+1,j+1) - &
     &               v2 * (kc(i,j) + kc(i,j) + kc_cross)

               elseif((i.eq.1).and.(j.eq.N_cells_y))then
               
                R = kc(i,j) * v(i+1,j) + kc(i,j) * v(i,j-1) + kc_cross * v(i+1,j-1) - &
     &               v2 * (kc(i,j) + kc(i,j) + kc_cross )
                
                elseif((i.eq.N_cells_x).and.(j.eq.1))then
               
                R = kc(i,j) * v(i-1,j) + kc(i,j) * v(i,j+1) + kc_cross * v(i-1,j+1) - &
     &               v2 * (kc(i,j) + kc(i,j) + kc_cross)
               
                elseif((i.eq.N_cells_x).and.(j.eq.N_cells_y))then
               
                R = kc(i,j) * v(i-1,j) + kc(i,j) * v(i,j-1) + kc_cross * v(i-1,j-1) - &
     &               v2 * (kc(i,j) + kc(i,j) + kc_cross)
               
                elseif(((i.eq.1).and.(j.gt.1)).and.(j.lt.N_cells_y))then
               
                R = kc(i,j) * v(i+1,j) + kc(i,j) * (v(i,j-1)+v(i,j+1)) + & 
     &             kc_cross * (v(i+1,j-1)+v(i+1,j+1)) - &
     &              v2 * (kc(i,j) + kc(i,j) * 2.0d0 + kc_cross * 2.0d0)
                
                elseif(((i.eq.N_cells_x).and.(j.gt.1)).and.(j.lt.N_cells_y))then
               
                R = kc(i,j) * v(i-1,j) + kc(i,j) * (v(i,j-1)+v(i,j+1)) +  &
     &              kc_cross * (v(i-1,j-1)+v(i-1,j+1)) - &
     &              v2 * (kc(i,j) + kc(i,j) * 2.0d0 + kc_cross * 2.0d0 )
                
                elseif(((j.eq.1).and.(i.gt.1)).and.(i.lt.N_cells_x))then
               
                R = kc(i,j) * (v(i-1,j)+v(i+1,j)) + kc(i,j) * (v(i,j+1)) + &
     &              kc_cross * (v(i-1,j+1)+v(i+1,j+1)) - &
     &              v2 * (kc(i,j) * 2.0d0 + kc(i,j) +kc_cross* 2.0d0)
                
                elseif(((j.eq.N_cells_y).and.(i.gt.1)).and.(i.lt.N_cells_x))then
               
                R = kc(i,j) * (v(i-1,j)+v(i+1,j)) + kc(i,j) * (v(i,j-1)) + &
     &              kc_cross * (v(i-1,j-1)+v(i+1,j-1)) - &
     &              v2 * (kc(i,j) * 2.0d0 + kc(i,j) + kc_cross * 2.0d0 )
                
                else
                
                R = kc(i,j) * (v(i-1,j)+v(i+1,j)) + kc(i,j) * (v(i,j-1)+v(i,j+1)) + &
     &              kc_cross * (v(i-1,j-1)+v(i+1,j-1)+v(i-1,j+1)+v(i+1,j+1)) - &
     &              v2 * (kc(i,j) * 2.0d0 + kc(i,j) * 2.0d0 + kc_cross * 4.0d0)
     
                endif
      
      
             
      
               g = u2 - v2 + R

             
      return
      END
      


 subroutine printing(t,numb,N_cells_x,N_cells_y,length_cell_x,length_cell_y,h,u,v)
 real(kind=8), allocatable, intent(in) :: u(:,:), v(:,:)
 real(kind=8), intent(in) :: h
 real(kind=8), intent(in) :: length_cell_x, length_cell_y
 integer, intent(in) :: N_cells_x, N_cells_y
 integer, intent(in) :: t
 integer, intent(inout) :: numb
 
 integer :: p1x, p1y
 integer :: p2x, p2y
 integer :: p3x, p3y
 integer :: p4x, p4y
 integer :: pcx, pcy
 character (len=30) :: filename1, filename2, filename3, filename4, filename5
 character (len=30) :: filename6, filename8, integ
 integer :: i, j
 
 
       filename1 = 'datap1.dat'
       filename2 = 'datap2.dat'
       filename3 = 'datap3.dat'
       filename4 = 'datap4.dat'
       filename5 = 'datac.dat'
       
       
       open(unit = 12, file = filename1, status = 'unknown')
       open(unit = 13, file = filename2, status = 'unknown')
       open(unit = 14, file = filename3, status = 'unknown')
       open(unit = 15, file = filename4, status = 'unknown')
       open(unit = 16, file = filename5, status = 'unknown')
       
       
  if(mod(t,int(1.0d-3/h)).eq.1)then  
  
   
      pcx = int(N_cells_x/2); pcy = int(N_cells_y/2)
      
      p1x = int(N_cells_x/4); p1y = pcy
      
      p2x = N_cells_x; p2y = pcy
   
      p3x = pcx; p3y = 1
      
      p4x = pcx; p4y = N_cells_y
        
          
   
     write(12,*)t * h,u(p1x,p1y) 
     write(13,*)t * h,u(p2x,p2y) 
     write(14,*)t * h,u(p3x,p3y) 
     write(15,*)t * h,u(p4x,p4y)
     write(16,*)t * h,u(pcx,pcy) 
     
  endif  
  
   if(mod(t,int(1.0d-1/h)).eq.0)then
   
       numb = numb + 1

       
       if(numb.lt.10)then
       write(integ,'(I1)')numb
       elseif(numb.ge.10 .and. numb.lt.100)then
       write(integ,'(I2)')numb
       elseif(numb.ge.100 .and. numb.lt.1000)then
       write(integ,'(I3)')numb
      else
       write(integ,'(I4)')numb
       endif
    
     filename6 = 'space/'//trim(integ)//'.dat'
     open(unit = 10, file = filename6, status = 'unknown')
     
     filename8 = 'space_2/'//trim(integ)//'.dat'
     open(unit = 11, file = filename8, status = 'unknown')
    
     do i = 1, N_cells_x 
      do j = 1, N_cells_y    
      
       write(10,*)i*length_cell_x,j*length_cell_y,u(i,j)
       
       if(j.eq.int(N_cells_y/2)) write(11,*)i*length_cell_x,u(i,j)
        
        
      end do
     end do
     
    endif
 
 
 end subroutine printing

end module global

program bz_cells
use global
implicit none

 real(kind=8), allocatable :: u(:,:), v(:,:)
 real(kind=8), allocatable :: kc(:,:)
 real(kind=8), allocatable :: ff(:,:)
 real(kind=8) :: ff0
 real(kind=8) :: h, length_cell_x, length_cell_y 
 real(kind=8) :: time, x, y
 real(kind=8) :: Du, Dv
 real(kind=8) :: kc_cross
 integer :: Ntime, N_cells_x, N_cells_y
 integer :: numb, numb1
 integer :: i, j, t
 
 
 ! Defining time and space
 
 time = 1000
 
 x = 100.0
 
 y = 10.0
  
 ! Defining the size of the time and space step
 
 h = 1.0d-4
  
 ! Supposing each cell is rectangular shaped. Its lengths are:  
  
 length_cell_x  = 0.5d0

 length_cell_y  = 0.5d0
 
 ! Number of cells 
 
 N_cells_x = int(x/length_cell_x)
 
 N_cells_y = int(y/length_cell_y)
 
 ! Number of steps
 
 Ntime = time/h
 
 allocate(u(N_cells_x,N_cells_y))
 allocate(v(N_cells_x,N_cells_y))
 allocate(kc(N_cells_x,N_cells_y))
 allocate(ff(N_cells_x,N_cells_y))
 
 ! Defining the diffusion coefficient
 
 Du = 1.0d0
 Dv = 1.0d0
 
 ff0 = 2.d0
 
 kc_cross = 0.0d0
 
 call in_cond(N_cells_x,N_cells_y,u,v)
 
 
 call coupling(N_cells_x,N_cells_y,kc)
 
 call epsilonp(N_cells_x,N_cells_y,ff0,ff)
 
 numb = 0
 numb1 = 0
do t = 1, Ntime

  print*,t
  
  call intg(h,length_cell_x,length_cell_y,N_cells_x,N_cells_y,kc,kc_cross,Du,Dv,ff,u,v)
 
 
  if(t .gt. int(time-200)/h)then
 
     call printing(t,numb,N_cells_x,N_cells_y,length_cell_x,length_cell_y,h,u,v)
   
  endif
   
enddo
 

end program bz_cells
