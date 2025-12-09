program rk4_3cuerpos
real :: t(0:30000), x(0:30000), y(0:30000), z(0:30000), vx(0:30000), vy(0:30000), vz(0:30000), g, m1, m2, m3, h, d, &
kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kz1, kz2, kz3, kz4, lx1, lx2, lx3, lx4,ly1, ly2, ly3, ly4,lz1, lz2, lz3, lz4

g = 1
h = 0.001
n = 25000 
d = 1
do i=0,n
 t(i)=i*h
enddo
! Condiciones iniciales:
write(*,*) "X(0)"
read(*,*) x(0)

write(*,*) "Y(0)"
read(*,*) y(0)

write(*,*) "Z(0)"
read(*,*) z(0)


write(*,*) "Vx(0)"
read(*,*) vx(0)

write(*,*) "Vy(0)"
read(*,*) vy(0)

write(*,*) "Vz(0)"
read(*,*) vz(0)

write(*,*) "Masa de la primera estrella (m2): "
read(*,*) m2

write(*,*) "Masa de la segunda estrella (m3): "
read(*,*) m3

!rk4 comienzo de iteracion
do i=0,n-1
kx1 = fx(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))
lx1 = fvx(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))

kx2 = fx(t(i)+(h/2),x(i)+(h/2)*kx1,y(i),z(i),vx(i)+(h/2)*lx1,vy(i),vz(i))
lx2 = fvx(t(i)+(h/2),x(i)+(h/2)*kx1,y(i),z(i),vx(i)+(h/2)*lx1,vy(i),vz(i))

kx3 = fx(t(i)+h/2,x(i)+(h/2)*kx2,y(i),z(i),vx(i)+(h/2)*lx2,vy(i),vz(i))
lx3 = fvx(t(i)+h/2,x(i)+(h/2)*kx2,y(i),z(i),vx(i)+(h/2)*lx2,vy(i),vz(i))

kx4 = fx(t(i)+h,x(i)+h*kx3,y(i),z(i),vx(i)+h*lx3,vy(i),vz(i))
lx4 = fvx(t(i)+h,x(i)+h*kx3,y(i),z(i),vx(i)+h*lx3,vy(i),vz(i))

x(i+1) = x(i) + (h/6)*(kx1+2*kx2+2*kx3+kx4)
vx(i+1) = vx(i) + (h/6)*(lx1+2*lx2+2*lx3+lx4)
! ------------------------------------------------
ky1 = fy(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))
ly1 = fvy(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))

ky2 = fy(t(i)+(h/2),x(i),y(i)+(h/2)*ky1,z(i),vx(i),vy(i)+(h/2)*ly1,vz(i))
ly2 = fvy(t(i)+(h/2),x(i),y(i)+(h/2)*ky1,z(i),vx(i),vy(i)+(h/2)*ly1,vz(i))

ky3 = fy(t(i)+h/2,x(i),y(i)+(h/2)*ky2,z(i),vx(i),vy(i)+(h/2)*ly2,vz(i))
ly3 = fvy(t(i)+h/2,x(i),y(i)+(h/2)*ky2,z(i),vx(i),vy(i)+(h/2)*ly2,vz(i))

ky4 = fy(t(i)+h,x(i),y(i)+h*ky3,z(i),vx(i),vy(i)+h*ly3,vz(i))
ly4 = fvy(t(i)+h,x(i),y(i)+h*ky3,z(i),vx(i),vy(i)+h*ly3,vz(i))

y(i+1) = y(i) + (h/6)*(ky1+2*ky2+2*ky3+ky4)
vy(i+1) = vy(i) + (h/6)*(ly1+2*ly2+2*ly3+ly4)
! ----------------------------------------------------------

kz1 = fz(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))
lz1 = fvz(t(i),x(i),y(i),z(i),vx(i),vy(i),vz(i))

kz2 = fz(t(i)+(h/2),x(i),y(i),z(i)+(h/2)*kz1,vx(i),vy(i),vz(i)+(h/2)*lz1)
lz2 = fvz(t(i)+(h/2),x(i),y(i),z(i)+(h/2)*kz1,vx(i),vy(i),vz(i)+(h/2)*lz1)

kz3 = fz(t(i)+h/2,x(i),y(i),z(i)+(h/2)*kz2,vx(i),vy(i),vz(i)+(h/2)*lz2)
lz3 = fvz(t(i)+h/2,x(i),y(i),z(i)+(h/2)*kz2,vx(i),vy(i),vz(i)+(h/2)*lz2)

kz4 = fz(t(i)+h,x(i),y(i),z(i)+h*kz3,vx(i),vy(i),vz(i)+h*lz3)
lz4 = fvz(t(i)+h,x(i),y(i),z(i)+h*kz3,vx(i),vy(i),vz(i)+h*lz3)

z(i+1) = z(i) + (h/6)*(kz1+2*kz2+2*kz3+kz4)
vz(i+1) = vz(i) + (h/6)*(lz1+2*lz2+2*lz3+lz4)

enddo
! guardar datos
open(unit=2, file="rk4_3cuerpos.dat")
do i=0,n
write(2,*) t(i), x(i), y(i), z(i), vx(i), vy(i), vz(i)
enddo
close(2)


contains
! Funciones 
real function fvx(t,x,y,z,vx,vy,vz)
fvx=-g*m2*x/(sqrt(x**2+y**2+(d-z)**2)**3) - g*(m3)*x/(sqrt(x**2+y**2+(d+z)**2)**3) 
end function

real function fvy(t,x,y,z,vx,vy,vz)
fvy=-g*m2*y/(sqrt(x**2+y**2+(d-z)**2)**3) - g*(m3)*y/(sqrt(x**2+y**2+(d+z)**2)**3) 
end function

real function fvz(t,x,y,z,vx,vy,vz)
fvz=-g*m2*(z-d)/(sqrt(x**2+y**2+(d-z)**2)**3) - g*(m3)*(z+d)/(sqrt(x**2+y**2+(d+z)**2)**3) 
end function

real function fx(t,x,y,z,vx,vy,vz)
fx=vx
end function

real function fy(t,x,y,z,vx,vy,vz)
fy=vy
end function

real function fz(t,x,y,z,vx,vy,vz)
fz = vz
end function

end 
