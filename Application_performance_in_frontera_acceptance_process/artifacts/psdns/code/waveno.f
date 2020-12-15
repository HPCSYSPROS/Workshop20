        subroutine waveno
c
	use com
	include 'intvars'
c
c  routine to form differential operators
c
c       bk1 = b(1,1)**2 * kx
c       bk2 = b(2,2)**2 * ky + b(1,2) * b(2,2) * kx
c       bk3 = b(3,3)**2 * kz
c
c       bkk12 + bkk3 = laplacian
c
c  loop over meshes
c
        do 10 m=1,2
c
c  form wave numbers for x and y differentiation
c
        do 1 x=1,nxhpad
        bk1(x,m)=b11(m)*kx(x)
 1	continue
c
#ifdef SHEAR
        do 2 y=1,nypad
        k2(x,y,m)=ky(y)+b12(m)*kx(x)
        bk2(x,y,m)=b22(m)*k2(x,y,m)
        bkk12(x,y,m)=bk1(x,m)*kx(x)+bk2(x,y,m)*k2(x,y,m)
 2	continue
#endif
c
c  form wave numbers for z differentiation
c
        do 3 z=1,nzpad
        bk3(z,m)=b33(m)*kz(z)
3       bkk3(z,m)=bk3(z,m)*kz(z)
c
 10     continue
c
        return
        end
