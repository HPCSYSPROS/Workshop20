      subroutine trfsp (uny,uy,m)
c
	use comsp
	implicit none
	include 'intvars'
c
	complex(b8) :: uny(ny,zjsz*xisz,3)
	complex(b8) :: uy(ny,zjsz*xisz,3)
	integer i,l
	
c
c icvc=1 not ported yet (check other routines)
c
c routine to calculate the transfer spectrum, called from sub. predic.
c between calls to proc3 and proc4a
c
c "un" contains the velocity field at time t,
c "u"  contains the convective term, divided by dt/2, at the same time
c (note the use of integrating factors)
c
c coded for isotropic turbulence only, 4/90
c
c	complex(b8), allocatable :: utmp1(:),utmp2(:),utmp3(:)
!	real(b8), allocatable :: ksq,kxsq,rksq
	real(b8) :: ksq,kxsq,rksq,der2,term

	real(b8), allocatable :: tij(:,:)
c
	real, allocatable :: tfk(:),tkij(:,:),tfky(:),tkijy(:,:)
 	save tfk,tkij,tfky,tkijy
c
	integer nkout,k,ik,lu,m
	real(b8) fmult,sum,deltak,s1,s2,s3,deltak_r
c
	complex(b8) vc(3),czero
	real(b8) kv(3),k3u,krdel(3,3)
c
	integer xst
	integer :: istat
c
      data krdel/1.,0.,0.,0.,1.,0.,0.,0.,1./
      data czero/(0.,0.)/
c
      if (icvc.ne.0) then
      write (6,*) 'trfsp: icvc.ne.0 not yet ported to MPL'
      stop 'aborts in sub. trfsp'
      end if
c
c      allocate (utmp1(nxh),utmp2(nxh),utmp3(nxh),stat=istat)
c      allocate (ksq(nxh),kxsq(nxh),rksq(nxh),stat=istat)
c      allocate (der2(nxh),term(nxh),tij(nxh,6),stat=istat)
c


	allocate (tfk(mxyz),tfky(mxyz))
	allocate (tkij(mxyz,6),tkijy(mxyz,6))
	tfky(:)=0.
	tkijy(:,:)=0.

	deltak=min(beta1,beta2,beta3)
        deltak_r = 1.0/deltak
c
c      do 2 x=1,nxh
c 2    kxsq(x)=b11(m)*kx(x)**2
c
      s1=sqrt(b11(m))
      s2=sqrt(b22(m))
      s3=sqrt(b33(m))
c
      if (icvc.eq.0) then

      xp = mystart_x
      z = mystart_z
      do l=1,num_al

         x = xp + xist-1

      kxsq = b11(m)*kx(x)**2
      k3u=b33(m)*kz(z)**2
c      kv(3)=kz(z)

      do 20 y=1,ny
      if (.not. mask(y,l)) go to 20

c
c
c      if (mask(xp,zp).eq.0.) go to 20
c
      der2=(b12(m)*kx(x)+ky(y))**2*b22(m)
      ksq=kxsq+der2+k3u
      ik=sqrt(ksq)*deltak_r+1.5
c
      term=real(conjg(uny(y,l,1))*uy(y,l,1))
      term=term+real(conjg(uny(y,l,2))*uy(y,l,2))
      term=term+real(conjg(uny(y,l,3))*uy(y,l,3))
c
      tfky(ik)=tfky(ik)+term*tfact(x)
c
 20   continue
c
      call next_xz(xp,z)
      enddo

#ifdef NOTYET

      else if (icvc.eq.1) then

      do 30 x=xst,nxh

      if (mask(x,y).eq.0.) go to 30

      kv(2)=(b12(m)*kx(x)+ky(y))*s2
      der2(x)=(b12(m)*kx(x)+ky(y))**2*b22(m)
      ksq(x)=kxsq(x)+der2(x)+k3u
      k=sqrt(ksq(x))+1.5

      kv(1)=kx(x)

 the fourier transform of u.grad(u)-grad(p) is the projection of the
 fourier transform of u.grad(u) normal to the wavevector
 (lesieur, p. 55)

 project convective term normal to the wave-vector

      rksq=1./ksq(x)
      do 33 i=1,3
      vc(i)=czero
      do 33 j=1,3
      vc(i)=vc(i)+uc(x,y,z,j)*(krdel(i,j)-kv(i)*kv(j)*rksq(x))
 33   continue

      utmp1(x)=conjg(un(x,y,z,1))
      utmp2(x)=conjg(un(x,y,z,2))
      utmp3(x)=conjg(un(x,y,z,3))
      tij(x,1)=real(utmp1(x)*vc(1))
      tij(x,2)=.5*real(utmp1(x)*vc(2)+utmp2(x)*vc(1))
      tij(x,3)=.5*real(utmp1(x)*vc(3)+utmp3(x)*vc(1))
      tij(x,4)=real(utmp2(x)*vc(2))
      tij(x,5)=.5*real(utmp2(x)*vc(3)+utmp3(x)*vc(2))
      tij(x,6)=real(conjg(un(x,y,z,3))*vc(3))
c
      term(x)=tij(x,1)+tij(x,4)+tij(x,6)
c
      do 35 ij=1,6
 35   tkijz(k,ij)=tkijz(k,ij)+tfact(x)*tij(x,ij)

 30   continue

 10   continue

 100  continue

#endif



      end if
c
c let task 0 collect and write the results
c
      call MPI_REDUCE (tfky,tfk,nxh,MPI_REAL,MPI_SUM,0,
     1                   MPI_COMM_WORLD,mpierr)
c
      if (taskid.ne.0) go to 90
c
	nkout=nxh
c
#ifdef NOTYET
      if (icvc.eq.1) then
      do 130 k=1,nkout
      tfk(k)=tkij(k,1)+tkij(k,4)+tkij(k,6)
 130  continue
      end if
      if (icvc.eq.1) then
      fmc=2.*fmult
      write (34,212) istep-1,time-dt
      do 220 ij=1,6
      do 225 k=1,nkout
      tkij(k,ij)=tkij(k,ij)*fmc
 225  continue
      write (34,202) ij,(tkij(k,ij),k=1,10)
      write (34,201) (tkij(k,ij),k=11,nkout)
 220  continue
      end if
#endif
c
      sum=0.
      fmult=2./dt
      do 210 k=1,nkout
      tfk(k)=tfk(k)*fmult
      sum=sum+tfk(k)
 210  continue
      write (6,*) 'total of t(k)=',sum
      lu=31
      write (lu,211) istep-1
      write (lu,201) (tfk(k),k=1,nkout)
c
c
 201  format ((3x,1p,10e13.5))
 202  format (i3,1p,10e13.5)
 211  format ('transfer spectrum at istep=',i5)
 212  format ('component transfer spectra at istep=',i5,
     1        4x,'time=',1p,e13.5)
c
 90   continue
c
c      deallocate (utmp1,utmp2,utmp3,stat=istat)
c      deallocate (ksq,kxsq,rksq,stat=istat)
c      deallocate (tij,stat=istat)
      deallocate (tfk,tfky,tkij,tkijy)
c,tkij,tkijy,stat=istat)
c


      return
      end
