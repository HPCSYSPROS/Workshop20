	subroutine bijdecomp (uny)
#ifdef MHD
c
	use comsp
	implicit none
	include 'intvars'
c
	complex :: uny(ny,zjsz*xisz,3)
	real, parameter :: pir = .25/atan(1.)
	real rk1,rk3,rk3sq,rk13sq,rk2,rkmag,rkmagr,rk23,emag
	real e1(3),e2(3),e3(3),s1,s2,s3,nvec(3)
	real factor,bije,bije0,bijz,bijz0,bijs,bijs0,bij0
	complex upol,utor
	real epol,etor
        real(b8) beta_min  
	integer ik,a,yst
        real et,ep,et0,ep0
c
#ifdef SHELL_DK
        beta_min=min(beta1,beta2,beta3)
#endif
c
        s1=sqrt(b11(2))
        s2=sqrt(b22(2))
        s3=sqrt(b33(2))
c broadcast ek to all MPI tasks for the formation of bije
        call MPI_BCAST (ek,mxyz,mpireal,0,MPI_COMM_WORLD,mpierr)
c
c unit vector along imaxis direction
c
        nvec=0.
        nvec(imaxis)=1.
c
	bijs = 0. !isotropic part
	bije = 0. !directional anisotropy part
	bijz = 0. !polarization anisotropy part
        et = 0. !toroidal energy
        ep = 0. !poloidal energy
c
	xp = mystart_x
	z = mystart_z
c
	do 100 a=1,num_al
c
	x = xp + xist - 1
c
         if(x .eq. 1) then
            tfact_x=1
         else
            tfact_x=2
         endif
c
        rk1=s1*kx(x)
	rk3=s3*kz(z)
	rk3sq=rk3**2
	rk13sq=rk1**2+rk3sq
c
      yst=1
c 
c special treatment for cross-product in e1, e2
c no need to compute Fourier modes parallel to imaxis
c since their contribution to Reynolds stress concerned is zero
c
      if (imaxis.eq.1.and.z.eq.1) yst=2
      if (imaxis.eq.3.and.x.eq.1) yst=2
      if (imaxis.eq.2.and.x.eq.1.and.z.eq.1) yst=ny+1
c
      do 110 y=yst,ny
         if(.not.mask(y,a)) goto 110
c
         rk2=s2*ky(y)
         rkmag=sqrt(rk13sq+rk2**2)
         rkmagr=1./rkmag
         rk23=sqrt(rk2*rk2+rk3sq)
c
c define the new unit vectors: e1 = k/|k| x n, e2 = k/|k| x e1
c
        e3(1)=rk1/rkmag
        e3(2)=rk2/rkmag
        e3(3)=rk3/rkmag
c 
c emag = cos(theta)
c where theta=angle between wavenumber vector k and imaxis vector
c factor = sin^2(theta)
c
	emag=nvec(1)*e3(1)+nvec(2)*e3(2)+nvec(3)*e3(3)
	factor=1.-emag*emag
c
        e1(1)=e3(2)*nvec(3)-e3(3)*nvec(2)
        e1(2)=e3(3)*nvec(1)-e3(1)*nvec(3)
        e1(3)=e3(1)*nvec(2)-e3(2)*nvec(1)
        emag=sqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
        e1(:)=e1(:)/emag
        e2(1)=e3(2)*e1(3)-e3(3)*e1(2)
        e2(2)=e3(3)*e1(1)-e3(1)*e1(3)
        e2(3)=e3(1)*e1(2)-e3(2)*e1(1)
        emag=sqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
        e2(:)=e2(:)/emag
c
#ifdef SHELL_DK
        ik=rkmag/beta_min+1.5
#else
        ik=rkmag+1.5
#endif
c
c \vec{uny}=utor \vec(e1) + upol \vec(e2)
c use dot product to obtain: 
c utor: toroidal (along latitudinal dir)
c upol: poloidal (along longitudinal dir)
c
	utor=s1*uny(y,a,1)*e1(1)+s2*uny(y,a,2)*e1(2)+s3*uny(y,a,3)*e1(3)
	upol=s1*uny(y,a,1)*e2(1)+s2*uny(y,a,2)*e2(2)+s3*uny(y,a,3)*e2(3)
c
c poloidal e & toroidal e
c 0.5 factor to be applied at the end by taskid 0
c
	epol = real(upol*conjg(upol))
	etor = real(utor*conjg(utor))
        et = et + tfact_x*etor
        ep = ep + tfact_x*epol
c
c bijz: polarization anisotropy
c       integral Z*sin^2(theta) over all modes
c       Z = epol - etor
c
c bije: directional anisotropy
c       integral (e - E(k)/(4pi k^2))*sin^2(theta) over all modes
c       e: kinetic energy in this mode
c       E(k): kinetic energy spectrum
c
c bijs: isotropic part
c       integral E(k)/(4pi k^2)*sin^2(theta) over all modes
c
c here 'i,j' correspond to index along imaxis (in mind)
c e.g. if imaxis.eq.1 i=1,j=1, we expect the relation
c       u1u1 = bijs + bije + bijz
c it is expected that in isotropic turbulence, 
c bijs = 2K/3
c thus bije and bijz characterize any anisotropy
c 
c u2u2, u3u3 can be expressed by u1u1 using formula
c see "Energy transfer in rotating turbulence"
c
	bijz = bijz + tfact_x*factor*(epol - etor)
c another 0.5 factor to be applied by taskid 0
	bije = bije + tfact_x*factor*
     1	(epol + etor - ek(ik)*.5*pir*rkmagr*rkmagr)
c 
c another 0.25 factor to be applied by taskid 0
	bijs = bijs + tfact_x*factor*(ek(ik)*pir*rkmagr*rkmagr)
c
 110	continue
c
c special treatment modes with k // imaxis
c they do not contribute to related anisotropy due to sin^2 theta
c but do contribute to toroidal and poloidal energy
        if (imaxis.eq.1.and.z.eq.1) then
        yst=1
        et = et + tfact_x*(b11(2)*uny(yst,a,1)*conjg(uny(yst,a,1)) +
     1          b22(2)*uny(yst,a,2)*conjg(uny(yst,a,2)) +
     1          b33(2)*uny(yst,a,3)*conjg(uny(yst,a,3)) )
        else if (imaxis.eq.3.and.x.eq.1) then
        yst=1
        else if (imaxis.eq.2.and.x.eq.1.and.z.eq.1) then

        end if
c
	call next_xz (xp,z)
 100	continue
c
	call MPI_REDUCE (bije,bije0,1,mpireal,
     1		MPI_SUM,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE (bijz,bijz0,1,mpireal,
     1		MPI_SUM,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE (bijs,bijs0,1,mpireal,
     1		MPI_SUM,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE (et,et0,1,mpireal,
     1		MPI_SUM,0,MPI_COMM_WORLD,mpierr)
	call MPI_REDUCE (ep,ep0,1,mpireal,
     1		MPI_SUM,0,MPI_COMM_WORLD,mpierr)
c
	if(taskid.eq.0) then
c
	if(imaxis.eq.1) then
	bij0=.5*corr(1)/tke - 1./3.
	else if(imaxis.eq.2) then
	bij0=.5*corr(4+nc)/tke - 1./3.
	else if(imaxis.eq.3) then
	bij0=.5*corr(6+2*nc)/tke - 1./3.
	end if
c
	bije0 = bije0*.25/tke
	bijz0 = bijz0*.25/tke
	bijs0 = bijs0*.125/tke
c
        if(istep.eq.0) then
        write(37,"(a6,5a13)") ' istep','    iso.part ','  dir.aniso ',
     1  '   pol.aniso ',' total aniso ','dir.ani.comp.'
        write(37,"(i6,1p,5e13.5)") istep,bijs0,bije0,bijz0,bij0,bij0-bijz0 
        else
        write(37,"(i6,1p,5e13.5)") istep,bijs0,bije0,bijz0,bij0,bij0-bijz0 
        end if
c 
c add et and ep spectra into axisymmetric spectra?
        et0 = .5*et0
        ep0 = .5*ep0
        write(6,"(a15,i6,1p,3e13.5)") "istep,et,ep,tke",istep,et0,ep0,tke
c
	end if
#endif
	return
c
	end subroutine
