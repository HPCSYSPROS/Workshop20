!
! output of Eulerian scalar statistiscs (calculated by sub. sptr)
!
      subroutine escout (lu,io,uny)
      
#ifndef NOSCALAR
      use comsp
      
      complex(b8) :: uny(ny,zjsz*xisz,3+nc)	  

        integer nshell
        real(b8) shell

      
      save ist,icgrad
      data ist/0/
      
      if (nc.eq.0) return
      
      if (shear.eq.0.) then
         write (lu,601) istep,time
      else
         write (lu,602) istep,time,shear*time
      end if
 
      q2=corr(kt(1))+corr(kt(2))+corr(kt(3))
      tau=q2/epslon
 
#ifndef SHELL_DK
        shell=1.0
        nshell=nxh*min(beta1,beta2,beta3)
#else
        shell=min(beta1,beta2,beta3)
        nshell=nxh
#endif

      do 10 is=1,nc
 
      i=3+is
 
! calculate fourth and sixth moments of the scalar spectrum
 
      spm4=0.
      spm6=0.
      do 15 l=1,nshell
      spm4=spm4+sijk(l,kt(i))
      spm6=spm6+gijk(l,kt(i))
 15   continue
      spm4=spm4/2./viscos
      spm6=spm6/2./viscos

#ifdef BOUSS
      write (lu,603) is,pr(is),0.,0., -1.
#else
      write (lu,603) is,pr(is),(grad(l,is),l=1,3)
#endif
 
      write (lu,201) real(uny(1,1,i)),corr(kt(i)),scdiss(is),spm4
 
      write (lu,202) corr(3+is),corr(5+nc+is),corr(6+2*nc+is),spm6
 
      write (lu,235) (lijk(kt(i),j),j=1,3)
      write (lu,236) (taylor(i,i,j),j=1,3)
 
! scalar time scale and mechanical-to-scalar time scale ratio
 
	sctime=0.
	ratio=0.
	if (scdiss(is).gt.0.) then
      sctime=corr(kt(i))/scdiss(is)
      ratio=tau/sctime
	end if
 
      write (lu,231) sctime,ratio
 
! include the Batchelor and Obukhov-Corrsin scales
! (e.g., M. Rogers p. 51)
 
      batch=klen/sqrt(pr(is))
      obukc=klen/pr(is)**.75
      write (lu,232) batch,obukc
 
      write (lu,223) min(beta1,beta2,beta3)
      write (lu,220) (eijk(l,kt(i)),l=1,nxh)
      write (lu,224)
      write (lu,220) (dijk(l,kt(i)),l=1,nxh)
 
 10   continue
 
! for scalars with mean gradients, write the scalar flux spectrum
! in the direction of the mean scalar gradient
 
      ist=ist+1
      if (ist.eq.1) then
      icgrad=0
      do 20 i=1,nc
      if (grad(1,i).ne.0..or.grad(2,i).ne.0..or.grad(3,i).ne.0.) then
      icgrad=1
      if (istart.eq.0) write (17,203) i,(grad(j,i),j=1,3)
      end if
 20   continue
      end if
 
      if (icgrad.eq.0) go to 90
 
      do 30 i=1,nc
 
      if (grad(1,i).ne.0.)  then
      ij=3+i
      write (17,204) i,'x',istep,time
      write (17,220) (eijk(k,ij),k=1,nshell)
      end if
 
      if (grad(2,i).ne.0.)  then
      ij=5+nc+i
      write (17,204) i,'y',istep,time
      write (17,220) (eijk(k,ij),k=1,nshell)
      end if
 
      if (grad(3,i).ne.0.)  then
      ij=6+2*nc+i
      write (17,204) i,'z',istep,time
      write (17,220) (eijk(k,ij),k=1,nshell)
      end if
 
 30   continue
 
 90   continue
 
      close (lu)
      open (lu,file='escout',form='formatted',iostat=ios,
     1      position='append')
 
! add output of scalar gradient spectrum, 5/8/03
 
      if (any(grad.ne.0.).or.rrate.gt.0.)  then
      do 40 ic=1,nc
      i=3+ic
      ij=kt(i)
      write (36,241) istep,time,ic
      write (36,220) (grijk(k,ij,1),k=1,nshell)
      write (36,242) istep,time,ic
      write (36,220) (grijk(k,ij,2),k=1,nshell)
      write (36,243) istep,time,ic
      write (36,220) (grijk(k,ij,3),k=1,nshell)
 40   continue
      end if
 
      return
 
601    format(/,'scalar output for step',i6,2x,'time=',1p,e14.5)
602     format(/,'scalar output for step',i6,2x,'time=',1p,e14.5,
     1         3x,'st=',1p,e14.5)
603     format('scalar #',i1,': Pr, mean grad=',1p,4e11.3)
!
 201  format (' mean, varce & dissip. :',1p,4e12.4)
 202  format (' scalar fluxes         :',1p,4e12.4)
 203  format ('scalar #',i1,'  mean gradients=',1p,3e12.4)
 204  format ('scalar #',i1,', ',a1,'-component scalar flux spectrum',
     1        ' at istep=',i5,'  time=',1p,e12.4)
!
 235  format (' scalar int. length    :',1p,3e12.4)
 236  format (' scalar taylor scale   :',1p,3e12.4)
 237  format ('time, scalar rms & diss:',1p,3e12.4)
 238  format ('number of scalars:',i2)
 231  format (' time scale and ratio  :',1p,2e12.4)
 232  format ('batchelor & obukhov-corrsin scales: ',1p,2e12.4)
 220  format ((3x,1p,10e13.5))
 223  format ('3-d scalar spectrum, shell thickness=',1p,e12.4)
 224  format ('3-d scalar dissipation spectrum')
 241  format ('spectrum of dphi/dx, istep=',i6,'  time=',1p,e12.4,
     1        '  scalar #',i1)
 242  format ('spectrum of dphi/dy, istep=',i6,'  time=',1p,e12.4,
     1        '  scalar #',i1)
 243  format ('spectrum of dphi/dz, istep=',i6,'  time=',1p,e12.4,
     1        '  scalar #',i1)
 
#endif
        end
