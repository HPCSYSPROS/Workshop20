      module precision
        
      include 'mpif.h'
!       integer, parameter:: b8 = selected_real_kind(6,70)
!       integer, parameter:: mpireal = MPI_REAL8

#ifndef DOUBLE_PREC	 
      integer, parameter:: b8 = 4
      integer, parameter:: mpireal = MPI_REAL
      integer, parameter:: mpicomplex = MPI_COMPLEX
#else
      integer, parameter:: b8 = 8
      integer, parameter:: mpireal = MPI_REAL8
      integer, parameter:: mpicomplex = MPI_COMPLEX16
#endif
 
      end module precision
 
!--------------------------------------------------
! global grid resolution in 3 dimensions 
 
      module param
          
      use precision
 
      integer  :: zero,one
      integer  :: nx,ny,nz,nxh,nxhp,nyh,nyhp,nzh,nzhp,nu,nxhp2pad,nyhpad,nyhppad
      integer  :: nxpad,nxhpad,nxhppad,nypad,nzpad,nyxhppad,nzhpad,nzhppad
      integer  :: nut
      integer  :: nc,ncd,ncop
      integer*8 Ntot
      integer  :: nxhp2,nyhp2,nzhp2
      real(b8) :: paramxx
      integer  :: kfor,k2fo
      integer  :: ncp
      integer  :: cache_l,NB1_X,NB1_Z,NB2_Z,NB2_Y 

      end module param
 
!----------------------------------------------------
! global variables related to MPI

      module mpicom
	
      use param
c
        integer, public :: num_thr
	integer, allocatable :: ia_st(:)
	integer, allocatable :: ia_sz(:),ia_xst(:),ia_zst(:)
	save ia_st,ia_sz,ia_xst,ia_zst

! mpi process info
      integer ierr, dims(3),  cartid(3)
      logical periodic(3),remain_dims(3)
      integer numtasks,taskid,mpierr
      integer iproc, jproc, ipid, jpid
      integer impid, ippid, jmpid, jppid

      real(b8) mpixx
 
      integer, dimension(:), allocatable :: iist,iien,iisz
      integer, dimension(:), allocatable :: jjst,jjen,jjsz
      integer, dimension(:), allocatable :: kist,kien,kisz
      integer, dimension(:), allocatable :: kistpad,kienpad,kiszpad
      integer, dimension(:), allocatable :: kjst,kjen,kjsz
      integer, save :: mystart_x,mystart_z
!array of one's
      integer, dimension(:), allocatable :: oneAr

      integer :: xist,xjst,xisz,xjsz,xien,xjen
      integer :: yist,yjst,yisz,yjsz,yien,yjen
      integer :: zist,zisz,zien
      integer :: zjst,zjsz,zjen
      integer :: num_al,max_al_x
      integer, allocatable :: max_al_z(:),cut_z(:),cut_z_i(:,:)
      integer, allocatable :: num_al_i(:)
	integer :: padx,pady,padz,szmax,szf

      integer mpi_comm_cart
      integer mpi_comm_row, mpi_comm_col
! mpi derived data types for implementing alltoallv using send-recvs
      integer IfCntMax,KfCntMax,JrCntMax,KrCntMax
      logical IfCntUneven,KfCntUneven,JrCntUneven,KrCntUneven
      integer,dimension(:),allocatable:: IfSndCnts,IfSndStrt
      integer,dimension(:),allocatable:: IfRcvCnts,IfRcvStrt
      integer,dimension(:),allocatable:: KfSndCnts,KfSndStrt
      integer,dimension(:),allocatable:: KfRcvCnts,KfRcvStrt
      integer,dimension(:),allocatable:: JrSndCnts,JrSndStrt
      integer,dimension(:),allocatable:: JrRcvCnts,JrRcvStrt
      integer,dimension(:),allocatable:: KrSndCnts,KrSndStrt
      integer,dimension(:),allocatable:: KrRcvCnts,KrRcvStrt
#ifdef INPUT_PEN
      integer,dimension(:),allocatable:: KfSndCnts_io,KfSndStrt_io
      integer,dimension(:),allocatable:: KfRcvCnts_io,KfRcvStrt_io
#endif
#ifdef OUTPUT_PEN
      integer,dimension(:),allocatable:: JrSndCnts_io,JrSndStrt_io
      integer,dimension(:),allocatable:: JrRcvCnts_io,JrRcvStrt_io
#endif
      integer,dimension(:,:),allocatable:: status
      integer,dimension(:), allocatable:: mymap,inverse_map
      integer, save :: bal_period

#ifdef ALLTOALLV
        integer, allocatable :: atac_start(:),atac_count(:)
        integer, allocatable :: atar_start(:),atar_count(:)
#endif

        integer, allocatable :: ipid_all(:),jpid_all(:)
        
      end module mpicom
 
!--------------------------------------------
 
! variables corresponding to old "com"
 
      module com
 
      use mpicom
 
c PKY, 7/10/11: rwall0 from main.f now in module
c
	real(8) rwall0
	real acccpu
	
 
! Jan 30, 2007: next 2 lines moved to the 'intvars' file
! integer coordinate indices
!     integer :: x,y,z,yg,zg,z1,z2,y1,y2,zp,yp,xp
 
! input/output logical units and filenames
      integer, allocatable :: luinit(:),kinit(:)
      character(8), allocatable ::  fninit(:)
	character(7), allocatable :: fnwu(:)
	integer, allocatable :: luwu(:)
	character*7 msrfn
 
! wavenumber components
      real(b8), allocatable :: kx(:),ky(:),kz(:)
      integer :: kmax

! integer to determine rkmethod
      integer(4) :: rkmethod

!variables that determine ratio for 3/2 dealiasing
	real(b8) :: ratio_pad_x, ratio_pad_y, ratio_pad_z

! grid metric factors and phase shifts
#ifdef RKFOUR
      real(b8) b11(2),b22(2),b33(2),b12(2),gsh(3,4)
#else
      real(b8) b11(2),b22(2),b33(2),b12(2),gsh(3,2)
#endif
      real(b8), dimension(:,:), allocatable :: bk1,bk3,bkk3
      complex(b8), dimension(:,:), allocatable :: sx,sy,sz
 
      complex(b8) imagi
c
c index for numerical sequence of checkpoints written by the code
	
	integer ichkpt

	integer mrcpin  ! size of batch for dwrtu1
	integer nbl_x   ! size of Xcoord batch for dwrtu1
 
!ccccccccccccccc old common blocks ccccccccccccccccccccccccc
 	integer  :: nsteps,iostep,istart,kstop,isave,ioflag,isflag,ictf
 	integer  :: iopflag, itsave(10),iocount,ioaxi(3)
 	real(b8) :: stend,entime,dtout,dtsave,tfout,tfsave,trmsh
 	logical  :: ctflag,schflag

	integer  :: irz,luran1,luran2,ksran,kranf1,kranf2
	real(b8) :: viscos 
	real(b8) :: a11,a22,a33,shear

	integer  :: ncps,nceach,lusc,luscn,luscd,luscxy,luscyp,scfreq
	real(b8) :: skc
	integer  :: idfdf,inorm,iscref
!	real(b8) :: dthist
	real(b8) :: beta1, beta2, beta3, bet12
	real(b8) :: cfl,dt,cnmax,dtmin
	integer  :: istep,jstep,kstep,istop,istep0
	real(b8) :: time0
	real(b8) :: velmax,pi,dt2
	integer  :: icase,iwfld,ifltr,icsptr,iovor,io1d
	
	integer  :: irmss,iwmss
 	character*20 msrdir,mswdir
 	character*25 scrdir
 	character*30 msrhead,mswhead

	integer  :: irot
	real(b8) :: rrate,rossby
	integer  :: iraxis

	integer  :: kshift,lustr,ioetd,npetd,luetd

	real(b8) :: time,s0
c
        character*256 :: indir_fn

!       common/run/nsteps,iostep,istart,kstop,isave,stend,entime,       
!    1             ioflag,isflag,dtout,dtsave,tfout,tfsave,ictf,ctflag,
!    1             iopflag,trmsh,itsave(10),iocount
!    1 /initc/irz,
!    1 /random1/luran1,luran2,ksran,kranf1,kranf2              
!    1 /fluid/viscos 
!    1 /strain/a11,a22,a33,shear              
!    1 /scalar/ncps,nceach,skc,
!    1         lusc,luscn,luscd,luscxy,luscyp,scfreq, 
!    1         idfdf,inorm,iscref,dthist,schflag
!    1 /scale/beta1,beta2,beta3,bet12         
!    1 /tstep/cfl,dt,cnmax,dtmin              
!    1 /steps/istep,jstep,kstep,istop,istep0,time0     
!    1 /misc/velmax,pi,dt2,icase,imagi,iwfld,ifltr,icsptr,iovor,io1d,
!    1       irmss,iwmss,msrfn,msrdir,mswdir,scrdir,msrhead,mswhead,
!    1       irot,rrate
!    1 /sft/kshift 
!    1 /chkpt/lustr            !,luwu,fnwu   
!    1 /etd/ioetd,npetd,luetd
!    1 /inv/time,s0
!ccccc end of old common blocks ccccccccccccccccccccccccc
 
!     integer nc,ncd
#ifndef NOSCALAR
      real(b8), allocatable :: pr(:),grad(:,:),iks(:),cb(:,:),svnorm(:)
      integer, allocatable :: kcps(:),kran(:)
	logical scgflag

	real(b8) :: gpdflim ! limits for PDF of scalar gradients
#endif
#ifdef BIF
	real bvfsq
#endif
 
      real(b8) :: tforce, epsfor, kforce, tfiuo
      real(b8), allocatable :: suo(:),svo(:)
      integer  nrproc
 
#ifdef SHEAR
      real(b8), dimension(:,:,:) allocatable: k2, bk2, bkk12
#endif

 
c      common/convec/icvc,flagc              
      integer icvc,npout
      logical flagc        

	integer seed_input

        integer*4, allocatable :: num_fft0(:)
        save  num_fft0
        integer*4, allocatable :: ixp1(:),ixpi(:)
        save ixp1,ixpi
        integer*4, allocatable :: iyp1(:),iypi(:)
        save iyp1,iypi
        integer*4, allocatable :: izp1(:),izpi(:)
        save izp1,izpi
#ifdef FFTW
        integer*8, save :: plan2_p1
        integer*8, allocatable  :: plan4_p1(:),plan4_p3(:)
	integer*8, allocatable  :: plan2_p2a(:),plan2_xk(:)
        integer*8, allocatable  :: plan1_p2b(:)
        integer*8, save :: plan4_p2a
        integer*8, save :: plan2_p2b, plan3_p2b
        integer*8, save :: plan2_p3
        integer*8, save :: plan1_kx,plan3_kx
        integer*8, save :: plan1_xk,plan3_xk
	save plan1_p2b
        save plan4_p1,plan4_p3
	save plan2_p2a,plan2_xk
        integer*8, allocatable  :: plan3_p2a(:),plan3_p2c(:),plan2_kx(:)
        save plan3_p2a,plan3_p2c,plan2_kx
#elif defined  ESSL
      integer :: cnaux,rnaux1,rnaux2
      real(8),allocatable :: caux1(:,:),caux2(:,:),raux1(:,:),raux2(:,:)
      real(8) :: raux3(1)
        integer cnaux1,cnaux2,rnaux3
#endif


 
c#ifdef FEK_FORC
        integer kf_shell
        real, allocatable ::  ek_nf1(:),ek_nf2(:),ekinf(:)
        real, allocatable ::  ekinf_input(:)
c#endif

#ifdef MHD
	integer  :: imaxis
        real conduc,stuart
        real, allocatable :: jdky(:,:),jdkc(:,:),jdk(:)
        real joule_diss(3)
#endif

#ifdef LVSTEP
c update velocity field only once every ivstep steps,
c with linear interpolation in time (done in physical space)
c
	integer ivstep
	real(b8), allocatable :: cvx(:,:,:,:,:,:)
#endif

! to test for unacceptable large t/s/p 
! the code inquires if the file tsp_max exists in com_set
! and reads the t/s/p beyond which the code stops
	real :: tsp_max

! to stop execution if something goes wrong. In main, the code
! inquires if teh file 'signal' exists. If it does, it reads 
! an integer: isignal.
! isignal=0: do a checkpoint and exit
! isignal>0: set nsteps=isignal 
! isignal=-1: read another line and set entime
	integer :: isignal

 
      end module com
 
!-----------------------------------------------------------
 
! variables corresponding to old "comp"
 
      module comp
 
      use com
 
! the main field arrays
 
      real(b8), allocatable :: un(:,:,:,:)
      real(b8), allocatable :: u(:,:,:,:)
#ifdef RKFOUR
      real(b8), allocatable :: u1(:,:,:,:)
#endif
c#ifdef CVC_PRESS
        complex, allocatable :: upy(:,:)
c#endif
c#ifdef VORPHY
        real, allocatable :: uom(:,:,:,:)
c#endif



#ifdef TEST_CYLIO
      complex(b8), allocatable :: ucomp(:,:,:,:)
#endif
 
! integrating factors in wavenumber space
 
#ifdef LINDBORG
      real(b8), allocatable :: xyfac(:,:),zfac(:)
	real(b8), allocatable :: xydif(:,:,:),zdif(:,:)
#else
      real(b8), allocatable :: xfac(:),yfac(:),zfac(:)
#ifndef NOSCALAR
	real(b8), allocatable :: xdif(:,:),ydif(:,:),zdif(:,:)
#endif
#endif
 
! mask for truncation of aliased modes
 
      logical, allocatable :: mask(:,:)

! energy input from forcing
      real(b8), allocatable :: efki(:,:),efkz(:,:),efk(:)
      real(b8) :: eirate(3),erate
      complex(b8), allocatable :: for(:,:,:,:)

	real(b8) :: tke,epslon
	real(b8), allocatable :: ek(:),dk(:)
 
! used in proc2a and proc2b and associated routines only
 
	complex(b8), allocatable :: ut(:,:,:)
 
!	complex(b8), allocatable :: bk2i(:)
 	complex(b8) bk2i

#ifdef BOUSS
	real(b8), allocatable :: ubuoy(:,:,:,:)
#endif

	real(8) clock_chkpt

      end module comp
 
!------------------------------------
 
      module comsp
 
      use comp
 
	integer :: mxyz
	real(b8) :: ekl,skew
	real(b8) :: klen,kvel,ktime
	real(b8), allocatable ::  laak(:,:),taak(:,:),raak(:)
	real(b8), allocatable ::  tmre(:), rms(:),ett(:)
	real(b8), allocatable ::  sk(:),tmij(:,:)
	real(b8), allocatable ::  lijk(:,:),corr(:)
	integer, allocatable :: kt(:)
 
	complex(b8), allocatable :: xrij(:,:,:),yrij(:,:,:),zrij(:,:,:)
!
#ifdef OUT1D
	complex(b8), allocatable :: eaak1(:,:),eaak2(:,:),eaak3(:,:)
	complex(b8), allocatable :: qijx(:,:),qijy(:,:),qijz(:,:)
#endif

        real(b8), allocatable :: eijky(:,:,:),eijk(:,:),tfact(:)
        real(b8), allocatable :: dijky(:,:,:),dijk(:,:)
        real(b8), allocatable :: sijky(:,:,:),sijk(:,:)
        real(b8), allocatable :: gijky(:,:,:),gijk(:,:)
        real(b8), allocatable :: grijky(:,:,:,:),grijk(:,:,:)
        real(b8), allocatable :: ekky(:,:),ekk(:),s(:)
        real(b8), allocatable :: taylor(:,:,:), der2y(:,:),der2ij(:)
 
	real(b8), allocatable :: kx2(:),ky2(:),kz2(:)
	real(b8), allocatable :: vijky(:,:,:),vijk(:,:),vk(:)
	real(b8) vij(6),vvij(3,3),tfact_x
 
	real(b8), allocatable :: umean(:),varce(:),skcof(:),flat(:)
      integer, allocatable :: imean(:)
 
#ifndef NOSCALAR
	integer ncgd
	real(b8), allocatable :: scdiss(:)
	real(b8), allocatable :: scgmsq(:,:),scgcor(:,:),scgcov(:,:)
	real(b8), allocatable :: scgvar(:,:)
#endif


      end module comsp

! ------------------------------------
 
      module force_mod
	use precision
	complex(b8), allocatable :: velinc(:,:,:,:)
	end module force_mod

! -------------------------------------	

	module timers_comm
      real(8) :: t1_comm,t2_comm,t3_comm,t4_comm,t4t_comm,tp1_comm
      real(8) :: gt1_comm,gt2_comm,gt3_comm,gt4_comm,gt4t_comm,gtp1_comm
	real(8) :: t_alltoall
	integer :: i1_comm, i2_comm, i3_comm, i4_comm, i4t_comm, ip_comm
	end module timers_comm

        module timers_comp
        real(8) :: tcpu_fft,tcpu_other
        end module timers_comp
c
	module timers_io
	integer :: iread_io,iwrite_io
	real(8) :: tread_io,twrite_io
	end module timers_io
c
        module timers_tran
        real(8) t_xk(4,2),t_kx(4,2)
        real(8) t_xkcomm1(4,3),t_xkcomm2(4,3)
        real(8) t_kxcomm1(4,3),t_kxcomm2t(4,3),t_kxcomm2(4,3)
        end module timers_tran
c
	module timers_rkstep
	integer ncpusteps,ncpusteps_io
	real(8) t_rks(10,2)
	real(8) t_itrans(4,3)
	real(8) t_trans(4,3)
	end module timers_rkstep

#ifdef MODEL_SPECTRUM
      ! Parameters and variables to be used when setting the initial energy
      ! spectrum to the model spectrum in Pope (pp. 232).
      MODULE spectrum
         ! Required modules.
         USE precision,ONLY: b8
         ! Parameters that the user will supply.
         REAL(b8) :: kmaxetaSpec, boxint, kol
         ! Model parameters.
         REAL(b8),PARAMETER :: p0 = 2.0
         REAL(b8),PARAMETER :: beta = 5.2
         REAL(b8),PARAMETER :: ceta = 0.4
         ! Parameters that are calculated after inputs are set.
         REAL(b8) :: cL, intlen, eta, eps
#ifdef SPECTRAL_TRUNCATION
         ! Cutoff wavenumbers used to truncate modes for non-cubic grids.
         REAL(b8) :: kLower, kUpper
#endif
      END MODULE spectrum
#endif

c --------------------------------------
#ifdef RANDOM_SEED
        module ranseed
        integer iseed,seed_size
        integer, allocatable :: rseed(:)
        end module ranseed
#endif
   
