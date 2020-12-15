	module mpilag
c
#ifdef LAG

	use mpicom
c
	integer nbx,nby,nbz
	integer ndim
	integer ndpart,ndprop,nbsdim
c
        real, allocatable :: p(:,:),q(:,:),t(:,:),denom(:),ss(:)
c
      integer,save:: bxistart,bxiend,bxisize
      integer,save:: bzistart,bziend,bzisize
      integer,save:: bzjstart,bzjend,bzjsize
      integer, save, allocatable ::
     & bxist(:),bxien(:),bxisz(:),
     &   bzist(:),bzien(:),bzisz(:),
     & bzjst(:),bzjen(:),bzjsz(:)
c
      integer, save, dimension(:),allocatable ::
     &   bcnts_xi,bcnts_zi,bcnts_zj,
     &   bcnts_yj,bdisp_xi,bdisp_zi,bdisp_zj,bdisp_yj
c
c Added by PK Yeung, 2/11/09, for choice of intbf/spcal implementation
        integer numbatch,ibflag
        integer pim_flag
c
#endif

#ifdef CF_LAG
        integer, allocatable :: rantasks(:)
#endif


#ifndef LAG_DOUBLE
	integer, parameter:: p8 = b8
	integer, parameter:: pptype = mpireal
#else
	integer, parameter:: p8 = 8
	integer, parameter:: pptype = MPI_REAL8
#endif

	end module mpilag
!###############################################
c
c--------------------------------------------------------------------
	module compart
c
#ifdef LAG
	use comsp
	use mpilag
c
	integer pstart,itetra
	integer nop,ngp,nsubset
	integer kipran
	integer nrec1,nrec2,nmrec
c
	real gxyz(3),xyzl(3)
	integer nxyz(3),nd(3),nplu,iod(3,10),ipvg(9),lupop,iovel,ioposn
c
	integer, allocatable :: lpfunc(:),lpgrad(:),lplapl(:),ipdis(:),
     1				lpsecd(:),lpsdis(:)
	integer, allocatable :: npo(:),ipj1(:),ipj2(:)
c
	real(p8), allocatable :: bs(:,:,:,:)
c
	real, allocatable :: gx(:),gy(:),gz(:), ps(:)
c
	integer ibxyz
c
	real(p8), allocatable :: pp(:,:)
	real(p8), allocatable :: bfpos(:,:)
c
	real, allocatable :: dxpr(:),dypr(:),dzpr(:)
	
	integer iout
c
#ifdef LGRAD
	real, allocatable :: udxzr(:,:,:,:)
#endif
c
	integer nop2, nom
c
c#ifdef LAGSC2
	integer ndpart2,ndprop2
	integer ibxyz2
        real(p8), allocatable :: pp2(:,:)
        real(p8), allocatable :: pp2h(:,:)
	real(p8), allocatable :: bfpos2(:,:)
c#endif
c
c#ifdef MOL
	integer nompp,ngm,nm,nmsubset,mstart,minitseed
        integer, allocatable :: mseed(:)

	real(p8), allocatable :: mpp(:,:,:)
	real(b8), allocatable :: molsc(:)
c#endif
c
#endif
	end module compart
!#########################################################

	module lag_timers
c
	real(8) cpu_lag
        real(8) cpu_partic
        real(8) cpu_partic2(2)
        real(8) cpu_lagout
        real(8) cpu_partsp(0:5,2)
        real(8) cpu_spxyz(0:10,2)
        real(8) cpu_intbf(4,3,2)
        real(8) cpu_spcal(4,3,2)
        real(8) cpu_rkpart(2,2)
        integer nspxyz,nintbf,nspcal,npartsp,itrack,icpvrms,nspcal_io,nintbf_io
c
	end module lag_timers
