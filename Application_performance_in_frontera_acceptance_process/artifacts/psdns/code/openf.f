	subroutine openf
	use com
	implicit none

	character*5 numer
	character*16 fn
	integer nchar,i,ii,numc

	character*11 form,unform
      data form,unform/'formatted','unformatted'/
	character*20 caux

c Lines made redundant upon introduction of the new inpen.f
c routine, PK Yeung Dec 2007 have been removed

c lustr: logical unit for checkpoint data	 


	if (taskid.eq.0) then

	caux='chkptin'
      if ((istart.ne.0).or.(istart.eq.0.and.tfiuo.lt.0.))
c    1                 call fopen1 (iabs(istart),'chkptin',form)
     1                 call fopen1 (1,caux,form)
	caux='ranin'
      if (istart.gt.0) call fopen1 (luran1,caux,form)

	endif

	if (taskid.ne.0) return

	call time_stamp ('openf: before log')

 	caux='log'
 	call fopen1 (7,caux,form)
	call time_stamp ('openf: log')

	caux='courant'
 	call fopen1 (8,caux,form)
	call time_stamp ('openf: courant')

	caux='sptr1d'
 	call fopen1 (13,caux,form)
	call time_stamp ('openf: sptr1d')

	caux='corr1d'
 	call fopen1 (14,caux,form)
	call time_stamp ('openf: corr1d')

	caux='eulstat'
 	call fopen1 (15,caux,form)
	call time_stamp ('openf: eulstat')

#ifdef MODEL_SPECTRUM
      caux='rstensor'
      CALL fopen1(996, caux, form)
      CALL time_stamp('openf: rstensor')
      caux='sptr3d'
      CALL fopen1(997,caux,form)
      CALL time_stamp('openf: sptr3d')
      caux='tmij'
      CALL fopen1(998,caux,form)
      CALL time_stamp('openf: tmij')
      caux='metrics'
      CALL fopen1(999,caux,form)
      CALL time_stamp('openf: metrics')
#endif
      caux='spectrum.free'
      CALL fopen1(991,caux,form)
      CALL time_stamp('openf: spectrum.free')


	if (nc.eq.0) then
 	caux='vgskew'
 	call fopen1 (18,caux,form)
	caux='vgmomt'
 	call fopen1 (190,caux,form)
	caux='vgmomtn'
 	call fopen1 (191,caux,form)
	end if
 	caux='corrmat'
 	call fopen1 (30,caux,form)
	caux='trfsptr'
 	call fopen1 (31,caux,form)
	caux='vorstat'
 	if (iovor.ge.1) call fopen1 (35,caux,form)
	caux='forcinp'
 	call fopen1 (50,caux,form)

#ifdef MHD
        caux='msqvelgrad'
        call fopen1 (24,caux,form)
        caux='mhd.stat'
        call fopen1 (25,caux,form)
        caux='tkcomp'
        call fopen1 (33,caux,form)
        caux='tkcomp2'
        call fopen1 (34,caux,form)
        caux='rsani'
        call fopen1 (37,caux,form)
#endif
c
#ifndef NOSCALAR
c
	if (nc.gt.0) then
c
	caux='escout'
	call fopen1 (16,caux,form)
	caux='scgrad'
	call fopen1 (38,caux,form)
	if (any(grad.ne.0.)) then
	caux='scflxsp'
	call fopen1 (17,caux,form)
	caux='scgmom'
	call fopen1 (94,caux,form)
	end if
	if (any(grad.ne.0.).or.rrate.gt.0.) then
	caux='scgsptr'
	call fopen1 (36,caux,form)
	end if
c
	caux='cosptr1'
      if (nc.eq.2) call fopen1 (27,caux,form)
      if (nc.eq.3) then
	caux='cosptr1'
      call fopen1 (27,caux,form)
	caux='cosptr2'
      call fopen1 (28,caux,form)
	caux='cosptr3'
      call fopen1 (29,caux,form)
      end if
c
	end if
c
#endif
c
c need a strategy for these files to be opened--probably subcall in bouss module
#ifdef BOUSS
	caux='hsptrz'
      call fopen1 (61,caux,form)
	caux='hsptr'
      call fopen1 (62,caux,form)
#endif
c
#ifdef BOXPHY
	caux='dns_pdf_epsenst'
      call fopen1 (751,caux,form)
	caux='dns_max_epsenst'
      call fopen1 (752,caux,form)
	caux='vort_skewflat'
      call fopen1 (753,caux,form)
#endif

	return
	end
