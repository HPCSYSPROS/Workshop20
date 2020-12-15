      subroutine pptoid(xloc,zloc,tid,nid,shif)

! map an x,z location to particular ipid, jpid, taskid
! and neighboring id
c
	use compart
        implicit none
c

        integer itid,jtid,tid,nid,jj,ib
        real(p8) xloc,zloc,xbox,zbox,shif
c
c


        xbox = xloc + shif*nx - gx(1)
        xbox = mod(xbox,1.*nx)
!        ib= floor(xbox/nx)
!        xbox = xbox - ib*float(nx)


        zbox = zloc + shif*nz - gz(1)
        zbox = mod(zbox,1.*nz)
!        ib = floor(zbox/nz)
!        zbox = zbox - ib*float(nz)

        itid =  xbox/xb
        jtid = zbox/zb

        tid = jtid*iproc + itid

        nid=-1

        if(tid.eq.taskid) nid=0

        if(tid.ne.taskid) then
        do jj=1,8
        if(tid.eq.nbtask(jj)) nid = jj
        enddo
        endif


        return

 
        end

