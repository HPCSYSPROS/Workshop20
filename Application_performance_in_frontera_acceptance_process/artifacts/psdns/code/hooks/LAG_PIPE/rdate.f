      subroutine rdate (dmy,hms,iwall)
c
c this routine computes wall time referred back to beginning of current year
c
c entry arguments dmy and hms must have been updated by a call datetime
c
      integer dmy(3),hms(3),iwall
c
      integer days (12)
      data days/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /
c
c the "days" array holds the total number of days in all previous
c months of the current year
c
      leap=0
      if (mod(dmy(3),4).eq.0.and.mod(dmy(3),100).ne.0.
     1   and.dmy(2).gt.2) leap=1
c
      iwall=(days(dmy(2))+leap+dmy(1)-1)*24*3600
     1     +hms(1)*3600+hms(2)*60+hms(3)
c
      return
      end
