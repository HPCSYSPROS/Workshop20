#define SPATIAL_DETERMINANT(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_) \
  (-(gxz_)**2*(gyy_) + 2.0d0*(gxy_)*(gxz_)*(gyz_) - (gxx_)*(gyz_)**2 - (gxy_)**2*(gzz_) \
   + (gxx_)*(gyy_)*(gzz_))

#define DOTP(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,x1_,y1_,z1_,x2_,y2_,z2_) \
 ( (gxx_)*(x1_)*(x2_)+(gyy_)*(y1_)*(y2_)+(gzz_)*(z1_)*(z2_)+ \
   (gxy_)*( (x1_)*(y2_)+(y1_)*(x2_) )+(gxz_)*( (x1_)*(z2_)+(z1_)*(x2_) )+\
   (gyz_)*( (y1_)*(z2_)+(z1_)*(y2_) ) )

#define DOTP2(gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,x_,y_,z_)        \
 ( (gxx_)*(x_)**2+(gyy_)*(y_)**2+(gzz_)*(z_)**2+ \
  2.0*( (gxy_)*(x_)*(y_)+(gxz_)*(x_)*(z_)+(gyz_)*(y_)*(z_) ) )


#define IF_BELOW_ATMO(rho, rho_min, rho_tol, r)  \
  dummy1 = atmo_tolerance_radius &&\
  dummy2 = atmo_falloff_radius &&\
  if (r .gt. atmo_tolerance_radius) then &&\
     dummy1 = r &&\
  endif &&\
  if (r .gt. atmo_falloff_radius) then &&\
     dummy2 = r &&\
  endif &&\
  if (rho .le. rho_min*(1.0d0 + rho_tol * (dummy1/atmo_tolerance_radius)**atmo_tolerance_power) * (atmo_falloff_radius/dummy2)**atmo_falloff_power)



!#define ATMOCHECK(rho, rho_min, rho_tol, r)  \
!  (atmo_type .eq. 0 .and. rho .le. rho_min*(1.0d0 + rho_tol)) .or. \
!  (atmo_type .eq. 1 .and. rho .le. rho_min*(1.0d0 + rho_tol * (atmo_tolerance_radius/r)**atmo_tolerance_power)) .or. \
!  (atmo_type .eq. 2 .and. rho .le. rho_min*(1.d00 + rho_tol) * (atmo_falloff_radius/r)**atmo_falloff_power) .or. \
!  (atmo_type .eq. 3 .and. rho .le. rho_min*(1.0d0 + rho_tol * (r/atmo_tolerance_radius)**atmo_tolerance_power) * (r/atmo_falloff_radius)**atmo_falloff_power)


#define SET_ATMO_MIN(rho_, rho_min, r) &&\
  dummy1 = atmo_falloff_radius &&\
  if (r .gt. atmo_falloff_radius) then &&\
     dummy1 = r &&\
  endif &&\
  rho_ = rho_min * (atmo_falloff_radius/dummy1)**atmo_falloff_power

