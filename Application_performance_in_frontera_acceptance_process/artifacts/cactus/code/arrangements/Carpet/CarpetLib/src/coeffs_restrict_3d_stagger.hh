// We are interpolating to point i-1/2 (a coarse gridpoint) using staggered fine
// grid _array index_ i.
//
// Interpolation coefficients are set by solving the following equation for
// (a,b,c)=(f11,f12,f13):
//  THE FOLLOWING DESCRIPTION IS ONLY CORRECT FOR INTERPOLATING TO STAGGERED
//  GRIDPOINTS.
//  PLEASE SEE THE MATHEMATICA WORKSHEETS (.math files) IN THIS DIRECTORY FOR
//  THE ACTUAL ALGORITHM.
// a f(x1) + b f(x2) + c f(x3) = f(x4), where (x1,x2,x3,x4)=(-2,-1,0,-1/2).
// If f(x) = k = constant, we have
// a + b + c = 1
// If f(x) = k x, we have
// a x1 + b x2 + c x3 = x4
// Similarly for f(x) = k x^2, so we have a quadratic polynomial fit.
/*
      double f0 = -0.125; // 2*i-2 on fine grid (there is a +1 offset below due
   to the staggering)
      double f1 = 0.75;   // 2*i-1 on fine grid
      double f2 = 0.375;  // 2*i   on fine grid
*/

/*
      double f0 = -1.0/16.0; // 2*i-2 on fine grid (there is a +1 offset below
   due to the staggering)
      double f1 = 9.0/16.0;  // 2*i-1 on fine grid
      double f2 = 9.0/16.0;  // 2*i   on fine grid
      double f3 = -1.0/16.0; // 2*i+1 on fine grid
*/

#define ORDER_STAG 3

T coeff[10];
T coeff_i[10];
// 2ND ORDER_STAG:
if (ORDER_STAG == 2) {
  // interpolate to point 2*i+1/2, using points 2*i-1,2*i,2*i+1
  coeff[1] = -0.125; // 2*i-1 on fine grid
  coeff[2] = 0.75;   // 2*i   on fine grid
  coeff[3] = 0.375;  // 2*i+1 on fine grid
}

// 3RD ORDER_STAG:
if (ORDER_STAG == 3) {
  // interpolate to point 2*i+1/2, using points 2*i-1,2*i,2*i+1,2*i+2
  coeff[1] = -1.0 / 16.0; // 2*i-1 on fine grid
  coeff[2] = 9.0 / 16.0;  // 2*i   on fine grid
  coeff[3] = 9.0 / 16.0;  // 2*i+1 on fine grid
  coeff[4] = -1.0 / 16.0; // 2*i+2 on fine grid

  coeff_i[1] = 1.0 / 8.0; // 2*i-1 on fine grid
  coeff_i[2] = 6.0 / 8.0; // 2*i   on fine grid
  coeff_i[3] = 1.0 / 8.0; // 2*i+1 on fine grid
}

// 4TH ORDER_STAG:
if (ORDER_STAG == 4) {
  // interpolate topoint 2*i+1/2, using points 2*i-2,2*i-1,2*i,2*i+1,2*i+2
  coeff[1] = 3.0 / 128.0;  // 2*i-2
  coeff[2] = -5.0 / 32.0;  // 2*i-1
  coeff[3] = 45.0 / 64.0;  // 2*i
  coeff[4] = 15.0 / 32.0;  // 2*i+1
  coeff[5] = -5.0 / 128.0; // 2*i+2

  coeff_i[1] = -1.0 / 128.0; // 2*i-2 on fine grid
  coeff_i[2] = 5.0 / 32.0;   // 2*i-1 on fine grid
  coeff_i[3] = 45.0 / 64.0;  // 2*i   on fine grid
  coeff_i[4] = 5.0 / 32.0;   // 2*i+1 on fine grid
  coeff_i[5] = -1.0 / 128.0; // 2*i+2 on fine grid
}
