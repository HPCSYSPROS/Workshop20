// Ax (the gridfunction we are prolongating) is defined on the semi-staggered
// grid (i,j+1/2,k+1/2)

// Here we set coefficients for interpolation to point i-1/4 (staggered coarse
// grid _array index_ i), using information i-2,i-1,i,i+1.
//  THE FOLLOWING DESCRIPTION IS ONLY CORRECT FOR INTERPOLATING TO STAGGERED
//  GRIDPOINTS.
//  PLEASE SEE THE MATHEMATICA WORKSHEETS (.math files) IN THIS DIRECTORY FOR
//  THE ACTUAL ALGORITHM.
// Interpolation coefficients are set by solving the following equation for
// (a,b,c,d)=(f11,f12,f13,f14):
// a f(x1) + b f(x2) + c f(x3) + d f(x4) = f(x5), where
// (x1,x2,x3,x4,x5)=(-2,-1,0,1,-1/4).
// If f(x) = k = constant, we have
// a + b + c + d = 1
// If f(x) = k x, we have
// a x1 + b x2 + c x3 + d x4 = x5
// Similarly for f(x) = k x^2 and f(x) = k x^3, so we have a cubic polynomial
// fit

T coeff[10][10];

if (ORDER == 2) {
  coeff[1][1] = 5.0 / 32.0;  // i-1
  coeff[1][2] = 15.0 / 16.0; // i
  coeff[1][3] = -3.0 / 32.0; // i+1
}

if (ORDER == 3) {
  /*
  coeff[1][1] = -45.0/256.0;  // i-2
  coeff[1][2] = 105.0/128.0;  // i-1
  coeff[1][3] =  0.0;         // i
  coeff[1][4] =  63.0/128.0;  // i+1
  coeff[1][5] = -35.0/256.0;  // i+2
  */
  coeff[1][1] = -5.0 / 128.0;  // i-2
  coeff[1][2] = 35.0 / 128.0;  // i-1
  coeff[1][3] = 105.0 / 128.0; // i
  coeff[1][4] = -7.0 / 128.0;  // i+1
}

if (ORDER == 4) {
  coeff[1][1] = -45.0 / 2048.0; // i-2
  coeff[1][2] = 105.0 / 512.0;  // i-1
  coeff[1][3] = 945.0 / 1024.0; // i
  coeff[1][4] = -63.0 / 512.0;  // i+1
  coeff[1][5] = 35.0 / 2048.0;  // i+2
  // coeff[1][1] = -5.0/128.0;  // i-2
  // coeff[1][2] = 15.0/32.0;   // i-1
  // coeff[1][3] = 45.0/64.0;   // i
  // coeff[1][4] = -5.0/32.0;   // i+1
  // coeff[1][5] =  3.0/128.0;  // i+2
}

// Here we set coefficients for interpolation to point i+1/4 (staggered coarse
// grid _array index_ i), using information i-1,i,i+1,i+2

if (ORDER == 2) {
  coeff[2][1] = -3.0 / 32.0; // i-1
  coeff[2][2] = 15.0 / 16.0; // i
  coeff[2][3] = 5.0 / 32.0;  // i+1
}

if (ORDER == 3) {
  /*
  coeff[2][1] = -35.0/256.0;  // i-2
  coeff[2][2] =  63.0/128.0;  // i-1
  coeff[2][3] =  0.0;         // i
  coeff[2][4] = 105.0/128.0;  // i+1
  coeff[2][5] = -45.0/256.0;  // i+2
  */
  coeff[2][1] = -7.0 / 128.0;  // i-1
  coeff[2][2] = 105.0 / 128.0; // i
  coeff[2][3] = 35.0 / 128.0;  // i+1
  coeff[2][4] = -5.0 / 128.0;  // i+2
}

if (ORDER == 4) {
  coeff[2][1] = 35.0 / 2048.0;  // i-2
  coeff[2][2] = -63.0 / 512.0;  // i-1
  coeff[2][3] = 945.0 / 1024.0; // i
  coeff[2][4] = 105.0 / 512.0;  // i+1
  coeff[2][5] = -45.0 / 2048.0; // i+2
}

// Here we set coefficients for interpolation to point i+1/2 (unstaggered coarse
// grid _array index_ i), using information i-1,i,i+1,i+2

if (ORDER == 2) {
  coeff[3][1] = -1.0 / 8.0;
  coeff[3][2] = 3.0 / 4.0;
  coeff[3][3] = 3.0 / 8.0;
}

if (ORDER == 3) {
  /*
  coeff[3][1] =-3.0/32.0;
  coeff[3][2] = 5.0/16.0;
  coeff[3][3] = 0.0;
  coeff[3][4] =15.0/16.0;
  coeff[3][5] =-5.0/32.0;
  */
  coeff[3][1] = -5.0 / 64.0;
  coeff[3][2] = 37.0 / 64.0;
  coeff[3][3] = 37.0 / 64.0;
  coeff[3][4] = -5.0 / 64.0;
  /*
  coeff[3][1] = -1.0/16.0;
  coeff[3][2] = 9.0/16.0;
  coeff[3][3] = 9.0/16.0;
  coeff[3][4] = -1.0/16.0;
  */
}

if (ORDER == 4) {
  coeff[3][1] = 63.0 / 2048.0;  // i-2
  coeff[3][2] = -103.0 / 512.0; // i-1
  coeff[3][3] = 781.0 / 1024.0; // i
  coeff[3][4] = 233.0 / 512.0;  // i+1
  coeff[3][5] = -97.0 / 2048.0; // i+2
}

// Here we set coefficients for interpolation to point i (unstaggered coarse
// grid _array index_ i), using information i-1,i,i+1

if (ORDER == 3) {
  coeff[4][1] = -1.0 / 32.0;
  coeff[4][2] = 34.0 / 32.0;
  coeff[4][3] = -1.0 / 32.0;
}

if (ORDER == 4) {
  coeff[4][1] = 7.0 / 2048.0;    // i-2
  coeff[4][2] = -23.0 / 512.0;   // i-1
  coeff[4][3] = 1109.0 / 1024.0; // i
  coeff[4][4] = -23.0 / 512.0;   // i+1
  coeff[4][5] = 7.0 / 2048.0;    // i+2
}
