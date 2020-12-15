void calc_prim(int EOS, double * EOS_const, double S_squared, double BdotS, double B_squared, double * con, double * prim, double * x, 
	       double B1_cov, double B2_cov, double B3_cov, struct nuc_eos_vars *struct_ptr){
  
  double W = x[0];
  double Z = x[1];
//  double P  = x[2];
  double T  = x[2];

  /* Recover the primitive variables from the scalars (W,Z) 
   * and conserved variables
   * 
   * Eq (23)-(25) in Cerdá-Durán et al.
   */
  
  prim[RHO] = con[D]/W;
  prim[v1_cov] = (con[S1_cov] + (BdotS)*B1_cov/Z)/(Z+B_squared);
  prim[v2_cov] = (con[S2_cov] + (BdotS)*B2_cov/Z)/(Z+B_squared);
  prim[v3_cov] = (con[S3_cov] + (BdotS)*B3_cov/Z)/(Z+B_squared);
//  if(EOS!=1)
//1  prim[EPS] = -1.0 + Z/(con[D]*W) - P*W/con[D];
//  prim[10] = P;
  prim[10] = P_iter_(EOS,EOS_const,prim,con,x,struct_ptr);
  prim[EPS] = -1.0 + Z/(con[D]*W) - prim[10]*W/con[D];
  prim[11] = W;
  prim[B1_con] = con[B1_con];
  prim[B2_con] = con[B2_con];
  prim[B3_con] = con[B3_con];
  if(EOS==3){
//    prim[YE] = con[YE]/con[D]; // THIS DOESN'T DO ANYTHING..... ***********
    prim[TEMP] = T;
  }
  
}

void NR_step_3D(int EOS, double * EOS_const, double S_squared, double BdotS, double B_squared, double * con, double * prim, 
		double * x, double dx[3], double f[3], double T1, double Ecalc1, double * Ecalc, double B1_cov, double B2_cov, double B3_cov,struct nuc_eos_vars *struct_ptr, int stepsize)
{
  /* Finding the roots of f(x):
   * 
   * x_{n+1} = x_{n} - f(x)/J = x_{n} + dx_{n}
   * 
   * where J is the Jacobian matrix J_{ij} = df_i/dx_j
   * 
   * Here, we compute dx = [dW, dZ, dP]
   */
  
  double J[3][3];
  double W = x[0];
  double Z = x[1];
// 1  double P = x[2];
  double T = x[2];
// 1  double Pcalc,dPdW,dPdZ;
  double E,P,dEdW,dEdZ,dEdT,dPdT;
  P = P_iter_(EOS,EOS_const,prim,con,x,struct_ptr);
  E = -1.0 + Z/(con[D]*W) - P*W/con[D];
//  printf("x = %15e,%15e,%15e\n",x[0],x[1],x[2]);
//  printf("E,P: %15e,%15e\n",E,P);
//  printf("NR.cl\n");

  // Need dEdW, and dEdZ
  if(EOS==3)
//    con[TEMP]=prim[TEMP];
// 1  P_dPdW_dPdZ(&Pcalc,&dPdW,&dPdZ,EOS,EOS_const,x,con,prim,struct_ptr);
  E_dEdW_dEdZ(Ecalc,&dEdW,&dEdZ,&dEdT,&dPdT,EOS,EOS_const,x,con,prim,struct_ptr,stepsize);
//  printf("Ecalc,dEdW,dEdZ,dEdT,dPdT: %15e,%15e,%15e,%15e,%15e\n",*Ecalc,dEdW,dEdZ,dEdT,dPdT);

  // If we are closer than 10^{-3} we switch to the secant method for dedt

//  if (f[2] < 1.0e-3) {
//    printf("Switching to secant method: Old dEdT = %15e\n", dEdT);
//    dEdT = -(*Ecalc-Ecalc1)/(T-T1);
//    printf("Switching to secant method: New dEdT = %15e\n", dEdT);
//  } 

  if(EOS==3)
  //  prim[TEMP]=con[TEMP];
  
  // d/dW (1)
  J[0][0] = 2.0*((Z + B_squared)*(Z + B_squared)-S_squared-(2.0*Z + B_squared)*((BdotS)*(BdotS))/(Z*Z))*W;
  double a = J[0][0];
  
  // d/dZ (1)
  J[0][1] = (2.0*(Z + B_squared) + (2.0/(Z*Z) + 2.0*B_squared/(Z*Z*Z))*(BdotS*BdotS))*W*W - 2.0*(Z + B_squared);
  double b = J[0][1];
  
  // d/dT (1)
  J[0][2] = 0;
  double c = J[0][2];
  
  // d/dW (2)
  J[1][0] = 2.0*(con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P)*W;
  double d = J[1][0];
  
  // d/dZ (2)
  J[1][1] = (-1.0-(BdotS*BdotS)/(Z*Z*Z))*W*W;
  double e = J[1][1];

  // d/dT (2)
  J[1][2] = W*W*dPdT;
  double fc = J[1][2];
  
  // d/dW (E-E(rho,T,Y_e))
  J[2][0] = dEdW;
  double g = J[2][0];
  
  // d/dZ (E-E(rho,T,Y_e))
  J[2][1] = dEdZ;
  double h = J[2][1];
  
  // d/dT (E-E(rho,T,Y_e))
  J[2][2] = dEdT;
  double k = J[2][2];
  
  // Compute f(x) from (84) and (85), and P-P(rho,eps)=0
  f[0] = ((Z+B_squared)*(Z+B_squared) - S_squared - (2.0*Z +B_squared)*(BdotS*BdotS)/(Z*Z))*W*W - ((Z + B_squared)*(Z + B_squared));
  f[1] = (con[TAU] + con[D] - Z - B_squared + (BdotS*BdotS)/(2.0*Z*Z) + P)*W*W + B_squared/2.0;
  f[2] = E - *Ecalc;

//  printf("f: %15e,%15e,%15e,%15e,%15e\n",f[0],f[1],f[2], E, *Ecalc);

//   prim[EPS] = fabs(prim[EPS]);
  
  // Compute the determinant
  double A = e*k-fc*h;
  double B = fc*g-d*k;
  double C = d*h-e*g;
  double detJ = a*(A) + b*(B) + c*(C);

  //Compute the matrix inverse
  double Ji[3][3];
  Ji[0][0] = A/detJ;
  Ji[1][0] = B/detJ;
  Ji[2][0] = C/detJ;
  Ji[0][1] = (c*h-b*k)/detJ;
  Ji[1][1] = (a*k-c*g)/detJ;
  Ji[2][1] = (g*b-a*h)/detJ;
  Ji[0][2] = (b*fc-c*e)/detJ;
  Ji[1][2] = (c*d-a*fc)/detJ;
  Ji[2][2] = (a*e-b*d)/detJ;

  // Compute the step size
  dx[0] = Ji[0][0]*f[0] + Ji[0][1]*f[1] + Ji[0][2]*f[2];
  dx[1] = Ji[1][0]*f[0] + Ji[1][1]*f[1] + Ji[1][2]*f[2];
  dx[2] = Ji[2][0]*f[0] + Ji[2][1]*f[1] + Ji[2][2]*f[2];
  
//  printf("dxf: %15e,%15e,%15e,%15e\n",dx[0],dx[1],dx[2], detJ);
//  printf("Ji: %15e,%15e,%15e\n",Ji[2][0],Ji[2][1],Ji[2][2]);
}




void calc_unknowns_from_prim(int EOS, double * EOS_const, double * prim, double * con, double * x, double g_con[4][4], struct nuc_eos_vars *struct_ptr){
   double v1_con = g_con[1][1]*prim[v1_cov]+g_con[1][2]*prim[v2_cov]+g_con[1][3]*prim[v3_cov];
   double v2_con = g_con[1][2]*prim[v1_cov]+g_con[2][2]*prim[v2_cov]+g_con[2][3]*prim[v3_cov];
   double v3_con = g_con[1][3]*prim[v1_cov]+g_con[2][3]*prim[v2_cov]+g_con[3][3]*prim[v3_cov];

  //v^2 = v^i * v_i
  double v_squared = prim[v1_cov]*v1_con + prim[v2_cov]*v2_con + prim[v3_cov]*v3_con;
//  double v_squared = prim[v1_cov]*prim[v1_cov] + prim[v2_cov]*prim[v2_cov] + prim[v3_cov]*prim[v3_cov];
  double W;
  if(v_squared > 1.0)
    W = 10000.0;
  else
    W = 1.0/sqrt(1.0-v_squared);
  
  x[0] = W;
  
  //Always calculate rho from D and W so that using D in EOS remains consistent
  prim[RHO] = con[D]/W;
   
  //DON'T HAVE EPSILON YET!!!!!!!
  
  //NOTE: done once, so prim can be reset
  
  double P = P_prim_(EOS,EOS_const,prim,struct_ptr);
//  x[2] = P;
  x[2] = prim[TEMP];

  double h;
  if(EOS==1){
    h = 1.0 + P/prim[RHO];
  }
  else{
    h = 1.0 + prim[EPS]+P/prim[RHO];
  }
  double Z = prim[RHO]*h*W*W;
  x[1] = Z;
//  printf("x = %15e,%15e,%15e\n",x[0],x[1],x[2]);
}

// EOS INDEPENDENT
int NR_3D(int EOS, double * EOS_const, double S_squared, 
          double BdotS, double B_squared, double * con, double * prim, 
	  double * x, double g_con[4][4], double B1_cov, double B2_cov, double B3_cov, 
          int CODETEST, struct nuc_eos_vars *struct_ptr, int stepsize){

  double dx[3];
  double f[3];
  double T1,Ecalc1,Ecalc;
  Ecalc1 = prim[EPS];
  T1 = prim[TEMP];
  int i;
  double TOL = 1.0e-10;
  if (stepsize == 5) TOL = 1.0e-8;
  int done = 0;
  double error[3];
  int count = 0; // count number of iterations
  int extra = 0;
  double P_guess = 0.0;
  double Z_guess = 0.0;
  double W_guess = 0.0;
  int maxcount_soft = 128;
  int maxcount_hard = 128;

  while(!done){
    count ++;
    
    // If codetest then start with safeguess values
    if(count==1){
      if(CODETEST==1){
        // Eq (41): W_guess
        W_guess = 10000.0;
  	x[0] = W_guess;
      
        // Eq (35): rho_max = D
        prim[RHO] = con[D];

        // Eq (37): eps_max = 1/D*[tau-B^2/2]
        if(EOS!=1)
	  prim[EPS] = 1.0/con[D] * (con[TAU] - B_squared/2.0);
      
        // NEED TO COMPUTE TMAX?????????????????????????
      
        // Eq (34): P_max = P(rho_max, eps_max)
        P_guess = T_prim_(EOS,EOS_const,prim,struct_ptr);
//        x[2] = P_guess;
	x[2] = prim[TEMP];

        // Eq (38): z_max = tau + P_max + D - B^2/2
        Z_guess = con[TAU] + P_guess + con[D] - B_squared/2.0;
        x[1] = Z_guess;
      }
      else if(CODETEST==0){
        calc_unknowns_from_prim(EOS,EOS_const,prim,con,x,g_con,struct_ptr);
      }
    }
//    printf("Count: %d\n",count);
    // Take step

    NR_step_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, dx, f, T1, Ecalc1, &Ecalc, B1_cov, B2_cov, B3_cov, struct_ptr, stepsize);

//    if (count % 5 == 0) {
//      printf("Reducing stepsize!\n");
//     dx[0] = 0.72*dx[0]; 
//     dx[1] = 0.822*dx[1]; 
//     dx[2] = 0.05*dx[2]; 
//    }

      // Save temperature and internal energy values from last iteration in case we switch to secant method
      T1 = x[2];
      Ecalc1 = Ecalc;

    // Update values
//    x[0] = x[0] + dx[0];
//    x[1] = x[1] + dx[1];
//    x[2] = x[2] + dx[2];
    x[0] = fabs(x[0] - dx[0]);
    x[1] = fabs(x[1] - dx[1]);
    x[2] = fabs(x[2] - dx[2]);
//    x[2] = x[2] - dx[2]

    // Calculate error
    double maxerror = 0.0;
    for(i=0;i<3;i++){
      error[i] = fabs(dx[i]/x[i]);
      if(i==0)
	maxerror = error[i];
      else if(error[i]>maxerror)
	maxerror = error[i];
    }  
  
    if (count > 127 && maxerror > 1.0e-10) {
      printf("Count, maxerror, error: %d, %15e, %15e,%15e,%15e \n",count, maxerror, error[0], error[1], error[2]);
      printf("x,dx2: %15e,%15e,%15e,%15e,%15e,%15e\n",x[0],x[1],x[2],dx[0],dx[1],dx[2]);
    }
    // if successful
    if(maxerror<TOL){	
      if(extra==0)
        done = 1;
      else
	extra--;
    }
    else if(count>=maxcount_soft){ // if not successful
//      calc_prim(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, B1_cov, B2_cov, B3_cov, struct_ptr);
      if(count>=maxcount_hard){ // if too many tries
	done = 1;
	break;
      }
//      if(prim[EPS]<0.0){  // if it's a solvable problem
//	prim[EPS] = fabs(prim[EPS]);
//	calc_unknowns_from_prim(EOS,EOS_const,prim,con,x,g_con,struct_ptr);
//
//      }
      
    }
  }
  return count;
}

double function(double W,double B_squared, double S_squared, double BdotS, int EOS, double * EOS_const,double * con, double * prim){
  
  double rho = con[D]/W;
  double gammapoly = EOS_const[0];
  double Kpoly = EOS_const[1];
  double P = Kpoly*pow(rho,gammapoly);
  double h = 1.0 + P/rho;
  double Z = rho*h*W*W;
  
  return (((Z+B_squared)*(Z+B_squared)) - S_squared - (2*Z+B_squared)/(Z*Z)*(BdotS*BdotS))*W*W - ((Z+B_squared)*(Z+B_squared));
}

int bisectionmethod(int EOS, double * EOS_const, double S_squared, double BdotS, 
                    double B_squared, double * con, double * prim, 
		    double * x, double B1_cov, double B2_cov, double B3_cov, 
                    int CODETEST, struct nuc_eos_vars *struct_ptr){
  double P,h,Z,W;
  
  double Wmin = 1.0 + 1.0e-5;
  double Wmax = 1.0 + 1.0e5;
  double Wmid = 0.5*(Wmin+Wmax);
  int Nmax = 64;
  int N = 0;
  
  W = 0.0;
  
  while(N<Nmax){
    
    Wmid = 0.5*(Wmin+Wmax);
    
    if(fabs(function(Wmid,B_squared,S_squared,BdotS,EOS,EOS_const,con,prim)) < 1.0e-22 || 0.5*(Wmax-Wmin)<1.0e-13)
      break;
        
    N = N + 1;
    
    if(function(Wmin,B_squared,S_squared,BdotS,EOS,EOS_const,con,prim)*function(Wmid,B_squared,S_squared,BdotS,EOS,EOS_const,con,prim)<0.0)
      Wmax = Wmid;
    else
      Wmin = Wmid;
  
  }
  W = Wmid;
  
  prim[RHO] = con[D]/W;
  P = P_prim_(EOS,EOS_const,prim,struct_ptr);
  h = 1.0 + P/prim[RHO];
  Z = prim[RHO]*h*W*W;
  
  x[0] = W;
  x[1] = Z;
  x[2] = P;
  
  return N;
}
 
