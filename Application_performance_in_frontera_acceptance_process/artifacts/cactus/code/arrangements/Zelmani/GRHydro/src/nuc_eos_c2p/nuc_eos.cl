// CalcEulate P from primitive variables
double P_iter_(int EOS, double * EOS_const, double * prim, double * con, double * x,struct nuc_eos_vars *struct_ptr){
  if(EOS==3){
    double W = x[0];
    double T = x[2];
    // outvars
    double xeps,xprs,xprs2;
    const double prec = 1.0e-10;
    int keytemp,keyerr,anyerr;
    double xrho = (con[D]/W)/rho_gf;
    double xrho2 = (con[D]/W);
    double xtemp = T;
    double xye = prim[YE];
    xeps = prim[EPS]/eps_gf;
    xprs = 0.0;
    keytemp = 1;

    nuc_eosK_P_only(xrho,&xtemp,xye,&xeps,&xprs,
		    keytemp,
		    &keyerr,prec,struct_ptr);

//    EOS_Omni_press(4,keytemp,prec,1,&xrho2,&xeps,&xtemp,&xye,&xprs2,&keyerr,&anyerr);
//    xprs2 = xprs2/press_gf;

//    if (xtemp > 1.5) {

//    if (anyerr != 0) printf("Anyerr1 %d Keyerr1: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",anyerr,keyerr,xtemp,xrho,xeps,xye);
//    if (fabs(xprs-xprs2) > 1.0e-10) {
//      printf("Press mismatch at xprs = %g, xprs2 = %g\n", xprs, xprs2);
//    }
//    }
//    prim[TEMP] = xtemp;

    if (keyerr != 0) printf("Keyerr2: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
//    printf("Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
    if(xprs<1.0e-10){

    }
    return xprs * press_gf;
//    return xprs2;
  }
  else{

  }
  
  return 0.0;
}

// Calculate P,T from primitive variables
double T_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr){

    // outvars
    double xeps,xprs;
    const double prec = 1.0e-10;
    int keytemp,keyerr;
    double xrho = prim[RHO]/rho_gf;
    double xtemp = prim[TEMP];
    double xye = prim[YE];
    xeps = prim[EPS]/eps_gf;
    xprs = 0.0;
    keytemp = 0;

    nuc_eosK_P_only(xrho,&xtemp,xye,&xeps,&xprs,
		    keytemp,
		    &keyerr,prec,struct_ptr);
    prim[TEMP] = xtemp;

    if (keyerr != 0) printf("Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
    printf("Keyerr2: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
    if(xprs<1.0e-10){

    }
    return xprs * press_gf;
}

// Calculate P from primitive variables
double P_prim_(int EOS, double * EOS_const, double * prim,struct nuc_eos_vars *struct_ptr){

  if(EOS==0){ // IDEAL GAS
    double gammaideal = EOS_const[0];
    return prim[RHO]*prim[EPS]*(gammaideal-1.0);
  }
  else if(EOS==1){ // POLYTROPIC
    double gammapoly = EOS_const[0];
    double Kpoly = EOS_const[1];
    return Kpoly*pow(prim[RHO],gammapoly);
  }
  // Hybrid EOS
  else if(EOS==2){

    double gamma1 = EOS_const[0];
    double K1 = EOS_const[1];
    double gamma2 = EOS_const[2];
    double gammath = EOS_const[3];
    double rho_nuc = EOS_const[4];
    
    // Equations 8-9 from Janka et al. 1993
    double E1 = K1/(gamma1-1.0);
    double E2 = (gamma1-1.0)/(gamma2-1.0)*E1*pow(rho_nuc,(gamma1-gamma2)); 
    double K2 = (gamma2-1.0)*E2;
    double E3 = (gamma2-gamma1)/(gamma2-1.0)*E1*pow(rho_nuc,(gamma1-1.0));
    
    double Kx, Ex, Gx, Ex3;
    double eps_p, eps_th,p_th,p_poly;
    
    if(prim[RHO]<rho_nuc){
      Kx = K1;
      Ex = E1;
      Gx = gamma1;
      Ex3 = 0.0;
    }
    else{
      Kx = K2;
      Ex = E2;
      Gx = gamma2;
      Ex3 = E3;
    }
    
    // Equation (7) from Janka et al. 1993
    eps_p = Ex*pow(prim[RHO],Gx) + Ex3*prim[RHO];
    
    eps_th = prim[EPS]*prim[RHO] - eps_p;
    
    // Equation (11) from Janka et al. 1993
    p_th = (gammath-1.0)*eps_th;
    
    p_th = p_th > 0.0 ? p_th : 0.0;
   
    p_poly = Kx*pow(prim[RHO],Gx);
    
    return p_poly + p_th;
  }
  else if(EOS==3){

    // outvars
    double xeps,xprs;
    const double prec = 1.0e-10;
    int keytemp,keyerr;
    double xrho = prim[RHO]/rho_gf;
    double xtemp = prim[TEMP];
    double xye = prim[YE];
    xeps = prim[EPS]/eps_gf;
    xprs = 0.0;
    keytemp = 1;

    nuc_eosK_P_only(xrho,&xtemp,xye,&xeps,&xprs,
		    keytemp,
		    &keyerr,prec,struct_ptr);
//    prim[TEMP] = xtemp;

    if (keyerr != 0) printf("Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
//    printf("Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",keyerr,xtemp,xrho,xeps,xye);
    if(xprs<1.0e-10){

    }
    return xprs * press_gf;
  }
  else{

  }
  
  return 0.0;
}

void P_dpdr_dpde(int EOS, double * EOS_const, double * x, double * con, double * prim, double * Pprim, double * dPdrho,double * dPde,struct nuc_eos_vars *struct_ptr){
  double W = x[0];
  double Z = x[1];
  double P = x[2];
  
  // Calculate rho, eps from (23), (25) of Cerdá-Durán
  double rho = con[D]/W;
  double eps = -1.0 + Z/(con[D]*W) - P*W/con[D];
  // DOESN'T MAKE SENSE TO RECALCULATE - SHOULD JUST USE PRIMS, OR NEVER CALCULATE PRIMS AT ALL!!!!!!!!!
  
  if(EOS==0){
    double gammaideal = EOS_const[0];
    *Pprim = rho*eps*(gammaideal-1.0);
    *dPdrho = (gammaideal-1.0)*eps;
    *dPde = (gammaideal-1.0)*rho;
  }
  // Polytropic EOS
  else if(EOS==1){
    double gammapoly = EOS_const[0];
    double Kpoly = EOS_const[1];
    *Pprim = Kpoly*pow(rho,gammapoly);
    *dPdrho = Kpoly*gammapoly*pow(rho,gammapoly-1.0);
    *dPde = 0.0;
  }
  // Hybrid EOS
  else if(EOS==2){
    double gamma1 = EOS_const[0];
    double K1 = EOS_const[1];
    double gamma2 = EOS_const[2];
    double gammath = EOS_const[3];
    double rho_nuc = EOS_const[4];
    
    double E1 = K1/(gamma1-1.0);
    double E2 = (gamma1-1.0)/(gamma2-1.0)*E1*pow(rho_nuc,(gamma1-gamma2));
    double K2 = (gamma2-1.0)*E2;
    double E3 = (gamma2-gamma1)/(gamma2-1.0)*E1*pow(rho_nuc,gamma1-1.0);
    
    double Kx, Ex, Gx, Ex3;
    double eps_p, eps_th,p_th,p_poly;
    
    if(rho<rho_nuc){
      Kx = K1;
      Ex = E1;
      Gx = gamma1;
      Ex3 = 0.0;
    }
    else{
      Kx = K2;
      Ex = E2;
      Gx = gamma2;
      Ex3 = E3;
    }
    
    //p_th = (gammath-1.0)*eps_th;
    double dp_th_drho = (gammath-1.0)*(eps - Ex*Gx*pow(rho,(Gx-1.0)) - Ex3);
   
    //p_poly = Kx*pow(rho,Gx);
    double dp_poly_drho = Gx*Kx*pow(rho,(Gx-1.0));
    *dPdrho = dp_poly_drho + dp_th_drho;
    
    
    double dp_th_deps = (gammath-1.0)*rho;
    double dp_poly_deps = 0.0;
   
    *dPde = dp_poly_deps + dp_th_deps;
    
    
    
    
    // Equation (7) from Janka et al. 1993
    eps_p = Ex*pow(rho,Gx) + Ex3*rho;
    
    eps_th = eps*rho - eps_p;
    
    // Equation (11) from Janka et al. 1993
    p_th = (gammath-1.0)*eps_th;
    p_th = p_th > 0.0 ? p_th : 0.0;
    p_poly = Kx*pow(rho,Gx);
    *Pprim = p_poly + p_th;
  }
  else if(EOS==3){

    // outvars
    double xeps,xprs,xdpderho,xdpdrhoe;
    const double prec = 1.0e-10;
    int keytemp,keyerr;
    double xrho = rho/rho_gf;
    double xtemp = prim[TEMP];
    double xye = prim[YE];
//    double xtemp = con[TEMP];
//    double xye = con[YE];
    xeps = eps/eps_gf;
    
    xprs = 0.0;
    xdpdrhoe = 0.0;
    xdpderho = 0.0;
    keytemp = 0;
    
    nuc_eosK_P_dpdr_dpde(xrho,&xtemp,xye,&xeps,&xprs,
		    &xdpderho,&xdpdrhoe,keytemp,
		    &keyerr,prec,struct_ptr);
    
//     con[TEMP] = xtemp;
    *dPdrho = xdpdrhoe * press_gf/rho_gf;
    *dPde = xdpderho * press_gf/eps_gf;
    *Pprim = xprs * press_gf;
    
  }
  else{

  }
  
  
  
}

void P_dPdW_dPdZ(double * Pprim, double * dPdW, double * dPdZ, int EOS, double * EOS_const, double * x, double * con, double * prim, struct nuc_eos_vars *struct_ptr){
  double W = x[0];
  double Z = x[1];
  double P = x[2];

  // Partial derivatives
  double drhodW = -con[D]/(W*W);
  double depsdW = -Z/(con[D]*W*W)-P/con[D];
  double drhodZ = 0.0;
  double depsdZ = 1.0/(con[D]*W);
  
  double xP,xdPdrho,xdPde;
  xP=0.0;
  xdPdrho=0.0;
  xdPde = 0.0;
  P_dpdr_dpde(EOS,EOS_const,x,con,prim,&xP,&xdPdrho,&xdPde,struct_ptr);
  *Pprim = xP;
  *dPdW = xdPdrho*drhodW + xdPde*depsdW;
  *dPdZ = xdPdrho*drhodZ + xdPde*depsdZ;
}

void nuc_eosK_dpdrho_dpdt_dedrho_dedt(double xrho, double *xtemp, double xye,
		     double *xenr, double *xenr2, double *xprs, double *xdedt2,
		     double *dlepsdlrho, double *dlepsdlt, double *dlPdlrho, double *dlPdlrho2, double *dlPdlt, int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr, int stepsize) 
{
  // This routine expects rho in g/cm^3 and T in MeV.
  // It will strictly return values in cgs or MeV

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)



  *keyerr = 0;

  if(xrho > struct_ptr->eos_rhomax) {
    *keyerr = 105;
    return;
  }
  if(xrho < struct_ptr->eos_rhomin) {
    *keyerr = 106;
    return;
  }
  if(xye > struct_ptr->eos_yemax) {
    *keyerr = 101;
    return;
  }
  if(xye < struct_ptr->eos_yemin) {
    *keyerr = 102;
    return;
  }
  if(keytemp == 1) {
    if(*xtemp > struct_ptr->eos_tempmax) {
      *keyerr = 103;
      return;
    }
    if(*xtemp < struct_ptr->eos_tempmin) {
      *keyerr = 104;
      return;
    }
  }

  // set up local vars
  double lr = log10(xrho);
  double lt = log10(*xtemp);
  double xeps = *xenr + struct_ptr->energy_shift;
  double leps = log10(MAX(xeps,1.0));



  // find temperature if need be
  if(keytemp == 0) {
    double nlt = 0.0;
    nuc_eosK_C_findtemp(lr,lt,xye,leps,rfeps,&nlt,keyerr,struct_ptr);

    if(*keyerr != 0) return;
    lt = nlt;
    *xtemp = pow(10.0,lt);
  } else if(keytemp == 2) {

  }

  double res[4];

  if (stepsize ==1) {
    nuc_eosK_C_linterp_some2(lr,lt,xye,res,struct_ptr->fourtables,struct_ptr->ivs_short,
  			   struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,4,struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,struct_ptr,
                           dlepsdlrho, dlepsdlt, dlPdlrho, dlPdlt);
  }
  else if (stepsize ==2) {
    nuc_eosK_C_linterp_some3(lr,lt,xye,res,struct_ptr->fourtables,struct_ptr->ivs_short,
  			   struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,4,struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,struct_ptr,
                           dlepsdlrho, dlepsdlt, dlPdlrho, dlPdlt);
  }
  else {
    nuc_eosK_C_linterp_some(lr,lt,xye,res,struct_ptr->fourtables,struct_ptr->ivs_short,
  			   struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,4,struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,struct_ptr);
  } 
  // assign results

  if(keytemp != 0) {
    // set up internal energy correctly, correcting for shift
    *xenr = pow(10.0,res[1]) - struct_ptr->energy_shift;
    *xenr2 = pow(10.0,res[1]);
  }

  *xprs     = pow(10.0,res[0]);
  *xdedt2 = res[2];
  *dlPdlrho2 = res[3];
//  *xdpderho = res[3];
  return;
}

void E_dEdr_dEdt_dPdr_dPdt(int EOS, double * EOS_const, double * x, double * con, double * prim, double * Eprim, double * dEdrho, double * dEdt, double * dPdrho, double * dPdt, struct nuc_eos_vars *struct_ptr, int stepsize){
  double W = x[0];
  double Z = x[1];
  double T = x[2];
  
  // Calculate rho, eps from (23), (25) of Cerdá-Durán
  double P = P_iter_(EOS,EOS_const,prim,con,x,struct_ptr);
  double rho = con[D]/W;
  double eps = -1.0 + Z/(con[D]*W) - P*W/con[D];
  // DOESN'T MAKE SENSE TO RECALCULATE - SHOULD JUST USE PRIMS, OR NEVER CALCULATE PRIMS AT ALL!!!!!!!!!
 
    // outvars
    double xeps,xprs,xdedrho,xdpdrho,xdedt,xdpdt;
    const double prec = 1.0e-10;
    int keytemp,keyerr;
    double xrho = rho/rho_gf;
    double xtemp = T;
    double xye = prim[YE];
    xeps = eps/eps_gf;
    double xeps2 = 0.0;
    
    xprs = 0.0;
    xdpdrho = 0.0;
    double xdpdrho2 = 0.0;
    xdedt = 0.0;
    double xdedt2 = 0.0;
    xdpdt = 0.0;
    xdedrho = 0.0;
    keytemp = 1;
    
    nuc_eosK_dpdrho_dpdt_dedrho_dedt(xrho,&xtemp,xye,&xeps,&xeps2,&xprs,&xdedt2,
		    &xdedrho,&xdedt,&xdpdrho,&xdpdrho2,&xdpdt,keytemp,
		    &keyerr,prec,struct_ptr,stepsize);
    
    *dEdrho = xeps2/xrho*xdedrho * eps_gf/rho_gf;
    *dEdt = xeps2/xtemp*xdedt*eps_gf;
    *dPdrho = xprs/xrho*xdpdrho*press_gf/rho_gf;
    *dPdt = xprs/xtemp*xdpdt*press_gf;
    *Eprim = xeps * eps_gf;
    if (stepsize == 3) {
//      *dEdt = xdedt2*eps_gf;
      *dEdt = xeps/xtemp*xdedt2*eps_gf;
    }
//     con[TEMP] = xtemp;
//    *dEdt = xdedt * eps_gf;
//    *dPdt = xdpdt * press_gf;
//    *dPdrho = xdpdrho2*press_gf/rho_gf;
//    printf("xeps,xrho,xprs,xtemp: %g,%g,%g,%g\n",xeps,xrho,xprs,xtemp);
//    printf("xdedrho,xdpxrho,xdedt,xdpdt,xdedt2: %g,%g,%g,%g,%g\n",xdedrho,xdpdrho,xdedt,xdpdt,xdedt2);
//    printf("dEdt2,dEdt,dEdt3,dPdrho,dPdrho2: %g,%g,%g,%g,%g\n",*dEdt,(xeps/xtemp*xdedt*eps_gf),(xeps2/xtemp*xdedt*eps_gf),*dPdrho,(xdpdrho2*press_gf/rho_gf));
}

void E_dEdW_dEdZ(double * Eprim, double * dEdW, double * dEdZ, double * dEdT, double * dPdT, int EOS, double * EOS_const, double * x, double * con, double * prim, struct nuc_eos_vars *struct_ptr, int stepsize){
  double W = x[0];
  double Z = x[1];
  double T = x[2];
  
  double P,xP,E,xE,xdEdrho,xdPdrho,xdEdt,xdPdt,dEdP,dPde,dEdt;
  P=0.0;
  xP=0.0;
  xdPdrho=0.0;

  P = P_iter_(EOS,EOS_const,prim,con,x,struct_ptr);
  
  E_dEdr_dEdt_dPdr_dPdt(EOS,EOS_const,x,con,prim,&xE,&xdEdrho,&xdEdt,&xdPdrho,&xdPdt,struct_ptr,stepsize);

// Partial derivatives
  double drhodW = -con[D]/(W*W);
  double depsdW = -Z/(con[D]*W*W)-P/con[D] + xdPdrho/W;
  double depsdP = -W/con[D];
  double drhodZ = 0.0;
  double depsdZ = 1.0/(con[D]*W);

  *Eprim = xE;
  *dEdW = depsdW - xdEdrho*drhodW;
  // *dEdZ = depsdZ - xdEdrho*drhodZ;
  *dEdZ = depsdZ;
  *dEdT = depsdP*xdPdt - xdEdt;
  *dPdT = xdPdt;
}



void nuc_eosK_P_dpdr_dpde(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     double* xdpderho,
		     double *xdpdrhoe, int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr) 

{
  // This routine expects rho in g/cm^3 and T in MeV.
  // It will strictly return values in cgs or MeV

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)



  *keyerr = 0;

  if(xrho > struct_ptr->eos_rhomax) {
    *keyerr = 105;
    return;
  }
  if(xrho < struct_ptr->eos_rhomin) {
    *keyerr = 106;
    return;
  }
  if(xye > struct_ptr->eos_yemax) {
    *keyerr = 101;
    return;
  }
  if(xye < struct_ptr->eos_yemin) {
    *keyerr = 102;
    return;
  }
  if(keytemp == 1) {
    if(*xtemp > struct_ptr->eos_tempmax) {
      *keyerr = 103;
      return;
    }
    if(*xtemp < struct_ptr->eos_tempmin) {
      *keyerr = 104;
      return;
    }
  }

  // set up local vars
  double lr = log10(xrho);
  double lt = log10(*xtemp);
  double xeps = *xenr + struct_ptr->energy_shift;
  double leps = log10(MAX(xeps,1.0));



  // find temperature if need be
  if(keytemp == 0) {
    double nlt = 0.0;
    nuc_eosK_C_findtemp(lr,lt,xye,leps,rfeps,&nlt,keyerr,struct_ptr);

    if(*keyerr != 0) return;
    lt = nlt;
    *xtemp = pow(10.0,lt);
  } else if(keytemp == 2) {

  }

  double res[4];
  nuc_eosK_C_linterp_some(lr,lt,xye,res,struct_ptr->fourtables,struct_ptr->ivs_short,
			 struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,4,struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,struct_ptr);
  
  // assign results

  if(keytemp != 0) {
    // set up internal energy correctly, correcting for shift
    *xenr = pow(10.0,res[1]) - struct_ptr->energy_shift;
  }

  *xprs     = pow(10.0,res[0]);
  *xdpdrhoe = res[2];
  *xdpderho = res[3];
  return;
}

void nuc_eosK_P_only(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs,
		     int keytemp,
		     int *keyerr,double rfeps,struct nuc_eos_vars *struct_ptr) 

{

  // This routine expects rho in g/cm^3 and T in MeV.
  // It will strictly return values in cgs or MeV

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)

  *keyerr = 0;

  if(xrho > struct_ptr->eos_rhomax) {
    *keyerr = 105;
    return;
  }
  if(xrho < struct_ptr->eos_rhomin) {
    *keyerr = 106;
    return;
  }
  if(xye > struct_ptr->eos_yemax) {
    *keyerr = 101;
    return;
  }
  if(xye < struct_ptr->eos_yemin) {
    *keyerr = 102;
    return;
  }
  if(keytemp == 1) {
    if(*xtemp > struct_ptr->eos_tempmax) {
      *keyerr = 103;
      return;
    }
    if(*xtemp < struct_ptr->eos_tempmin) {
      *keyerr = 104;
      return;
    }
  }

  // set up local vars
  double lr = log10(xrho);
  double lt = log10(*xtemp);
  double xeps = *xenr + struct_ptr->energy_shift;
  double leps = log10(MAX(xeps,1.0));



  // find temperature if need be
  if(keytemp == 0) {
    double nlt = 0.0;

    nuc_eosK_C_findtemp(lr,lt,xye,leps,rfeps,&nlt,keyerr,struct_ptr);

    if(*keyerr != 0) return;
    lt = nlt;

    *xtemp = pow(10.0,lt);
  } else if(keytemp == 1) {
    
  }

  double res[4];
  nuc_eosK_C_linterp_some(lr,lt,xye,res,struct_ptr->fourtables,struct_ptr->ivs_short,
			 struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,4,struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,struct_ptr);
  
  // assign results

  if(keytemp != 0) {
    // set up internal energy correctly, correcting for shift
    *xenr = pow(10.0,res[1]) - struct_ptr->energy_shift;
  }

  *xprs     = pow(10.0,res[0]);
  return;
}


