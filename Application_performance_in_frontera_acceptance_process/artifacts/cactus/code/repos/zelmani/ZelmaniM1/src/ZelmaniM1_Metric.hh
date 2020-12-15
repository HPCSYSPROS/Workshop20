#ifndef _THREEMETRICCLASS_
#define _THREEMETRICCLASS_

#include <math.h>
using namespace std;


// This is a storage capsule for local properties of the spacetime that are useful for radiative transfer
class ThreeMetric {

protected:	
	// Input quantities
	double detg; 
	 
public:
	double alpha,betax,betay,betaz; //Lapse and shift (raised index)
	double gxx,gxy,gxz,gyy,gyz,gzz; //gamma_{ij}
	double guxx,guxy,guxz,guyy,guyz,guzz; //gamma^{ij}
        double Kxx,Kxy,Kxz,Kyy,Kyz,Kzz; //K_{ij}
		
	// Load the three metric
	void loadG(double gxxi, double gxyi, double gxzi,
		   double gyyi, double gyzi, double gzzi) {
		gxx = gxxi; gxy = gxyi; gxz = gxzi; 
		gyy = gyyi; gyz = gyzi; gzz = gzzi;
		
		detg = (gxx*gyy*gzz + 2.0*gxy*gxz*gyz 
		  - gyy*gxz*gxz - gxx*gyz*gyz - gzz*gxy*gxy);
		
		double idenom = 1.0/detg;
		
		guxx = (gyy*gzz - gyz*gyz)*idenom;
		guxy = (gxz*gyz - gxy*gzz)*idenom;
		guxz = (gxy*gyz - gxz*gyy)*idenom;
		guyy = (gxx*gzz - gxz*gxz)*idenom;
		guyz = (gxy*gxz - gxx*gyz)*idenom;
		guzz = (gxx*gyy - gxy*gxy)*idenom;
	
	}
 
	void loadLapseShift(double alphai, double betaxi){
		alpha = alphai; betax = betaxi;
	}
        
	void loadLapseShiftAll(double alphai, double betaxi, 
	                       double betayi, double betazi){
		alpha = alphai; betax = betaxi;
                betay = betayi; betaz = betazi;	
	}
	void loadExtCurv(double Kxxi, double Kxyi, double Kxzi, 
	                 double Kyyi, double Kyzi, double Kzzi) {
	 	Kxx = Kxxi; Kxy = Kxyi; Kxz = Kxzi;
		Kyy = Kyyi; Kyz = Kyzi;
		Kzz = Kzzi;
		return;
	}
        
	double contractKllTuuSym(double Txx, double Txy, double Txz,
	                         double Tyy, double Tyz, double Tzz) {
		return        Txx*Kxx + Tyy*Kyy + Tzz*Kzz 
		       + 2.0*(Txy*Kxy + Txz*Kxz + Tyz*Kyz); 
	}

	double getDetg()  { return detg;  } 
	double getAlpha() { return alpha; } 
	double getBetax() { return betax; } 
	double getBetay() { return betay; } 
	double getBetaz() { return betaz; } 

	// Contract the three metric with three vectors 
        double contractAlBu(double Ax, double Ay, double Az,
			    double Bx, double By, double Bz) {
        	return Ax*Bx + Ay*By + Az*Bz;
	}
	double contractGllAuBu(double Ax, double Ay, double Az,
			       double Bx, double By, double Bz) {
		return  gxx*Ax*Bx + gxy*Ax*By + gxz*Ax*Bz 
		      + gxy*Ay*Bx + gyy*Ay*By + gyz*Ay*Bz 
		      + gxz*Az*Bx + gyz*Az*By + gzz*Az*Bz;
	}

        double contractGuuAlBl(double Ax, double Ay, double Az,
			       double Bx, double By, double Bz) {
		return  guxx*Ax*Bx + guxy*Ax*By + guxz*Ax*Bz 
		      + guxy*Ay*Bx + guyy*Ay*By + guyz*Ay*Bz 
		      + guxz*Az*Bx + guyz*Az*By + guzz*Az*Bz;
	}

        // Raise and lower three vector indices	
	void lowerAu(double Axu, double Ayu, double Azu,
		     double *Axl, double *Ayl, double *Azl) {
		*Axl = gxx*Axu + gxy*Ayu + gxz*Azu;
		*Ayl = gxy*Axu + gyy*Ayu + gyz*Azu;
		*Azl = gxz*Axu + gyz*Ayu + gzz*Azu;
		
		return;
	}
	
	void raiseSAll(double  Axxl, double  Axyl, double  Axzl, 
	               double  Ayyl, double  Ayzl, double  Azzl,
                       double *Axxu, double *Axyu, double *Axzu, 
		       double *Ayyu, double *Ayzu, double *Azzu) {
		
		*Axxu =   guxx*guxx*Axxl + 2.0*guxx*(guxy*Axyl + guxz*Axzl)
		        + guxy*guxy*Ayyl + 2.0*guxy* guxz*Ayzl
		        + guxz*guxz*Azzl; 

		*Axyu =   guxy*guxx*Axxl + 2.0*guxy*(guxy*Axyl + guxz*Axzl)
		        + guyy*guxy*Ayyl + 2.0*guyy* guxz*Ayzl
		        + guyz*guxz*Azzl; 
		
		*Axzu =   guxz*guxx*Axxl + 2.0*guxz*(guxy*Axyl + guxz*Axzl)
		        + guxz*guxy*Ayyl + 2.0*guyz* guxz*Ayzl
		        + guzz*guxz*Azzl; 

		*Ayyu =   guxy*guxy*Axxl + 2.0*guxy*(guyy*Axyl + guyz*Axzl)
		        + guyy*guyy*Ayyl + 2.0*guyy* guyz*Ayzl
		        + guyz*guyz*Azzl; 

		*Ayzu =   guxz*guxy*Axxl + 2.0*guxz*(guyy*Axyl + guyz*Axzl)
		        + guyz*guyy*Ayyl + 2.0*guyz* guyz*Ayzl
		        + guzz*guyz*Azzl; 

		*Azzu =   guxz*guxz*Axxl + 2.0*guxz*(guyz*Axyl + guzz*Axzl)
		        + guyz*guyz*Ayyl + 2.0*guyz* guzz*Ayzl
			+ guzz*guzz*Azzl; 

	        return;
	}
	
	void lowerSAuu(double  Axxu, double  Axyu, double  Axzu, 
	               double  Ayyu, double  Ayzu, double  Azzu,
                       double *Axxl, double *Axyl, double *Axzl, 
		       double *Ayyl, double *Ayzl, double *Azzl) {
		
		*Axxl =   gxx*gxx*Axxu + 2.0*gxx*(gxy*Axyu + gxz*Axzu)
		        + gxy*gxy*Ayyu + 2.0*gxy*gxz*Ayzu
			+ gxz*gxz*Azzu; 

		*Axyl =   gxy*gxx*Axxu + 2.0*gxy*(gxy*Axyu + gxz*Axzu)
		        + gyy*gxy*Ayyu + 2.0*gyy*gxz*Ayzu
			+ gyz*gxz*Azzu; 
		
		*Axzl =   gxz*gxx*Axxu + 2.0*gxz*(gxy*Axyu + gxz*Axzu)
		        + gxz*gxy*Ayyu + 2.0*gyz*gxz*Ayzu
			+ gzz*gxz*Azzu; 

		*Ayyl =   gxy*gxy*Axxu + 2.0*gxy*(gyy*Axyu + gyz*Axzu)
		        + gyy*gyy*Ayyu + 2.0*gyy*gyz*Ayzu
			+ gyz*gyz*Azzu; 

		*Ayzl =   gxz*gxy*Axxu + 2.0*gxz*(gyy*Axyu + gyz*Axzu)
		        + gyz*gyy*Ayyu + 2.0*gyz*gyz*Ayzu
			+ gzz*gyz*Azzu; 

		*Azzl =   gxz*gxz*Axxu + 2.0*gxz*(gyz*Axyu + gzz*Axzu)
		        + gyz*gyz*Ayyu + 2.0*gyz*gzz*Ayzu
			+ gzz*gzz*Azzu; 

	        return;
	}
	
	void raiseAl(double Axl, double Ayl, double Azl,
		     double *Axu, double *Ayu, double *Azu) {
		*Axu = guxx*Axl + guxy*Ayl + guxz*Azl;
		*Ayu = guxy*Axl + guyy*Ayl + guyz*Azl;
		*Azu = guxz*Axl + guyz*Ayl + guzz*Azl;
		
		return;
	}
};

#endif
