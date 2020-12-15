#ifndef PARM_H
#define PARM_H

#include <stdio.h>
 
//#define DOUBLE

#ifdef DOUBLE
#define _REAL		double
#else
#define _REAL		float
#endif


typedef struct parm {
	char	ititl[81];
	int 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap;
	char 	*AtomNames, *ResNames, *AtomSym, *AtomTree;
	_REAL	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB6;
	_REAL	Box[3], Cutcap, Xcap, Ycap, Zcap;
	int 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		*AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		*BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		*AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		*DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		*DihAt3, *DihAt4, *DihNum, *Boundary;
	int	popn;	// A flag for whether the file is compressed or not
	int	data_read;	// A flag to note whether the arrays are filled

	parm();
	~parm();
	FILE  *genopen (const char *name);
	void  genclose (FILE *);
	char  *get (int);
	void  preadln (FILE *, const char *, char *);
	int  readparm (char *);
	int  firstwat();
	int  read_fortran_12I6 (FILE *, int *, int);  // Read FORTRAN 12I6 format data
	int  moveto (FILE *, char *);  // Move to a section (AMBER 7)
} Ambertoppar;
#endif
