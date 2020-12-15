// Code written by Lei Huang for wham post-process of FEP data in NAMD

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define KBOLTZ	(1.987191E-03)
#define NCYCLE	(100000)
//#define TOL		(1.0E-7)
#define TOL           (1.0E-6)
#define max_data	(1000000)
#define max_wind	(20)

#define CHECK_POINTER(p)	{if((p)==NULL) {printf("Invalid pointer.\nQuit\n");	exit(1);}}

int *NPoints;
double **epprtd;
double *lambda;
double *F;
double *F2;
int nb_wind;
double kbt;
char szData[256];
char szFile_FE[256];

// max_data - max number of data points
// max_wind - max number of windows

void wham_FEP(void);
void wham(void);
void wham_E_Histogram(void);
double** Allocate_Mem_2D_Array(int Ny, int Nx); 
void Deallocate_Mem_2D_Array(double **p, int Ny, int Nx);
int** Allocate_Mem_2D_Array_Int(int Ny, int Nx); 
void Deallocate_Mem_2D_Array_Int(int **p, int Ny, int Nx);
void Free_Memory(void);


int main(int argc, char *argv[])
{
	if(argc <3)	{
		printf("Usage: wham_fep data_file output_wham_file\n");
		exit(1);
	}
	strcpy(szData, argv[1]);
	strcpy(szFile_FE, argv[2]);

	printf("----------------------------------\n");
	printf("Weigthed Histogram Analysis Method\n");
	printf("Author:  Benoit Roux (1995).\n");
	printf("Ported to c by Lei Huang (2012).\n\n");
	
	
	NPoints = new int[max_wind];
	epprtd = Allocate_Mem_2D_Array(max_wind, max_data);
	lambda = new double[max_wind];
	
	F = new double[max_wind];
	F2 = new double[max_wind];

	wham_FEP();

	Free_Memory();

	return 0;
}


//  call WHAM1(maxwind, maxtime, nb_wind, Ntime, epprtd, &
//     Lambda, F, F2)

void wham_FEP(void)
{
	int i, Time_1, Time_2, Time_3;
	double Temp, lambx, eppx;
	int Ntot1, Ntot2, ReadItem;
	int ioffset=0;	// the window used as the reference
	FILE *fIn, *fOut;
	
	Ntot1 = Ntot2 = 0;
	Temp = 300.0;	
	kbt=KBOLTZ*Temp;
	
	nb_wind = 0;
	NPoints[0] = 0;
	lambda[0] = 0.0;
	
	fIn = fopen(szData, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szData);
		Free_Memory();
		exit(1);
	}

	while(1)	{
		ReadItem = fscanf(fIn, "%lf%lf", &lambx, &eppx);
		if(ReadItem != 2)	{
			break;
		}
		
		if( lambx != lambda[nb_wind] )	{
			nb_wind++;
			if(nb_wind > max_wind)	{
				printf("Number of windows exceeds limit.\nQuit\n");
				exit(1);
			}
			NPoints[nb_wind] = 0;
			if(Ntot1 > Ntot2) Ntot2=Ntot1;
			Ntot1 = 0;
			lambda[nb_wind] = lambx;
		}
		
		
		Ntot1 ++;
		
		if(NPoints[nb_wind] < max_data)	{
			epprtd[nb_wind][NPoints[nb_wind]] = eppx;
			NPoints[nb_wind] = NPoints[nb_wind] + 1;
		}
	}
	fclose(fIn);
	
	if(Ntot2 > max_data)	{
		printf("Some data points have been left out.\n");
	}
	else	{
		printf("All data points have been read.\n");
	}
	
	printf("Window       #_of_Data       Lambda\n");
	
	for(i=0; i<=nb_wind; i++)	{
		printf("%6d       %10d      %.5lf\n", i+1, NPoints[i], lambda[i]);
        F[i] = 0.0;
        F2[i] = 0.0;
	}
	
//	Time_1 = time(NULL);
//	wham();		// the original wham by Benoit Roux
//	Time_2 = time(NULL);
//	printf("Ellapsed time %d seconds in wham().\n", Time_2-Time_1);

//	printf("\n\nWindow       #_of_Data     Lambda    F\n");
//	for(i=0; i<=nb_wind; i++)	{
//		printf("    %-d   %10d        %.5lf  %.5lf\n", i+1, NPoints[i], lambda[i], F2[i]-F2[ioffset]);
//	}

	Time_2 = time(NULL);
	wham_E_Histogram();	// loop over the histogram
	Time_3 = time(NULL);
	printf("\nEllapsed time %d seconds in wham().\n", Time_3-Time_2);

	printf("\n\nWindow       #_of_Data     Lambda    F\n");
	for(i=0; i<=nb_wind; i++)	{
		printf(" %4d   %10d        %.5lf  %.5lf\n", i+1, NPoints[i], lambda[i], F2[i]-F2[ioffset]);
	}

	fOut = fopen(szFile_FE, "w");
	fprintf(fOut, " 100  %.5lf\n", F2[nb_wind]-F2[0]);
	fclose(fOut);
	

}

void wham(void)
{
	int i, j, k, t, Converged=0, icycle;
	double arg1, arg2, argmax=0.0, Bot_i, Top_i, diff, *pDataTmp, kbt_Inv;
	
	kbt_Inv = 1.0/kbt;
	pDataTmp = new double[nb_wind+1];

	for(icycle=1; icycle<=NCYCLE; icycle++)	{
		for(k=0; k<=nb_wind; k++)	{
			F[k]=0.0;
			
			for(i=0; i<=nb_wind; i++)	{
				for(t=0; t<NPoints[i]; t++)	{
					arg1   = (-lambda[k]*epprtd[i][t])*kbt_Inv;
					argmax = arg1;
					
					for(j=0; j<=nb_wind; j++)	{
						arg2  = (F2[j]-lambda[j]*epprtd[i][t])*kbt_Inv;
						
						if(arg2 > argmax)	argmax = arg2;
						pDataTmp[j] = arg2;
					}
					Top_i = exp(arg1-argmax);
					
					Bot_i  = 0.0;
					for(j=0; j<=nb_wind; j++)	{
						arg2  = pDataTmp[j]-argmax;
						Bot_i  = Bot_i  + NPoints[j] * exp(arg2);
					}
					F[k] = F[k] + Top_i / Bot_i;
					
				}
			}
		}
		
		Converged = 1;
		for(k=0; k<=nb_wind; k++)	{
			F[k] = -kbt*log( F[k] );
//			printf("window %d: F(k) = %.5lf\n", k+1, F[k]);
			diff = fabs( F2[k] - F[k] );
			
			F2[k] = F[k] - F[0];
			
			if(diff > TOL) Converged = 0;
		}
		
		if(Converged) {
			printf("\nConverged at step: %d\n", icycle);
			delete []pDataTmp;
			return;
		}
		
	}
	
	printf("\nFail to converge in %d steps.\n", NCYCLE);
	delete []pDataTmp;
}

#define E_BIN_SIZE	(0.0004)
//#define E_BIN_SIZE      (0.0001)

void wham_E_Histogram(void)
{
	int i, j, k, t, Converged=0, icycle, Idx;
	double arg1, arg2, argmax=0.0, Bot_i, Top_i, diff, *pDataTmp, kbt_Inv, E_Local;
	double *pEmin=NULL, *pEmax=NULL, **pE_Hist=NULL;
	int *pN_Bin=NULL, **nHist_E=NULL;

	//start	to build the histogram of E for each window and free the memory for all data points, epprtd[][]
	pEmin = new double[nb_wind+1];
	CHECK_POINTER(pEmin)
	pEmax = new double[nb_wind+1];
	CHECK_POINTER(pEmax)
	pN_Bin = new int[nb_wind+1];
	CHECK_POINTER(pN_Bin)
	nHist_E = new int*[nb_wind+1];
	CHECK_POINTER(nHist_E)
	pE_Hist = new double*[nb_wind+1];
	CHECK_POINTER(pE_Hist)
	for(i=0; i<=nb_wind; i++)	{
		pEmin[i] = 1.0E100;
		pEmax[i] = -1.0E100;
		for(j=0; j<NPoints[i]; j++)	{
			if(epprtd[i][j] > pEmax[i])	{
				pEmax[i] = epprtd[i][j];
			}
			if(epprtd[i][j] < pEmin[i])	{
				pEmin[i] = epprtd[i][j];
			}
		}
		pN_Bin[i] = (int)((pEmax[i]-pEmin[i])/E_BIN_SIZE)+1;
		nHist_E[i] = new int[pN_Bin[i]];
		CHECK_POINTER(nHist_E[i])
		pE_Hist[i] = new double[pN_Bin[i]];
		CHECK_POINTER(pE_Hist[i])
		memset(nHist_E[i], 0, sizeof(int)*pN_Bin[i]);

		for(j=0; j<pN_Bin[i]; j++)	{	// assign the list of energies histogram
			pE_Hist[i][j] = pEmin[i] + j*E_BIN_SIZE;
		}

		for(j=0; j<NPoints[i]; j++)	{
			Idx = (int)((epprtd[i][j] - pEmin[i])/E_BIN_SIZE);
			nHist_E[i][Idx]++;
		}
	}
	Deallocate_Mem_2D_Array(epprtd, max_wind, max_data);
	epprtd = NULL;
	//end	to build the histogram of E for each window and free the memory for all data points, epprtd[][]


	kbt_Inv = 1.0/kbt;
	pDataTmp = new double[nb_wind+1];

	for(icycle=1; icycle<=NCYCLE; icycle++)	{
		for(k=0; k<=nb_wind; k++)	{
			F[k]=0.0;
			
			for(i=0; i<=nb_wind; i++)	{
				for(t=0; t<pN_Bin[i]; t++)	{
					if(nHist_E[i][t] == 0)	{
						continue;
					}
					E_Local = pE_Hist[i][t];
					arg1   = (-lambda[k]*E_Local)*kbt_Inv;
					argmax = arg1;
					
					for(j=0; j<=nb_wind; j++)	{
						arg2  = (F2[j]-lambda[j]*E_Local)*kbt_Inv;
						
						if(arg2 > argmax)	argmax = arg2;
						pDataTmp[j] = arg2;
					}
					Top_i = exp(arg1-argmax) * nHist_E[i][t];
					
					Bot_i  = 0.0;
					for(j=0; j<=nb_wind; j++)	{
						arg2  = pDataTmp[j]-argmax;
						Bot_i  = Bot_i  + NPoints[j] * exp(arg2);
					}
					F[k] = F[k] + Top_i / Bot_i;
				}
			}
		}
		
		Converged = 1;
		for(k=0; k<=nb_wind; k++)	{
			F[k] = -kbt*log( F[k] );
//			printf("window %d: F(k) = %.5lf\n", k+1, F[k]);
			diff = fabs( F2[k] - F[k] );
			
			F2[k] = F[k] - F[0];
			
			if(diff > TOL) Converged = 0;
		}
		
		if(Converged) {
			printf("\nConverged at step: %d\n", icycle);
			delete []pDataTmp;

			for(i=0; i<=nb_wind; i++)	{
				delete [](nHist_E[i]);
				delete [](pE_Hist[i]);
			}
			delete []nHist_E;
			delete []pE_Hist;

			delete []pEmin;
			delete []pEmax;
			delete []pN_Bin;

			return;
		}
		
	}
	
	printf("\nFail to converge in %d steps.\n", NCYCLE);
	delete []pDataTmp;

	for(i=0; i<=nb_wind; i++)	{
		delete [](nHist_E[i]);
		delete [](pE_Hist[i]);
	}
	delete []nHist_E;
	delete []pE_Hist;
	
	delete []pEmin;
	delete []pEmax;
	delete []pN_Bin;

	return;
}
#undef	E_BIN_SIZE

double** Allocate_Mem_2D_Array(int Ny, int Nx)
{
	double **p;
	int i;
	
	p = new double*[Ny];
	if(p == NULL)	{
		printf("Fail to allocate memory in Allocate_Mem_2D_Array().\nQuit\n");
		exit(1);
	}

	for(i=0; i<Ny; i++)	{
		p[i] = new double[Nx];
		if(p[i] == NULL)	{
			printf("Fail to allocate memory in Allocate_Mem_2D_Array().\nQuit\n");
			exit(1);
		}
	}

	return p;
}

void Deallocate_Mem_2D_Array(double **p, int Ny, int Nx)
{
	int i;

	for(i=0; i<Ny; i++)	{
		delete [](p[i]);
	}
	delete []p;
	p=NULL;
}


int** Allocate_Mem_2D_Array_Int(int Ny, int Nx)
{
	int **p;
	int i;
	
	p = new int*[Ny];
	if(p == NULL)	{
		printf("Fail to allocate memory in Allocate_Mem_2D_Array_Int().\nQuit\n");
		exit(1);
	}

	for(i=0; i<Ny; i++)	{
		p[i] = new int[Nx];
		if(p[i] == NULL)	{
			printf("Fail to allocate memory in Allocate_Mem_2D_Array_Int().\nQuit\n");
			exit(1);
		}
	}

	return p;
}

void Deallocate_Mem_2D_Array_Int(int **p, int Ny, int Nx)
{
	int i;

	for(i=0; i<Ny; i++)	{
		delete [](p[i]);
	}
	delete []p;
	p=NULL;
}

void Free_Memory(void)
{
	if(epprtd)	{
		Deallocate_Mem_2D_Array(epprtd, max_wind, max_data);
	}
	if(NPoints)	delete []NPoints;
	if(lambda)	delete []lambda;
	if(F)	delete []F;
	if(F2)	delete []F2;
}

