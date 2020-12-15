// Code written by Lei Huang for wham post-process of FEP data in NAMD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N_WIN	(128)

int n_Win, nLine[MAX_N_WIN];
double rcut2[MAX_N_WIN];

inline int Get_Index_rcut(double rcut);

int main(int argc, char *argv[])
{
	FILE *fIn, *fOut[MAX_N_WIN];
	int ReadItem, i, Idx;
	double fTmp=0.0, rCut_1, rCut_2, lambda, dE;
	char szEType[256], szOutput[256], szLine[256], *ReadLine;

	if(argc != 3)	{
		printf("Usage: rcut_2_list data_file\nThe output will be repu_*.wham.\n");
		exit(1);
	}

	fIn = fopen(argv[1], "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", argv[1]);
		exit(1);
	}

	n_Win=0;
	while(1)	{
		ReadItem = fscanf(fIn, "%lf", &(rcut2[n_Win]));
		if(ReadItem == 1)	{
			nLine[n_Win] = 0;
			n_Win++;
		}
		else	{
			break;
		}
	}
	fclose(fIn);



	fIn = fopen(argv[2], "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", argv[2]);
		exit(1);
	}

	for(i=0; i<n_Win; i++)	{
		sprintf(szOutput, "repu_%d.wham", i);
		fOut[i] = fopen(szOutput, "w");
		if(fOut[i] == NULL)	{
			printf("Fail to open file: %s for write.\nQuit\n", szOutput);
			exit(1);
		}
	}

//	d_Cut = rcut2[1]-rcut2[0];
//	d_Cut_Inv = 1.0/d_Cut;

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 128, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%s%lf%lf%lf%lf", szEType, &rCut_1, &rCut_2, &lambda, &dE);
		if(ReadItem == 5)	{
			if(strcmp(szEType, "FEP_WCA_REP")==0)	{
//				Idx = (int)(rCut2*d_Cut_Inv + 0.1) - 1;
				Idx = Get_Index_rcut(rCut_2);
				if(Idx < 0)	{
					printf("Error: unrecognized rcut2  %lf\n", rCut_2);
				}
				else	{
					fprintf(fOut[Idx], "%.4lf  %.4lf\n", lambda, dE);
					nLine[Idx]++;
				}
			}
		}
//		else	{
//			break;
//		}
	}
	fclose(fIn);

	for(i=0; i<n_Win; i++)	{
		fclose(fOut[i]);
	}


	printf("Window   rcut2    #_data\n");
	for(i=0; i<n_Win; i++)	{
		printf(" %6d  %6.4lf  %7d\n", i+1, rcut2[i], nLine[i]);
	}


	return 0;
}

inline int Get_Index_rcut(double rcut)
{
	for(int i=0; i<n_Win; i++)	{
//		if(fabs(rcut2[i]-rcut) < 2.0E-4)	{
		if(rcut2[i] == rcut)	{
			return i;
		}
	}
	return -1;
}

