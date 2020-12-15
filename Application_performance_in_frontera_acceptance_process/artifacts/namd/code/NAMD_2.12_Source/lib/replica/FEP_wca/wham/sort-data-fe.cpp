// Code written by Lei Huang for wham post-process of FEP data in NAMD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN_COL_1	(0.0)
#define MAX_COL_1	(1.0)
#define N_BIN		(200)

#ifndef min(a,b)
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

#define MAX_LINE (10000000)
double *pcol_1, *pcol_2, BinSize;
int ID[N_BIN], BinCount[N_BIN], *pIDList;

int nData=0;
char szInput[512], szOutput[512];

int ReadData(int nLine);
void Assign_Bin_ID(void);
void Write_Sorted_Data(void);

int main(int argc, char *argv[])
{
	if(argc<3)	{
		printf("Usage: sort-data you-data new-data #-lines\n");
		exit(1);
	}
	else if(argc==3)	{
		nData = MAX_LINE;
	}
	else if(argc==4)	{
		nData = atoi(argv[3]);
	}

	strcpy(szInput, argv[1]);
	strcpy(szOutput, argv[2]);

	pcol_1 = new double[nData];
	if(pcol_1 == NULL)	{
		printf("Fail to allocate memory for pcol_1 in sort-data-wham.\nQuit\n");
		exit(1);
	}

	pcol_2 = new double[nData];
	if(pcol_2 == NULL)	{
		printf("Fail to allocate memory for pcol_2 in sort-data-wham.\nQuit\n");
		delete []pcol_1;
		exit(1);
	}

	pIDList = new int[nData];

	if(pIDList == NULL)	{
		printf("Fail to allocate memory for pIDList in sort-data-wham.\nQuit\n");
		delete []pcol_1;
		delete []pcol_2;
		exit(1);
	}

	ReadData(nData);
	Assign_Bin_ID();

	Write_Sorted_Data();


	delete	[]pcol_1;
	delete	[]pcol_2;
	delete	[]pIDList;

	return 0;
}

int ReadData(int nLine)
{
	FILE *fIn;
	int nCount = 0, ReadItem;
	double fTmp=0.0;

	fIn = fopen(szInput, "r");
	if(fIn == NULL)	{
		printf("Fail to open file %s for read.\nQuit\n", szInput);
		exit(1);
	}

	while(1)	{
		ReadItem = fscanf(fIn, "%lf %lf", &(pcol_1[nCount]), &(pcol_2[nCount]));
		if(ReadItem == 2)	{
			nCount++;
			if(nCount >= nData)	{
				break;
			}
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	nData = nCount;

	return nCount;
}

void Assign_Bin_ID(void)
{
	int i, Idx, MaxID;
	double BinSize_Inv;

	MaxID = N_BIN - 1;
	BinSize = (MAX_COL_1-MIN_COL_1)/N_BIN;
	BinSize_Inv = 1.0/BinSize;
	memset(BinCount, 0, sizeof(int)*N_BIN);

	for(i=0; i<nData; i++)	{
		Idx = int((pcol_1[i] - MIN_COL_1)*BinSize_Inv);
		Idx = min(MaxID, max(0, Idx));

		pIDList[i] = Idx;

		BinCount[Idx] ++;
	}
	
}

void Write_Sorted_Data(void)
{
	FILE *fOut;
	int ID, i;

	fOut = fopen(szOutput, "w");
	if(fOut == NULL)	{
		printf("Fail to open file %s for read.\nQuit\n", szOutput);
		exit(1);
	}

	for(ID=0; ID<N_BIN; ID++)	{
		if(BinCount[ID] > 0)	{
			for(i=0; i<nData; i++)	{
				if(pIDList[i] == ID)	{
					fprintf(fOut, "%8.5lf %12.5lf\n", pcol_1[i], pcol_2[i]);
				}
			}
		}
	}

	fclose(fOut);
}

