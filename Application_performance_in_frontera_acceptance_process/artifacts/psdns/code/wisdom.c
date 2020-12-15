#include<fftw3.h>
#include<stdio.h>

#ifdef NUS_XCOMP
void epfftw_wisdom(int *iopt,int *taskid, int *iprec)
#else
void epfftw_wisdom_(int *iopt,int *taskid, int *iprec)
#endif
{
FILE *fp;
int ierr;

// iopt==1: try to read from wisdom
// iopt==2: write wisdom

// iprec==1: single precision 
// iprec==2: double precision 

if (*iopt==0) // read wisdom if exists 
{
  if (*iprec==1){ fp=fopen("wisdoms","r");};
  if (*iprec==2){ fp=fopen("wisdomd","r");};
  if (fp) {
    if (*iprec==1){ ierr=fftwf_import_wisdom_from_file(fp);};
    if (*iprec==2){ ierr=fftw_import_wisdom_from_file(fp);};
    if (ierr!=1){
//      printf("Cannot import wisdom in taskid: %d\n",*taskid);
//      printf("Measured plans are to be created\n");
    } else {
      fclose(fp);
    };	
  } else { 
//   printf("Cannot open wisdom in taskid: %d\n",*taskid);
//   printf("Measured plans are to be created\n");
  };
}

// Save the wisdom into the file : wisdom 
if (*iopt==1) {
  if (*taskid==0) {
    if (*iprec==1){ fp=fopen("wisdoms","w");};
    if (*iprec==2){ fp=fopen("wisdomd","w");};
    if (fp) {
      if (*iprec==1){ fftwf_export_wisdom_to_file(fp);};
      if (*iprec==2){ fftw_export_wisdom_to_file(fp);};
      fclose(fp);
	printf ("Wisdom saved in file wisdom\n");
    } else {
      printf ("Cannot save wisdom\n");
    };
  };
};

}



