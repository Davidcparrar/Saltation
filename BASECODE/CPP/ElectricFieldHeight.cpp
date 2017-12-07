#include"../SRC/Collider.h"
#include<fstream> 
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>
#include<time.h>
#include<dirent.h>
#include"../SRC/Random.h"

using namespace std;

//---------------Main Function---------------------------

void usage()
{
  printf("\n");
  fprintf(stderr, "usage:  Saltation -i (name) -n (int) -w (int) [rtcfxyzRsTCMh] \n\
\n\
-n (int) \t number of particles (required)    \n\
-b (int) \t Initial file number (required)    \n\
-e (int) \t Last file number (required)    \n\
\n OPTIONAL arguments\n\
-h   \t show this message  \n\n");

}				
int main(int argc, char** argv){
  FILE *fin=NULL;
  char outfile1[256];
  char ginfile[256];
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int numberOfFiles = 0;
  int numberOfFiles_0 = 0;
  int numberOfFiles_02 = 0;
  int numberOfFiles2 = 0;
  double runningTime = 100.;
  int counter=0;
  while((option_letter = getopt(argc, argv, "b:n:e:w:i:g:v:ph")) != -1 ){
  
    if (option_letter == 'e') {    
      numberOfFiles = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'g') {    
      strcpy(ginfile, optarg);
      usage_flag++;      
    }
    else if (option_letter == 'b') {    
      numberOfFiles_0 = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'w') {    
      numberOfFiles2 = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'i') {    
      numberOfFiles_02 = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'n') {    
      N = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'v') {    
      WINDVEL = atof(optarg);
      usage_flag++;      
    }    
    else if (option_letter == 'h') {    
      usage();
      exit(1);      
    } 	    
  }   

  if (usage_flag < 3) {   
    usage();
    exit(1);      
  } 
  Collider Aeolian;
  Body *Sand;
  Sand = new Body[N];
  int i,j;
  ofstream fout;
  double x,y,z,Vx,Vy,Vz,R,q;
  double lowE, highE, lowET, highET;
  int filei=static_cast<int>(100.0*WINDVEL);
  char *base1;
  base1="FILES/file";
  char input[15];
  clock_t timeCheck;
  timeCheck=clock();
  Aeolian.Init(N);
  double *EFTM;
  int *CountEF;
  double ***Carga;
  Carga= new double**[Aeolian.Lxp];
  for(i=0;i<Aeolian.Lxp;i++){
      Carga[i]= new double*[Aeolian.Lyp];
    for(j=0;j<Aeolian.Lyp;j++)
    Carga[i][j]= new double[Aeolian.Lzp];
  }

  int k;

  Aeolian.ReadGreenFunction(ginfile);
  for(k=0;k<Aeolian.Lzp;k++)
  for(i=0;i<Aeolian.Lxp;i++)
    for(j=0;j<Aeolian.Lyp;j++){

      Carga[i][j][k]=0.;
    }


  clog<<(numberOfFiles2-numberOfFiles_02)<<endl;

  sprintf(outfile1, "ElectricFieldHeights%d.txt",filei);
  fout.open(outfile1);

  for(k=0;k<Aeolian.Lzp;k++)
  for(i=0;i<Aeolian.Lxp;i++)
  for(j=0;j<Aeolian.Lyp;j++)
    Carga[i][j][k]=0.;
     
for(int t=numberOfFiles_02;t<=numberOfFiles2;t++){
    counter++;
    sprintf (input, "%s%05d", "PROFILES/ChargeMap", t);
    clog<<"Procesing file .... "<<t<<endl;
    if((fin=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"File "<<input<<endl;
      exit(1);
    }
  for(k=0;k<Aeolian.Lzp;k++)
      for(j=0;j<Aeolian.Lyp;j++)
        for(i=0;i<Aeolian.Lxp;i++){
        fscanf(fin, "%lf	%lf	%lf	%lf\n", &x, &y, &z, &q);
        Carga[i][j][k]+=q/ELEMENTARYCHARGE;
    }
    fclose(fin);

  Aeolian.PotentialSolverTotal(Carga);
  Aeolian.ElectricFieldTotal();
    fout<<t*0.2<<"\t"<<Aeolian.Ey[5][5][1]*ELEMENTARYCHARGE*KE/1000<<"\t"<<Aeolian.Ey[5][75][1]*ELEMENTARYCHARGE*KE/1000<<"\t"<<Aeolian.Ey[5][175][1]*ELEMENTARYCHARGE*KE/1000<<"\t"<<Aeolian.Ey[5][275][1]*ELEMENTARYCHARGE*KE/1000<<endl;
  for(k=0;k<Aeolian.Lzp;k++)
      for(j=0;j<Aeolian.Lyp;j++)
        for(i=0;i<Aeolian.Lxp;i++){
        Carga[i][j][k]=0;
        Aeolian.potential[i][j][k]=0;
    }
  }
     
  fout.close();
  timeCheck = clock() - timeCheck;
  clog<<"Running time "<<runningTime<<endl;
  clog<<"Elapsed time "<<(static_cast<double>(timeCheck))/CLOCKS_PER_SEC<<endl;


  return 0;
}
