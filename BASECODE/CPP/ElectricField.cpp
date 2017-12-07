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
  } 

  for(int m=0;m<Aeolian.Lzp;m++)
    for(int h=0;h<Aeolian.Lxp;h++)
      for(int l=0;l<Aeolian.Lyp;l++)
        Carga[h][l][m]=Carga[h][l][m]/(counter);

  Aeolian.PotentialSolverTotal(Carga);
  Aeolian.ElectricFieldTotal();

  clog<<(numberOfFiles2-numberOfFiles_02)<<endl;
  double HighEF=0;
  double HeightMaxEF=0;
  sprintf(outfile1, "ElectricFieldRunTime%d.txt",filei);
  fout.open(outfile1);
  double AverageChargeHeight=0;
  for(int l=1;l<Aeolian.Lyp;l++){
    double AverageFieldx=0;
    double AverageFieldy=0;
    double AverageFieldz=0;
    AverageChargeHeight=0;
    for(int m=0;m<Aeolian.Lzp;m++)
    for(int h=0;h<Aeolian.Lxp;h++){              
      AverageFieldx+=Aeolian.Ex[h][l][m];
      AverageFieldy+=Aeolian.Ey[h][l][m];
      AverageFieldz+=Aeolian.Ez[h][l][m];
      AverageChargeHeight+=Carga[h][l][m];
      }
  fout<<l*Aeolian.D_CELL<<"\t"<<sqrt(AverageFieldx*AverageFieldx+AverageFieldy*AverageFieldy+AverageFieldz*AverageFieldz)*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldx*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldz*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageChargeHeight*ELEMENTARYCHARGE/(Aeolian.Lxp*Aeolian.Lzp*Aeolian.D_CELL*Aeolian.D_CELL*Aeolian.D_CELL)<<"\t"<<AverageFieldx*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)*AverageChargeHeight*ELEMENTARYCHARGE<<"\t"<<-AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)*AverageChargeHeight*ELEMENTARYCHARGE<<"\t"<<AverageFieldz*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)*AverageChargeHeight*ELEMENTARYCHARGE<<endl;
    if(HighEF==0) HighEF=AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000;
    if(AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<0 && HeightMaxEF==0) HeightMaxEF=l*Aeolian.D_CELL;
    }

  fout.close();
 sprintf(outfile1, "ElectricFieldRunTimeScaled%d.txt",filei);
  fout.open(outfile1);
  for(int l=1;l<Aeolian.Lyp;l++){
    double AverageFieldx=0;
    double AverageFieldy=0;
    double AverageFieldz=0;
    AverageChargeHeight=0;
    for(int m=0;m<Aeolian.Lzp;m++)
    for(int h=0;h<Aeolian.Lxp;h++){              
      AverageFieldy+=Aeolian.Ey[h][l][m];
      }
        if((AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/(1000))>0){
	fout<<l*Aeolian.D_CELL/HeightMaxEF<<"\t"<<AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/(1000*HighEF)<<endl;}
   }
  fout.close();

  for(k=0;k<Aeolian.Lzp;k++)
  for(i=0;i<Aeolian.Lxp;i++)
  for(j=0;j<Aeolian.Lyp;j++)
    Carga[i][j][k]=0.;
  counter=0;
     
  for(j=numberOfFiles_0;j<=numberOfFiles;j++){
    counter++;
    sprintf (input, "%s%05d", base1, j);

    if((fin=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"File "<<input<<endl;
      exit(1);
    }
    for(i=0;i<N;i++){
        fscanf(fin, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R, &lowE, &highE, &lowET, &highET, &q);
      Sand[i].Init(x,y,z,Vx,Vy,Vz,R,highET, highET, q/ELEMENTARYCHARGE);
    }
    
    Aeolian.PrintElectric(Sand);
    for(int m=0;m<Aeolian.Lzp;m++)
    for(int h=0;h<Aeolian.Lxp;h++)
      for(int l=0;l<Aeolian.Lyp;l++){
        Carga[h][l][m]+=Aeolian.ChargeGrid[h][l][m];
      }
    fclose(fin);
  }
      for(int m=0;m<Aeolian.Lzp;m++)
  for(int h=0;h<Aeolian.Lxp;h++)
      for(int l=0;l<Aeolian.Lyp;l++)
        Carga[h][l][m]=Carga[h][l][m]/(counter);

  Aeolian.PotentialSolverTotal(Carga);
  Aeolian.ElectricFieldTotal();

  sprintf(outfile1, "ElectricFieldMap%d.txt",filei);
  fout.open(outfile1);
  double MaxE=0;
  counter=0;
  for(int h=0;h<Aeolian.Lxp;h++)
      for(int l=0;l<Aeolian.Lyp;l++){
        if(sqrt(Aeolian.Ex[h][l][Aeolian.Lzp/2]*Aeolian.Ex[h][l][Aeolian.Lzp/2]+Aeolian.Ey[h][l][Aeolian.Lzp/2]*Aeolian.Ey[h][l][Aeolian.Lzp/2])>MaxE)MaxE=sqrt(Aeolian.Ex[h][l][Aeolian.Lzp/2]*Aeolian.Ex[h][l][Aeolian.Lzp/2]+Aeolian.Ey[h][l][Aeolian.Lzp/2]*Aeolian.Ey[h][l][Aeolian.Lzp/2]) ;
  }
  clog<<MaxE*ELEMENTARYCHARGE*KE/1000<<" kV/m"<<endl;


    for(int h=0;h<Aeolian.Lxp;h++)
      for(int l=0;l<Aeolian.Lyp;l++){
        fout<<scientific<<h<<"\t"<<l<<"\t"<<sqrt(Aeolian.Ex[h][l][Aeolian.Lzp/2]*Aeolian.Ex[h][l][Aeolian.Lzp/2]+Aeolian.Ey[h][l][Aeolian.Lzp/2]*Aeolian.Ey[h][l][Aeolian.Lzp/2])*ELEMENTARYCHARGE*KE<<"\t"<<Aeolian.Ex[h][l][Aeolian.Lzp/2]*ELEMENTARYCHARGE*KE<<"\t"<<Aeolian.Ey[h][l][Aeolian.Lzp/2]*ELEMENTARYCHARGE*KE<<"\t"<<Carga[h][l][Aeolian.Lzp/2]*ELEMENTARYCHARGE<<"\t"<<Aeolian.Ex[h][l][Aeolian.Lzp/2]/MaxE<<"\t"<<Aeolian.Ey[h][l][Aeolian.Lzp/2]/MaxE<<endl;
      }

  fout.close();
  sprintf(outfile1, "ElectricFieldVsHeight%d.txt",filei);
  fout.open(outfile1);
 
  AverageChargeHeight=0;
  for(int l=0;l<Aeolian.Lyp;l++){
    double AverageFieldx=0;
    double AverageFieldy=0;
    double AverageFieldz=0;
    AverageChargeHeight=0;
    for(int m=0;m<Aeolian.Lzp;m++)
    for(int h=0;h<Aeolian.Lxp;h++){              
      AverageFieldx+=Aeolian.Ex[h][l][m];
      AverageFieldy+=Aeolian.Ey[h][l][m];
      AverageFieldz+=Aeolian.Ez[h][l][m];
      AverageChargeHeight+=Carga[h][l][m];
      }
  fout<<l*Aeolian.D_CELL<<"\t"<<sqrt(AverageFieldx*AverageFieldx+AverageFieldy*AverageFieldy+AverageFieldz*AverageFieldz)*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldx*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldy*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageFieldz*ELEMENTARYCHARGE*KE/(Aeolian.Lxp*Aeolian.Lzp)/1000<<"\t"<<AverageChargeHeight*ELEMENTARYCHARGE/(Aeolian.Lxp*Aeolian.Lzp*Aeolian.D_CELL*Aeolian.D_CELL*Aeolian.D_CELL)<<endl;
    }

  fout.close();

  timeCheck = clock() - timeCheck;
  clog<<"Running time "<<runningTime<<endl;
  clog<<"Elapsed time "<<(static_cast<double>(timeCheck))/CLOCKS_PER_SEC<<endl;


  return 0;
}
