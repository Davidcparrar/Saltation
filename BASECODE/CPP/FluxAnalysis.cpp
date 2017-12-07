#include<fstream> 
#include<iostream>
#include<cmath>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>
#include"../SRC/Parameters.h"

using namespace std;

class ChargeInfo{
public:
  long int *Small;
  long int *Big;
  long int *Total;
  long int *SmallCharged;
  long int *BigCharged;
  long int *TotalCharged;
  int size;
  ChargeInfo(int n);
  ~ChargeInfo(void);
  void addSmall(int n);
  void addBig(int n);
  void addSmallCharged(int n);
  void addBigCharged(int n);
};

ChargeInfo::ChargeInfo(int n){
  size = n;
  int i;
  Small = new long int [size];
  Big = new long int [size];
  Total = new long int [size];
  SmallCharged = new long int [size];
  BigCharged = new long int [size];
  TotalCharged = new long int [size];

  for(i=0;i<size;i++){
    Small[i]=0;
    Big[i]=0;
    Total[i]=0;
    SmallCharged[i]=0;
    BigCharged[i]=0;
    TotalCharged[i]=0;
  }

};

ChargeInfo::~ChargeInfo(void){
  delete [] Small;
  delete [] Big;
  delete [] Total;
  delete [] SmallCharged;
  delete [] BigCharged;
  delete [] TotalCharged;
};

void ChargeInfo::addSmall(int n){
  Small[n]=Small[n]+1;
  Total[n]=Total[n]+1;
}

void ChargeInfo::addBig(int n){
  Big[n]=Big[n]+1;
  Total[n]=Total[n]+1;
}
void ChargeInfo::addSmallCharged(int n){
  SmallCharged[n]=SmallCharged[n]+1;
  TotalCharged[n]=TotalCharged[n]+1;
}

void ChargeInfo::addBigCharged(int n){
  BigCharged[n]=BigCharged[n]+1;
  TotalCharged[n]=TotalCharged[n]+1;
}


void usage()
{
  printf("\n");
  fprintf(stderr, "usage:  Saltation -i (name) -n (int) -h (int) [rtcfxyzRsTCMh] \n\
\n\
-n (int) \t number of Particles (required)    \n\
-i (int)\t first file (required)\n\
-w (int)\t last file (required)\n\
\n OPTIONAL arguments\n\
-h   \t show this message  \n\n");

}				
int main(int argc, char** argv){

  FILE *fin=NULL;
  char *base1;
  char input[15];
  char ginfile[256];
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int ini,last,N;
  ofstream fpout;
  ofstream fout;
  int SIZE = static_cast<int>(HEIGHT/D_MEAN);
  double timeSteady=1.;
  char outfile[256];
  while((option_letter = getopt(argc, argv, "i:n:w:e:v:g:ph")) != -1 ){
  
    if (option_letter == 'i') {    
      ini = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'n') {    
      N = atoi(optarg);
      usage_flag++;      
    }
    else if (option_letter == 'w') {    
      last = atoi(optarg);
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
    else if (option_letter == 'g') {    
      strcpy(ginfile, optarg);
      usage_flag++;      
    }	    
  }   

  if (usage_flag < 3) {   
    usage();
    exit(1);      
  } 

  int i,j,runningTime;


  FILE *fp=NULL;
  runningTime = last-ini;
  if(runningTime<0){
    printf("Cannot Open File\n");
    exit(1);	
  }    
  clog<<last<<endl;

  base1="FILES/file";
  double Area=DEEP*WIDTH;
  double m; 
  double rho=DENSITY;
  double rhoair=AIRDENSITY;
  double QTotal=0;
  double QPos=0;
  double QNeg=0;
  double QNeu=0;
  double x, y, z, Vx, Vy, Vz, R;
  double time, energy;
  int iy;

  //BI SIZE DIST

  double flux;
  int nPart=0;
  int iyf=0;
  int maxY=0;
  int counter=0;
  double q;
  double Low,High,qcharge,qchargeDensity;
  int filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "Flux%d.txt",filei);
  fpout.open(outfile);
  double AverageCharge=0;

    sprintf(outfile, "Energy%d.txt",filei);
    if((fin=fopen(outfile, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"Energia"<<endl;
      exit(1);
    }
  double *dataF;
  double *dataFP;
  double *dataFN;

  dataF=new double [runningTime+1];
  dataFP=new double [runningTime+1];
  dataFN=new double [runningTime+1];

  double media=0;
  for(j=0;j<=runningTime;j++){ 

    fscanf(fin,"%lf	%lf	%lf	%lf	%lf\n",&time,&energy,&flux,&qcharge,&qchargeDensity);

    sprintf (input, "%s%05d", base1, j);

    if((fp=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"File"<<endl;
      exit(1);
    }

    for(i=0;i<N;i++){  
      fscanf(fp, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R, &Low,&Low,&High,&High,&q);

      m=rho*R*R*R*M_PI*4.0/3.0;
      QTotal+=m*Vx;
      if(q>0)
        QPos+=m*Vx;
      else if(q<0)
        QNeg+=m*Vx;
      else
        QNeu+=m*Vx;
      

    }
    QTotal=(QTotal/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));      
    QPos=(QPos/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));      
    QNeg=(QNeg/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));      
    QNeu=(QNeu/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));     
    if(j>15){
    dataF[j]=QTotal;
    dataFP[j]=QPos;
    dataFN[j]=QNeg; 
    media+=(QPos+QNeg)/QTotal;
    }


    fpout<<scientific<<j*0.2<<"\t"<<QTotal<<"\t"<<QPos<<"\t"<<QNeg<<"\t"<<QNeu<<"\t"<<(QPos+QNeg)/QTotal<<endl;
    time=0;
    QTotal=0;
    QPos=0;
    QNeg=0;
    QNeu=0;

    fclose(fp);
  }
  fpout.close();

  media/=(runningTime-15);
  double standard_dev=0;
  QTotal=0;
  QPos=0;
  QNeg=0;

  for(j=16;j<=runningTime;j++){ 
    standard_dev+=((dataFN[j]+dataFP[j])/dataF[j]-media)*((dataFN[j]+dataFP[j])/dataF[j]-media);
    QTotal+=dataF[j];
    QPos+=dataFP[j];
    QNeg+=dataFN[j]; 
  }

  ofstream fcout;

  fcout.open(ginfile,ios::app);

  fcout<<WINDVEL<<"\t"<<media<<"\t"<<standard_dev/(runningTime-15)<<"\t"<<sqrt(standard_dev/(runningTime-15))<<"\t"<<(runningTime-15)<<"\t"<<QTotal/(runningTime-15)<<"\t"<<QNeg/(runningTime-15)<<"\t"<<QPos/(runningTime-15)<<endl;

  fcout.close();

  return 0;
}
