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
  char *base1,*base2;
  char input[15];
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int ini,last,N;
  ofstream fpout;
  ofstream fout;
  int SIZE = static_cast<int>(HEIGHT/D_MEAN);
  double timeSteady=1.;
  char outfile[256];
  while((option_letter = getopt(argc, argv, "i:n:w:e:v:ph")) != -1 ){
  
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
  }   

  if (usage_flag < 3) {   
    usage();
    exit(1);      
  } 

  int i,j,runningTime;


  FILE *fp=NULL;
  FILE *fw=NULL;
  runningTime = last-ini;
  if(runningTime<0){
    printf("Cannot Open File\n");
    exit(1);	
  }    
  clog<<last<<endl;

  base1="FILES/file";
  base2="PROFILES/Profile";
  double Area=DEEP*WIDTH;
  double m; 
  double rho=DENSITY;
  double rhoair=AIRDENSITY;
  double Q=0;
  double Qflux=0;
  double x, y, z, Vx, Vy, Vz, R;
  double time, energy;
  int iy;
  long double windVelocity[SIZE];
  long double dwindVelocity[SIZE];
  long double Tau[SIZE];
  long double previousWindVelocity[SIZE];
  long double concentration[SIZE];
  long double Carga[SIZE];
  long double Masa[SIZE];
  long double VelocityX[SIZE];
  long int nParticles[SIZE];
  long int Total[SIZE];
  long int TotalCh[SIZE];
  //ChargeInfo Dist(SIZE);

  //BI SIZE DIST

  double AuxWindVelocity;
  double AuxdwindVelocity;
  double AuxTau;
  double AuxPreviousWindVelocity;
  double AuxHeight;
  double AuxTauSingle;
  double AuxnTau;
  double AuxForce;
  double AuxI;
  double flux;
  int nPart=0;
  int iyf=0;
  int maxY=0;
  int counter=0;
  double q;
  double Low,High,qcharge,qchargeDensity;
  int filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "Data%d.txt",filei);
  fpout.open(outfile);
  double AverageCharge=0;

   for(i=0;i<SIZE;i++){  
     windVelocity[i]=0;
     dwindVelocity[i]=0;
     Tau[i]=0;
     previousWindVelocity[i]=0;
     concentration[i]=0;
     VelocityX[i]=0;
     nParticles[i]=0;
     Carga[i]=0;
     Masa[i]=0; 
     Total[i]=0; ; 
     TotalCh[i]=0; 
   }
    sprintf(outfile, "Energy%d.txt",filei);
    if((fin=fopen(outfile, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"Energia"<<endl;
      exit(1);
    }

  for(j=0;j<=runningTime;j++){ 

    fscanf(fin,"%lf	%lf	%lf	%lf	%lf\n",&time,&energy,&flux,&qcharge,&qchargeDensity);
    
    sprintf (input, "%s%05d", base2, j);  
    if((fw=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"Profile"<<endl;
      exit(1);
    }

    for(i=0;i<SIZE;i++){  

      fscanf(fw,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",&AuxWindVelocity,&AuxPreviousWindVelocity,&AuxdwindVelocity,&AuxTau,&AuxTauSingle,&AuxnTau,&AuxI,&AuxHeight,&AuxForce);
        if(time>timeSteady){
            windVelocity[i]+=AuxWindVelocity;
  	    previousWindVelocity[i]+=AuxPreviousWindVelocity;
            Tau[i]+=AuxTau;
            dwindVelocity[i]+=AuxdwindVelocity;
        }
    }
    sprintf (input, "%s%05d", base1, j);

    if((fp=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"File"<<endl;
      exit(1);
    }

    for(i=0;i<N;i++){  
      fscanf(fp, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R, &Low,&Low,&High,&High,&q);
      iy =static_cast<int>(y/0.0001);
      iyf =static_cast<int>(y/0.0002);

      if(iy > maxY)
        maxY=iy;

      m=rho*R*R*R*M_PI*4.0/3.0;
      if(time>timeSteady){
        concentration[iy]+=R*R*R*M_PI*4.0/3.0;
        VelocityX[iy]+=Vx;
        nParticles[iy]++;
        Carga[iy]+=q;
        Masa[iy]+=m;

      }
      Q+=m*Vx;
      if(iyf>=AuxHeight*2.0){
      Qflux+=m*Vx;
      nPart++;
      AverageCharge+=abs(q);
      }
    }
    Q=(Q/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));      
    Qflux=(Qflux/Area)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN));  
    if(time>=timeSteady) counter ++;    
    fpout<<scientific<<time<<"\t"<<Q<<"\t"<<Qflux<<"\t"<<nPart<<"\t"<<maxY*0.0001<<"\t"<<AverageCharge/nPart<<endl;
    time=0;
    Q=0;
    Qflux=0;
    nPart=0;
    maxY=0;
    AverageCharge=0;
    iyf=0;

    fclose(fp);
    fclose(fw);
  }
  sprintf(outfile, "MeanProfile%d.txt",filei);
  fout.open(outfile);
  int NCELL=5;
  double auxconcentration = 0;
  double auxforce = 0;
  double auxVelocity = 0;
  double auxWindVel=0;
  double auxq=0;
  double auxMass=0;
  int auxnPar=0;
  int conter=0;
  int TotalTotal=0;
  int TotalTotalCharged=0;
  for(i=0;i<SIZE;i++){
    auxconcentration+=concentration[i];
    auxforce+=0;
    auxVelocity+=VelocityX[i];
    auxWindVel+=windVelocity[i];
    auxnPar+=nParticles[i];
    auxq+=Carga[i];
    auxMass+=Masa[i];
    conter++;
    if((conter)%NCELL==0){
      if(auxnPar!=0){
        auxVelocity=auxVelocity/((double)auxnPar);
        auxMass=auxq/auxMass;
        auxq=auxq/((double)auxnPar);;
      }

      fout<<scientific<<(i-(NCELL-1))*D_MEAN<<"\t"<<auxconcentration/(counter*NCELL*0.0001*DEEP*WIDTH)<<"\t"<<auxforce/(counter*NCELL)<<"\t"<<auxVelocity<<"\t"<<auxWindVel/(counter*NCELL)<<"\t"<<auxVelocity*auxconcentration/(counter*NCELL*0.0001*DEEP*WIDTH)<<"\t"<<auxq<<"\t"<<auxMass<<endl;
      auxconcentration = 0;
      auxforce = 0;
      auxVelocity = 0;
      auxWindVel=0;
      auxnPar=0;
      auxq=0;
      auxMass=0;
    }
    TotalTotal+=Total[i];
    TotalTotalCharged+=TotalCh[i];
  }
  fout.close();

  return 0;
}
