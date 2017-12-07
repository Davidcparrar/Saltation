#include"SRC/Collider.h"
#include<fstream> 
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>
#include<time.h>
#include"SRC/Random.h"

using namespace std;

void printSaltationInfo(ofstream &fout, double runningTime){
  time_t tim; 
  time(&tim); 
  fout<<endl;
  fout<<"--------Saltation Program 3D --------"<<endl;
  fout << ctime(&tim) << endl;
  fout<<endl;
  fout<<scientific<<"  Running Time: "<<runningTime<<" (s)"<<endl;
  fout<<"  Number of Particles: "<<N<<endl;
  #ifdef BCY
  fout<<"  Upper and Lower Walls"<<endl;
  #endif
  #ifdef WIND
  fout<<"  Wind: "<<WINDVEL<<" (m/s)"<<endl;
  #endif
  fout<<scientific<<"  Surface Charge Density: "<<RHO_H<<" (# electrons/m^2)"<<endl;
  fout<<scientific<<"  Height of Energy Barrier: "<<Eb0/ELECTRONVOLT<<" (eV)"<<endl;
  #ifdef ELECF
  fout<<scientific<<"  Charge rate from low energy level to high energy level : "<<LIGHT<<" (s^-1)"<<endl;
  #endif
  fout<<endl;
  fout<<"--------Parameters--------"<<endl;
  fout<<endl;
  fout<<scientific<<"  Step time: "<<DELTAT<<" (s)"<<endl;
  fout<<scientific<<"  Running Time: "<<runningTime<<" (s)"<<endl;
  #ifdef WIND
  fout<<scientific<<"  Air density: "<<AIRDENSITY<<" (Kg/m^3)"<<endl;
  fout<<scientific<<"  Air dynamic viscosity: "<<AIRVISCOSITY<<" (Kg/m s)"<<endl;
  #endif
  fout<<scientific<<"  Sand density: "<<DENSITY<<" (Kg/m^3)"<<endl;
  fout<<scientific<<"  Sand stiffness: "<<K<<endl;
  fout<<scientific<<"  Walls stiffness: "<<KWALL<<endl;
  fout<<scientific<<"  Sand energy dissipation coefficient: "<<GAMMA<<endl;
  fout<<scientific<<"  Walls energy dissipation coefficient: "<<GAMMAWALL<<endl;
  fout<<scientific<<"  System height: "<<HEIGHT*1000<<" (mm)"<<endl;
  fout<<scientific<<"  System width: "<<WIDTH*1000<<" (mm)"<<endl;
  fout<<scientific<<"  System depth: "<<DEEP*1000<<" (mm)"<<endl;
  fout<<endl;
  fout<<"--------------------------"<<endl;

}


bool checkExists(string file){
    ifstream file_to_check (file.c_str());
    if(file_to_check.is_open())
      return true;
    return false;
    file_to_check.close();
}
void LMT(double DM){
  LENGTH  = 1.0/10000.0;
  MASS    = 1.0;
  TIME    = 1.0;
  AMPERE  = 1.0;//ELEMENTARYCHARGE/TIME;
}

void DT(void){
  WIDTH = WIDTH/LENGTH;
  HEIGHT = HEIGHT/LENGTH;
  K=K*TIME*TIME;
  GAMMA=GAMMA*TIME/MASS;
  GAMMAWALL=GAMMAWALL*TIME/MASS;
  RADIUS_MAX=RADIUS_MAX/LENGTH;
  AIRVISCOSITY=AIRVISCOSITY*LENGTH*TIME/MASS;
  ROUGHNESS=ROUGHNESS/LENGTH;
  AIRDENSITY=AIRDENSITY*LENGTH*LENGTH*LENGTH/MASS;
  DEEP=DEEP/LENGTH;
  WINDVEL=WINDVEL*TIME/LENGTH;
  DENSITY=DENSITY*LENGTH*LENGTH*LENGTH/MASS;
  PLANCK=PLANCK*TIME/(LENGTH*LENGTH*MASS);
  REDUCEDPLANCK=REDUCEDPLANCK*TIME/(LENGTH*LENGTH*MASS);
  ELECTRONMASS=ELECTRONMASS/MASS;
  ELEMENTARYCHARGE = ELEMENTARYCHARGE/(AMPERE*TIME);
  ELECTRONVOLT=ELECTRONVOLT*TIME*TIME/(MASS*LENGTH*LENGTH);
  Eb0=Eb0*TIME*TIME/(MASS*LENGTH*LENGTH);
  EbMAX=EbMAX*TIME*TIME/(MASS*LENGTH*LENGTH);
  RHO_H=RHO_H*LENGTH*LENGTH;
  TCOLL=TCOLL/TIME;
  A2=A2/(LENGTH*LENGTH);
  E0=E0*LENGTH*LENGTH*LENGTH*MASS/(TIME*TIME*TIME*TIME*AMPERE*AMPERE); //AMPERE NOT INCLUDED
  KE=KE/(LENGTH*LENGTH*LENGTH*MASS/(TIME*TIME*TIME*TIME*AMPERE*AMPERE)); //AMPERE NOT INCLUDED
  D_BRIDGE=D_BRIDGE/LENGTH;
  D_BRIDGEK=D_BRIDGEK/LENGTH;
  DELTAT = DELTAT/TIME;
  GRAVITY=GRAVITY/LENGTH;
  D_MEAN=D_MEAN/LENGTH;
  MOBILITY=MOBILITY*LENGTH*LENGTH;
/*  K=0.5;
  KWALL=1;
  GAMMA=0.000012;
  GAMMAWALL=0.000033;*/

  K=500;
  GAMMA=0.00017;
  KWALL=800;
  GAMMAWALL=0.0005;
  clog<<endl;
  clog<<"--------Parameters--------"<<endl;
  clog<<scientific<<"  Step time: "<<DELTAT*TIME<<endl;
  clog<<scientific<<"  Time: "<<TIME<<endl;
  clog<<scientific<<"  Mass: "<<MASS<<endl;
  clog<<scientific<<"  Length: "<<LENGTH<<endl;
  clog<<scientific<<"  Current: "<<AMPERE<<endl;
  clog<<scientific<<"  Height: "<<HEIGHT<<endl;
  clog<<scientific<<"  Width: "<<WIDTH<<endl;
  clog<<"--------------------------"<<endl;
  clog<<endl; 
}


//---------------Main Function---------------------------

void usage()
{
  printf("\n");
  fprintf(stderr, "usage:  Saltation -i (name) -n (int) -w (int) [rtcfxyzRsTCMh] \n\
\n\
-n (int) \t number of particles (required)    \n\
-i (file)\t input file (required)\n(format: x y z vx vy vz angvx angvy angvz r temp color_int)\n\
\n OPTIONAL arguments\n\
-h   \t show this message  \n\n");

}				
int main(int argc, char** argv){
  FILE *fin=NULL;
  char infile[256];
  char ginfile[256];
  char outfile1[256];
  #ifdef WIND
  char outfile2[256];
  char outfile3[256];
  #endif
  #ifdef ELECF
  char outfile4[256];
  #endif
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int numberOfFiles = 0;
  double runningTime = 100;

  while((option_letter = getopt(argc, argv, "i:n:w:v:r:g:b:eph")) != -1 ){
  
    if (option_letter == 'i') {    
      strcpy(infile, optarg);
      usage_flag++;      
    }
    else if (option_letter == 'g') {    
      strcpy(ginfile, optarg);
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
    else if (option_letter == 'r') {    
      RHO_H = atof(optarg);
      usage_flag++;      
    }   
    else if (option_letter == 'b') {    
      numberOfFiles = atoi(optarg);
      usage_flag++;      
    }   
    else if (option_letter == 'h') {    
      usage();
      exit(1);      
    } 	    
  }   
     
  if (usage_flag < 4) {   
    usage();
    exit(1);      
  } 
  Body *Sand;
  Sand = new Body[N];
  #ifdef WIND
  Wind Hurricane;
  #endif
  double RateFile=1.0/5.0;
  double t=(numberOfFiles-1)*RateFile;
  double writeTime;
  if(t<0)t=0;
  if(t==0) writeTime = 0; else writeTime = t;//((numberOfFiles)*RateFile);

  Collider Aeolian;
  Aeolian.Init(1);
//Print Green Function
  #ifdef ELECF
  #ifdef PGF
  ofstream greenF;
  greenF.open("GreenFunction_500.txt");
  Aeolian.PrintGreenFunction(greenF);
  greenF.close();
  exit(0);
  #else
  Aeolian.ReadGreenFunction(ginfile);
  #endif
  #endif


  #ifdef WIND
  Hurricane.Init(Sand);
  #endif
//------------------Simulation------------------------

  Aeolian.ChargeGrid[0][167][0]=100000;
  Aeolian.ChargeGrid[1][167][0]=100000;
  Aeolian.ChargeGrid[2][167][0]=100000;
  Aeolian.ChargeGrid[3][167][0]=100000;
  Aeolian.ChargeGrid[4][167][0]=100000;
  Aeolian.ChargeGrid[5][167][0]=100000;
  Aeolian.ChargeGrid[6][167][0]=100000;
  Aeolian.ChargeGrid[7][167][0]=100000;
  Aeolian.ChargeGrid[0][167][1]=100000;
  Aeolian.ChargeGrid[1][167][1]=100000;
  Aeolian.ChargeGrid[2][167][1]=100000;
  Aeolian.ChargeGrid[3][167][1]=100000;
  Aeolian.ChargeGrid[4][167][1]=100000;
  Aeolian.ChargeGrid[5][167][1]=100000;
  Aeolian.ChargeGrid[6][167][1]=100000;
  Aeolian.ChargeGrid[7][167][1]=100000;
  int num=0;
   for(int j=0;j<Aeolian.Lyp;j++)
      for(int i=0;i<Aeolian.Lxp;i++)
        for(int k=0;k<Aeolian.Lzp;k++){
              Aeolian.smallArray[num]=i+Aeolian.Lxp*j+Aeolian.Lxp*Aeolian.Lyp*k;
              num++;
        }
  Aeolian.sizeSA=num;

  Aeolian.PotentialSolverParallel();
  ofstream fout;
  fout.open("AllesEF.txt");
   for(int j=0;j<Aeolian.Lyp;j++)
      for(int i=0;i<Aeolian.Lxp;i++)
        for(int k=0;k<Aeolian.Lzp;k++){
        Aeolian.Done[i][j][k];
        Aeolian.ElectricFieldReal(i,j,k);
    fout<<scientific<<i<<"\t"<<j<<"\t"<<k<<"\t"<<Aeolian.Ex[i][j][k]<<"\t"<<Aeolian.Ey[i][j][k]<<"\t"<<Aeolian.Ez[i][j][k]<<"\t"<<Aeolian.ChargeGrid[i][j][k]<<endl;
  }
  fout.close();
  return 0;
}
