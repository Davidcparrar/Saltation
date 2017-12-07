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
  int i,Nair=0;
  ofstream fout;
  ofstream fpout;
  ofstream fcout;
  ofstream flout;
  ofstream fmout;
  double x,y,z,Vx,Vy,Vz,R;
  double lowE, highE, lowET, highET,q=0;
  double En=0;
  int filei=static_cast<int>(100.0*WINDVEL);
  fout.open("INFO");
  printSaltationInfo(fout, runningTime);
  fout.close();
//-------------Adimentionalization-------------------
  LMT(D_MEAN);
//----------------Inicialization---------------------

  if(numberOfFiles!=0)
    sprintf (infile, "FILES/file%05d", (numberOfFiles-1));
    clog<<infile<<endl;
  if((fin=fopen(infile, "r"))==NULL){
    printf("Cannot Open File\n");
    exit(1);	
  }    
  long long int densityEle;
  RADIUS_MAX=0;
  Crandom r(0);
  R=4.e-4;
  densityEle = M_PI*R*R*RHO_H;
  double *altura;
  altura= new double [1000];
  int contAlt=0;


    
  //double dx=0.0002; double dy =0.006; double dz=0.0002; 
  for(i=0;i<N;i++){
    if(numberOfFiles!=0)
      fscanf(fin, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R, &lowE, &highE, &lowET, &highET, &q);
    else
      fscanf(fin, "%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R);
    densityEle = 4*M_PI*R*R*RHO_H;
    Sand[i].Init(x,y,z,Vx,Vy,Vz,R,densityEle,densityEle,q/ELEMENTARYCHARGE);
    /*if(dx>=0.9*WIDTH) dx = 0.0002;
    else dx+= 0.0002;
    if(dx==0.0002)dz+=0.0002;
    if(dz>=0.9*DEEP) dz = 0.0002;
    if(dz==0.0002 && dx==0.0002)dy+=0.0002;
     //clog<<scientific<<x<<"\t"<<y<<"\t"<<z<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<R<<endl;*/
    
  }
//Dimentionalization*/
  for(i=0;i<N;i++){
    Sand[i].Ad();
  }
  DT();
  Collider Aeolian;
  Aeolian.Init(N);
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


  runningTime=runningTime/TIME;
  sprintf(outfile1, "Energy%d.txt",filei);
  fout.open(outfile1,ios::app);

  //-----------Start Leap Frog--------  
  clog<<"Running ..."<<endl;
  Aeolian.Collide(Sand,HEIGHT);
    double sumAlt;
  #ifdef WIND
    Hurricane.ComputingWindForces(Sand,0);
  #endif

  for(i=0;i<N;i++) {Sand[i].UpdatePos(DELTAT);}   
  Aeolian.Collide(Sand,HEIGHT);
  #ifdef WIND
    Hurricane.ComputingWindForces(Sand,0);
  #endif
  clog<<"Running ..."<<endl;
  for(i=0;i<N;i++) {Sand[i].UpdateVel(DELTAT,0);}   

  int mia; 
  for(mia = 0;mia<1000;mia++)      
    altura[mia]=HEIGHT;
  clock_t timeCheck;
  timeCheck=clock();
  double qflux=0,qcharge=0,qchargeNonAbs=0;

  clog<<t<<"\t"<<writeTime<<"\t"<<numberOfFiles<<endl;

  while(t<=runningTime+DELTAT){
   if(t>=writeTime){ 
      sprintf (outfile1, "FILES/file%05d", numberOfFiles);
      fpout.open(outfile1);
      #ifdef WIND
        sprintf (outfile2, "PROFILES/Profile%05d", numberOfFiles);    
        fcout.open(outfile2);
        sprintf (outfile3, "CONCENTRATION/Concentration%05d", numberOfFiles);
        flout.open(outfile3); 
      #endif
      #ifdef ELECF
        sprintf (outfile4, "PROFILES/ChargeMap%05d", numberOfFiles);
        fmout.open(outfile4);             
      #endif
      for(i=0;i<N;i++){
        fpout<<scientific<<Sand[i].x*LENGTH<<"\t"<<Sand[i].y*LENGTH<<"\t"<<Sand[i].z*LENGTH<<"\t"<<Sand[i].Vx*LENGTH/TIME<<"\t"<<Sand[i].Vy*LENGTH/TIME<<"\t"<<Sand[i].Vz*LENGTH/TIME<<"\t"<<Sand[i].R*LENGTH<<"\t"<<0<<"\t"<<0<<"\t"<<Sand[i].H<<"\t"<<Sand[i].OH<<"\t"<<Sand[i].q*ELEMENTARYCHARGE*AMPERE*TIME<<endl;
        qflux+=Sand[i].Vx*Sand[i].m;
        En+=Sand[i].getE();
        qcharge+=abs(Sand[i].q);
        qchargeNonAbs+=Sand[i].q;
        #ifdef WIND
        if(Sand[i].y>=D_MEAN*Hurricane.height){
          Nair++;
        }
        #endif
      }
        
      #ifdef WIND
        Hurricane.PrintVelProfile(fcout);
        fcout.close();
        Hurricane.PrintConcentrationProfile(flout);
        flout.close();
        if(Nair==0) Nair=1;
      #else
        Nair=N;
      #endif
      #ifdef ELECF
        Aeolian.PrintCharge(fmout);
        fmout.close();
      #endif
      fout<<scientific<<t*TIME<<"\t"<<En*MASS*LENGTH*LENGTH/(TIME*TIME)<<"\t"<<qflux/(DEEP*WIDTH*DENSITY*sqrt((DENSITY/AIRDENSITY-1.)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))<<"\t"<<qcharge*AMPERE*TIME*ELEMENTARYCHARGE/Nair<<"\t"<<qchargeNonAbs*AMPERE*TIME*ELEMENTARYCHARGE<<"\t"<<sumAlt<<endl;
      clog<<scientific<<t*TIME<<"\t"<<En*MASS*LENGTH*LENGTH/(TIME*TIME)<<"\t"<<qflux/(DEEP*WIDTH*DENSITY*sqrt((DENSITY/AIRDENSITY-1.)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))<<"\t"<<qcharge*AMPERE*TIME*ELEMENTARYCHARGE/Nair<<"\t"<<qchargeNonAbs*AMPERE*TIME*ELEMENTARYCHARGE<<"\t"<<sumAlt<<endl;
      
      En=0;
      qflux=0;
      qcharge=0;
      qchargeNonAbs=0;
      Nair=0;
      writeTime+=RateFile/(TIME);
      fpout.close();
      numberOfFiles++;
    }
    //-----------Velocity Stormer Verlet--------

    #ifdef WIND
    sumAlt=0;
    for(i=0;i<N;i++) {Sand[i].UpdatePos(DELTAT);}   
    for(mia = 0;mia<1000;mia++)      
      sumAlt+=altura[mia];
    sumAlt/=1000.;  
    Aeolian.Collide(Sand,sumAlt);
    Hurricane.ComputingWindForces(Sand,t);
    altura[contAlt]=2.*D_MEAN*Hurricane.height;
    #else
    for(i=0;i<N;i++) {Sand[i].UpdatePos(DELTAT);}   
      Aeolian.Collide(Sand,0);
    #endif
    for(i=0;i<N;i++) {Sand[i].UpdateVel(DELTAT,sumAlt);}   
    t+=DELTAT;

    contAlt++;
     if(contAlt>=1000)contAlt=0;
  }

  fout.close();
  timeCheck = clock() - timeCheck;

  clog<<"Running time "<<runningTime<<endl;
  clog<<"Elapsed time "<<(static_cast<double>(timeCheck))/CLOCKS_PER_SEC<<endl;
  return 0;
}
