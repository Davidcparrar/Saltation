#include <iostream>
#include <cmath> 
#include<fstream> 
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include "Particle.h"

using namespace std;
double LIMIT;
class WindInfo{
public:
  double *ForceX;
  double *ForceY;
  double *Concentration;
  double *VelocityX;
  double *VelocityY;
  double *Charge;
  double *ChargeToArea;
  double *ChargeToMass;
  double *ChargeNonAbs;
  double *Mass;
  int *nParticles;
  int Box;
  int time2Average;
  int size;
  void init(int n);
  ~WindInfo(void);
  void reset(void);
  friend class Wind;
};

void WindInfo::init(int n){
  size = n;
  int i;
  ForceX = new double [size];
  ForceY = new double [size];
  Concentration = new double [size];
  VelocityX = new double [size];
  VelocityY = new double [size];
  nParticles = new int [size];
  Charge= new double [size];
  ChargeToArea= new double [size];
  ChargeToMass= new double [size];
  Mass= new double [size];
  ChargeNonAbs= new double [size];
  time2Average = 0;
  Box=5*D_MEAN/2;
  for(i=0;i<size;i++){
    ForceX[i]=0;
    ForceY[i]=0;
    Concentration[i]=0;
    VelocityX[i]=0;
    VelocityY[i]=0;
    nParticles[i]=0;
    Charge[i]=0;
    ChargeToArea[i]=0;
    ChargeToMass[i]=0;
    Mass[i]=0;
    ChargeNonAbs[i]=0;
  }
}

WindInfo::~WindInfo(){
  delete [] ForceX;
  delete [] ForceY;
  delete [] Concentration;
  delete [] VelocityX;
  delete [] VelocityY;
  delete [] nParticles;
  delete [] Charge;
  delete [] ChargeToArea;
  delete [] ChargeToMass;
  delete [] Mass;
  delete [] ChargeNonAbs;
}

void WindInfo::reset(void){
  int i=0;
  time2Average = 0;
  for(i=0;i<size;i++){
    ForceX[i]=0;
    ForceY[i]=0;
    Concentration[i]=0;
    VelocityX[i]=0;
    VelocityY[i]=0;
    nParticles[i]=0;
    Charge[i]=0;
    ChargeToArea[i]=0;
    ChargeToMass[i]=0;
    Mass[i]=0;
    ChargeNonAbs[i]=0;
  }
}
class List{
public:
  int done, head, nParticles;
public:
  List();
};

List::List(){
  head=-1;
  done=0;
}

class Wind{
  public:
  int SIZE;
  double *windVelocity;
  double *dwindVelocity;
  double *Tau;
  double *nTau;
  double *TauSingle;
  double *nTauSingle;
  double *histogram;
  double *previousWindVelocity;
  double *previousWindVelocity2;
  double *momentumTransfer;
  double SHEARAREA;
  double *shearGrain;
  int contShear;
  int *inty;
  int maxHeight;
  List *heads;
  WindInfo viento;
  double *Forcex;
  double *Forcey;
  double *LiftForce;
  int height;
  double Average;
  double OldAverage;
  int iteration;
  ~Wind(void);
  double VelocityProfile(double y, int height);
  double Cd(double Re);
  double dVelocityProfile(double y, double &Stress);
  void Init(Body *Particles);
  void ComputingWindForces(Body *Particles, double time);
  void PrintConcentrationProfile(ofstream &fout);
  void PrintVelProfile(ofstream &fout);
  void PrintForceProfile(ofstream &fout, Body *Particles);
  friend class Collider;
};

Wind::~Wind(void){
  delete [] windVelocity;
  delete [] dwindVelocity;
  delete [] Tau;
  delete [] nTau;
  delete [] TauSingle;
  delete [] nTauSingle;
  delete [] histogram;
  delete [] previousWindVelocity; 
  delete [] previousWindVelocity2; 
  delete [] momentumTransfer; 
  delete [] Forcex;
  delete [] Forcey;  
  delete [] LiftForce; 
  delete [] inty;
  delete [] heads;
  delete [] shearGrain;
}

double Wind::VelocityProfile(double y, int height){
  return WINDVEL/0.4*log((y-height*1.0)*30);
}

double Wind::Cd(double Re){
  return pow(pow(32.0/Re,2.0/3.0)+1,1.5);
}

double Wind::dVelocityProfile(double y, double &Stress){
  if(Stress<LIMIT){
    return WINDVEL*sqrt(1-Stress/(LIMIT))/(y*0.4);
  }
  else{
    return 0;
  }
}

void Wind::Init(Body *Particles){
  int i,j,k,n,set,setted,l;
  ofstream finfo,faverage,fvel;
  double ReAux;
  LIMIT=AIRDENSITY*WINDVEL*WINDVEL;
  SIZE = HEIGHT/(0.5*D_MEAN);
  SHEARAREA=1.0/(DEEP*WIDTH);
  viento.init(SIZE);
  windVelocity = new double [SIZE];
  dwindVelocity = new double [SIZE];
  Tau = new double [SIZE];
  nTau = new double [SIZE];
  TauSingle = new double [SIZE];
  nTauSingle = new double [SIZE];
  histogram = new double [SIZE];
  previousWindVelocity = new double [SIZE];
  previousWindVelocity2 = new double [SIZE];
  momentumTransfer = new double [SIZE];
  inty = new int [N];
  heads = new List [SIZE];
  Forcex = new double [N];
  Forcey = new double [N];
  LiftForce = new double [N];
  shearGrain = new double [SIZE];
  contShear = 0;

  for(i=0;i<SIZE;i++)
    histogram[i]=0;
  //Initializacion
  int iy=0;
  //Making histogram	
  for(i=0;i<N;i++){
    iy = (int) (Particles[i].y/(D_MEAN*0.5));
    histogram[iy]+=1;
  }
   set = 0;
   setted= 0;
   height=0;
  //Initial Velocity Profile
  finfo.open("InitialVelocityProfile.txt");
  for(i=0;i<SIZE;i++){
    windVelocity[i]=WINDVEL*log((i*0.5+0.5+1.0/30)/(1.0/30))/0.4;
    previousWindVelocity[i]=0;
    finfo<<fixed<<windVelocity[i]*LENGTH/TIME<<"\t"<<((i+1)+1.0/30)*LENGTH<<"\t"<<histogram[i]<<endl;
    Tau[i]=0.0; 
    nTau[i]=0.0; 
    heads[i].head = -1;
    heads[i].done = 0;
    heads[i].nParticles = 0;
    TauSingle[i]=0;
    nTauSingle[i]=0;
    shearGrain[i]=0;
  }

  finfo.close();
  // Calculating Tau
  Average=0;
  OldAverage=10000;
  double Waux=0;
  int newHeight=0;
  double OldCriteria=-1000;
  double NewCriteria=0;
  double Taux=0;
  clog<<"Size: "<<SIZE<<endl;
  clog<<"Stress Limit: "<<WINDVEL*WINDVEL*AIRDENSITY<<endl;
  clog<<"Speed Limit: "<<WINDVEL*0.1*LENGTH/TIME<<endl;
  faverage.open("VelocityAverage.txt");
  iteration=0;
  int UpLimit;
  int cell;
  int downLimit;
  double VelocityAux,VelocityAuy,VelAir,VelocityMag;
  OldAverage=0;
  Average=0;
  height=0;
  iteration = 0;
  maxHeight = 0;
  double ParamV=0;
  char outfile1[20];
  int OldHeight;
  double c=0.3;
  for(i=0;i<SIZE;i++){
    Tau[i]=0;
    TauSingle[i]=0;
    nTau[i]=1;
    heads[i].nParticles=0;
    heads[i].head=-1;
    windVelocity[i]=WINDVEL*log((i*0.5+0.5+1.0/30)/(1.0/30))/0.4;
    previousWindVelocity[i]=0;
    shearGrain[i]=0;
  }

   for(i=0;i<N;i++){
    inty[i] = (int) (Particles[i].y/(D_MEAN*0.5));

    if(inty[i]>maxHeight) maxHeight = inty[i];

    Particles[i].nextW = heads[inty[i]].head;
    heads[inty[i]].head = i; 
    heads[inty[i]].nParticles++;
  }
  set=0; 
  newHeight=0;
 
  for(j=0;j<10000;j++){ 

    set=0; 

    newHeight=0;
    for(i=0;i<SIZE;i++){
      previousWindVelocity[i] = windVelocity[i];
      windVelocity[i]=0;
      if(nTau[i]<0)
        dwindVelocity[i]= -1.0*WINDVEL*log(((i+1-newHeight)*0.5+1/30.0)/((i-newHeight)*0.5+1/30.0))*sqrt(-nTau[i])/0.4;
      else
        dwindVelocity[i]= WINDVEL*log(((i+1-newHeight)*0.5+1/30.0)/((i-newHeight)*0.5+1/30.0))*sqrt(nTau[i])/0.4;

      if(i==0)
        windVelocity[0]=dwindVelocity[0];
      else 
        windVelocity[i]=c*(dwindVelocity[i]+windVelocity[i-1])+(1-c)*previousWindVelocity[i];

      if(windVelocity[i]<0)windVelocity[i]=0;

      if(i<=maxHeight)
      Average+=windVelocity[i];

      if(windVelocity[i]<=0.1*WINDVEL && set==0){
        newHeight++;
      }
      else{
        set = 1;
      }

    }
    height=newHeight;


    for(i=0;i<SIZE;i++){
      Tau[i]=0;
      TauSingle[i]=0;
      nTau[i]=0;
    }
    for(k=0;k<N;k++){
      Forcex[k]=0;
      VelAir=log(((Particles[k].y/(0.5*D_MEAN)-newHeight)*0.5+1/30.0)/((inty[k]-newHeight)*0.5+1/30.0))/log(((inty[k]+1-newHeight)*0.5+1/30.0)/((inty[k]-newHeight)*0.5+1/30.0))*(windVelocity[inty[k]+1]-windVelocity[inty[k]])+windVelocity[inty[k]];
     if ((VelAir != VelAir) || (VelAir>fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]))){ 
       VelAir=fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]);
     }
      VelocityAux = Particles[k].Vx - VelAir;
      VelocityMag= sqrt(VelocityAux*VelocityAux+Particles[k].Vy*Particles[k].Vy);
      ReAux = AIRDENSITY*VelocityMag*Particles[k].R*2/AIRVISCOSITY;
      if(VelocityMag!=0)
        Forcex[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAux*VelocityMag;
      TauSingle[inty[k]]+=Forcex[k]*SHEARAREA;
    }
 
    for(k=SIZE;k>=0;k--){
      if(k<SIZE-1)
        Tau[k]+=TauSingle[k]+Tau[k+1];
      else
        Tau[k]+=TauSingle[k];
      nTau[k]=(1-Tau[k]/LIMIT);
    }

    faverage<<fixed<<j<<"\t"<<Average/WINDVEL<<"\t"<<OldAverage/WINDVEL<<"\t"<<abs(Average-OldAverage)/WINDVEL<<"\t"<<newHeight<<endl;
    if(j<10){
    sprintf(outfile1, "VelProfile%05d", iteration);
    fvel.open(outfile1);
    PrintVelProfile(fvel);
    fvel.close();
    }
    
    OldCriteria=NewCriteria;
    OldAverage=Average;
    Average=0;
    iteration ++; 

  }

  sprintf(outfile1, "VelProfile%05d", iteration);
  fvel.open(outfile1);
  PrintVelProfile(fvel);
  fvel.close();
  height=newHeight;
  finfo.open("InitialStressProfile.txt");
  for(i=0;i<SIZE;i++)
    finfo<<fixed<<(i)<<"\t"<<Tau[i]*MASS/(TIME*TIME*LENGTH)<<"\t"<<windVelocity[i]*LENGTH/TIME<<endl;   
  finfo.close();

  finfo.open("InitialVelocityStressProfile.txt");
  for(i=0;i<SIZE;i++){
    finfo<<fixed<<(i+1.0/30)*LENGTH<<"\t"<<windVelocity[i]*LENGTH/TIME<<"\t"<<previousWindVelocity[i]*LENGTH/TIME<<endl;
  }
  finfo.close();
}

void Wind::ComputingWindForces(Body *Particles, double time){
  int i,j,k,n,set=0;
  int  newHeight=-1;
  double Waux=0;
  double OldCriteria=-1000;
  double NewCriteria=0;
  double setted = -1;
  double ReAux;
  double Taux=0;
  int OldHeight=newHeight;
  double VelocityAux,VelocityAuy,VelAir,VelocityMag;
  iteration = 0;
  char outfile1[20];
  int downlimit=0;
  int altu=0;
  ofstream faverage,fvel;
  OldAverage=-10000;
  //faverage.open("TestVelocityAverage.txt");
  double c=0.99;
  for(i=0;i<SIZE;i++){
    Tau[i]=0.0; 
    heads[i].head = -1;
    heads[i].done = 0;
    heads[i].nParticles = 0;
    TauSingle[i]=0;
    nTauSingle[i]=0;
  }

   for(i=0;i<N;i++){
    inty[i] = (int) (Particles[i].y/(D_MEAN*0.5));

    if(inty[i]>maxHeight) maxHeight = inty[i];

    Particles[i].nextW = heads[inty[i]].head;
    heads[inty[i]].head = i; 
    heads[inty[i]].nParticles++;
    
  }
 
  for(j=0;j<4;j++){ 

    set=0; 
    newHeight=0;
    for(i=0;i<SIZE;i++){
      previousWindVelocity[i] = windVelocity[i];
      windVelocity[i]=0;

      if(i==0)
        windVelocity[0]=dwindVelocity[0];
      else{
        if(nTau[i]<0)
          windVelocity[i]=(1-c)*(-1.0*WINDVEL*log(((i+1-newHeight)*0.5+1/30.0)/((i-newHeight)*0.5+1/30.0))*sqrt(-nTau[i])/0.4+windVelocity[i-1])+(c)*previousWindVelocity[i];
        else
          windVelocity[i]=(1-c)*(WINDVEL*log(((i+1-newHeight)*0.5+1/30.0)/((i-newHeight)*0.5+1/30.0))*sqrt(nTau[i])/0.4+windVelocity[i-1])+(c)*previousWindVelocity[i];
      }

      if(windVelocity[i]<0)windVelocity[i]=0;

      if(i<=maxHeight)
      Average+=windVelocity[i];

      if(windVelocity[i]<=0.1*WINDVEL && set==0){
        newHeight++;
      }
      else{
        set = 1;
      }

    }
    height=newHeight;

    for(i=0;i<SIZE;i++){
      Tau[i]=0;
      TauSingle[i]=0;
      nTau[i]=0;
    }
     for(k=0;k<N;k++){
      Forcex[k]=0;
      VelAir=log(((Particles[k].y/(0.5*D_MEAN)-newHeight)*0.5+1/30.0)/((inty[k]-newHeight)*0.5+1/30.0))/log(((inty[k]+1-newHeight)*0.5+1/30.0)/((inty[k]-newHeight)*0.5+1/30.0))*(windVelocity[inty[k]+1]-windVelocity[inty[k]])+windVelocity[inty[k]];
     if ((VelAir != VelAir) || (VelAir>fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]))){ 
       VelAir=fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]);
     }
      VelocityAux = Particles[k].Vx - VelAir;
      VelocityAuy = Particles[k].Vy;
      VelocityMag= sqrt(VelocityAux*VelocityAux+VelocityAuy*VelocityAuy);
      ReAux = AIRDENSITY*VelocityMag*Particles[k].R*2/AIRVISCOSITY;
      if(VelocityMag!=0){
        Forcex[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAux*VelocityMag;
        Forcey[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAuy*VelocityMag;
      }
      TauSingle[inty[k]]+=Forcex[k]*SHEARAREA;
    }

    for(k=SIZE;k>=0;k--){
      if(k<SIZE-1)
        Tau[k]+=TauSingle[k]+Tau[k+1];
      else
        Tau[k]+=TauSingle[k];
      nTau[k]=(1-Tau[k]/LIMIT);
    }
    if(abs(Average-OldAverage)/WINDVEL<0.0001)
    j=100003;

    NewCriteria=abs(Average-OldAverage)/WINDVEL;
    OldAverage=Average;
    Average=0;
    iteration ++; 
    

  }
  height=newHeight;

  for(i=0;i<N;i++){
    LiftForce[i]=0;
    if(inty[i]>=height){
      Particles[i].Fy+=Forcey[i];
      Particles[i].Fx+=Forcex[i];
      shearGrain[inty[i]]+=Forcex[i];
      viento.ForceX[inty[i]]+=Forcex[i];
      viento.ForceY[inty[i]]+=Forcey[i];
      viento.Concentration[inty[i]]+=Particles[i].Vol;
      viento.VelocityX[inty[i]]+=Particles[i].Vx;
      viento.VelocityY[inty[i]]+=Particles[i].Vy;
    }
    viento.nParticles[inty[i]]++; 
    viento.Charge[inty[i]]+=abs(Particles[i].q*ELEMENTARYCHARGE);
    viento.ChargeToArea[inty[i]]+=abs(Particles[i].q*ELEMENTARYCHARGE/(Particles[i].R*Particles[i].R));
    viento.ChargeToMass[inty[i]]+=abs(Particles[i].q*ELEMENTARYCHARGE/(Particles[i].m));
  }
  height=newHeight;
  contShear++;
  viento.time2Average++;

}

void Wind::PrintConcentrationProfile(ofstream &fout){
  int i;
  double auxC=0,auxShearX=0,auxShearY=0,auxVelX=0,auxVelY=0,auxCharge=0,auxChargeToMass=0,auxChargeToArea=0;
  int cont = 0;
  int part=0;
  for(i=0;i<viento.size;i++){
    cont++;
    auxC+=viento.Concentration[i];
    auxShearX+=viento.ForceX[i];
    auxShearY+=viento.ForceY[i];
    auxVelX+=viento.VelocityX[i];
    auxVelY+=viento.VelocityY[i];
    part+=viento.nParticles[i];
    auxCharge+=viento.Charge[i];
    auxChargeToMass+=viento.ChargeToMass[i];
    auxChargeToArea+=viento.ChargeToArea[i];
    if(cont%5==0){
      if(part!=0){
        auxVelY=auxVelY/part;
        auxVelX=auxVelX/part;
        auxCharge=auxCharge/part;
        auxChargeToMass=auxChargeToMass/part;
        auxChargeToArea=auxChargeToArea/part;

      }
     fout<<scientific<<i-4<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)<<"\t"<<auxShearX*MASS*LENGTH/(TIME*TIME*viento.time2Average*viento.Box)<<"\t"<<auxShearY*MASS*LENGTH/(viento.time2Average*viento.Box*TIME*TIME)<<"\t"<<auxVelX*LENGTH/TIME<<"\t"<<auxVelY*LENGTH/TIME<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)*auxVelX*LENGTH/TIME<<"\t"<<auxCharge*AMPERE*TIME<<"\t"<<auxChargeToMass/(MASS)*AMPERE*TIME<<"\t"<<auxChargeToArea/(LENGTH*LENGTH)*AMPERE*TIME<<endl;
      part=0;
      auxC=0;auxShearX=0;auxShearY=0;auxVelX=0;auxVelY=0;auxCharge=0;auxChargeToMass=0;auxChargeToArea=0;
     }
     
  }
  viento.reset();
 
}
void Wind::PrintVelProfile(ofstream &fout){
  int i;
  for(i=0;i<SIZE;i++)
    fout<<scientific<<windVelocity[i]*LENGTH/TIME<<"\t"<<previousWindVelocity[i]*LENGTH/TIME<<"\t"<<dwindVelocity[i]*LENGTH/TIME<<"\t"<<Tau[i]*MASS/(LENGTH*TIME*TIME)<<"\t"<<TauSingle[i]*MASS/(LENGTH*TIME*TIME)<<"\t"<<Tau[i]*MASS/(LENGTH*TIME*TIME)/SHEARAREA<<"\t"<<i<<"\t"<<height<<"\t"<<shearGrain[i]*MASS*LENGTH/(TIME*TIME)/(contShear*D_MEAN*0.5)<<endl;
}

void Wind::PrintForceProfile(ofstream &fout, Body *Particles){
  int i;
  for(i=0;i<N;i++)
    fout<<scientific<<i<<"\t"<<Forcex[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<Forcey[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<LiftForce[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<Particles[i].m*GRAVITY*MASS*LENGTH/(TIME*TIME)<<endl;
}
