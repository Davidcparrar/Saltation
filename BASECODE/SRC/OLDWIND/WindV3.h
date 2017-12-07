#include <iostream>
#include <cmath> 
#include<fstream> 
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include"Particle.h"

using namespace std;
double LIMIT;
class WindInfo{
public:
  double *ForceX;
  double *ForceY;
  double *ForceZ;
  double *Concentration;
  double *VelocityX;
  double *VelocityY;
  double *Charge;
  double *ChargeToArea;
  double *ChargeToMass;
  double *ChargeNonAbs;
  double *Mass;
  double *VelProf;
  #ifdef BINDIST
  long int *Small;
  long int *Big;
  long int *Total;
  long int *SmallCh;
  long int *BigCh;
  long int *TotalCh;
  double *ChargeBig;
  double *ChargeSmall;
  #endif
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
  ForceZ = new double [size];
  Concentration = new double [size];
  VelocityX = new double [size];
  VelocityY = new double [size];
  nParticles = new int [size];
  Charge= new double [size];
  ChargeToArea= new double [size];
  ChargeToMass= new double [size];
  Mass= new double [size];
  ChargeNonAbs= new double [size];
  VelProf= new double [size];
  #ifdef BINDIST
  Small= new long int [size];
  Big= new long int [size];
  Total= new long int [size];
  SmallCh= new long int [size];
  BigCh= new long int [size];
  TotalCh= new long int [size];
  ChargeSmall= new double [size];
  ChargeBig= new double [size];
  #endif
  time2Average = 0;
  Box=10*D_MEAN/2;
  for(i=0;i<size;i++){
    ForceX[i]=0;
    ForceY[i]=0;
    ForceZ[i]=0;
    Concentration[i]=0;
    VelocityX[i]=0;
    VelocityY[i]=0;
    nParticles[i]=0;
    Charge[i]=0;
    ChargeToArea[i]=0;
    ChargeToMass[i]=0;
    Mass[i]=0;
    ChargeNonAbs[i]=0;
    VelProf[i]=0;
    #ifdef BINDIST
    Small[i]=0;
    Big[i]=0;
    Total[i]=0;
    SmallCh[i]=0;
    BigCh[i]=0;
    TotalCh[i]=0;
    ChargeSmall[i]=0;
    ChargeBig[i]=0;
    #endif
  }
}

WindInfo::~WindInfo(){
  delete [] ForceX;
  delete [] ForceY;
  delete [] ForceZ;
  delete [] Concentration;
  delete [] VelocityX;
  delete [] VelocityY;
  delete [] nParticles;
  delete [] Charge;
  delete [] ChargeToArea;
  delete [] ChargeToMass;
  delete [] Mass;
  delete [] ChargeNonAbs;
  delete [] VelProf;
  #ifdef BINDIST
  delete [] Small;
  delete [] Big;
  delete [] SmallCh;
  delete [] BigCh;
  delete [] ChargeSmall;
  delete [] ChargeBig;
  #endif
}

void WindInfo::reset(void){
  int i=0;
  time2Average = 0;
  for(i=0;i<size;i++){
    ForceX[i]=0;
    ForceY[i]=0;
    ForceZ[i]=0;
    Concentration[i]=0;
    VelocityX[i]=0;
    VelocityY[i]=0;
    nParticles[i]=0;
    Charge[i]=0;
    ChargeToArea[i]=0;
    ChargeToMass[i]=0;
    Mass[i]=0;
    ChargeNonAbs[i]=0;
    VelProf[i]=0;
    #ifdef BINDIST
    Small[i]=0;
    Big[i]=0;
    SmallCh[i]=0;
    BigCh[i]=0;
    ChargeSmall[i]=0;
    ChargeBig[i]=0;
    #endif
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
  double SHEARAREA;
  double *shearGrain;
  int contShear;
  int *inty;
  int maxHeight;
  List *heads;
  WindInfo viento;
  double *Forcex;
  double *Forcey;
  double *Forcez;
  double *LiftForce;
  int height;
  double Average;
  double OldAverage;
  int iteration;
  ~Wind(void);
  #ifdef BINDIST
  double SIZE1;
  double SIZE2;
  #endif
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
  delete [] Forcex;
  delete [] Forcey;  
  delete [] Forcez;  
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
  int i,j,k,set;
  ofstream finfo,faverage,fvel;
  double ReAux;
  LIMIT=AIRDENSITY*WINDVEL*WINDVEL;
  SIZE = HEIGHT/(D_MEAN);
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
  inty = new int [N];
  heads = new List [SIZE];
  Forcex = new double [N];
  Forcey = new double [N];
  Forcez = new double [N];
  LiftForce = new double [N];
  shearGrain = new double [SIZE];
  contShear = 0;

  #ifdef BINDIST
  SIZE1=0.5;
  SIZE2=1;
  #endif

  for(i=0;i<SIZE;i++)
    histogram[i]=0;
  //Initialization
  int iy=0;
  //Making histogram	
  for(i=0;i<N;i++){
    iy = static_cast<int>(Particles[i].y/(D_MEAN));
    histogram[iy]+=1;
  }
   set = 0;
   height=0;
  //Initial Velocity Profile
  #ifdef DATAW
  finfo.open("InitialVelocityProfile.txt");
  #endif
  for(i=0;i<SIZE;i++){
    windVelocity[i]=WINDVEL*log((i+1+1.0/30)/(1.0/30))/0.4;
    previousWindVelocity[i]=0;
    #ifdef DATAW
    finfo<<fixed<<windVelocity[i]*LENGTH/TIME<<"\t"<<((i+1)+1.0/30)*LENGTH<<"\t"<<histogram[i]<<endl;
    #endif
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
  int newHeight=0;
  clog<<"Size: "<<SIZE<<endl;
  clog<<"Stress Limit: "<<WINDVEL*WINDVEL*AIRDENSITY<<endl;
  clog<<"Speed Limit: "<<WINDVEL*0.1*LENGTH/TIME<<endl;
  #ifdef DATAW
  faverage.open("VelocityAverage.txt");
  #endif
  iteration=0;
  double VelocityAux,VelocityAuy,VelocityAuz,VelAir,VelocityMag;
  OldAverage=0;
  Average=0;
  height=0;
  iteration = 0;
  maxHeight = 0;
  #ifdef DATAW
  char outfile1[20];
  #endif
  double c=0.3;
  for(i=0;i<SIZE;i++){
    Tau[i]=0;
    TauSingle[i]=0;
    nTau[i]=1;
    heads[i].nParticles=0;
    heads[i].head=-1;
    windVelocity[i]=WINDVEL*log((i+1+1.0/30)/(1.0/30))/0.4;
    previousWindVelocity[i]=0;
    shearGrain[i]=0;
  }

   for(i=0;i<N;i++){
    inty[i] =  static_cast<int>(Particles[i].y/(D_MEAN));

    if(inty[i]>maxHeight) maxHeight = inty[i];

    Particles[i].nextW = heads[inty[i]].head;
    heads[inty[i]].head = i; 
    heads[inty[i]].nParticles++;
    Forcex[i] = 0.;
    Forcey[i] = 0.;
    Forcez[i] = 0.;
    LiftForce[i] = 0.;
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
        dwindVelocity[i]= -1.0*WINDVEL*log(((i+1-newHeight)+1/30.0)/((i-newHeight)+1/30.0))*sqrt(-nTau[i])/0.4;
      else
        dwindVelocity[i]= WINDVEL*log(((i+1-newHeight)+1/30.0)/((i-newHeight)+1/30.0))*sqrt(nTau[i])/0.4;

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
      VelAir=log(((Particles[k].y/(D_MEAN)-newHeight)+1/30.0)/((inty[k]-newHeight)+1/30.0))/log(((inty[k]+1-newHeight)+1/30.0)/((inty[k]-newHeight)+1/30.0))*(windVelocity[inty[k]+1]-windVelocity[inty[k]])+windVelocity[inty[k]];
     if ((VelAir != VelAir) || (VelAir>fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]))){ 
       VelAir=fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]);
     }
      VelocityAux = Particles[k].Vx - VelAir;
      VelocityAuy = Particles[k].Vy;
      VelocityAuz = Particles[k].Vz;
      VelocityMag= sqrt(VelocityAux*VelocityAux+VelocityAuy*VelocityAuy+VelocityAuz*VelocityAuz);
      ReAux = AIRDENSITY*VelocityMag*Particles[k].R*2/AIRVISCOSITY;
      if(VelocityMag!=0){
        Forcex[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAux*VelocityMag;
        Forcey[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAuy*VelocityMag;
        Forcez[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAuz*VelocityMag;
      }
      TauSingle[inty[k]]+=Forcex[k]*SHEARAREA;
    }
 
    for(k=SIZE-1;k>=0;k--){
      if(k<SIZE-1)
        Tau[k]+=TauSingle[k]+Tau[k+1];
      else
        Tau[k]+=TauSingle[k];
      nTau[k]=(1-Tau[k]/LIMIT);
    }
    #ifdef DATAW
    faverage<<fixed<<j<<"\t"<<Average/WINDVEL<<"\t"<<OldAverage/WINDVEL<<"\t"<<abs(Average-OldAverage)/WINDVEL<<"\t"<<newHeight<<endl;
    if(j<10){
    sprintf(outfile1, "VelProfile%05d", iteration);
    fvel.open(outfile1);
    PrintVelProfile(fvel);
    fvel.close();
    }
    #endif  
    OldAverage=Average;
    Average=0;
    iteration ++; 

  }
  #ifdef DATAW
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
  faverage.close();
  #endif
}

void Wind::ComputingWindForces(Body *Particles, double time){
  int i,j,k,set=0;
  int  newHeight=-1;
  double ReAux;
  double VelocityAux,VelocityAuy,VelocityAuz,VelAir,VelocityMag;
  iteration = 0;
  #ifdef DATAW
  char outfile1[20];
  #endif
  ofstream faverage,fvel;
  OldAverage=-10000;
  #ifdef DATAW
  faverage.open("TestVelocityAverage.txt");
  #endif
  double c=0.99;
  for(i=0;i<SIZE;i++){
    Tau[i]=0.0; 
    heads[i].head = -1;
    heads[i].done = 0;
    heads[i].nParticles = 0;
    TauSingle[i]=0;
    nTauSingle[i]=0;
    viento.VelProf[i]+=windVelocity[i];
  }

   for(i=0;i<N;i++){
    inty[i] = static_cast<int>(Particles[i].y/(D_MEAN));

    if(inty[i]>maxHeight) maxHeight = inty[i];

    Particles[i].nextW = heads[inty[i]].head;
    heads[inty[i]].head = i; 
    heads[inty[i]].nParticles++;
    Forcex[i] = 0.;
    Forcey[i] = 0.;
    Forcez[i] = 0.;
    LiftForce[i] = 0.;
  }
 
  for(j=0;j<20;j++){ 

    set=0; 
    newHeight=0;
    for(i=0;i<SIZE;i++){
      previousWindVelocity[i] = windVelocity[i];
      windVelocity[i]=0;

      if(i==0)
        windVelocity[0]=dwindVelocity[0];
      else{
        if(nTau[i]<0)
          windVelocity[i]=(1-c)*(-1.0*WINDVEL*log(((i+1-newHeight)+1/30.0)/((i-newHeight)+1/30.0))*sqrt(-nTau[i])/0.4+windVelocity[i-1])+(c)*previousWindVelocity[i];
        else
          windVelocity[i]=(1-c)*(WINDVEL*log(((i+1-newHeight)+1/30.0)/((i-newHeight)+1/30.0))*sqrt(nTau[i])/0.4+windVelocity[i-1])+(c)*previousWindVelocity[i];
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
      VelAir=log(((Particles[k].y/(D_MEAN)-newHeight)+1/30.0)/((inty[k]-newHeight)+1/30.0))/log(((inty[k]+1-newHeight)+1/30.0)/((inty[k]-newHeight)+1/30.0))*(windVelocity[inty[k]+1]-windVelocity[inty[k]])+windVelocity[inty[k]];
     if ((VelAir != VelAir) || (VelAir>fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]))){ 
       VelAir=fmax(windVelocity[inty[k]+1],windVelocity[inty[k]]);
     }
      VelocityAux = Particles[k].Vx - VelAir;
      VelocityAuy = Particles[k].Vy;
      VelocityAuz = Particles[k].Vz;
      VelocityMag= sqrt(VelocityAux*VelocityAux+VelocityAuy*VelocityAuy+VelocityAuz*VelocityAuz);
      ReAux = AIRDENSITY*VelocityMag*Particles[k].R*2/AIRVISCOSITY;
      if(VelocityMag!=0){
        Forcex[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAux*VelocityMag;
        Forcey[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAuy*VelocityMag;
        Forcez[k]=-0.5*M_PI*(Particles[k].R)*(Particles[k].R)*AIRDENSITY*Cd(ReAux)*VelocityAuz*VelocityMag;
      }
      TauSingle[inty[k]]+=Forcex[k]*SHEARAREA;
    }

    for(k=SIZE-1;k>=0;k--){
      if(k<SIZE-1)
        Tau[k]+=TauSingle[k]+Tau[k+1];
      else
        Tau[k]+=TauSingle[k];
      nTau[k]=(1-Tau[k]/LIMIT);
    }
    if(abs(Average-OldAverage)/WINDVEL<0.0001)
    j=100003;

    OldAverage=Average;
    Average=0;
    iteration ++; 
    

  }
  height=newHeight;

  for(i=0;i<N;i++){
    LiftForce[i]=0;
    if(inty[i]>=height){
      Particles[i].Fz+=Forcez[i];
      Particles[i].Fy+=Forcey[i];
      Particles[i].Fx+=Forcex[i];
      shearGrain[inty[i]]+=Forcex[i];
      viento.ForceX[inty[i]]+=Forcex[i];
      viento.ForceY[inty[i]]+=Forcey[i];
      viento.ForceZ[inty[i]]+=Forcez[i];
      viento.Concentration[inty[i]]+=Particles[i].Vol;
      viento.VelocityX[inty[i]]+=Particles[i].Vx;
      viento.VelocityY[inty[i]]+=Particles[i].Vy;
    }
      #ifdef BINDIST
      if(Particles[i].R<0.75){
        viento.Small[inty[i]]++;
        viento.ChargeSmall[inty[i]]+=(Particles[i].q);
        if(Particles[i].q!=0){
          viento.SmallCh[inty[i]]++;
        }
      }
      else if(Particles[i].R>0.75){
        viento.Big[inty[i]]++;
        viento.ChargeBig[inty[i]]+=(Particles[i].q);
        if(Particles[i].q!=0){
          viento.BigCh[inty[i]]++;
        }

      }
      #endif
    viento.nParticles[inty[i]]++; 
    viento.Charge[inty[i]]+=abs(Particles[i].q);
    viento.ChargeToArea[inty[i]]+=abs(Particles[i].q/(Particles[i].R*Particles[i].R));
    viento.ChargeToMass[inty[i]]+=abs(Particles[i].q/(Particles[i].m));
    viento.Mass[inty[i]]+=Particles[i].m;
    viento.ChargeNonAbs[inty[i]]+=Particles[i].q;
  }
  height=newHeight;
  contShear++;
  viento.time2Average=viento.time2Average+1;

}

void Wind::PrintConcentrationProfile(ofstream &fout){
  int i;
  double auxC=0,auxShearX=0,auxShearY=0,auxVelX=0,auxVelY=0,auxCharge=0,auxChargeToMass=0,auxChargeToArea=0,auxMass=0,auxChargeNonAbs=0,auxChargeNonAbs2=0,auxVelProf=0;
  #ifdef BINDIST
  double small=0, big=0, smallCh=0, bigCh=0,ChargeB=0,ChargeS=0;
  #endif
  int cont = 0;
  double part=0;
  for(i=0;i<viento.size;i++){
    cont++;
    auxC+=viento.Concentration[i];
    auxShearX+=viento.ForceX[i];
    auxShearY+=viento.ForceY[i];
    auxVelX+=viento.VelocityX[i];
    auxVelY+=viento.VelocityY[i];
    part+=static_cast<double>(viento.nParticles[i]);
    auxCharge+=viento.Charge[i];
    auxChargeToMass+=viento.ChargeToMass[i];
    auxMass+=viento.Mass[i];
    auxChargeToArea+=viento.ChargeToArea[i];
    auxChargeNonAbs+=viento.ChargeNonAbs[i];
    auxVelProf+=viento.VelProf[i];
    #ifdef BINDIST
    small+=static_cast<double>(viento.Small[i]);
    big+=static_cast<double>(viento.Big[i]);
    smallCh+=static_cast<double>(viento.SmallCh[i]);
    bigCh+=static_cast<double>(viento.BigCh[i]);
    ChargeB+=viento.ChargeBig[i];
    ChargeS+=viento.ChargeSmall[i];
    #endif
    if(cont%5==0){
      if(part!=0){
        auxVelY=auxVelY/part;
        auxVelX=auxVelX/part;
        auxChargeNonAbs2=auxChargeNonAbs/auxMass;
        //auxMass=auxCharge/auxMass;
        //auxCharge=auxCharge/part;
        auxChargeToMass=auxChargeToMass/part;
        auxChargeToArea=auxChargeToArea/part;
      }
    #ifdef BINDIST
    fout<<scientific<<(i-4)*D_MEAN*LENGTH<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)<<"\t"<<auxShearX*MASS*LENGTH/(TIME*TIME*viento.time2Average*viento.Box)<<"\t"<<auxShearY*MASS*LENGTH/(viento.time2Average*viento.Box*TIME*TIME)<<"\t"<<auxVelX*LENGTH/TIME<<"\t"<<auxVelY*LENGTH/TIME<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)*auxVelX*LENGTH/TIME<<"\t"<<auxChargeNonAbs*ELEMENTARYCHARGE*AMPERE*TIME<<"\t"<<auxMass*MASS<<"\t"<<auxChargeToArea*ELEMENTARYCHARGE/(LENGTH*LENGTH)*AMPERE*TIME<<"\t"<<part/viento.time2Average<<"\t"<<auxChargeNonAbs2*ELEMENTARYCHARGE/(MASS)*AMPERE*TIME<<"\t"<<auxVelProf*LENGTH/(TIME*viento.time2Average*viento.Box)<<"\t"<<small/viento.time2Average<<"\t"<<big/viento.time2Average<<"\t"<<smallCh/viento.time2Average<<"\t"<<bigCh/viento.time2Average<<"\t"<<ChargeS*ELEMENTARYCHARGE/viento.time2Average<<"\t"<<ChargeB*ELEMENTARYCHARGE/viento.time2Average<<endl;   
    #else
    fout<<scientific<<(i-4)*D_MEAN*LENGTH<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)<<"\t"<<auxShearX*MASS*LENGTH/(TIME*TIME*viento.time2Average*viento.Box)<<"\t"<<auxShearY*MASS*LENGTH/(viento.time2Average*viento.Box*TIME*TIME)<<"\t"<<auxVelX*LENGTH/TIME<<"\t"<<auxVelY*LENGTH/TIME<<"\t"<<auxC/(viento.time2Average*viento.Box*DEEP*WIDTH)*auxVelX*LENGTH/TIME<<"\t"<<auxChargeNonAbs*ELEMENTARYCHARGE*AMPERE*TIME<<"\t"<<auxMass*MASS<<"\t"<<auxChargeToArea*ELEMENTARYCHARGE/(LENGTH*LENGTH)*AMPERE*TIME<<"\t"<<part/viento.time2Average<<"\t"<<auxChargeNonAbs2*ELEMENTARYCHARGE/(MASS)*AMPERE*TIME<<"\t"<<auxVelProf*LENGTH/(TIME*viento.time2Average*viento.Box)<<endl;
      #endif
      part=0;
      auxC=0;auxShearX=0;auxShearY=0;auxVelX=0;auxVelY=0;auxCharge=0;auxChargeToMass=0;auxChargeToArea=0,auxMass=0,auxChargeNonAbs=0,auxChargeNonAbs2=0;auxVelProf=0;
      #ifdef BINDIST
      small=0; big=0; smallCh=0; bigCh=0; ChargeB=0; ChargeS=0;
      #endif
     }
     
  }
  viento.reset();
 
}
void Wind::PrintVelProfile(ofstream &fout){
  int i;
  for(i=0;i<SIZE;i++)
    fout<<scientific<<windVelocity[i]*LENGTH/TIME<<"\t"<<previousWindVelocity[i]*LENGTH/TIME<<"\t"<<dwindVelocity[i]*LENGTH/TIME<<"\t"<<Tau[i]*MASS/(LENGTH*TIME*TIME)<<"\t"<<TauSingle[i]*MASS/(LENGTH*TIME*TIME)<<"\t"<<Tau[i]*MASS/(LENGTH*TIME*TIME)/SHEARAREA<<"\t"<<i<<"\t"<<height<<"\t"<<shearGrain[i]*MASS*LENGTH/(TIME*TIME)/(contShear*D_MEAN)<<endl;
}

void Wind::PrintForceProfile(ofstream &fout, Body *Particles){
  int i;
  for(i=0;i<N;i++)
    fout<<scientific<<i<<"\t"<<Forcex[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<Forcey[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<LiftForce[i]*MASS*LENGTH/(TIME*TIME)<<"\t"<<Particles[i].m*GRAVITY*MASS*LENGTH/(TIME*TIME)<<endl;
}
