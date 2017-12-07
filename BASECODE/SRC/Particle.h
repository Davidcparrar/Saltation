#include"Parameters.h"
#include<stdlib.h>
#include<stdio.h>
#include<fstream>

using namespace std;

class Body{
public:
  double m,R; 
  double x,y,z,Vx,Vy,Vz,Fx,Fy,Fz,Foldx,Foldy,Foldz,Vol;
  long int q;
  long long int lowE, highE, lowET, highET;
  long long int H, OH, H_0, OH_0;
  int nextx,nexty,nextz,nextxCh,nextyCh,nextzCh,nextW;
  double W1,W2,W3,W4,W5,W6;
  double a,b,c,d,e,f,g,h;
  double At,WallFlux;
  double Vel[100];
  double VelAverage;
  int counter;
  int Checked;
  int Charged;
  int Bed;
  double EFx, EFy;
public:
  void Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0, long long int lowE0, long long int highE0, long long int lowET0, long long int highET0);
  void Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0);
  void Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0, long long int H_0, long long int OH_0, long long int q0);
  void Ad(void);
  void UpdatePos(double dt);
  void UpdateVel(double dt, double alt);
  void Plot(ofstream &fout,double t);
  void Draw(void);
  void setCharge(long long int deltaH, long long int deltatL);
  void IonTransfer(long long int deltaH, long long int deltatL);
  void IonTransferVolume(double deltaH);
  void neutralization(long long int deltaL);
  void ElectronTransfer(double dt);
  void setChargeWall(void);
  void neutralizationIon(long long int Ti);
  void setChargeWallIon(void);
  double getR(){return R;};
  double getM(){return m;};
  double getE(){return m*((Vx*Vx+Vy*Vy+Vz*Vz)/2.0/*+GRAVITY*y*/);};
  double getAverageVel(void);
  friend class Collider;
  friend class Wind;
};

void Body::Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0, long long int lowE0, long long int highE0, long long int lowET0, long long int highET0){
  x=x0;y=y0;z=z0;  R=R0;  Vx=Vx0;Vy=Vy0;Vz=Vz0; m=(DENSITY)*(4.0/3.0*M_PI*R*R*R);

  lowE=lowE0;  highE=highE0;
  lowET=lowET0;  highET=highET0;

  q=(highE-highET+lowE-lowET);

  Fx=0;Fy=0;Fz=0;Foldx=0;Foldy=0;Foldz=0;
  nextx=nexty=nextz=-1;
  nextW=-1;
  Vol=4.0/3.0*M_PI*R*R*R;

  W1=0;W2=0;W3=0;W4=0;  W5=0;W6=0;
  a=0;b=0;c=0;d=0;  e=0;f=0;g=0;h=0;

  counter=0;
  for(int i =0;i<100;i++)
    Vel[i]=0;
  Checked=0;
  VelAverage=0;
  Charged=0;
  Bed=0;
  EFx=0; EFy=0;
}

void Body::Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0){
  x=x0;y=y0;z=z0;  R=R0;  Vx=Vx0;Vy=Vy0;Vz=Vz0; m=(DENSITY)*(4.0/3.0*M_PI*R*R*R);
  lowE=0;  highE=0;
  lowET=0;  highET=0;
  q=0;
  Fx=0;Fy=0;Fz=0;Foldx=0;Foldy=0;Foldz=0;
  nextx=nexty=nextz=-1;
  Vol=4.0/3.0*M_PI*R*R*R;

  W1=0;W2=0;W3=0;W4=0;  W5=0;W6=0;
  a=0;b=0;c=0;d=0;  e=0;f=0;g=0;h=0;
  H=0; OH=0;
  H_0=0; OH_0=0;
  counter=0;
  for(int i =0;i<100;i++)
    Vel[i]=0;
  Checked=0;
  VelAverage=0;
  Charged=0;
  Bed=0;
  EFx=0; EFy=0;
}

void Body::Init(double x0,double y0, double z0, double Vx0, double Vy0, double Vz0, double R0, long long int H_0, long long int OH_0, long long int q0){
  x=x0;y=y0;z=z0;  R=R0;  Vx=Vx0;Vy=Vy0;Vz=Vz0; m=(DENSITY)*(4.0/3.0*M_PI*R*R*R);
  H=H_0;
  OH=OH_0;
  q=q0;
  At=4*M_PI*R*R*ALPHA;
  Fx=0;Fy=0;Fz=0;Foldx=0;Foldy=0;Foldz=0;
  nextx=nexty=nextz=-1;
  Vol=4.0/3.0*M_PI*R*R*R;

  W1=0;W2=0;W3=0;W4=0;  W5=0;W6=0;
  a=0;b=0;c=0;d=0;  e=0;f=0;g=0;h=0;
  WallFlux=0;
  counter=0;
  for(int i =0;i<100;i++)
    Vel[i]=0;
  Checked=0;
  VelAverage=0;
  Charged=0;
  Bed=0;
  EFx=0; EFy=0;
}

void Body::Ad(void){
  x=x/LENGTH;  Vx=Vx*TIME/LENGTH; 
  y=y/LENGTH;  Vy=Vy*TIME/LENGTH; 
  z=z/LENGTH;  Vz=Vz*TIME/LENGTH; 
  R=R/LENGTH;  m=m/MASS;
  Vol=Vol/(LENGTH*LENGTH*LENGTH);
  q=q/AMPERE;
  At=At/(LENGTH*LENGTH);
}


void Body::UpdatePos(double dt){

  x+= dt*Vx+Fx*dt*dt/(m*2.0);
  y+= dt*Vy+Fy*dt*dt/(m*2.0);
  z+= dt*Vz+Fz*dt*dt/(m*2.0);
  Foldx = Fx;
  Foldy = Fy;
  Foldz = Fz;

  if(x<0)x+=WIDTH;
  if(x>WIDTH)x-=WIDTH;
  if(y<0)y+=HEIGHT;
  if(y>HEIGHT)y-=HEIGHT;
  if(z<0)z+=DEEP;
  if(z>DEEP)z-=DEEP;
}
void Body::UpdateVel(double dt, double alt){
  Vx += (Fx+Foldx)*dt/(2.0*m);
  Vy += (Fy+Foldy)*dt/(2.0*m);
  Vz += (Fz+Foldz)*dt/(2.0*m); 
  Vel[counter]=sqrt(Vx*Vx+Vy*Vy+Vz*Vz);
  if(Vel[counter]<0.01/LENGTH && y <= alt*1.0 )  q=0;//{
    //q=0.50*q; if(abs(q)<=1) q=0;}

  counter++;
  if(counter==100)counter=0;
}
void Body::setCharge(long long int deltaH, long long int deltaL){
  lowET+=deltaL;  highET-=deltaH;
  q=(highE-highET+lowE-lowET);
}


void Body::IonTransfer(long long int DCi, long long int DCj){
  q+=-(DCj-DCi);
}

void Body::IonTransferVolume(double DCi){
  //clog<<"Entro al menos"<<endl;
  q+=DCi;
}
void Body::neutralization(long long int deltaL){
  OH+=deltaL;  
  q=(H-OH);
}
void Body::neutralizationIon(long long int Ti){
  q+=Ti;
}
void Body::ElectronTransfer(double dt){
  if(lowET>0){
    long long int transferE=static_cast<long long int>(LIGHT*lowET*dt);
    if(transferE>lowET)transferE=lowET;
    highET+= transferE;
    lowET-= transferE;
  }
}
void Body::setChargeWallIon(void){
  if(q<0){
    WallFlux+=(-q*MOBILITY*ALPHA);
    if(WallFlux>0)
    q+=static_cast<long long int>(WallFlux);
    WallFlux-=static_cast<long int>(WallFlux);
  }
  else
    WallFlux+=(q*MOBILITY*ALPHA);
    if(WallFlux>0)
    q-=static_cast<long long int>(WallFlux);
    WallFlux-=static_cast<long long int>(WallFlux);
}

void Body::setChargeWall(void){
  //q=0.95*q;
  //if(abs(q)<=1) q=0;
  //q=0.50*q; if(abs(q)<=1) q=0;
  q=0;
}

double Body::getAverageVel(void){
  //clog<<"Entro al menos"<<Checked<<endl;
  if(Checked==0){ 
  VelAverage=0;
  for(int i=0;i<100;i++){
    VelAverage+=Vel[i];
    }
  VelAverage/=100.;
  }
  Checked=1;
  return VelAverage;
}
