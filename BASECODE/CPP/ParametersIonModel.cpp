#include<fstream> 
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>
#include<time.h>
#include<dirent.h>
#include<math.h>


using namespace std;


int main(int argc, char** argv){
  double Ri=1.3e-4;
  double Rj=8e-5;
  double ELEMENTARYCHARGE = 1.60217646e-19;
  double RHO_H = 1e+20;
  double E0=8.85418782e-12;
  double KE = 1.0/(4.0*M_PI*E0);
  double MOBILITY = 1e15;
  double DIFFUSIVITY=0;
  double ALPHA=1e-10;
  double Ai=4*M_PI*Ri*Ri;
  double Aj=4*M_PI*Rj*Rj;
  double Ati=Ai*ALPHA;
  double Atj=Aj*ALPHA;

  long int Ci=RHO_H*Ai;
  long int Ci0=RHO_H*Ai;
  long int Cj=RHO_H*Aj;
  long int Cj0=RHO_H*Aj;

  long int qi=0;
  long int qj=0;
  long int DCi=0;
  long int DCj=0;
  double Ei=0;
  double Ej=0;
  double Ti=0;
  double Tj=0;
  double Di=0;
  double Dj=0;
  double Ti2=0;
  double Tj2=0;
  cout<<MOBILITY*KE<<"\t"<<Ati<<"\t"<<Ci<<"\t"<<MOBILITY<<"\t"<<1/KE<<endl;
  double sumContacti=0;
  double sumContactj=0;
  for(int t=0;t<100015;t++){
     if(t==0){
       DCi=static_cast<long int>(Ati*Ci/Ai);
       DCj=static_cast<long int>(Atj*Cj/Aj);
     }
     else{
       DCi=0;
       DCj=0;
     }

     Ei=-ELEMENTARYCHARGE*(qi/Ai-qj/Aj);
     Ei=-ELEMENTARYCHARGE*(qj/Aj-qi/Ai);

     if(Ei<0) sumContacti+=(-Ei*MOBILITY*KE*Ati);
     else sumContacti=0;
     if(Ej<0) sumContactj+=(-Ej*MOBILITY*KE*Atj);
     else sumContactj=0;

     if(sumContacti>0){
       Ti=static_cast<long int>(sumContacti);
       sumContacti-=static_cast<long int>(sumContacti);
     }
     else
       Ti=0;
     if(sumContactj>0){
       Tj=static_cast<long int>(sumContactj);
       sumContactj-=static_cast<long int>(sumContactj);
     }
     else
       Tj=0;
     /*Di=(Ci/Ai-Cj/Aj);
     Dj=(Cj/Aj-Ci/Ai);
     if(Di>0)
     Ti2=(Di*DIFFUSIVITY*Ati);
     else
     Ti2=0;
     if(Dj>0)
     Tj2=Dj*DIFFUSIVITY*Atj;
     else
     Tj2=0;*/


     //if(t>100000)
     clog<<scientific<<t<<"\t"<<Ci<<"\t"<<Cj<<"\t"<<(qi)<<"\t"<<(qj)<<"\t"<<DCi<<"\t"<<DCj<<"\t"<<Ei<<"\t"<<Ej<<"\t"<<Ti<<"\t"<<Tj<<"\t"<<sumContacti<<"\t"<<sumContactj<<endl;


     qi+=(DCj-DCi+Ti-Tj);
     qj+=(DCi-DCj-Ti+Tj);

     if(t>10000)
     qj=0;

     Ci=Ci0+qi;
     Ci=Cj0+qj;

  }

 

  return 0;
}
