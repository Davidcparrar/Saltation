#include<fstream> 
#include<iostream>
#include<cmath>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>
#include"../SRC/Parameters.h"

using namespace std;

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

  char *base1,*base2;
  string file32("CONCENTRATION/Concentration");

  char input[256];
  char input2[256];
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int ini,last;
  ofstream fpout;
  ofstream fout;
  char outfile[256];
  char ginfile[256];
  while((option_letter = getopt(argc, argv, "i:n:w:v:e:g:ph")) != -1 ){
  
    if (option_letter == 'i') {    
      ini = atoi(optarg);
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
    else if (option_letter == 'g') {    
      strcpy(ginfile, optarg);
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

  int i,j,runningTime,maxHeight=0;


  FILE *fw=NULL;
  FILE *fe=NULL;
  runningTime = last-ini;
  if(runningTime<0){
    printf("Cannot Open File\n");
    exit(1);	
  }    
  clog<<last<<endl;

  base2="CONCENTRATION/Concentration";
  int SIZESYSTEM = 200;
  double Q=0;
  long double HauxC[SIZESYSTEM];
  long double HAuxShearX[SIZESYSTEM];
  long double HAuxShearY[SIZESYSTEM];
  long double HauxVelX[SIZESYSTEM];
  long double HauxVelY[SIZESYSTEM];
  long double HauxCVelX[SIZESYSTEM];
  long double HauxChargeNonAbs[SIZESYSTEM];
  long double HauxMass[SIZESYSTEM];
  long double HauxChargeToArea[SIZESYSTEM];
  long double Hpart[SIZESYSTEM];
  long double HauxChargeNonAbs2[SIZESYSTEM];
  long double HauxVelProf[SIZESYSTEM];
  #ifdef BINDIST
  long double HauxSmall[SIZESYSTEM];
  long double HauxBig[SIZESYSTEM];
  long double HauxSmallCh[SIZESYSTEM];
  long double HauxBigCh[SIZESYSTEM];
  long double HauxTotalCh[SIZESYSTEM];
  long double HauxChargeB[SIZESYSTEM];
  long double HauxChargeS[SIZESYSTEM];
  long double MassFlux[SIZESYSTEM];
  long double CargeFlux[SIZESYSTEM];
  long double Total=0, TotalCh=0;
  #endif
  double EF[374];
  for(i=0;i<374;i++) EF[i]=0;

  double rho=DENSITY;
  double rhoair=AIRDENSITY;
  int counter=0;
  int filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "Concentration%d.txt",filei);
  fpout.open(outfile);

  filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "Flux%d.txt",filei);
  fout.open(outfile);
  sprintf(input2, "ElectricFieldRunTime%d.txt",filei);

  if((fe=fopen(input2, "r"))==NULL) {
      printf("Cannot Open File\n");
      exit(1);
  }
  double a1,a2,a3,a4,a5,a6,a7,a8,a9;
  for(i=0;i<374;i++){
      fscanf(fe,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9);
     EF[i]=a4*1e6;
  }
  fclose(fe);

  double auxC=0,AuxShearX=0,AuxShearY=0,auxVelX=0,auxVelY=0,auxCVelX=0,auxChargeNonAbs=0,auxMass=0,auxChargeToArea=0,part=0,auxChargeNonAbs2=0,auxVelProf=0,auxBCharge=0,auxSCharge=0,ausMassF=0,ausChargeF=0;
  #ifdef BINDIST
  double small=0, big=0, smallCh=0, bigCh=0;
  #endif
   for(i=0;i<SIZESYSTEM;i++){  
     HauxC[i]=0;
     HAuxShearX[i]=0;
     HAuxShearY[i]=0;
     HauxVelX[i]=0;
     HauxVelY[i]=0;
     HauxCVelX[i]=0;
     HauxChargeNonAbs[i]=0;
     HauxMass[i]=0;
     HauxChargeToArea[i]=0;
     Hpart[i]=0;
     HauxChargeNonAbs2[i]=0;
     HauxVelProf[i]=0;
     HauxSmall[i]=0;
     HauxBig[i]=0;
     HauxSmallCh[i]=0;
     HauxBigCh[i]=0;
     HauxTotalCh[i]=0;
     HauxChargeB[i]=0;
     HauxChargeS[i]=0;
     MassFlux[i]=0;
     CargeFlux[i]=0;
    }

  for(j=ini;j<=last;j++){ 

    sprintf (input, "%s%05d", base2, j);  
    if((fw=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      exit(1);
    }

  double sumCharge=0, sumMass = 0;
  double sumFlux = 0;
  int set = 0;
    for(i=0;i<SIZESYSTEM;i++){  
    #ifdef BINDIST
      fscanf(fw,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",&Q, &auxC, &AuxShearX, &AuxShearY, &auxVelX, &auxVelY, &auxCVelX, &auxChargeNonAbs, &auxMass, &auxChargeToArea, &part, &auxChargeNonAbs2, &auxVelProf, &small, &big, &smallCh, &bigCh, &auxSCharge, &auxBCharge, &ausChargeF, &ausMassF);
    #else
      fscanf(fw,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",&Q,&auxC,&AuxShearX,&AuxShearY,&auxVelX,&auxVelY,&auxCVelX,&auxChargeNonAbs,&auxMass,&auxChargeToArea,&part,&auxChargeNonAbs2,&auxVelProf);
    #endif
       //clog<<i<<"\t"<<auxC<<"\t"<<AuxShearX<<"\t"<<AuxShearY<<"\t"<<auxVelX<<"\t"<<auxVelY<<"\t"<<auxCVelX<<"\t"<<auxChargeNonAbs<<"\t"<<auxMass<<"\t"<<auxChargeToArea<<"\t"<<part<<"\t"<<auxChargeNonAbs2<<"\t"<<auxVelProf<<endl;
     HauxC[i]+=auxC;
     HAuxShearX[i]+=AuxShearX;
     HAuxShearY[i]+=AuxShearY;
     HauxVelX[i]+=auxVelX*part*200000;
     HauxVelY[i]+=auxVelY*part*200000;
     HauxCVelX[i]+=auxCVelX;
     HauxChargeNonAbs[i]+=auxChargeNonAbs;
     HauxMass[i]+=auxMass;
     HauxChargeToArea[i]+=auxChargeToArea;
     Hpart[i]+=part;
     HauxChargeNonAbs2[i]+=auxChargeNonAbs2;
     HauxVelProf[i]+=auxVelProf;
     #ifdef BINDIST 
     HauxSmall[i]+=small;
     HauxBig[i]+=big;
     HauxSmallCh[i]+=smallCh;
     HauxBigCh[i]+=bigCh;
     HauxTotalCh[i]+=bigCh+smallCh;
     HauxChargeB[i]+=auxBCharge/ELEMENTARYCHARGE;
     HauxChargeS[i]+=auxSCharge/ELEMENTARYCHARGE;
     Total+=(small+big);
     TotalCh+=(smallCh+bigCh);
     MassFlux[i]+=ausMassF-sumMass;
     sumFlux+=ausMassF-sumMass;
     CargeFlux[i]+=ausChargeF-sumCharge;
     //sumMass+=(ausMassF-sumMass);
     //sumCharge+=(ausChargeF-sumCharge);

     #endif
    if(Hpart[i]==0 && set == 0){
      maxHeight=i;
      set = 1;
    }
     //clog<<i<<"\t"<<HauxC[i]/counter<<"\t"<<HAuxShearX[i]/counter<<"\t"<<HAuxShearY[i]/counter<<"\t"<<HauxVelX[i]/counter<<"\t"<<HauxVelY[i]/counter<<"\t"<<HauxCVelX[i]/counter<<"\t"<<HauxChargeNonAbs[i]/counter<<"\t"<<HauxMass[i]/counter<<"\t"<<HauxChargeToArea[i]/counter<<"\t"<<Hpart[i]/counter<<"\t"<<HauxChargeNonAbs2[i]/counter<<"\t"<<HauxVelProf[i]/counter<<"\t"<<endl;

    }
    set = 0;
    counter ++;    
    Q=0;
    fout<<0.2*j<<"\t"<<sumFlux/(DEEP*WIDTH)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))<<endl;
    sumFlux=0;
    fclose(fw);
  } 
  fout.close();
  double EfForce=0; int index1,index2;
  double massOhyeah=0,  chargeOhyeah=0;
  for(i=1;i<maxHeight;i++){
    if(i%2==0){
     index1=static_cast<int>(2.5*i);
     EfForce = EF[index1];
    }
    else{
     index1=static_cast<int>(2.5*i);
     index2=static_cast<int>(2.5*i)+1;
     EfForce = (EF[index1]+EF[index2])/2.0;
    }
    if(Hpart[i]>0)
    fpout<<scientific<<(i-1)*0.0002*10<<"\t"<<HauxC[i]/counter<<"\t"<<HAuxShearX[i]/(counter*DEEP*WIDTH)<<"\t"<<-HAuxShearY[i]/counter<<"\t"<<HauxVelX[i]/(Hpart[i]*200000)<<"\t"<<HauxVelY[i]/(Hpart[i]*200000)<<"\t"<<HauxVelX[i]/(Hpart[i]*200000)*HauxC[i]/counter<<"\t"<<HauxChargeNonAbs[i]*1e6/HauxMass[i]<<"\t"<<HauxChargeToArea[i]/counter<<"\t"<<HauxChargeNonAbs[i]/(Hpart[i]*200000)<<"\t"<<HauxChargeNonAbs2[i]/counter<<"\t"<<HauxVelProf[i]/counter<<"\t"<<MassFlux[i]/(counter*DEEP*WIDTH)<<"\t"<<1e6*CargeFlux[i]/(counter*DEEP*WIDTH)<<"\t"<<HauxMass[i]*GRAVITY/(Hpart[i]*200000)<<"\t"<<abs(HauxChargeNonAbs[i]/(Hpart[i]*200000)*EfForce)<<"\t"<<abs(HauxChargeNonAbs[i]/abs((Hpart[i]*200000)*EfForce)/(HauxMass[i]*GRAVITY/(Hpart[i]*200000)))<<"\t"<<HauxChargeNonAbs[i]/HauxMass[i]*(EfForce/GRAVITY)<<"\t"<<HauxChargeNonAbs[i]*1e6/(i*0.0002*10*0.0002*12*0.0002*48)<<"\t"<<HauxChargeNonAbs[i]*1e6<<"\t"<<HauxMass[i]<<endl;
    else
    fpout<<scientific<<(i-1)*0.0002*10<<"\t"<<HauxC[i]/counter<<"\t"<<HAuxShearX[i]/counter<<"\t"<<HAuxShearY[i]/counter<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<HauxChargeToArea[i]/counter<<"\t"<<0<<"\t"<<HauxChargeNonAbs2[i]/counter<<"\t"<<HauxVelProf[i]/counter<<"\t"<<MassFlux[i]/(counter*DEEP*WIDTH)<<"\t"<<1e6*CargeFlux[i]/(counter*DEEP*WIDTH)<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
    massOhyeah+=MassFlux[i];
    chargeOhyeah+=CargeFlux[i];
  }

  fpout.close();

  #ifdef BINDIST
  filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "SizeDist%d.txt",filei);
  fout.open(outfile);
  for(i=1;i<maxHeight;i++){
    if(TotalCh!=0)
      if(HauxTotalCh[i]!=0)
    fout<<scientific<<(i-1)*0.0002*10<<"\t"<<1.0*HauxSmall[i]/(1.0*Total)<<"\t"<<1.0*HauxBig[i]/(1.0*Total)<<"\t"<<1.0*HauxSmallCh[i]/(1.0*TotalCh)<<"\t"<<1.0*HauxBigCh[i]/(1.0*TotalCh)<<"\t"<<1.0*HauxSmallCh[i]/(1.0*Total)<<"\t"<<1.0*HauxBigCh[i]/(1.0*Total)<<"\t"<<HauxSmallCh[i]/HauxTotalCh[i]<<"\t"<<-HauxChargeS[i]/(0.0002*10*WIDTH*DEEP)*ELEMENTARYCHARGE/counter<<"\t"<<HauxChargeB[i]/(0.0002*10*WIDTH*DEEP)*ELEMENTARYCHARGE/counter<<"\t"<<(HauxChargeB[i]/counter+HauxChargeS[i]/counter)*ELEMENTARYCHARGE/(0.0002*10*WIDTH*DEEP)<<"\t"<<(HauxChargeB[i]/counter+HauxChargeS[i]/counter)*ELEMENTARYCHARGE/(Hpart[i]*200000)<<"\t"<<HauxSmallCh[i]/counter<<"\t"<<HauxBigCh[i]/counter<<"\t"<<(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i])<<"\t"<<(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i])*MassFlux[i]<<"\t"<<(1-(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i]))*MassFlux[i]<<"\t"<<MassFlux[i]/(DEEP*WIDTH)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))/counter<<"\t"<<HauxSmall[i]/Hpart[i]<<"\t"<<HauxBig[i]/Hpart[i]<<"\t"<<Hpart[i]/counter<<endl;
      else
    fout<<scientific<<(i-1)*0.0002*10<<"\t"<<1.0*HauxSmall[i]/(1.0*Total)<<"\t"<<1.0*HauxBig[i]/(1.0*Total)<<"\t"<<1.0*HauxSmallCh[i]/(1.0*TotalCh)<<"\t"<<1.0*HauxBigCh[i]/(1.0*TotalCh)<<"\t"<<1.0*HauxSmallCh[i]/(1.0*Total)<<"\t"<<1.0*HauxBigCh[i]/(1.0*Total)<<"\t"<<0<<"\t"<<-HauxChargeS[i]/(0.0002*10*WIDTH*DEEP)*ELEMENTARYCHARGE/counter<<"\t"<<HauxChargeB[i]/(0.0002*10*WIDTH*DEEP)*ELEMENTARYCHARGE/counter<<"\t"<<(HauxChargeB[i]/counter+HauxChargeS[i]/counter)*ELEMENTARYCHARGE/(0.0002*10*WIDTH*DEEP)<<"\t"<<(HauxChargeB[i]/counter+HauxChargeS[i]/counter)*ELEMENTARYCHARGE/(Hpart[i]*200000)<<"\t"<<HauxSmallCh[i]/counter<<"\t"<<HauxBigCh[i]/counter<<"\t"<<(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i])<<"\t"<<(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i])*MassFlux[i]<<"\t"<<(1-(HauxSmallCh[i]+HauxBigCh[i])/(Hpart[i]))*MassFlux[i]<<"\t"<<MassFlux[i]/(DEEP*WIDTH)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))/counter<<"\t"<<HauxSmall[i]/Hpart[i]<<"\t"<<HauxBig[i]/Hpart[i]<<"\t"<<Hpart[i]/counter<<endl;
    else
    fout<<scientific<<(i-1)*0.0002*10<<"\t"<<1.0*HauxSmall[i]/(1.0*Total)<<"\t"<<1.0*HauxBig[i]/(1.0*Total)<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<"\t"<<0<<endl;
  }
  fout.close();
  #endif

  ofstream fcout;

  fcout.open(ginfile,ios::app);

  fcout<<WINDVEL<<"\t"<<massOhyeah/(DEEP*WIDTH*counter)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))<<"\t"<<chargeOhyeah/(DEEP*WIDTH*counter)/(rho*sqrt((rho/rhoair-1)*GRAVITY*D_MEAN*D_MEAN*D_MEAN))<<endl;

  fcout.close();

  #ifdef BINDIST
  filei=static_cast<int>(100.0*WINDVEL);
  sprintf(outfile, "HeightProb%d.txt",filei);
  fout.open(outfile);

  double ProbH[SIZESYSTEM];
  double ProbHPos[SIZESYSTEM];
  double ProbHNeg[SIZESYSTEM];
  double tH=0,tHPos=0,tHNeg=0;
  for(i=1;i<SIZESYSTEM;i++){
    ProbH[i]=0;
    ProbHPos[i]=0;
    ProbHNeg[i]=0;
  }

  for(i=SIZESYSTEM-2;i>=0;i--){
    ProbH[i]=Hpart[i]-Hpart[i+1];
    ProbHPos[i]=HauxBig[i]-HauxBig[i+1];
    ProbHNeg[i]=HauxSmall[i]-HauxSmall[i+1];
    tH+=Hpart[i];
    tHPos+=HauxBig[i];
    tHNeg+=HauxSmall[i];
  }
  for(i=1;i<SIZESYSTEM;i++)
    fout<<(i-1)*0.0002*10<<"\t"<<ProbH[i]<<"\t"<<ProbHPos[i]<<"\t"<<ProbHNeg[i]<<endl;
  fout.close();
  #endif

  return 0;
}
