#include<fstream> 
#include<iostream>
#include<cmath>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<getopt.h>

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

  FILE *fin=NULL;
  char *base1;
  char input[15];
  extern char *optarg;
  int option_letter=0, usage_flag=0;
  int ini,last,N;
  ofstream fpout;
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
  runningTime = last-ini;
  if(runningTime<0){
    printf("Cannot Open File\n");
    exit(1);	
  }    

  base1="FILES/file";
  double Q=0;
  double Qflux=0;
  double x, y, z, Vx, Vy, Vz, R;
  double ix,iy,iz,iR;
  double *Radii;
  double *Charge;

  Radii = new double [N];
  Charge = new double [N];

  //BI SIZE DIST

  int nPart=0;
  int iyf=0;
  double q;
  double Low,High,qcharge,qchargeDensity;

  for(j=ini;j<=last;j++){ 

    sprintf (input, "%s%05d", base1, j);

    if((fp=fopen(input, "r"))==NULL) {
      printf("Cannot Open File\n");
      clog<<"File"<<endl;
      exit(1);
    }
  sprintf(outfile, "VTU/Visual_%d.vtu",j);
  fpout.open(outfile);

    fpout<<"<?xml version=\"1.0\"?>"<<endl;
    fpout<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"<<endl;
    fpout<<" <UnstructuredGrid GhostLevel=\"0\">"<<endl;
    fpout<<"  <Piece NumberOfPoints=\""<<N<<"\" NumberOfCells=\""<<N<<"\">"<<endl;
    fpout<<"   <Points>"<<endl;
    fpout<<"    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>"<<endl;


    for(i=0;i<N;i++){  
      fscanf(fp, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", &x, &y, &z, &Vx, &Vy, &Vz, &R, &Low,&Low,&High,&High,&q);
      ix =(x/0.0001);
      iy =(y/0.0001);
      iz =(z/0.0001);
      iR =(R/0.0001);

    fpout<<ix<<" "<<iy<<" "<<iz<<endl;
    Radii[i]=iR; 
    Charge[i]=(q/1.6e-19);
    }
    fpout<<"   </Points>"<<endl;
    fpout<<"   <PointData>"<<endl;
    fpout<<"    <DataArray type=\"Float64\" Name=\"radius\"/>"<<endl;  
    for(i=0;i<N;i++){
      fpout<<Radii[i]<<endl;
    }
    fpout<<"    <DataArray type=\"Float64\" Name=\"charge\"/>"<<endl;

    for(i=0;i<N;i++){
      fpout<<Charge[i]<<endl;
    }

    fpout<<"   </PointData>"<<endl;
    fpout<<"   <CellData>"<<endl;
    fpout<<"   </CellData>"<<endl;
    fpout<<"   <Cells>"<<endl;
    fpout<<"    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
    for(i=0;i<N;i++){
      fpout<<0<<endl;
    }
    fpout<<"    </DataArray>"<<endl;
    fpout<<"    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;

    for(i=0;i<N;i++){
      fpout<<i+1<<endl;
    }
    fpout<<"    </DataArray>"<<endl;
    fpout<<"    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">"<<endl;
    for(i=0;i<N;i++){
      fpout<<1<<endl;
    }
    fpout<<"    </DataArray>"<<endl;
    fpout<<"   </Cells>"<<endl;
    fpout<<"  </Piece>"<<endl;
    fpout<<" </UnstructuredGrid>"<<endl;
    fpout<<"</VTKFile>"<<endl;
    fclose(fp);
    fpout.close();
  }

  return 0;
}
