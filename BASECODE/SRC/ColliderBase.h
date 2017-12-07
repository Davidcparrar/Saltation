#include"Wind.h"
#include<stdlib.h>
#include<complex>
//#include<omp.h>
double maxQQ=0,Q1=0,Q2=0;
class Collision{
  public:
  double hold;
  double Vel;
  long int Collisions;
  long int time;
  long int Charge;
  double ChargeCol;
  void init(void);
  void reset(void);
  friend class Collider;
};

void Collision::init(void){
  hold=0;
  Vel=0;
  Collisions=0;
  time=0;
  Charge=0;
  ChargeCol=0;
}

void Collision::reset(void){
  hold=0;
  Vel=0;
  Collisions=0;
  time=0;
  Charge=0;
  ChargeCol=0;
}

class Collider{
public:
  int XSIZE,YSIZE,ZSIZE,Lx,Ly,Lz,Lxp,Lyp,Lzp,RUNTIME;
  double D_MAX,CUTOFF,WSOR,D_CELL;
  List *listz,**listy,***listx;
  List *listzCh,**listyCh,***listxCh;
  Collision *CollitionsGrain;
  double ***potential,***ChargeGrid,***Ex,***Ey,***Ez;
  short ***SolveCell;
  short ***Done;
  double ***Carga;
  double ***GreenF;
  int *smallArray;
  unsigned sizeSA;
  double GreenRF[4][4][4][3];
  double D0;
public:
  ~Collider(void);
  void Init(int N1);
  void BodyVsBody(Body &Bodyi,Body &Bodyj,int i, int j, double Height);
  void ChargeVsCharge(Body &Bodyi,Body &Bodyj,int i, int j, double Height);
  void Collide(Body *Body, double Height);
  void Wall(Body *Body, int i);
  void WallUp(Body *Body, int i);
  void ContactDetection(Body *Body, double Height);
  void ChargeDetection(Body *Body, double Height);
  void ChargeAssignment(Body *Body);
  void ChargeAssignmentDirectEF(Body *Body);
  void GreenFunctionReal(void);
  void PotentialSolver(void);
  void PotentialSolverParallel(void);
  void ElectricFieldReal(int i, int j, int k);
  void ElectricField(void);
  void ElectricForceReal(Body *Body);
  void ReferenceForce(void);
  void PrintPotential(ofstream &fout,Body *Body);
  void PrintCollisionMap(ofstream &fout);
  void PrintElectric(Body *Body);
  void PrintCharge(ofstream &fout);
  void PotentialSolverTotal(double ***Carga);
  void ElectricFieldTotal(void);
  void PrintGreenFunction(ofstream &fout);
  void ReadGreenFunction(char *infile);
  friend class Body;
  friend class Wind;
};

void Collider::Init(int N1){
  int i,coll;
  D_MAX = D_MEAN*1.*2.0;
  D_CELL = D_MEAN*1.*6.0;
  CUTOFF = D_CELL*16.;
  XSIZE = static_cast<int>(WIDTH/D_MAX);
  YSIZE = static_cast<int>(HEIGHT/D_MAX);
  ZSIZE = static_cast<int>(DEEP/D_MAX);

  #ifdef BCY
  Lx=static_cast<int>(WIDTH/D_CELL);
  Ly=static_cast<int>(2.*HEIGHT/D_CELL);
  Lz=static_cast<int>(DEEP/D_CELL);

  Lxp=Lx;
  Lyp=Ly/2.;
  Lzp=Lz;
  #else
  Lx=static_cast<int>(WIDTH/D_CELL);
  Ly=static_cast<int>(HEIGHT/D_CELL);
  Lz=static_cast<int>(DEEP/D_CELL);

  Lxp=Lx;
  Lyp=Ly;
  Lzp=Lz;
  #endif

  coll = N1*(N1-1)/2;
  CollitionsGrain = new Collision [coll];
  for(i=0;i<coll;i++)
    CollitionsGrain[i].init();

  listz = new List [ZSIZE];
  listy = new List* [2];
  listy[0] = new List [YSIZE];
  listy[1] = new List [YSIZE]; 
 
  listx = new List** [2];
  listx[0] = new List* [3];
  listx[1] = new List* [3];  
  listx[0][0] = new List [XSIZE];
  listx[0][1] = new List [XSIZE];  
  listx[0][2] = new List [XSIZE];  
  listx[1][0] = new List [XSIZE];  
  listx[1][1] = new List [XSIZE];  
  listx[1][2] = new List [XSIZE];  

  for(i=0;i<ZSIZE;i++){
    listz[i].head = -1; 
    listz[i].done = 0; 
  }
  for(i=0;i<YSIZE;i++){
    listy[0][i].head = -1; 
    listy[0][i].done = 0; 
    listy[1][i].head = -1; 
    listy[1][i].done = 0; 
  }
  for(i=0;i<XSIZE;i++){
    listx[0][0][i].head = -1;
    listx[0][1][i].head = -1;  
    listx[0][2][i].head = -1;  
    listx[1][0][i].head = -1;
    listx[1][1][i].head = -1;  
    listx[1][2][i].head = -1;
    listx[0][0][i].done = 0;
    listx[0][1][i].done = 0;  
    listx[0][2][i].done = 0; 
    listx[1][0][i].done = 0;
    listx[1][1][i].done = 0;  
    listx[1][2][i].done = 0;  
  }

  listzCh = new List [Lzp];
  listyCh = new List* [2];
  listyCh[0] = new List [Lyp];
  listyCh[1] = new List [Lyp]; 
 
  listxCh = new List** [2];
  listxCh[0] = new List* [3];
  listxCh[1] = new List* [3];  
  listxCh[0][0] = new List [Lxp];
  listxCh[0][1] = new List [Lxp];  
  listxCh[0][2] = new List [Lxp];  
  listxCh[1][0] = new List [Lxp];  
  listxCh[1][1] = new List [Lxp];  
  listxCh[1][2] = new List [Lxp];  

  for(i=0;i<Lzp;i++){
    listzCh[i].head = -1; 
    listzCh[i].done = 0; 
  }
  for(i=0;i<Lyp;i++){
    listyCh[0][i].head = -1; 
    listyCh[0][i].done = 0; 
    listyCh[1][i].head = -1; 
    listyCh[1][i].done = 0; 
  }
  for(i=0;i<Lxp;i++){
    listxCh[0][0][i].head = -1;
    listxCh[0][1][i].head = -1;  
    listxCh[0][2][i].head = -1;  
    listxCh[1][0][i].head = -1;
    listxCh[1][1][i].head = -1;  
    listxCh[1][2][i].head = -1;
    listxCh[0][0][i].done = 0;
    listxCh[0][1][i].done = 0;  
    listxCh[0][2][i].done = 0; 
    listxCh[1][0][i].done = 0;
    listxCh[1][1][i].done = 0;  
    listxCh[1][2][i].done = 0;  
  }

  #ifdef ELECF
  int j,k;
  sizeSA=0;
  RUNTIME=0;
  smallArray=new int [Lxp*Lyp*Lzp];
  potential = new double**[Lxp];
  Ex = new double**[Lxp];
  Ey = new double**[Lxp];
  Ez = new double**[Lxp];
  Done = new short**[Lxp];
  SolveCell = new short**[Lxp];
  ChargeGrid = new double**[Lxp];
  Carga = new double**[Lxp];
  GreenF = new double**[Lxp];
  for(i=0;i<Lxp;i++){

    potential[i]= new double*[Lyp];
    Ex[i]= new double*[Lyp];
    Ey[i]= new double*[Lyp];
    Ez[i]= new double*[Lyp];
    Done[i] = new short*[Lyp];
    SolveCell[i] = new short*[Lyp];
    ChargeGrid[i] = new double*[Lyp];
    Carga[i] = new double*[Lyp];
    GreenF[i] = new double*[Ly];

    for(j=0;j<Lyp;j++){
      potential[i][j]= new double [Lzp];
      Ex[i][j]= new double [Lzp];
      Ey[i][j]= new double [Lzp];
      Ez[i][j]= new double [Lzp];
      Done[i][j] = new short[Lzp];
      SolveCell[i][j] = new short[Lzp];
      ChargeGrid[i][j] = new double[Lzp];
      Carga[i][j] = new double[Lzp];

    }
    for(j=0;j<Ly;j++)
      GreenF[i][j] = new double[Lzp];
  }

  for(i=0;i<Lxp;i++)
    for(j=0;j<Lyp;j++)
      for(k=0;k<Lzp;k++){
        potential[i][j][k]=0.;
 	Done[i][j][k]=0.;
        ChargeGrid[i][j][k]=0.;
        Carga[i][j][k]=0.;

        SolveCell[i][j][k]=0;
  }
  for(i=0;i<Lxp;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Lzp;k++)
        GreenF[i][j][k]=0.;

  ReferenceForce();
  #ifdef PGF
  GreenFunctionReal();
  #endif
  #endif
  D0=2e-9/LENGTH;//REDUCEDPLANCK/sqrt(8.0*ELECTRONMASS*Eb0)*log(REDUCEDPLANCK*TCOLL/(ETA*ELECTRONMASS*A2));
}

Collider::~Collider(void){
  delete [] listz;
  delete [] listy[0];
  delete [] listy[1];
  delete [] listy; 
  delete [] listx[0][0];
  delete [] listx[0][1];
  delete [] listx[1][0];
  delete [] listx[1][1];
  delete [] listx[0][2];
  delete [] listx[1][2];
  delete [] listx[0];
  delete [] listx[1];
  delete [] listx;
  delete [] CollitionsGrain;
  #ifdef ELECF
  int i,j;
  for(i=0;i<Lxp;i++){
   for(j=0;j<Lyp;j++){
    delete [] Ex[i][j];
    delete [] Ey[i][j];
    delete [] Ez[i][j];
    delete [] potential[i][j];
    delete [] Done[i][j];
    delete [] SolveCell[i][j];
    delete [] ChargeGrid[i][j];
    delete [] Carga[i][j];
    }

   for(j=0;j<Ly;j++)
      delete [] GreenF[i][j];

    delete [] Ex[i];
    delete [] Ey[i];
    delete [] Ez[i];
    delete [] potential[i];
    delete [] Done[i];
    delete [] SolveCell[i];
    delete [] ChargeGrid[i];
    delete [] Carga[i];
    delete [] GreenF[i];
  }
    delete [] Ex;
    delete [] Ey;
    delete [] Ez;
    delete [] potential;
    delete [] Done;
    delete [] SolveCell;
    delete [] ChargeGrid;
    delete [] Carga;
    delete [] GreenF;
  #endif
}


void Collider::ContactDetection(Body *Body, double Height){
  int i, j, k, l, m;
  int ix,iy,iz,ipx,ipy,ipz;

/* zero all headers and done flags for iz-list */
  for(i=0;i<ZSIZE;i++){
    listz[i].done=0;
    listz[i].head=-1;
  }

/* zero all headers and done flags for iy lists and iz and iz-1 columns */
  for(i=0;i<YSIZE;i++){
    listy[0][i].done=0;
    listy[0][i].head=-1;

    listy[1][i].done=0;
    listy[1][i].head=-1;
  }

/* zero headers and done flags for ix lists in the iy, iy+1, and iy-1 cols */
  for(i=0;i<XSIZE;i++){
    listx[0][0][i].done=0;
    listx[0][0][i].head=-1;

    listx[0][1][i].done=0;
    listx[0][1][i].head=-1;

    listx[0][2][i].done=0;
    listx[0][2][i].head=-1;

    listx[1][0][i].done=0;
    listx[1][0][i].head=-1;

    listx[1][1][i].done=0;
    listx[1][1][i].head=-1;

    listx[1][2][i].done=0;
    listx[1][2][i].head=-1;	  
  }



/* sorts Bodys into iz cells before we begin */
  for(i=0;i<N;i++){
    iz=static_cast<int>(Body[i].z/D_MAX);
    Body[i].nextz=listz[iz].head;
    listz[iz].head=i;
  }

/* loop over all Bodys */
  for(i=N-1;i>=0;i--){
/* check which iz-list Body i is in */
    iz=static_cast<int>(Body[i].z/D_MAX);
/* if iy lists have not yet been generated for this iz list, do it */
    if(listz[iz].done!=1){
/* mark it as done, so we do not do it again */
      listz[iz].done=1;			    			
/* we need to start at the top (head) of iz-list (already generated) */
      j=listz[iz].head;
      while(j!=-1){
/* loop over all Bodys in this iz list to place them in iy lists */
        iy=static_cast<int>(Body[j].y/D_MAX);		    	
        Body[j].nexty=listy[0][iy].head;
        listy[0][iy].head=j;
        j=Body[j].nextz;
      }
/* if iz == 0 we must fix to loop through list iz-1 */
      if(iz==0) ipz=ZSIZE-1;
      else ipz=iz-1;
      j=listz[ipz].head;
      while(j!=-1){
/* loop over all Bodys in iz-1 list to place them in iy lists */
        iy=static_cast<int>(Body[j].y/D_MAX);	    			
        Body[j].nexty=listy[1][iy].head;
        listy[1][iy].head=j;
        j=Body[j].nextz;
      }
/* now that iy lists are generated for iz and iz-1, start again to get the ix lists (for iz, iz-1, and all the related iys) -- resetting iz now.... */
      j=listz[iz].head;
      while(j!=-1){
/* check which iy-list Body j is in */
        iy=static_cast<int>(Body[j].y/D_MAX);		    	
/* if ix lists have not yet been generated for this iy list, do it */
        if(listy[0][iy].done!=1){
/* mark it as done, so we do not do it again */
          listy[0][iy].done=1;		    			
/* we need to start at the top (head) of iy-list (already generated)*/
          k=listy[0][iy].head;
	  while(k!=-1){
/*loop over all Bodys in this iy list, place them in ix lists*/
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    Body[k].nextx=listx[0][0][ix].head;		
            listx[0][0][ix].head=k;
	    k=Body[k].nexty;
	  }
	  k=listy[1][iy].head;
	  while(k!=-1){
/*loop over Bodys in iy (iz-1) list, place them in ix lists*/
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    Body[k].nextx=listx[1][0][ix].head;	     
	    listx[1][0][ix].head=k;
	    k=Body[k].nexty;
	  }
	  if(iy==0) ipy=YSIZE-1;
	  else ipy=iy-1;
	  k=listy[0][ipy].head;
	  while(k!=-1){
/* loop over Bodys in iy-1 (iz) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    Body[k].nextx=listx[0][1][ix].head; 	
	    listx[0][1][ix].head=k;
	    k=Body[k].nexty;
	  }
	  k=listy[1][ipy].head;
	  while(k!=-1){
/*loop over Bodys in iy-1 (iz-1) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    Body[k].nextx=listx[1][1][ix].head;       
	    listx[1][1][ix].head=k;
	    k=Body[k].nexty;
	  }
	  if(iy==(YSIZE-1)) ipy=0;
	  else ipy=iy+1;
	  k=listy[1][ipy].head;
	  while (k!=-1){
/*loop over Bodys in iy+1 (iz-1) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    Body[k].nextx=listx[1][2][ix].head;       
	    listx[1][2][ix].head=k;
	    k=Body[k].nexty;
	  }
/* lists are now complete, let's check if anyone is touching! */
/* starting with iz ([0]) list and loop through iy */
	  k=listy[0][iy].head;             			
	  while(k!=-1){
/* check which ix-list Body i is in */
	    ix=static_cast<int>(Body[k].x/D_MAX);
/* is this the first time checking collisions for this iz-iy-ix */
	    if(listx[0][0][ix].done!=1){
/* mark it as done, so we do not do it again */
	      listx[0][0][ix].done=1;
/* loop through all ixs in ix-iy-iz list */
	      l=listx[0][0][ix].head;

	      while(l!=-1){
                m=Body[l].nextx;
/* check collisions between ix-iy-iz and ix-iy-iz */
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[0][1][ix].head;
/* check collisions between ix-iy-iz and ix-iy-1-iz */
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][0][ix].head;
/* check collisions between ix-iy-iz and ix-iy-iz-1 */
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][1][ix].head;
/* check collisions between ix-iy-iz and ix-iy-1-iz-1 */
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][2][ix].head;
/* check collisions between ix-iy-iz and ix-iy+1-iz-1 */
                while (m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        if(ix==0) ipx=XSIZE-1;
	        else ipx=ix-1;
	        m=listx[0][0][ipx].head;
/* check collisions between ix-iy-iz and ix-1-iy-iz */
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[0][1][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][0][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][1][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][2][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        if(ix==(XSIZE-1)) ipx=0;
	        else ipx=ix+1;
	        m=listx[0][1][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][0][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][1][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        m=listx[1][2][ipx].head;
	        while(m!=-1){
		  BodyVsBody(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextx;
	        }
	        l=Body[l].nextx;
	      }
	    }
	    k=Body[k].nexty;
	  }
	  k=listy[0][iy].head;
/* loop over all ix-iy-iz lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    listx[0][0][ix].head=-1;
	    listx[0][0][ix].done=0;
	    k=Body[k].nexty;
	  }
	  k=listy[1][iy].head;
/* loop over ix-iy-iz1 lists to discard them*/
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    listx[1][0][ix].head=-1;
	    listx[1][0][ix].done=0;
	    k=Body[k].nexty;
	  }
	  if(iy==0) ipy=YSIZE-1;
	  else ipy=iy-1;
	  k=listy[0][ipy].head;
/* loop over all ix-iy1-iz lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    listx[0][1][ix].head=-1;
	    listx[0][1][ix].done=0;
	    k=Body[k].nexty;
	  }
	  k=listy[1][ipy].head;
/* loop over all ix-iy1-iz1 lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    listx[1][1][ix].head=-1;
	    listx[1][1][ix].done=0;
	    k=Body[k].nexty;
	  }
	  if(iy==(YSIZE-1)) ipy=0;
	  else ipy=iy+1;
	  k=listy[1][ipy].head;
/* loop over all ix-iy2-iz1 lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_MAX);
	    listx[1][2][ix].head=-1;
	    listx[1][2][ix].done=0;
	    k=Body[k].nexty;
	  }
        } 
        j=Body[j].nextz;
      }
      j=listz[iz].head;
/* loop over iz-list  Bodys to discard iz-iy lists */
      while (j!=-1){
        iy=static_cast<int>(Body[j].y/D_MAX);
        listy[0][iy].head=-1;
        listy[0][iy].done=0;
        j=Body[j].nextz;
      }
      if(iz==0) ipz=ZSIZE-1;
      else ipz=iz-1;
      j=listz[ipz].head;
      while(j!=-1){
        iy=static_cast<int>(Body[j].y/D_MAX);
        listy[1][iy].head=-1;
        listy[1][iy].done=0;
        j=Body[j].nextz;
      }
    }		
  }
}

void Collider::ChargeDetection(Body *Body, double Height){
  int i, j, k, l, m;
  int ix,iy,iz,ipx,ipy,ipz;

/* zero all headers and done flags for iz-list */
  for(i=0;i<Lzp;i++){
    listzCh[i].done=0;
    listzCh[i].head=-1;
  }

/* zero all headers and done flags for iy lists and iz and iz-1 columns */
  for(i=0;i<Lyp;i++){
    listyCh[0][i].done=0;
    listyCh[0][i].head=-1;

    listyCh[1][i].done=0;
    listyCh[1][i].head=-1;
  }

/* zero headers and done flags for ix lists in the iy, iy+1, and iy-1 cols */
  for(i=0;i<Lxp;i++){
    listxCh[0][0][i].done=0;
    listxCh[0][0][i].head=-1;

    listxCh[0][1][i].done=0;
    listxCh[0][1][i].head=-1;

    listxCh[0][2][i].done=0;
    listxCh[0][2][i].head=-1;

    listxCh[1][0][i].done=0;
    listxCh[1][0][i].head=-1;

    listxCh[1][1][i].done=0;
    listxCh[1][1][i].head=-1;

    listxCh[1][2][i].done=0;
    listxCh[1][2][i].head=-1;	  
  }



/* sorts Bodys into iz cells before we begin */
  for(i=0;i<N;i++){
    if(Body[i].q!=0){
      iz=static_cast<int>(Body[i].z/D_CELL);
      Body[i].nextzCh=listzCh[iz].head;
      listzCh[iz].head=i;
    }
  }

/* loop over all Bodys */
  for(i=N-1;i>=0;i--){
/* check which iz-list Body i is in */
    if(Body[i].q!=0){
    iz=static_cast<int>(Body[i].z/D_CELL);
/* if iy lists have not yet been generated for this iz list, do it */
    if(listzCh[iz].done!=1){
/* mark it as done, so we do not do it again */
      listzCh[iz].done=1;			    			
/* we need to start at the top (head) of iz-list (already generated) */
      j=listzCh[iz].head;
      while(j!=-1){
/* loop over all Bodys in this iz list to place them in iy lists */
        iy=static_cast<int>(Body[j].y/D_CELL);		    	
        Body[j].nextyCh=listyCh[0][iy].head;
        listyCh[0][iy].head=j;
        j=Body[j].nextzCh;
      }
/* if iz == 0 we must fix to loop through list iz-1 */
      if(iz==0) ipz=Lzp-1;
      else ipz=iz-1;
      j=listzCh[ipz].head;
      while(j!=-1){
/* loop over all Bodys in iz-1 list to place them in iy lists */
        iy=static_cast<int>(Body[j].y/D_CELL);	    			
        Body[j].nextyCh=listyCh[1][iy].head;
        listyCh[1][iy].head=j;
        j=Body[j].nextzCh;
      }
/* now that iy lists are generated for iz and iz-1, start again to get the ix lists (for iz, iz-1, and all the related iys) -- resetting iz now.... */
      j=listzCh[iz].head;
      while(j!=-1){
/* check which iy-list Body j is in */
        iy=static_cast<int>(Body[j].y/D_CELL);		    	
/* if ix lists have not yet been generated for this iy list, do it */
        if(listyCh[0][iy].done!=1){
/* mark it as done, so we do not do it again */
          listyCh[0][iy].done=1;		    			
/* we need to start at the top (head) of iy-list (already generated)*/
          k=listyCh[0][iy].head;
	  while(k!=-1){
/*loop over all Bodys in this iy list, place them in ix lists*/
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    Body[k].nextxCh=listxCh[0][0][ix].head;		
            listxCh[0][0][ix].head=k;
	    k=Body[k].nextyCh;
	  }
	  k=listyCh[1][iy].head;
	  while(k!=-1){
/*loop over Bodys in iy (iz-1) list, place them in ix lists*/
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    Body[k].nextxCh=listxCh[1][0][ix].head;	     
	    listxCh[1][0][ix].head=k;
	    k=Body[k].nextyCh;
	  }
	  if(iy==0) ipy=Lyp-1;
	  else ipy=iy-1;
	  k=listyCh[0][ipy].head;
	  while(k!=-1){
/* loop over Bodys in iy-1 (iz) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    Body[k].nextxCh=listxCh[0][1][ix].head; 	
	    listxCh[0][1][ix].head=k;
	    k=Body[k].nextyCh;
	  }
	  k=listyCh[1][ipy].head;
	  while(k!=-1){
/*loop over Bodys in iy-1 (iz-1) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    Body[k].nextxCh=listxCh[1][1][ix].head;       
	    listxCh[1][1][ix].head=k;
	    k=Body[k].nextyCh;
	  }
	  if(iy==(Lyp-1)) ipy=0;
	  else ipy=iy+1;
	  k=listyCh[1][ipy].head;
	  while (k!=-1){
/*loop over Bodys in iy+1 (iz-1) list, place them in ix lists */
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    Body[k].nextxCh=listxCh[1][2][ix].head;       
	    listxCh[1][2][ix].head=k;
	    k=Body[k].nextyCh;
	  }
/* lists are now complete, let's check if anyone is touching! */
/* starting with iz ([0]) list and loop through iy */
	  k=listyCh[0][iy].head;             			
	  while(k!=-1){
/* check which ix-list Body i is in */
	    ix=static_cast<int>(Body[k].x/D_CELL);
/* is this the first time checking collisions for this iz-iy-ix */
	    if(listxCh[0][0][ix].done!=1){
/* mark it as done, so we do not do it again */
	      listx[0][0][ix].done=1;
/* loop through all ixs in ix-iy-iz list */
	      l=listxCh[0][0][ix].head;

	      while(l!=-1){
                m=Body[l].nextxCh;
/* check collisions between ix-iy-iz and ix-iy-iz */
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[0][1][ix].head;
/* check collisions between ix-iy-iz and ix-iy-1-iz */
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][0][ix].head;
/* check collisions between ix-iy-iz and ix-iy-iz-1 */
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][1][ix].head;
/* check collisions between ix-iy-iz and ix-iy-1-iz-1 */
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][2][ix].head;
/* check collisions between ix-iy-iz and ix-iy+1-iz-1 */
                while (m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        if(ix==0) ipx=Lxp-1;
	        else ipx=ix-1;
	        m=listxCh[0][0][ipx].head;
/* check collisions between ix-iy-iz and ix-1-iy-iz */
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[0][1][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][0][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][1][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][2][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        if(ix==(Lxp-1)) ipx=0;
	        else ipx=ix+1;
	        m=listxCh[0][1][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][0][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][1][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        m=listxCh[1][2][ipx].head;
	        while(m!=-1){
		  ChargeVsCharge(Body[l],Body[m],l,m,Height);
		  m=Body[m].nextxCh;
	        }
	        l=Body[l].nextxCh;
	      }
	    }
	    k=Body[k].nextyCh;
	  }
	  k=listyCh[0][iy].head;
/* loop over all ix-iy-iz lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    listxCh[0][0][ix].head=-1;
	    listxCh[0][0][ix].done=0;
	    k=Body[k].nextyCh;
	  }
	  k=listyCh[1][iy].head;
/* loop over ix-iy-iz1 lists to discard them*/
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    listxCh[1][0][ix].head=-1;
	    listxCh[1][0][ix].done=0;
	    k=Body[k].nextyCh;
	  }
	  if(iy==0) ipy=Lyp-1;
	  else ipy=iy-1;
	  k=listyCh[0][ipy].head;
/* loop over all ix-iy1-iz lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    listxCh[0][1][ix].head=-1;
	    listxCh[0][1][ix].done=0;
	    k=Body[k].nextyCh;
	  }
	  k=listyCh[1][ipy].head;
/* loop over all ix-iy1-iz1 lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    listxCh[1][1][ix].head=-1;
	    listxCh[1][1][ix].done=0;
	    k=Body[k].nextyCh;
	  }
	  if(iy==(Lyp-1)) ipy=0;
	  else ipy=iy+1;
	  k=listyCh[1][ipy].head;
/* loop over all ix-iy2-iz1 lists to discard them */
	  while(k!=-1){
	    ix=static_cast<int>(Body[k].x/D_CELL);
	    listxCh[1][2][ix].head=-1;
	    listxCh[1][2][ix].done=0;
	    k=Body[k].nextyCh;
	  }
        } 
        j=Body[j].nextzCh;
      }
      j=listzCh[iz].head;
/* loop over iz-list  Bodys to discard iz-iy lists */
      while (j!=-1){
        iy=static_cast<int>(Body[j].y/D_CELL);
        listyCh[0][iy].head=-1;
        listyCh[0][iy].done=0;
        j=Body[j].nextzCh;
      }
      if(iz==0) ipz=Lzp-1;
      else ipz=iz-1;
      j=listzCh[ipz].head;
      while(j!=-1){
        iy=static_cast<int>(Body[j].y/D_CELL);
        listyCh[1][iy].head=-1;
        listyCh[1][iy].done=0;
        j=Body[j].nextzCh;
      }
    }
    }		
  }
}

void Collider::ReferenceForce(void){
  int i,j,k;
  GreenRF[0][0][0][0]=0;
  GreenRF[0][0][0][1]=0;
  GreenRF[0][0][0][2]=0;
  double ix,jy,kz;
  double denomPlus,denomLess;
  double h=1*D_CELL;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++){
        if(i==0 && j ==0 && k == 0)
        {
          GreenRF[i][j][k][0]=0;
          GreenRF[i][j][k][1]=0;
          GreenRF[i][j][k][2]=0;
        }
        else
        {
            ix=i*D_CELL;  jy=j*D_CELL;  kz=k*D_CELL;
            
            denomPlus=sqrt((ix+h)*(ix+h)+jy*jy+kz*kz);
            denomLess=sqrt((ix-h)*(ix-h)+jy*jy+kz*kz);

            if(denomPlus==0)
              GreenRF[i][j][k][0]=1/denomLess;
            else if(denomLess==0)
              GreenRF[i][j][k][0]=-1/denomPlus;
            else
              GreenRF[i][j][k][0]=1/denomLess-1/denomPlus;

            denomPlus=sqrt(ix*ix+(jy+h)*(jy+h)+kz*kz);
            denomLess=sqrt(ix*ix+(jy-h)*(jy-h)+kz*kz);

            if(denomPlus==0)
              GreenRF[i][j][k][1]=1/denomLess;
            else if(denomLess==0)
              GreenRF[i][j][k][1]=-1/denomPlus;
            else
              GreenRF[i][j][k][1]=1/denomLess-1/denomPlus;

            denomPlus=sqrt(ix*ix+jy*jy+(kz+h)*(kz+h));
            denomLess=sqrt(ix*ix+jy*jy+(kz-h)*(kz-h));

            if(denomPlus==0)
              GreenRF[i][j][k][2]=1/denomLess;
            else if(denomLess==0)
              GreenRF[i][j][k][2]=-1/denomPlus;
            else
              GreenRF[i][j][k][2]=1/denomLess-1/denomPlus;

              GreenRF[i][j][k][0]=0.5*GreenRF[i][j][k][0]/(h);
              GreenRF[i][j][k][1]=0.5*GreenRF[i][j][k][1]/(h);
              GreenRF[i][j][k][2]=0.5*GreenRF[i][j][k][2]/(h);
        }
      }
}

void Collider::BodyVsBody(Body &Bodyi,Body &Bodyj,int i, int j, double Height){
  double Ri,Rj,radij,separation,h;	
  double rxij,ryij,rzij;
  int min=i, max=j;
  int col,alt;
  double fx,fy,fz;
  double qxy=Bodyi.q*ELEMENTARYCHARGE*Bodyj.q*ELEMENTARYCHARGE;
  if(min>max){
    min=j;
    max=i;
  }
  col=(2*N-min-1)*min/2+max-min-1;

  Ri=Bodyi.R;
  Rj=Bodyj.R;
  radij=Ri+Rj;

  rxij=Bodyi.x-Bodyj.x;
  if(rxij>0.5*WIDTH)rxij-=WIDTH;
  else if(rxij<-0.5*WIDTH)rxij+=WIDTH;

  ryij=Bodyi.y-Bodyj.y;
  #ifndef BCY
  if(ryij>0.5*HEIGHT)ryij-=HEIGHT;
  else if(ryij<-0.5*HEIGHT)ryij+=HEIGHT;
  #endif

  rzij=Bodyi.z-Bodyj.z;
  if(rzij>0.5*DEEP)rzij-=DEEP;
  else if(rzij<-0.5*DEEP)rzij+=DEEP;

  separation=sqrt(rxij*rxij+ryij*ryij+rzij*rzij);
  h=radij-separation;

  alt=Bodyi.y-Ri;
  if(alt>Bodyj.y-Rj)  
    alt=Bodyj.y-Rj;
  int altint = static_cast<int>(alt);

  if(h>0){

     double Vxij, Vyij, Vzij;
     Vxij= Bodyi.Vx-Bodyj.Vx;
     Vyij= Bodyi.Vy-Bodyj.Vy;
     Vzij= Bodyi.Vz-Bodyj.Vz;
    
    fx=K*(h*rxij)/separation-GAMMA*(Vxij);
    fy=K*(h*ryij)/separation-GAMMA*(Vyij);
    fz=K*(h*rzij)/separation-GAMMA*(Vzij);

    Bodyi.Fx+=fx;	     Bodyj.Fx-=fx;
    Bodyi.Fy+=fy;	     Bodyj.Fy-=fy;
    Bodyi.Fz+=fz;	     Bodyj.Fz-=fz;
    
    if(CollitionsGrain[col].hold==0 && altint*1.0 > 0.4*Height){

      double Ai,Aj;
      long long int ions_i=0,ions_j=0;
      int Diff_Transfer=0;
      Ai=2*M_PI*Bodyi.R*Bodyi.R*(-D_BRIDGE*(2*Bodyj.R-D_BRIDGE)/(2*Bodyi.R*(D_BRIDGE-Bodyi.R-Bodyj.R)));
      Aj=2*M_PI*Bodyj.R*Bodyj.R*(-D_BRIDGE*(2*Bodyi.R-D_BRIDGE)/(2*Bodyj.R*(D_BRIDGE-Bodyi.R-Bodyj.R)));
      ions_i=static_cast<long long int>(Ai*((Bodyi.H+Bodyi.q)/(4*M_PI*Bodyi.R*Bodyi.R)));
      ions_j=static_cast<long long int>(Aj*((Bodyj.H+Bodyj.q)/(4*M_PI*Bodyj.R*Bodyj.R)));
      Diff_Transfer=sqrt((ions_i-ions_j)*(ions_i-ions_j));
      

      if(ions_i>ions_j){

        if(Bodyi.getAverageVel()>1e-1/LENGTH)
          Bodyi.IonTransferVolume(-Diff_Transfer/2.0);
        if(Bodyj.getAverageVel()>1e-1/LENGTH)
          Bodyj.IonTransferVolume(Diff_Transfer/2.0);

      }
      else{
        if(Bodyi.getAverageVel()>1e-1/LENGTH)
          Bodyi.IonTransferVolume(Diff_Transfer/2.0);
        if(Bodyj.getAverageVel()>1e-1/LENGTH)
          Bodyj.IonTransferVolume(-Diff_Transfer/2.0);
      }

      if(qxy<0){
        Bodyi.q=Bodyi.q*0.95;
        Bodyj.q=Bodyj.q*0.95;
      }
    }

    /*if(Bodyi.q!=0 && altint*1.0 < 1.2*Height)
      if(Bodyi.getAverageVel()<1e-1/LENGTH)
        Bodyi.setChargeWall();
    if(Bodyj.q!=0 && altint*1.0 < 1.2*Height)
      if(Bodyj.getAverageVel()<1e-1/LENGTH)
        Bodyj.setChargeWall();*/

    CollitionsGrain[col].hold=h;
  }
  else{
    if(CollitionsGrain[col].hold>0){

      CollitionsGrain[col].hold=0;        

    }
  } 
}

void Collider::ChargeVsCharge(Body &Bodyi,Body &Bodyj,int i, int j, double Height){
  double Ri,Rj,radij,separation,h;	
  double rxij,ryij,rzij;
  double qxy=Bodyi.q*ELEMENTARYCHARGE*Bodyj.q*ELEMENTARYCHARGE;
  double fx,fy,fz;

  Ri=Bodyi.R;
  Rj=Bodyj.R;
  radij=Ri+Rj;

  rxij=Bodyi.x-Bodyj.x;
  if(rxij>0.5*WIDTH)rxij-=WIDTH;
  else if(rxij<-0.5*WIDTH)rxij+=WIDTH;

  ryij=Bodyi.y-Bodyj.y;
  #ifndef BCY
  if(ryij>0.5*HEIGHT)ryij-=HEIGHT;
  else if(ryij<-0.5*HEIGHT)ryij+=HEIGHT;
  #endif

  rzij=Bodyi.z-Bodyj.z;
  if(rzij>0.5*DEEP)rzij-=DEEP;
  else if(rzij<-0.5*DEEP)rzij+=DEEP;

  separation=sqrt(rxij*rxij+ryij*ryij+rzij*rzij);
  h=radij-separation;

  if(Bodyi.Charged!=0 && Bodyj.Charged!=0){
    double RefForceijx=0;
    double RefForceijy=0;
    double RefForceijz=0;
    int axi,byi,czi;
    int axj,byj,czj;
    int rx,ry,rz;
    double sigx=1,sigy=1,sigz=1;

    axi=(Bodyi.a);
    byi=(Bodyi.b);
    czi=(Bodyi.c);
 
    axj=(Bodyj.a);
    byj=(Bodyj.b);
    czj=(Bodyj.c);
 
    //Reference Force Lee-Wonsok
    //all to j a,b,e to Manual Switch off signs
    ry=(byi-byj);
    if(ry<0){ry=-ry; sigy=-1;} else {sigy=1;}
    rz=(czi-czj);
    if(rz>0.5*Lzp)rz-=Lzp;
    else if(rz<-0.5*Lzp)rz+=Lzp;
    if(rz<0){rz=-rz; sigz=-1;} else {sigz=1;}

    if(axi!=axj || byi!=byj || czi!=czj){

      rx=(axi-axj);
      if(rx>0.5*Lxp)rx-=Lxp;
      else if(rx<-0.5*Lxp)rx+=Lxp;
      if(rx<0){rx=-rx; sigx=-1;} else {sigx=1;}

      RefForceijx+=sigx*GreenRF[rx][ry][rz][0];
      RefForceijy+=sigy*GreenRF[rx][ry][rz][1];
      RefForceijz+=sigz*GreenRF[rx][ry][rz][2];
    }

    Bodyi.Fx-=RefForceijx*qxy*KE;	     Bodyj.Fx+=RefForceijx*qxy*KE;
    Bodyi.Fy-=RefForceijy*qxy*KE;	     Bodyj.Fy+=RefForceijy*qxy*KE;
    Bodyi.Fz-=RefForceijz*qxy*KE;	     Bodyj.Fz+=RefForceijz*qxy*KE;

//    if(separation<radij/2.0)
  //    clog<<"Warning: Penetrating particles: "<<Bodyi.x<<"\t"<<Bodyi.y<<"\t"<<Bodyi.z<<"\t"<<Bodyj.x<<"\t"<<Bodyj.y<<"\t"<<Bodyj.z<<endl;
    if(h>0){
      separation=radij;
      //if(qxy<-MaxQ) qxy=-MaxQ;
      //if(qxy>MaxQ) qxy=MaxQ;
    }

    fx=(KE*qxy*rxij)/(separation*separation*separation);
    fy=(KE*qxy*ryij)/(separation*separation*separation);
    fz=(KE*qxy*rzij)/(separation*separation*separation);
    if(maxQQ<sqrt(fx*fx+fy*fy+fz*fz)){maxQQ=sqrt(fx*fx+fy*fy+fz*fz);Q1=Bodyi.q;Q2=Bodyj.q;}

     
    Bodyi.Fx+=fx;	     Bodyj.Fx-=fx;
    Bodyi.Fy+=fy;	     Bodyj.Fy-=fy;
    Bodyi.Fz+=fz;	     Bodyj.Fz-=fz;
  }
}
void Collider::Wall(Body *Body, int i){
  double h;
  double Vxij, Vyij, Vzij;
  double fx,fy,fz;
  double ryij;

  ryij=Body[i].y;
  Vxij=Body[i].Vx;
  Vyij=Body[i].Vy;
  Vzij=Body[i].Vz;

  h=Body[i].R-ryij;

  fx=-GAMMAWALL*(Vxij);
  fy=KWALL*h-GAMMAWALL*(Vyij);
  fz=-GAMMAWALL*(Vzij);

  Body[i].Fx+=fx;
  Body[i].Fy+=fy;
  Body[i].Fz+=fz;
  
}
void Collider::WallUp(Body *Body, int i){
  double h;
  double Vxij, Vyij, Vzij;
  double fx,fy,fz;
  double ryij;

  ryij=Body[i].y-HEIGHT;
  Vxij=Body[i].Vx;
  Vyij=Body[i].Vy;
  Vzij=Body[i].Vz;

  h=Body[i].R-sqrt(ryij*ryij);
  fx=-GAMMAWALL*(Vxij);
  fy=-KWALL*h-GAMMAWALL*(Vyij);
  fz=-GAMMAWALL*(Vzij);

  Body[i].Fx+=fx;
  Body[i].Fy+=fy;
  Body[i].Fz+=fz;
}

void Collider::ChargeAssignmentDirectEF(Body *Body){
  int i,j,k,n,a,b,c,num=0;
  double af,bf,cf;
  //Cloud In Cell Algorithm  
   for(j=0;j<Lyp;j++)
      for(i=0;i<Lxp;i++)
        for(k=0;k<Lzp;k++){
          Done[i][j][k]=0;
          SolveCell[i][j][k]=0;
          potential[i][j][k]=0.;
          ChargeGrid[i][j][k]=0.;
          smallArray[num]=-1;
          num++;
       }
  RUNTIME++;
  for(n=0;n<N;n++){
    if(Body[n].q!=0){
      a=static_cast<int>(Body[n].x/D_CELL);
      af = (Body[n].x/D_CELL - 1.0*a);
      if(af>0.5)
      a=(a+1)%Lxp;
      b=static_cast<int>(Body[n].y/D_CELL);
      bf = (Body[n].y/D_CELL - 1.0*b);
      if(bf>0.5){
      #ifdef BCY
        if(b+1<=Lyp-1)
          b=(b+1);
      #else
        b=(b+1)%Lyp;
      #endif
      }
      c=static_cast<int>(Body[n].z/D_CELL);
      cf = (Body[n].z/D_CELL - 1.0*c);
      if(cf>0.5)
      c=(c+1)%Lzp;
      ChargeGrid[a][b][c]+=Body[n].q;
      Carga[a][b][c]+=Body[n].q;
      Body[n].a=a;
      Body[n].b=b;
      Body[n].c=c;

      SolveCell[a][b][c]=1;
      #ifdef BCY

      if(b==Lyp-1){
        SolveCell[a][b-1][c]=1;
      }
      else if(b==0){
        SolveCell[a][b+1][c]=1;
      }
      else{
        SolveCell[a][b-1][c]=1;
        SolveCell[a][b+1][c]=1;
     }
      if(a==Lxp-1){
        SolveCell[a-1][b][c]=1;
        SolveCell[0][b][c]=1;
      }
      else if(a==0){
        SolveCell[Lxp-1][b][c]=1;
        SolveCell[a+1][b][c]=1;
      }
      else{
        SolveCell[a+1][b][c]=1;
        SolveCell[a-1][b][c]=1;
     }
      if(c==Lzp-1){
        SolveCell[a][b][c-1]=1;
        SolveCell[a][b][0]=1;
      }
      else if(c==0){
        SolveCell[a][b][Lzp-1]=1;
        SolveCell[a][b][c+1]=1;
      }
      else{
        SolveCell[a][b][c+1]=1;
        SolveCell[a][b][c-1]=1;
     }
     #else
     if(b==Lyp-1){
        SolveCell[a][b-1][c]=1;
        SolveCell[a][0][c]=1;
      }
      else if(b==0){
        SolveCell[a][b+1][c]=1;
        SolveCell[a][Lyp-1][c]=1;
      }
      else{
        SolveCell[a][b-1][c]=1;
        SolveCell[a][b+1][c]=1;
     }
      if(a==Lxp-1){
        SolveCell[a-1][b][c]=1;
        SolveCell[0][b][c]=1;
      }
      else if(a==0){
        SolveCell[Lxp-1][b][c]=1;
        SolveCell[a+1][b][c]=1;
      }
      else{
        SolveCell[a+1][b][c]=1;
        SolveCell[a-1][b][c]=1;
     }
      if(c==Lzp-1){
        SolveCell[a][b][c-1]=1;
        SolveCell[a][b][0]=1;
      }
      else if(c==0){
        SolveCell[a][b][Lzp-1]=1;
        SolveCell[a][b][c+1]=1;
      }
      else{
        SolveCell[a][b][c+1]=1;
        SolveCell[a][b][c-1]=1;
     }
     #endif
    }
  }

   num=0;
   sizeSA=0;
   for(j=0;j<Lyp;j++)
      for(i=0;i<Lxp;i++)
        for(k=0;k<Lzp;k++){
            if(SolveCell[i][j][k]==1){
              smallArray[num]=i+Lxp*j+Lxp*Lyp*k;
              num++;
            }
        }
    sizeSA=num;
}
void Collider::GreenFunctionReal(void){
  int j,i,k;
  int x,z;
  int images=1000;
  double dx,dz;
  for(k=0;k<Lzp;k++)
    for(j=0;j<Ly;j++)
      for(i=0;i<Lxp;i++){
         if(i==0 && j==0 && k==0){
           GreenF[i][j][k]=0;
           for(x=-images;x<=images;x++)
             for(z=-6*images;z<=6*images;z++){
                  dx=Lxp*x;
                  dz=Lzp*z;  
                  if(dx==0 && dz==0) 
                    GreenF[i][j][k]+=0.0;
		  else
                    GreenF[i][j][k]+=1.0/sqrt(dx*dx+dz*dz);
             }
         }
         else{
           GreenF[i][j][k]=0;
           for(x=-images;x<=images;x++)
             for(z=-6*images;z<=6*images;z++){
                  dx=i+Lxp*x;
                  dz=k+Lzp*z;  
                  GreenF[i][j][k]+=1.0/sqrt(dx*dx+j*j+dz*dz);
           }
         }
      }
}
void Collider::PotentialSolver(void){
  int l,m,n,j,i,k;
  int GreenX,GreenY,GreenZ;
  for(k=0;k<Lzp;k++)
    for(j=0;j<Lyp;j++)
      for(i=0;i<Lxp;i++){
        if(SolveCell[i][j][k]==1)
          for(n=0;n<Lzp;n++)
            for(m=0;m<Lyp;m++)
              for(l=0;l<Lxp;l++){
                if(SolveCell[l][m][n]==1){
                  double PotAux=0;
                  GreenX=static_cast<int>(sqrt(1.0*(l-i)*(l-i))); GreenY=static_cast<int>(sqrt(1.0*(m-j)*(m-j))); GreenZ=static_cast<int>(sqrt(1.0*(n-k)*(n-k)));
                  if(GreenX>0.5*Lxp)GreenX=Lxp-GreenX;
                  if(GreenZ>0.5*Lzp)GreenZ=Lzp-GreenZ;
                  PotAux=GreenF[GreenX][GreenY][GreenZ];
                  potential[i][j][k]+=PotAux*ChargeGrid[l][m][n]/D_CELL;
               }
        }
    }
   
}
void Collider::PotentialSolverParallel(void){
  int Area=Lxp*Lyp;

  for(unsigned int ijk=0;ijk<sizeSA;ijk++){
    int num = smallArray[ijk];
    int i=num%Lxp;
    int j=static_cast<int>((num%(Area))/Lxp);
    int k=static_cast<int>(num/Area);
      //if(j!=0)							/*Bed Charge ---- No need to solve (No "Charged Particles" there)*/
      for(unsigned int lmn=0;lmn<sizeSA;lmn++){
        int num2 = smallArray[lmn];
        int l=num2%Lxp;
        int m=static_cast<int>((num2%(Area))/Lxp);
        int n=static_cast<int>(num2/Area);
          int GreenX=static_cast<int>(sqrt(1.0*(l-i)*(l-i))); int GreenY=static_cast<int>(sqrt(1.0*(m-j)*(m-j))); int GreenZ=static_cast<int>(sqrt(1.0*(n-k)*(n-k)));
          //int GreenM=static_cast<int>(sqrt(1.0*(-m-j)*(-m-j)));
          if(GreenX>0.5*Lxp)GreenX=Lxp-GreenX;
          if(GreenZ>0.5*Lzp)GreenZ=Lzp-GreenZ;
          potential[i][j][k]+=(GreenF[GreenX][GreenY][GreenZ]/*-GreenF[GreenX][GreenM][GreenZ]*/)*ChargeGrid[l][m][n]/D_CELL;
    }
  }
   
}

void Collider::ElectricFieldReal(int i, int j, int k){
  if(Done[i][j][k]==0){
  double a,b,c,d,e,f;
  double deltaF=-0.5/(1.*D_CELL);
  #ifdef BCY
    a=potential[(i+1)%Lxp][j][k];
    b=potential[(i-1+Lxp)%Lxp][j][k];
    e=potential[i][j][(k+1)%Lzp];
    f=potential[i][j][(k-1+Lzp)%Lzp];

    if(j==Lyp-1)
      c=0;
    else
      c=potential[i][(j+1)%Lyp][k];

    if(j==0)
      d=0;
    else
      d=potential[i][(j-1+Lyp)%Lyp][k];
     
  #else
    a=potential[(i+1)%Lxp][j][k];
    b=potential[(i-1+Lxp)%Lxp][j][k];
    c=potential[i][(j+1)%Lyp][k];
    d=potential[i][(j-1+Lyp)%Lyp][k];
    e=potential[i][j][(k+1)%Lzp];
    f=potential[i][j][(k-1+Lzp)%Lzp];
  #endif
  Ex[i][j][k]=(a-b)*deltaF;
  #ifdef BCY
  if(j==0 )
    Ey[i][j][k]=(c)*deltaF;
  else if( j == Lyp-1)
    Ey[i][j][k]=0;
  else
  Ey[i][j][k]=(c-d)*deltaF;
  #else
  Ey[i][j][k]=(c-d)*deltaF;
  #endif
  Ez[i][j][k]=(e-f)*deltaF;
  Done[i][j][k]=1;
  }
}


void Collider::ElectricForceReal(Body *Body){
  int a,b,c,n;
  double C1;
  //NGP Algorithm  
  for(n=0;n<N;n++){
    if(Body[n].q!=0){

      a=Body[n].a;
      b=Body[n].b;
      c=Body[n].c;

      C1=ELEMENTARYCHARGE*ELEMENTARYCHARGE*Body[n].q*KE;

      ElectricFieldReal(a,b,c);

      Body[n].Fx+=Ex[a][b][c]*C1;

      Body[n].Fy+=Ey[a][b][c]*C1;

      Body[n].Fz+=Ez[a][b][c]*C1;
    }

  }

}

void Collider::Collide(Body *Body, double Height){
  int i;
  
  for(i=0;i<N;i++){
    Body[i].Fx=0;
    Body[i].Fy=-GRAVITY*Body[i].m;
    Body[i].Fz=0;
    Body[i].nextx=-1;
    Body[i].nexty=-1;
    Body[i].nextz=-1;
    Body[i].nextxCh=-1;
    Body[i].nextyCh=-1;
    Body[i].nextzCh=-1;
    Body[i].Checked=0;
    if(Body[i].q==0)
      Body[i].Charged=0;
    else
      Body[i].Charged=1;
    #ifdef BCY
    if(Body[i].y-Body[i].R < 0.0){
       Wall(Body, i);
    }
    else if(Body[i].R+Body[i].y-HEIGHT > 0.0){
       WallUp(Body, i);
    }
    #endif
  }

  #ifdef ELECF
    ChargeAssignmentDirectEF(Body);
    PotentialSolverParallel();
    ElectricForceReal(Body);
    ChargeDetection(Body,Height);
  #endif
  ContactDetection(Body,Height);

  
}
void Collider::PrintPotential(ofstream &fout,Body *Body){
  int i,j,k;
  j=Lyp/2;
  for(k=0;k<Lzp;k++)
    for(i=0;i<Lxp;i++){
      fout<<scientific<<i<<"\t"<<k<<"\t"<<potential[i][j][k]*ELEMENTARYCHARGE*LENGTH*LENGTH/(4.*M_PI*E0)<<"\t"<<KE*ELEMENTARYCHARGE*LENGTH*LENGTH*Body[0].q/sqrt((i*D_CELL-Body[0].x)*(i*D_CELL-Body[0].x)+(j*D_CELL-Body[0].y)*(j*D_CELL-Body[0].y)+(k*D_CELL-Body[0].z)*(k*D_CELL-Body[0].z))+KE*ELEMENTARYCHARGE*LENGTH*LENGTH*Body[1].q/sqrt((i*D_CELL-Body[1].x)*(i*D_CELL-Body[1].x)+(j*D_CELL-Body[1].y)*(j*D_CELL-Body[1].y)+(k*D_CELL-Body[1].z)*(k*D_CELL-Body[1].z))<<"\t"<<GreenF[i][j][k]<<endl;
    }
}

#ifdef ELECF
void Collider::PrintGreenFunction(ofstream &fout){
  for(int i=0;i<Lxp;i++)
    for(int j=0;j<Ly;j++)
      for(int k=0;k<Lzp;k++)
	fout<<scientific<<i<<"\t"<<j<<"\t"<<k<<"\t"<<GreenF[i][j][k]<<endl;
}
void Collider::ReadGreenFunction(char *infile){
  FILE *fin=NULL;
  if((fin=fopen(infile, "r"))==NULL){
    printf("Cannot Open Green Function File\n");
    exit(1);	
  }
  double GF=-1;
  double ii,jj,kk;
  for(int i=0;i<Lxp;i++)
    for(int j=0;j<Ly;j++)
      for(int k=0;k<Lzp;k++){
        fscanf(fin, "%lf	%lf	%lf	%lf\n", &ii, &jj, &kk, &GF);
	GreenF[i][j][k]=GF;
      }
}

#endif
void Collider::PrintCharge(ofstream &fout){
  int i,j,k;
  for(k=0;k<Lzp;k++)
  for(j=0;j<Lyp;j++)
    for(i=0;i<Lxp;i++){
      fout<<scientific<<i<<"\t"<<j<<"\t"<<k<<"\t"<<Carga[i][j][k]*ELEMENTARYCHARGE/RUNTIME<<endl;
      Carga[i][j][k]=0;
    }
  RUNTIME=0;
}

void Collider::PrintElectric(Body *Body){
    ChargeAssignmentDirectEF(Body);
}

void Collider::PotentialSolverTotal(double ***Carga){
  int l,m,n,j,i,k;
  int GreenX,GreenY,GreenZ,GreenM;
  for(k=0;k<Lzp;k++)
    for(j=0;j<Lyp;j++)
      for(i=0;i<Lxp;i++){
          for(n=0;n<Lzp;n++)
            for(m=0;m<Lyp;m++)
              for(l=0;l<Lxp;l++){
                  GreenX=static_cast<int>(sqrt(1.0*(l-i)*(l-i))); GreenY=static_cast<int>(sqrt(1.0*(m-j)*(m-j))); GreenZ=static_cast<int>(sqrt(1.0*(n-k)*(n-k)));
		  GreenM=static_cast<int>(sqrt(1.0*(-m-j)*(-m-j))); 
                  potential[i][j][k]+=(GreenF[GreenX][GreenY][GreenZ]-GreenF[GreenX][GreenM][GreenZ])*Carga[l][m][n]/D_CELL;
        }
    }

}

void Collider::ElectricFieldTotal(){
  for(int k=0;k<Lzp;k++)
  for(int i=0;i<Lxp;i++)
  for(int j=0;j<Lyp;j++){
  double a,b,c,d,e,f;
  double deltaF=-0.5/(1.*D_CELL);
  #ifdef BCY
    a=potential[(i+1)%Lxp][j][k];
    b=potential[(i-1+Lxp)%Lxp][j][k];
    e=potential[i][j][(k+1)%Lzp];
    f=potential[i][j][(k-1+Lzp)%Lzp];

    if(j==Lyp-1)
      c=0;
    else
      c=potential[i][(j+1)%Lyp][k];

    if(j==0)
      d=0;
    else
      d=potential[i][(j-1+Lyp)%Lyp][k];
     
  #else
    a=potential[(i+1)%Lxp][j][k];
    b=potential[(i-1+Lxp)%Lxp][j][k];
    c=potential[i][(j+1)%Lyp][k];
    d=potential[i][(j-1+Lyp)%Lyp][k];
    e=potential[i][j][(k+1)%Lzp];
    f=potential[i][j][(k-1+Lzp)%Lzp];
  #endif
  Ex[i][j][k]=(a-b)*deltaF;
  #ifdef BCY
  if(j==0 )
    Ey[i][j][k]=(c)*deltaF;
  else if( j == Lyp-1)
    Ey[i][j][k]=0;
  else
    Ey[i][j][k]=(c-d)*deltaF;
  #else
  Ey[i][j][k]=(c-d)*deltaF;
  #endif
  Ez[i][j][k]=(e-f)*deltaF;
  }
}
