#include <iostream>
#include <cmath>

using namespace std;
int counterofCollisions=0;

const double MU=0.3;
int N=500;				/* m */
double D_MEAN=0.0002;					/* m */
double ROUGHNESS=D_MEAN/30.0;					/* m */
double RADIUS_MAX;
int NUMBEROFIMAGES=100;

double DENSITY=2650;						/* kg/m^3	1000				*/
double AIRDENSITY=1.174;
double AIRVISCOSITY=1.8702e-05;			/* kg/m^3	1000				*/
double HEIGHT=2004*D_MEAN,WIDTH=48*D_MEAN,DEEP=12*D_MEAN;
double WINDVEL=0.9;
double DELTAT=1e-6;
double GRAVITY=9.81;
double LENGTH,MASS,TIME,AMPERE;
double GAMMA =0.00012;
double GAMMAWALL = 0.000033;
double K=500;
double KWALL=1.0;
double KARMAN = 0.4;
double PLANCK = 6.626068e-34;
double REDUCEDPLANCK = PLANCK/(2.0*M_PI);
double ELECTRONMASS = 9.10938188e-31;
double ELEMENTARYCHARGE = 1.60217646e-19;
double ELECTRONVOLT = 1.602176565e-19;
double Eb0 =0.05*ELECTRONVOLT;
double EbMAX =10.*ELECTRONVOLT;
double ETA = (1.12/(M_PI*M_PI))*(2.0+M_PI);
double RHO_H = 1e+23;
double TCOLL=1e-8;
double A2=1e-20;
double E0=8.85418782e-12;
double KE = 1.0/(4.0*M_PI*E0);
double LIGHT = 1.0/(D_MEAN*2.)*DELTAT;
double RATE_NEUTRAL = 10;
double MOBILITY = 1e-10;//1e-12;
double ALPHA =1e-10;
double D_BRIDGE=2e-9;
double D_BRIDGEK=D_BRIDGE/2.0;
double MaxQ=1e12*ELEMENTARYCHARGE*ELEMENTARYCHARGE;
