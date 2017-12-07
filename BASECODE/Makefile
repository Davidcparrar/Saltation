#Compiler
CC=g++
# Output -fopenmp -g -pg -lfftw3 -lm 
CFLAGS= -o Saltation -DBCY -Wall -O3 -DELECF -DBINDIST -DWIND -DPGF#-DPGF -O3
CFLAGSB=-o Process -Wall
CFLAGSC=-o EF -DELECF -DBCY -g -Wall -DRFOR
CFLAGSD=-o Concentration -DBINDIST -g -Wall
CFLAGSF=-o Flux
CFLAGSE=-o HEF -DELECF -DBCY -g -Wall -DRFOR

OUTPUT=> Sim.sim   2> DatosSaltation.txt

all:
	$(CC) main.cpp $(CFLAGS)  
	$(CC) $(CFLAGSB) CPP/ProcesingData.cpp 
	$(CC) $(CFLAGSC) CPP/ElectricField.cpp 
	$(CC) $(CFLAGSD) CPP/ConcentrationProfiles.cpp 
	$(CC) $(CFLAGSF) CPP/FluxAnalysis.cpp 
	$(CC) $(CFLAGSE) CPP/ElectricFieldHeight.cpp 
clean:
	-$(RM)  *.txt FILES/fil* FORCES/Force* PROFILES/Profile* PROFILES/ChargeMap* CONCENTRATION/Concentration* *~ FRAMES/frame* VelProfile* Saltation Process EF Concentration

