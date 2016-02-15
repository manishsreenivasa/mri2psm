/* MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany

Licensed under the zlib license. See LICENSE for more details.

This program computes inverse dynamics torques using the RBDL dynamics library. Inputs are the lua model, inverse kinematics (joint angles) and external forces (ground reaction forces in c3d file). Outputs are saved to text files, and may be plotted using the matlab scipts available in the subfolder "matlab_plot_functions".

In addition to RBDL, you will need BTK, Eigen3 and Boost to compile and run this code.
*/

#include <fstream>
#include "csvInfo.h"
#include "ID.h"

using namespace std;

bool app_analyze_options(int argc, char** argv, char** filenames);  

int main (int argc, char* argv[]) {
	char *filenames[argc];
	bool ForcesOn = app_analyze_options(argc, argv, filenames);
	ID modelID(filenames[1]);
	modelID.ReadQValuesFromFile(filenames[2]);
	if (ForcesOn)
		modelID.ReadForcesAndStepsFromFile(filenames[3],0);     
	else
		cerr<<"The c3d file is missing, the forces will be set to zero"<<endl;
	modelID.CalculateTauWithInverseDynamics();
	modelID.WriteTauValuesToFile("./id_res.txt");
	return 0;
}

bool app_analyze_options(int argc, char** argv, char** filenames) {
	bool ConsideringForces = 0;
	if (argc <3) {
		cerr << "There have to be at least 2 arguments with the formats .lua and .csv"<<endl;
		cerr << "One extra argument with the format .c3d containing the forces can be included"<<endl;
	}    
	for(int c=1; c<argc; c++) {
		if(strstr(argv[c], ".lua")>0) {filenames[1] = argv[c];}
		else if(strstr(argv[c], ".csv")>0) {filenames[2] = argv[c];}
		else if(strstr(argv[c], ".c3d")>0) {
			filenames[3] = argv[c];
			ConsideringForces = 1;
		}
		else {	
			cerr << "*** WARNING: undefined option: " << argv[c] << endl;
			cerr << "There have to be 3 arguments with the following formats: .lua .csv .c3d "<<endl;
		}
	}
	return ConsideringForces;
}
