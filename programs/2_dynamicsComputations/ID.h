/* MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany

Licensed under the zlib license.
*/
#include <btkConfigure.h>
#include <btkAcquisitionFileReader.h>	
#include <btkForcePlatformsExtractor.h>
#include <btkGroundReactionWrenchFilter.h>

#include <rbdl/rbdl.h>
#include <rbdl/addons/luamodel/luamodel.h>
#include <rbdl/addons/luamodel/luatables.h>

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Utils;

class ID {
friend class ForcePlatformsExtractor;
	vector<VectorNd> Q;                      
	vector<VectorNd> QDot;			
	vector<VectorNd> QDDot;		
	vector<VectorNd> Tau;   	
	double T0;
	double Tf;
	double delta_t;
	int NumT;
 	Vector3d Footstartposition_r;
        Vector3d Footstartposition_l;
	int NumQ;
	int FootPrintLeftInterval;
	int FootPrintRightInterval;		
	int InitialLeftStep;
	int InitialRightStep;	
	int EndRightStep;
	int EndLeftStep;
	unsigned int Id_Left;
	unsigned int Id_Right;                   
	vector<VectorNd> LForcePlate;
	vector<VectorNd> RForcePlate;	         
	vector<VectorNd> LInverseForceMomentum;
	vector<VectorNd> RInverseForceMomentum;	
	Model* model;

	public:
		ID(const char *FilenameLuamodel);
		void ReadQValuesFromFile(const char *filename);
		void ReadForcesAndStepsFromFile(const string& filename, bool ExtraSteps); 		
		void CalculateTauWithInverseDynamics();

		void WriteQValuesToFile(const char *filenameOfOutputFile);
		void WriteQDotValuesToFile(const char *filenameOfOutputFile);
		void WriteQDDotValuesToFile(const char *filenameOfOutputFile);
		void WriteTauValuesToFile(const char *filenameOfOutputFile);
		void WriteLForcePlateValuesToFile(const char *filenameOfOutputFile);
		void WriteRForcePlateValuesToFile(const char *filenameOfOutputFile);	

		void PrintLabelsFromc3dFile(const string& filename);   
   		~ID();	

	private:
		void GetDerivativeswithFiniteDifferences(); 
		void PrintObjectofTypevectorVectorNDtoFile(vector<VectorNd>* PointerToValues, const char *filename);
		
		btk::Acquisition::Pointer readAcquisition(const string& filename);  
		btk::WrenchCollection::Pointer readbtkForces(const string& filename);  
		std::string GetPointTypeAsString(btk::Point::Type t);		  
		int ReadStepInterval(btk::Acquisition::Pointer Acq); 
		VectorNd SetupApplicationPoint(VectorNd ForceToApply, VectorNd InstantPosition, Vector3d ApplicationPoint, unsigned int BodyPartId); 
		vector<SpatialVector>  SetupExternalForces(VectorNd ForcesMomentumsLeft,VectorNd ForcesMomentumsRight);
}; // End of class ID

ID::ID(const char *FilenameLuamodel) {
	model = new Model();
	Addons::LuaModelReadFromFile(FilenameLuamodel, model, false);
	NumQ=model->dof_count;
	Id_Left =model->GetBodyId("foot_l");
	Id_Right=model->GetBodyId("foot_r");
}

void ID::ReadQValuesFromFile(const char *filenameAnimationFile) { 
	csvInfo Animation (filenameAnimationFile); 
	Animation.loadCSV();
	T0 = Animation.GetStartTime();
	Tf = Animation.GetEndTime();
	NumT = Animation.GetTimeSteps();
	Q = Animation.GetQ();	

	if (NumQ !=Animation.GetDegreesOfFreedom()) {   
		cerr << "*** WARNING: Animation Degrees Of Freedom: " << Animation.GetDegreesOfFreedom() << endl;
		cerr << "*** Number of Q in the model: " << NumQ << endl;
	}
	delta_t = (Tf-T0)/(NumT-1);
	
	for (int i=0;i<NumT;i++) {
		QDot.push_back(VectorNd::Zero(NumQ));
		QDDot.push_back(VectorNd::Zero(NumQ));
		LForcePlate.push_back(VectorNd::Zero(6));
		RForcePlate.push_back(VectorNd::Zero(6));
		LInverseForceMomentum.push_back(VectorNd::Zero(6));		
		RInverseForceMomentum.push_back(VectorNd::Zero(6));
	}
	GetDerivativeswithFiniteDifferences();
}                   

void ID::WriteQValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Q, filenameOfOutputFile);
}
void ID::WriteQDotValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&QDot, filenameOfOutputFile);
}
void ID::WriteQDDotValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&QDDot, filenameOfOutputFile);
}
void ID::WriteTauValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Tau, filenameOfOutputFile);
}
void ID::WriteLForcePlateValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Tau, filenameOfOutputFile);
}
void ID::WriteRForcePlateValuesToFile(const char *filenameOfOutputFile) {
	PrintObjectofTypevectorVectorNDtoFile(&Tau, filenameOfOutputFile);
}

void ID::CalculateTauWithInverseDynamics() {
  
	vector<Math::VectorNd> f_left_mat;
	vector<Math::VectorNd> f_right_mat;

	for(int i=0; i<NumT; i++) {
	  Tau.push_back(VectorNd::Zero(Q[0].size()));
	  vector<Math::SpatialVector> f = SetupExternalForces(LInverseForceMomentum[i],RInverseForceMomentum[i]);   
	  f_left_mat.push_back(f[Id_Left]);
	  f_right_mat.push_back(f[Id_Right]);
	  InverseDynamics(*model, Q[i], QDot[i], QDDot[i], Tau[i],&f);
	} 
}    

ID::~ID(){}

void ID::ReadForcesAndStepsFromFile(const string& filename, bool ExtraSteps) { 

	btk::Acquisition::Pointer acq = readAcquisition(filename);		
	ReadStepInterval(acq);

	int Nframes = (acq->GetPoint("SACR"))->GetFrameNumber();
	if (Nframes != NumT) {
		cerr << "*** WARNING: The number of frames in the c3d file doesn`t match the ones in the animation" << endl;
	}

	Eigen::Matrix<double,Eigen::Dynamic,3> LForce;
	Eigen::Matrix<double,Eigen::Dynamic,3> RForce;

	btk::WrenchCollection::Pointer fpws = readbtkForces(filename);
	btk::WrenchCollection::ConstIterator itWrench = fpws->Begin();
	int numberOfForcePlates = fpws->GetItemNumber();
	int  numberOfFrames = (*itWrench)->GetPosition()->GetFrameNumber();
	int SamplePerFrame = static_cast<double>(acq->GetNumberAnalogSamplePerFrame()); 	
	if (numberOfFrames != Nframes*SamplePerFrame) {
		cerr << "*** WARNING: The number of frames in the c3d file doesn`t match the ones in the resample forces" << endl;
	}
	Eigen::MatrixXd force_array[numberOfForcePlates];
	Eigen::MatrixXd moment_array[numberOfForcePlates];
	Eigen::MatrixXd origin_array[numberOfForcePlates];
	Eigen::MatrixXd resize_force_array[numberOfForcePlates];
	Eigen::MatrixXd resize_moment_array[numberOfForcePlates];
	Eigen::MatrixXd resize_origin_array[numberOfForcePlates];
	Eigen::MatrixXd position_array[numberOfForcePlates];
	Eigen::MatrixXd Rstep_PlateValues[numberOfForcePlates];
	Eigen::MatrixXd Lstep_PlateValues[numberOfForcePlates];
	int R_plate_id = 0;
	int L_plate_id = 0;
	double Rstep_Norm[numberOfForcePlates];
	double Lstep_Norm[numberOfForcePlates];	
	for (int i = 0; i < numberOfForcePlates; ++i) {
		resize_force_array[i]  = Eigen::MatrixXd::Zero(Nframes,3);
		resize_moment_array[i] = Eigen::MatrixXd::Zero(Nframes,3);
		resize_origin_array[i] = Eigen::MatrixXd::Zero(Nframes,3);
		position_array[i]      = Eigen::MatrixXd::Zero(Nframes,3);
		Rstep_PlateValues[i]   = Eigen::MatrixXd::Zero(FootPrintRightInterval,3);
		Lstep_PlateValues[i]   = Eigen::MatrixXd::Zero(FootPrintLeftInterval,3);
	}
	for (int i = 0; i < numberOfForcePlates; ++i) {
		if (itWrench->get() != 0) {
			force_array[i]  = Eigen::MatrixXd::Zero(numberOfFrames,3);
			force_array[i]  = (*itWrench)->GetForce()->GetValues();
			moment_array[i] = Eigen::MatrixXd::Zero(numberOfFrames,3);
			moment_array[i] = (*itWrench)->GetMoment()->GetValues();
			origin_array[i]  = Eigen::MatrixXd::Zero(numberOfFrames,3);
			origin_array[i]  = (*itWrench)->GetPosition()->GetValues();
		}
		else {
	  		cerr << "Force platform wrench is empty." << endl;
		}
		++itWrench;
	}
		
	for (int i = 0; i < numberOfForcePlates; ++i) {
		for (int k=0;k<3; ++k) {	
			for (int j=1;j<(Nframes-1); ++j) {	
				for (int n=-(SamplePerFrame/2);n< (SamplePerFrame-SamplePerFrame/2); ++n) {
					resize_force_array[i](j,k)    = resize_force_array[i](j,k)    +  (force_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
					resize_moment_array[i](j,k)   = resize_moment_array[i](j,k)   + (moment_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
					resize_origin_array[i](j,k) = resize_origin_array[i](j,k) + (origin_array[i](SamplePerFrame*j+n,k))/SamplePerFrame;
				}
			}
			for (int n=0;n<SamplePerFrame/2; ++n) {
				resize_force_array[i](0,k)            = resize_force_array[i](0,k)            + (force_array[i](n,k))/(SamplePerFrame/2);
				resize_moment_array[i](0,k)           = resize_moment_array[i](0,k)           + (moment_array[i](n,k))/(SamplePerFrame/2);
				resize_origin_array[i](0,k)         = resize_origin_array[i](0,k)         + (origin_array[i](n,k))/(SamplePerFrame/2);
				resize_force_array[i](Nframes-1,k)    = resize_force_array[i](Nframes-1,k)    + (force_array[i](Nframes-1-n,k))/(SamplePerFrame/2);
				resize_moment_array[i](Nframes-1,k)   = resize_moment_array[i](Nframes-1,k)   + (moment_array[i](Nframes-1-n,k))/(SamplePerFrame/2);
				resize_origin_array[i](Nframes-1,k) = resize_origin_array[i](Nframes-1,k) + (origin_array[i](Nframes-1-n,k))/(SamplePerFrame/2);
			}
		}		
		for (int m=0; m<Nframes; ++m) {
			Vector3d FF = resize_force_array[i].row(m);
			Vector3d MM = resize_moment_array[i].row(m);
			Vector3d OO = resize_origin_array[i].row(m);
			position_array[i].row(m)=FF.cross(MM)/FF.squaredNorm() + OO;
			position_array[i](m,2)=0;
		}
	}
	  
	for (int i = 0 ; i < numberOfForcePlates ; ++i) {
		for (int k=0;k<3; ++k) {		
			for (int j=0;j<FootPrintRightInterval; ++j) {
					Rstep_PlateValues[i](j,k) = resize_force_array[i](j+InitialRightStep,k);
			}
			for (int j=0;j<FootPrintLeftInterval; ++j) {
					Lstep_PlateValues[i](j,k) = resize_force_array[i](j+InitialLeftStep,k);
			}
		}
		Rstep_Norm[i] = Rstep_PlateValues[i].squaredNorm();
		Lstep_Norm[i] = Lstep_PlateValues[i].squaredNorm();
		if (Rstep_Norm[i] > Rstep_Norm[R_plate_id]) {
			R_plate_id=i;
		}
		if (Lstep_Norm[i] > Lstep_Norm[L_plate_id]) {
			L_plate_id=i;
		}
	}	
	if (R_plate_id==L_plate_id) {
		cerr << "ERROR: Force platform wrench identified as left and right at the same time. Probably more than one plate was pressed by the same foot" << endl;
	}

	Vector3d ApplicationPointRight = Vector3d::Zero();
	Vector3d ApplicationPointLeft  = Vector3d::Zero();
	
	for(int i=0; i<Nframes; i++) {  	
		Vector3d FFR = resize_force_array[R_plate_id].row(i);
		Vector3d FFL = resize_force_array[L_plate_id].row(i);
		Vector3d PPR = position_array[R_plate_id].row(i);
		Vector3d PPL = position_array[L_plate_id].row(i);
		Vector3d MMR = PPR.cross(FFR)/1000;
		Vector3d MML = PPL.cross(FFL)/1000;
		for (int j=0;j<3;j++) {
			RInverseForceMomentum[i](j)= MMR(j);
			LInverseForceMomentum[i](j)= MML(j);
			RInverseForceMomentum[i](j+3)=FFR(j);
			LInverseForceMomentum[i](j+3)=FFL(j);
		}
	}
}

VectorNd ID::SetupApplicationPoint (VectorNd ForceToApply, VectorNd InstantPosition, Vector3d ApplicationPoint, unsigned int BodyPartId) {
	SpatialVector ForcesMomentums; 
	Vector3d Basepoint = CalcBodyToBaseCoordinates(*model,InstantPosition,BodyPartId,ApplicationPoint,1); 
	ForcesMomentums[0]=Basepoint[1]*ForceToApply[2]-Basepoint[2]*ForceToApply[1];            
	ForcesMomentums[1]=Basepoint[2]*ForceToApply[0]-Basepoint[0]*ForceToApply[2];               
	ForcesMomentums[2]=Basepoint[0]*ForceToApply[1]-Basepoint[1]*ForceToApply[0];  
	ForcesMomentums[3]=ForceToApply[0];
	ForcesMomentums[4]=ForceToApply[1];
	ForcesMomentums[5]=ForceToApply[2];
	VectorNd OutputForcesMomentums = ForcesMomentums; 
	return OutputForcesMomentums;
}

vector<Math::SpatialVector> ID::SetupExternalForces(VectorNd ForcesMomentumsLeft,VectorNd ForcesMomentumsRight) {   
	vector<SpatialVector> fext;
	for(int i=0; i<model->mBodies.size(); i++) {    
		fext.push_back(SpatialVector::Zero(6));
	}
	for (int i=0;i<6;i++) {
		fext[Id_Left][i]=ForcesMomentumsLeft[i];                   
		fext[Id_Right][i]=ForcesMomentumsRight[i];
	}     
	return fext;  
}

int ID::ReadStepInterval(btk::Acquisition::Pointer Acq) { 
	InitialRightStep = -1;
	InitialLeftStep  = -1;
	EndRightStep = -1;
	EndLeftStep  = -1;

	for(int i=0; i<Acq->GetEventNumber(); i++) {
		const std::string & label = (Acq->GetEvent(i))->GetLabel();
		const std::string & side = (Acq->GetEvent(i))->GetContext();
		int Frame = (Acq->GetEvent(i))->GetFrame();
		if(label=="Foot Strike") {
			if(side=="Right") {
				if((InitialRightStep == -1) or (Frame < InitialRightStep)) {
					InitialRightStep = Frame;
				}
			}
			else if(side=="Left") {
				if((InitialLeftStep  == -1) or (Frame < InitialLeftStep)) {
					InitialLeftStep = Frame;
				}
			}
		}
	}
	for(int i=0; i<Acq->GetEventNumber(); i++) {
		const std::string & label = (Acq->GetEvent(i))->GetLabel();
		const std::string & side = (Acq->GetEvent(i))->GetContext();
		int Frame = (Acq->GetEvent(i))->GetFrame();
		if(label=="Foot Off") {
			if(side=="Right") {
				if (((EndRightStep ==-1) or (Frame < EndRightStep)) and (Frame > InitialRightStep)) {
					EndRightStep = Frame;
				}
			}
			else if(side=="Left") {
				if (((EndLeftStep  ==-1) or (Frame < EndLeftStep)) and (Frame > InitialLeftStep)) {
					EndLeftStep = Frame;
				}
			}
		}
	}	
	int FirstFrame         = Acq->GetFirstFrame();
	EndRightStep           = EndRightStep     - FirstFrame;
	EndLeftStep            = EndLeftStep      - FirstFrame;
	InitialRightStep       = InitialRightStep - FirstFrame;
	InitialLeftStep        = InitialLeftStep  - FirstFrame;
	FootPrintRightInterval = EndRightStep     - InitialRightStep;
	FootPrintLeftInterval  = EndLeftStep      - InitialLeftStep;
}

void ID::GetDerivativeswithFiniteDifferences(){
       	for(int i=0; i<NumQ; i++) {		
		for (int j=0; j<NumT-1; j++) {
			(QDot[j])[i]=((Q[j+1])[i]-(Q[j])[i])/(delta_t);
		}
		QDot[NumT-1][i]=QDot[NumT-2][i];
	}
  	for(int i=0; i<NumQ; i++) {		
		for (int j=0; j<NumT-1; j++) {
			(QDDot[j])[i]=((QDot[j+1])[i]-(QDot[j])[i])/(delta_t);
		}
		QDDot[NumT-1][i]=QDDot[NumT-2][i];
	}
}

void ID::PrintObjectofTypevectorVectorNDtoFile(vector<VectorNd>* PointerToValues, const char * filename) {
	ofstream output_file (filename, ios::out);
	vector<VectorNd> Values=*PointerToValues;
	if (!output_file) {
		cerr << "Error: could not open file " << filename << "." << endl;abort();
	}
	for (int i = 0; i < Values.size(); i++) {    
		output_file << i << ", ";
		for(int j = 0; j<Values[i].size(); j++) {
			output_file << (Values[i])[j];
			if (j != Values[i].size() - 1) {
				output_file << ", ";
			}
		}
		output_file << endl;
	}
	output_file.close();
}

/** **************************************************************** **/
/** ** The following functions are derived from the BTK tutorials ** **/
/** **************************************************************** **/

btk::Acquisition::Pointer ID::readAcquisition(const std::string& filename) {
  	btk::AcquisitionFileReader::Pointer reader = btk::AcquisitionFileReader::New();
	reader->SetFilename(filename);
	reader->Update();
	return reader->GetOutput();
}

btk::WrenchCollection::Pointer ID::readbtkForces(const std::string& filename) {
	btk::Acquisition::Pointer reader = readAcquisition(filename);
	btk::ForcePlatformsExtractor::Pointer fpExtractor = btk::ForcePlatformsExtractor::New();
	fpExtractor->SetInput(reader);
	btk::ForcePlatformWrenchFilter::Pointer fpwFilter = btk::ForcePlatformWrenchFilter::New();
	fpwFilter->SetInput(fpExtractor->GetOutput());
	btk::WrenchCollection::Pointer fpwrs = fpwFilter->GetOutput();	
	fpwFilter->SetTransformToGlobalFrame(true);
	fpwrs->Update();
	return fpwrs;
}

std::string ID::GetPointTypeAsString(btk::Point::Type t) {
	if (t == btk::Point::Marker) return "Marker";
	else if (t == btk::Point::Angle) return "Angle";
	else if (t == btk::Point::Force) return "Force";
	else if (t == btk::Point::Moment) return "Moment";
	else if (t == btk::Point::Power) return "Power";
	else if (t == btk::Point::Scalar) return "Scalar";
}

void ID::PrintLabelsFromc3dFile(const string& filename) {  
	btk::Acquisition::Pointer acq = readAcquisition(filename);
	for (btk::Acquisition::PointConstIterator it = acq->BeginPoint() ; it != acq->EndPoint() ; ++it) {
		std::cout << (*it)->GetLabel() << " (" << (*it)->GetDescription() << "): " << GetPointTypeAsString((*it)->GetType())<<endl;
	}
}