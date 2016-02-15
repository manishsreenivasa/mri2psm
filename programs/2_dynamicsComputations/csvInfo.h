/* MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
*/
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <btkConfigure.h>
#include <rbdl/rbdl.h>

using namespace std;
using namespace RigidBodyDynamics::Math;

class csvInfo {                  
	vector<VectorNd> Q;
	const char* FileName;
	int DegreesOfFreedom;
	double StartTime;
	double EndTime;
	int TimeSteps;
	int LineWith_DATA;
	bool LineWith_DATAFOUND;
	string Line;	
 
	public:
		csvInfo(const char* NameOfFile);
		void loadCSV();
		int GetTimeSteps();
		int GetDegreesOfFreedom();
		double GetStartTime();
		double GetEndTime();
		vector<VectorNd> GetQ();   

	private:	
		void Initiate();        
		template<class string> void tokenizeV(const std::string &s,std::vector<string> &o);
};
csvInfo::csvInfo(const char* NameOfFile) {FileName=NameOfFile;}
void csvInfo::loadCSV() { 	
	Initiate();
	ifstream Animationfile(FileName);
	if(LineWith_DATAFOUND) {
		for (int i=0;i<=LineWith_DATA;i++)
			getline(Animationfile, Line);
	}
	getline(Animationfile, Line);
	std::vector<double> VectorFirstRow;
	tokenizeV(Line, VectorFirstRow);
	StartTime = VectorFirstRow[0];
	DegreesOfFreedom = VectorFirstRow.size()-1;
	Q.push_back(VectorNd::Zero(DegreesOfFreedom));
	for (int i=0;i<DegreesOfFreedom;i++)
		Q[0][i] = VectorFirstRow[i+1];
	TimeSteps = 1; 

	while(getline(Animationfile, Line)) {   
		std::vector<double> Row;
		tokenizeV(Line, Row);
		EndTime = Row[0];
		Q.push_back(VectorNd::Zero(DegreesOfFreedom));
		for(int i=0;i<DegreesOfFreedom;i++) 
			Q[TimeSteps][i] = Row[i+1];		
		TimeSteps++;
	}
}

void csvInfo::Initiate() {
	ifstream Animationfile(FileName);
	LineWith_DATA = 0;
	LineWith_DATAFOUND = false;
	while(getline(Animationfile, Line)) {  
		if(!strstr(&Line[0],"DATA")) 
			LineWith_DATA++;
		else {
			LineWith_DATAFOUND = true;
			break;
		}
	}
}

template<class string> void csvInfo::tokenizeV(const std::string &s,std::vector<string> &o) { 
	typedef boost::tokenizer<boost::escaped_list_separator<char> >  tok_t;
	tok_t tok(s);
	for(tok_t::iterator j (tok.begin());j != tok.end();++j) {
		std::string f(*j);
		boost::trim(f);
		o.push_back(boost::lexical_cast<string>(f));
	}
}
	
int csvInfo::GetTimeSteps() {return TimeSteps;}
int csvInfo::GetDegreesOfFreedom() {return DegreesOfFreedom;}
double csvInfo::GetStartTime() {return StartTime;}
double csvInfo::GetEndTime() {return EndTime;}
vector<VectorNd> csvInfo::GetQ() {return Q;}
