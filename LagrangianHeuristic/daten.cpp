using namespace std;

#include "daten.h"

Block::Block(){
	_id = -1;
	_xyz.x=0;
	_xyz.y=0;
	_xyz.z=0;
}
Block::~Block(){}
//enquiry
int Block::getID(){ return _id;}
coordinate Block::getXYZ(){ return _xyz;}

int Block::getX(){ return _xyz.x;}
int Block::getY(){ return _xyz.y;}
int Block::getZ(){ return _xyz.z;}

std::vector<int> * Block::getPreds(){ return &_pred;}
int Block::getNumPred(){ return _pred.size();}
bool Block::isPred(int b){
	int m = _pred.size();
	bool YesNo = 0;
	for(int i =0; i<m; i++){
		if(_pred[i] == b){
			YesNo = 1;
			break;
		} 
	}
	return YesNo;
}

std::vector<double> * Block::getAttr(){return &_attributes;}
double Block::getAttr(int i){return _attributes[i];}
int Block::getnAttr(){return _attributes.size();}	

double Block::getProfit(int d){return _profit[d];}
std::vector<double> * Block::getProfit(){return &_profit;}

double Block::getRCoef(int d, int r){return _rCoef[d][r];}
std::vector< std::vector<double> > * Block::getRCoef(){return &_rCoef;}

//methods
void Block::setID(int b){_id = b;}
void Block::setXYZ(coordinate xyz){_xyz.x= xyz.x; _xyz.y =xyz.y; _xyz.z = xyz.z;}
void Block::setX(int x){_xyz.x = x;}
void Block::setY(int y){_xyz.y = y;}
void Block::setZ(int z){_xyz.z = z;}

void Block::addPred(int b) {_pred.push_back(b);}

void Block::setAttr(int i, double v){_attributes[i]=v;}
void Block::addAttr(double v){_attributes.push_back(v);}
	
void Block::setProfit(int d, double p){_profit[d]=p;}
void Block::setProfit(std::vector<double> p){_profit=p;}

void Block::setRCoef(int d, int r, double c){ _rCoef[d][r]=c; }
void Block::initialize_rCoef( int nDestinations, int nResource){
	_rCoef = std::vector< std::vector<double> > (nDestinations, std::vector<double>(nResource, 0.0));
}

Daten::Daten(const char * filename, char prob = 'p'){
	// Initialize
	_prob_type = prob; 									// UPIT = 'u',  CPIT = 'c', PCPSP = 'p';
	_nBlock = 0;										
	_t_max = 0;												// # periods
	_d_max = 1;												// number of destinations
	_r_max = 0;										// # resources = number of destinations
	_rate = 1.0;
	readData(filename);
}
Daten::~Daten(){}
//enquiry
int Daten::getNBlock(){return _nBlock;}
char Daten::getProbType(){return _prob_type;}
int Daten::getNPeriod(){return _t_max;}
int Daten::getnDestination(){return _d_max;}
int Daten::getnResources(){return _r_max;}
double Daten::getDiscountRate(){return _rate;}
std::vector< std::vector< std::vector<int> > > * Daten::getLimit(){return &_rLimits;}
int Daten::getLimit(int r, int t, int i){return _rLimits[r][t][i];}

std::vector< std::vector<char> > * Daten::getResConstrType(){return &_rConstrType;}
char Daten::getResConstrType(int r, int t){return _rConstrType[r][t];}

std::vector<Block> * Daten::getBlock(){return &_blocks;}
Block Daten::getBlock(int i){return _blocks[i];}

void Daten::readData(const char * filename){
	
	// First read the problem
	 readProb(filename);
	// Read Block Description
	readBlocks(filename);
	// Read Precedence 
	readPrec(filename);
}
void Daten::initialize_rLimits(){
	_rLimits = std::vector< std::vector< std::vector<int> > > (_r_max, std::vector< std::vector<int> > (_t_max, std::vector<int> (2,0)));
	_rConstrType = std::vector< std::vector<char> >  (_r_max, std::vector<char> (_t_max,'L'));
}
void Daten::readProb(const char * filename){
	std::string datafile(filename);
	std::string line, strng;
	std::string::size_type foundPos, foundEnd;
	foundPos = datafile.find_last_of("/");
	if(foundPos==std::string::npos) strng="";
	else strng = datafile.substr(foundPos);
	datafile += strng;
	if(_prob_type == 'u') datafile += ".upit";
	else if(_prob_type == 'c'){ 
			datafile += ".cpit";
	}else datafile += ".pcpsp";
	std::ifstream inputfile;
	inputfile.open(datafile.c_str());
	if (inputfile.is_open()){
		
		// skip lines: NAME, TYPE, NGENERAL_SIDE_CONSTRAINTS	
		while(!inputfile.eof()){ 
			getline(inputfile,line);
			if(line == "") continue; 
			foundPos = line.find_first_not_of(" \t");
			if(line[foundPos] == '%') continue;  // skip comments
			if(line.find("NBLOCKS", foundPos) != std::string::npos){ 
				foundPos = line.find(":",foundPos);
				foundPos = line.find_first_not_of(" \t,", foundPos+1);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				_nBlock = atoi(strng.c_str()); 
				_blocks = std::vector<Block>(_nBlock);
				for(int i=0; i<_nBlock; i++) _blocks[i].setID(i);
			}else if(line.find("NPERIODS", foundPos) != std::string::npos){ 
				foundPos = line.find(":",foundPos);
				foundPos = line.find_first_not_of(" \t,", foundPos+1);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				_t_max = atoi(strng.c_str()); 
			}else if(line.find("NDESTINATIONS", foundPos) != std::string::npos){ 
				foundPos = line.find(":",foundPos);
				foundPos = line.find_first_not_of(" \t,", foundPos+1);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				_d_max = atoi(strng.c_str()); 
			}else if(line.find("NRESOURCE_SIDE_CONSTRAINTS", foundPos) != std::string::npos){ 
				foundPos = line.find(":",foundPos);
				foundPos = line.find_first_not_of(" \t,", foundPos+1);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				_r_max = atoi(strng.c_str()); 
			}else if(line.find("DISCOUNT_RATE", foundPos) != std::string::npos){ 
				foundPos = line.find(":",foundPos);
				foundPos = line.find_first_not_of(" \t,", foundPos+1);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				_rate = atof(strng.c_str()); 
			}else if(line.find("RESOURCE_CONSTRAINT_LIMITS", foundPos) != std::string::npos){ 
				// initialize the container
				initialize_rLimits();
				for(int i=0; i< _r_max*_t_max; i++){
					getline(inputfile,line); // get next line
					// get the resource index
					foundPos = line.find_first_not_of(" \t");
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos);
					int r = atoi(strng.c_str());
					//get the period index
					foundPos = line.find_first_not_of(" \t", foundEnd);
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos);
					int p = atoi(strng.c_str());
					//get constraint type
					foundPos = line.find_first_not_of(" \t", foundEnd);
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos); 
					bool interval = 0;
					if(strng.compare("R")==0) _rConstrType[r][p]='R';
					else if(strng.compare("I")==0){ 
											_rConstrType[r][p]='I';
											interval = 1;
					}	
					// get the limit
					foundPos = line.find_first_not_of(" \t", foundEnd);
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos); 
					_rLimits[r][p][0]=atoi(strng.c_str()); // value 1
					if(interval){ // get value 2
						foundPos = line.find_first_not_of(" \t", foundEnd);
						foundEnd = line.find_first_of(" \t", foundPos+1);
						strng = line.substr(foundPos, foundEnd-foundPos); 
						_rLimits[r][p][1]=atoi(strng.c_str());
						}
				}
			}else if(line.find("OBJECTIVE_FUNCTION", foundPos) != std::string::npos){ 
				for(int i =0; i< _nBlock; i++) _blocks[i].setProfit(std::vector<double> (_d_max, 0));
				for(int i =0; i< _nBlock; i++){
					getline(inputfile, line);
					foundPos = line.find_first_not_of(" \t");
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos); 
					int b=atoi(strng.c_str()); // block id
					// get profits
					for(int j = 0; j<_d_max; j++){
						foundPos = line.find_first_not_of(" \t", foundEnd);
						foundEnd = line.find_first_of(" \t", foundPos+1);
						strng = line.substr(foundPos, foundEnd-foundPos); 
						_blocks[b].setProfit(j, atof(strng.c_str()));
					}
				}
			}else if(line.find("RESOURCE_CONSTRAINT_COEFFICIENTS:", foundPos) != std::string::npos){ 
				for(int i =0; i< _nBlock; i++) _blocks[i].initialize_rCoef(_d_max, _r_max);
				while(!inputfile.eof()){
					getline(inputfile, line);
					if(line == "") continue; 
					foundPos = line.find_first_not_of(" \t");
					if(foundPos != std::string::npos){
							if(line[foundPos] == '%') continue;  // skip comments
							foundPos = line.find_first_not_of(" \t");
							foundEnd = line.find_first_of(" \t", foundPos+1);
							strng = line.substr(foundPos, foundEnd-foundPos); 
							int b=atoi(strng.c_str());	// must be block id
							int d = 0;
							if(_prob_type=='p'){
								foundPos = line.find_first_not_of(" \t", foundEnd);
								foundEnd = line.find_first_of(" \t", foundPos+1);
								strng = line.substr(foundPos, foundEnd-foundPos); 
								d = atoi(strng.c_str());	// must be destination index
							}
							foundPos = line.find_first_not_of(" \t", foundEnd);
							foundEnd = line.find_first_of(" \t", foundPos+1);
							strng = line.substr(foundPos, foundEnd-foundPos); 
							int r = atoi(strng.c_str());	// must be resource index
							
							foundPos = line.find_first_not_of(" \t", foundEnd);
							foundEnd = line.find_first_of(" \t", foundPos+1);
							strng = line.substr(foundPos, foundEnd-foundPos); 
							_blocks[b].setRCoef(d, r, atof(strng.c_str()));
							
//							std::cout<<"value = "<< atof(strng.c_str())<<"  \t string "<<strng<<std::endl;
						}
					}
				}
		}
		inputfile.close();
	}else{
			std::cerr<<"\nCould NOT find the file "<<datafile<<std::endl<<std::endl;
			usage();
			exit(-1);
	}
	
	
}

void Daten::readPrec(const char * filename){
	std::string datafile(filename);
	std::string line, strng;
	std::string::size_type foundPos, foundEnd;
	foundPos = datafile.find_last_of("/");
	if(foundPos==std::string::npos) strng="";
	else strng = datafile.substr(foundPos);
	datafile += strng;
	datafile += ".prec";	
	std::ifstream inputfile;
	inputfile.open(datafile.c_str());
	if (inputfile.is_open()){
		while(!inputfile.eof()){ // read from the file
			getline(inputfile,line);
			if(line == "") continue; 
			foundPos = line.find_first_not_of(" \t");
			if(foundPos != std::string::npos){
				if(line[foundPos] == '%') continue; 
				// get block id
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int b=atoi(strng.c_str()); // block id
				// get number of predecessors 
				foundPos = line.find_first_not_of(" \t", foundEnd);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int n = atoi(strng.c_str());
				// get predecessors
				for(int i=0; i<n; i++){
					foundPos = line.find_first_not_of(" \t", foundEnd);
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos); 
					int p = atoi(strng.c_str());	
					_blocks[b].addPred(p);
				}
			}
		}
		inputfile.close();
	}else{
			std::cerr<<"\nCould NOT find the file "<<datafile<<std::endl<<std::endl;
			usage();
			exit(-1);
	}
}
void Daten::readBlocks(const char * filename){
	std::string datafile(filename);
	std::string line, strng;
	std::string::size_type foundPos, foundEnd;
	foundPos = datafile.find_last_of("/");
	if(foundPos==std::string::npos) strng="";
	else strng = datafile.substr(foundPos);
	datafile += strng;
	datafile += ".blocks";	
	std::ifstream inputfile;
	inputfile.open(datafile.c_str());
	if (inputfile.is_open()){
		while(!inputfile.eof()){ // read from the file
			getline(inputfile,line);
			if(line == "") continue; 
			foundPos = line.find_first_not_of(" \t");
			if(foundPos != std::string::npos){
				if(line[foundPos] == '%') continue; 
				// get block id
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int b=atoi(strng.c_str()); // block id
				// get coordinates
				foundPos = line.find_first_not_of(" \t", foundEnd);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int x = atoi(strng.c_str());	
				foundPos = line.find_first_not_of(" \t", foundEnd);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int y = atoi(strng.c_str());	
				foundPos = line.find_first_not_of(" \t", foundEnd);
				foundEnd = line.find_first_of(" \t", foundPos+1);
				strng = line.substr(foundPos, foundEnd-foundPos); 
				int z = atoi(strng.c_str());	
				
				_blocks[b].setX(x);
				_blocks[b].setY(y);
				_blocks[b].setZ(z);
				
				// get attributes
				foundPos = line.find_first_not_of(" \t", foundEnd);
				while(foundPos != std::string::npos){
					foundEnd = line.find_first_of(" \t", foundPos+1);
					strng = line.substr(foundPos, foundEnd-foundPos); 
					double a = atof(strng.c_str());	
					_blocks[b].addAttr(a);
					foundPos = line.find_first_not_of(" \t", foundEnd);
				}
			}
		}
		inputfile.close();
	}else{
			std::cerr<<"\nCould NOT find the file "<<datafile<<std::endl<<std::endl;
			usage();
			exit(-1);
	}
}
void Daten::usage(){
	std::cout<<"Usage : ./exe [OPTION] PATH "<<std::endl;
	std::cout<<"       Options"<<std::endl;
	std::cout<<"          -u    UPIT"<<std::endl;
	std::cout<<"          -c    CPIT"<<std::endl;
	std::cout<<"          -p    PCPSP"<<std::endl;
	std::cout<<"Sample "<<std::endl;
	std::cout<<"      ./mine_planning -p Data/zuck_small"<<std::endl;
}

