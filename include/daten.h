// Class data
/* ******************************************************
 *
 * NOTE: General side constraints are not considered(add if necessary)
 *
 * ******************************************************/
#ifndef __DATEN_H__
#define __DATEN_H__
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

typedef struct {
	int x;
	int y;
	int z;
} coordinate;

class Block{
public:
	Block();
	~Block();
	//enquiry
	int getID() const;

	coordinate getXYZ() const;
	int getX() const;
	int getY() const;
	int getZ() const;

        std::vector<int> * getPreds(); // should probably get rid of this version
	const std::vector<int> &getPreds() const {return _pred;}
	std::vector<int> getPred() {return _pred;}

	int getNumPred() const;
	bool isPred(int b) const;

	std::vector<double> * getAttr();
	const std::vector<double> &getAttr() const {return _attributes;}
	double getAttr(int i) const;
	int getnAttr() const;

	double getProfit(int d)const;
	std::vector<double> * getProfit();
	const std::vector<double> &getProfit() const {return _profit;}

	double getRCoef(int d, int r) const;
	std::vector< std::vector<double> >  * getRCoef();
	const std::vector< std::vector<double> >  &getRCoef() const {return _rCoef;}

	//methods
	void setID(int b);
	void setXYZ(coordinate xyz);
	void setX(int x);
	void setY(int y);
	void setZ(int z);

	void addPred(int b);

	// methods for block groups
	void addMember(int b){_members.push_back(b);}
	std::vector<int> getMembers(){return _members;}

	void setAttr(int i, double v);
	void addAttr(double v);

	void setProfit(int d, double p);
	void setProfit(std::vector<double> p);

	void setRCoef(int d, int r, double c);
	void initialize_rCoef( int nDestinations, int nResources);

private:
	int _id;									// block id
	coordinate _xyz;							// x, y, z coordinates
	std::vector<int> _pred;						// list of immediate predecessors

	std::vector<int> _members; 		// list of members (for block groups)

	std::vector<double> _attributes;			// atribute => value

	/* Obj coefficients*/
	std::vector<double> _profit;  				// destination => profit

	/* resource requirements */
	std::vector< std::vector<double> >  _rCoef;  // destination -> resource => coefficient

};
//---------------------------------------------------------------------------------------
class Daten{
public:
	Daten(const char * filename, char prob = 'c'); 		// prob = 'u'(UPIT)), 'c'(CPIT), 'p'(PCPSP)
	Daten(){}
	~Daten();
	//enquiry
	std::string getName() const {return _name;}
	int getNBlock() const       {return _nBlock;}
	char getProbType() const    {return _prob_type;}
	int getNPeriod() const      {return _t_max;}
	int getnDestination() const {return _d_max;}
	int getnResources() const   {return _r_max;}
  	int getnArcs() const 	    {return _nArcs;}
	double getDiscountRate() const {return _rate;}
	std::vector< std::vector< std::vector<int> > > * getLimit();
	const std::vector< std::vector< std::vector<int> > > &getLimit() const {return _rLimits;}
	int getLimit(int r, int t, int i=0) const ;
	std::vector< std::vector<char> >  * getResConstrType();
	const std::vector< std::vector<char> >  &getResConstrType() const {return _rConstrType;}
	char getResConstrType(int r, int t) const;
	std::vector<Block> * getBlock();
	const std::vector<Block> &getBlock() const {return _blocks;} // just get reference, don't copy
    const Block &getBlock(int i) const {return _blocks[i];}
	std::vector< std::vector< std::vector<int> > > getResLimits(){return _rLimits;}

protected:
//members
	std::string _name;
	char _prob_type;
	int _nBlock;
	int _t_max;
	int _d_max; // # destinations (default = 1)
	int _r_max;
	double _rate;
	int _nArcs;
	std::vector< std::vector< std::vector<int> > >	_rLimits;  	// resource limits : resource -> period => (value 1; value 2)  value 2 = 0   if Type is L or G
	std::vector< std::vector<char> >  	_rConstrType; 			// resource -> period => constraint type( 'L' , 'G' or 'I')
	std::vector<Block> _blocks;
//methods
	void initialize_rLimits();
	void readData(const char * filename);
	void readProb(const char * filename);
	void readBlocks(const char * filename);
	void readPrec(const char * filename);

  void usage();
};

#endif
