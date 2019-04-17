//
//  4Q.h
//  stap++
//
//  Created by six on 2017/11/9.
//

#ifndef _Q_h
#define _Q_h

#include "Element.h"

using namespace std;

//! Bar element class
class C4Q : public CElement
{
private:
    double BmatTotal[32];
    double JacobiDet[4];
    double Cttmat[3];
public:

//!	Constructor
	C4Q();

//!	Desconstructor
	~C4Q();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream OutputFile
	virtual void Write(COutputter& output, unsigned int Ele);

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element gravity force
	virtual void ElementGravity(double* bodyforce, double Gravity);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
    
    void CalCttmat();

    void CalBmatTotal();

    void CalBmatSingle(double* Bmat,double* JacobiDet,double*coor,double yita,double psi);

    void ElementStiffnessSingle(double* StiffMatrix,double* Bmat,double* Cttmat,double JacobiDet, double weight);
};



#endif /* _Q_h */
