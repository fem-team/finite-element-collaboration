//Edited by Yunzhe

#pragma once

#include "Element.h"

using namespace std;

//! 4Q element class
class CQ4 : public CElement
{
public:

//Constructor
	CQ4();
//Deconstructor
	~CQ4();

//Read Data
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//Write Data
	virtual void Write(COutputter& output, unsigned int Ele);

//Generation LM
	virtual void GenerateLocationMatrix();

//Stiffness Matrix
	virtual void ElementStiffness(double* Matrix);

//Stress at Gauss Point
	virtual void ElementStress(double* stress, double* Displacement);

//
	virtual unsigned int SizeOfStiffnessMatrix();

	virtual double Gravity();

//
	virtual void ElementStrainMatrix(double* Be, double eta, double psi);

// 
	virtual void GaussPointLocalAssembly(double* CM, double* Matrix, double x_gauss, double y_gauss);
//
	virtual void ElementCoordinates (double * coord);

private:

	double GP;

	
};