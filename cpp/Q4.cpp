// Edited by Yunzhe


#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;


// Constructor
CQ4::CQ4()
{
	NEN_ = 4; //  4Q has 4 nodes.
	nodes_ = new CNode*[NEN_];

	ND_ = 8;  //  4 X 2 = 8 degree of freedom  
	LocationMatrix_ = new unsigned int [ND_];

	ElementMaterial_ = nullptr;

}

//Deconstructor
CQ4::~CQ4()
{
}


bool CQ4::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N ; // Element Number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;  // Material Set
	unsigned int N1, N2, N3, N4; // Node Number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[0]->bcode[2]=1;
	nodes_[1]->bcode[2]=1;
	nodes_[2]->bcode[2]=1;
	nodes_[3]->bcode[2]=1;

	return true;
}

void CQ4::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
	<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber
		<< setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

// Waiting for Certification
void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];


		
}

unsigned int CQ4::SizeOfStiffnessMatrix() { return 36; }
	
void CQ4::ElementStrainMatrix(double* Be,double eta, double psi)
	
{
	//Get Inversion of Jacobi
	double GN[4]; //
	double  J[4]; //J=[1 2;3 4]
	double  detJ;
	double  InvJ[4];
	GN[0] = 0.25*(eta-1);
	GN[1] = 0.25*(eta+1);
	GN[2] = 0.25*(psi-1);
	GN[3] = 0.25*(psi+1);

	J[0] = GN[0]*nodes_[0]->XYZ[0] - GN[0]*nodes_[1]->XYZ[0] + GN[1]*nodes_[2]->XYZ[0] - GN[1]*nodes_[3]->XYZ[0];
	J[2] = GN[0]*nodes_[0]->XYZ[1] - GN[0]*nodes_[1]->XYZ[1] + GN[1]*nodes_[2]->XYZ[1] - GN[1]*nodes_[3]->XYZ[1];

	J[1] = GN[2]*nodes_[0]->XYZ[0] - GN[3]*nodes_[1]->XYZ[0] + GN[3]*nodes_[2]->XYZ[0] - GN[2]*nodes_[3]->XYZ[0];
	J[3] = GN[2]*nodes_[0]->XYZ[1] - GN[3]*nodes_[1]->XYZ[1] + GN[3]*nodes_[2]->XYZ[1] - GN[2]*nodes_[3]->XYZ[1];

	detJ = J[0]*J[3]-J[1]*J[2];

	InvJ[0] = J[3]/detJ;
	InvJ[1] = -J[1]/detJ;
	InvJ[2] = -J[2]/detJ;
	InvJ[3] = J[0]/detJ;

	//Get StrainMatrix Be = [N1/x N1/y N2/x N2/y N3/x N3/y N4/x N4/y]
	Be[0] = detJ;
	Be[1] = InvJ[0]*GN[0] + InvJ[2]*GN[2];
	Be[2] = InvJ[1]*GN[0] + InvJ[3]*GN[2];
	Be[3] = -InvJ[0]*GN[0] - InvJ[2]*GN[3];
	Be[4] = -InvJ[1]*GN[0] - InvJ[3]*GN[3];
	Be[5] = InvJ[0]*GN[1] + InvJ[2]*GN[3];
	Be[6] = InvJ[1]*GN[1] + InvJ[3]*GN[3];
	Be[7] = -InvJ[0]*GN[1] - InvJ[2]*GN[2];
	Be[8] = -InvJ[1]*GN[1] - InvJ[3]*GN[2];
	

	
}

void CQ4::GaussPointLocalAssembly(double* CM, double* Matrix, double x_gauss, double y_gauss)
{
	

	double Be[9];
	ElementStrainMatrix(Be,x_gauss, y_gauss);
	Matrix[0] = Matrix[0] + ( Be[1]*CM[0]*Be[1] + Be[2]*CM[2]*Be[2] )*Be[0]; //b0*c0*conj(b0) + b1*c2*conj(b1)
	Matrix[1] = Matrix[1] + ( Be[1]*CM[2]*Be[1] + Be[2]*CM[0]*Be[2] )*Be[0]; //b0*c2*conj(b0) + b1*c0*conj(b1)
	Matrix[2] = Matrix[2] + ( Be[2]*CM[1]*Be[1] + Be[1]*CM[2]*Be[2] )*Be[0]; //b1*c1*conj(b0) + b0*c2*conj(b1)
	Matrix[3] = Matrix[3] + ( Be[3]*CM[0]*Be[3] + Be[4]*CM[2]*Be[4] )*Be[0]; //b2*c0*conj(b2) + b3*c2*conj(b3)
	Matrix[4] = Matrix[4] + ( Be[3]*CM[1]*Be[2] + Be[4]*CM[2]*Be[1] )*Be[0]; //b2*c1*conj(b1) + b3*c2*conj(b0)
	Matrix[5] = Matrix[5] + ( Be[3]*CM[0]*Be[1] + Be[4]*CM[2]*Be[2] )*Be[0]; //b2*c0*conj(b0) + b3*c2*conj(b1)
	Matrix[6] = Matrix[6] + ( Be[3]*CM[2]*Be[3] + Be[4]*CM[0]*Be[4] )*Be[0]; //b2*c2*conj(b2) + b3*c0*conj(b3)
	Matrix[7] = Matrix[7] + ( Be[4]*CM[1]*Be[3] + Be[3]*CM[2]*Be[4] )*Be[0]; //b3*c1*conj(b2) + b2*c2*conj(b3)
	Matrix[8] = Matrix[8] + ( Be[3]*CM[2]*Be[1] + Be[4]*CM[0]*Be[2] )*Be[0]; //b2*c2*conj(b0) + b3*c0*conj(b1)
	Matrix[9] = Matrix[9] + ( Be[4]*CM[1]*Be[1] + Be[3]*CM[2]*Be[2] )*Be[0]; //b3*c1*conj(b0) + b2*c2*conj(b1)
	Matrix[10] = Matrix[10] + ( Be[5]*CM[0]*Be[5] + Be[6]*CM[2]*Be[6] )*Be[0]; //b4*c0*conj(b4) + b5*c2*conj(b5)
	Matrix[11] = Matrix[11] + ( Be[5]*CM[1]*Be[4] + Be[6]*CM[2]*Be[3] )*Be[0]; //b4*c1*conj(b3) + b5*c2*conj(b2)
	Matrix[12] = Matrix[12] + ( Be[5]*CM[0]*Be[3] + Be[6]*CM[2]*Be[4] )*Be[0]; //b4*c0*conj(b2) + b5*c2*conj(b3)
	Matrix[13] = Matrix[13] + ( Be[5]*CM[1]*Be[2] + Be[6]*CM[2]*Be[1] )*Be[0]; //b4*c1*conj(b1) + b5*c2*conj(b0)
	Matrix[14] = Matrix[14] + ( Be[5]*CM[0]*Be[1] + Be[6]*CM[2]*Be[2] )*Be[0]; //b4*c0*conj(b0) + b5*c2*conj(b1)
	Matrix[15] = Matrix[15] + ( Be[5]*CM[2]*Be[5] + Be[6]*CM[0]*Be[6] )*Be[0]; //b4*c2*conj(b4) + b5*c0*conj(b5)
	Matrix[16] = Matrix[16] + ( Be[6]*CM[1]*Be[5] + Be[5]*CM[2]*Be[6] )*Be[0]; //b5*c1*conj(b4) + b4*c2*conj(b5)
	Matrix[17] = Matrix[17] + ( Be[5]*CM[2]*Be[3] + Be[6]*CM[0]*Be[4] )*Be[0]; //b4*c2*conj(b2) + b5*c0*conj(b3)
	Matrix[18] = Matrix[18] + ( Be[6]*CM[1]*Be[3] + Be[5]*CM[2]*Be[4] )*Be[0]; //b5*c1*conj(b2) + b4*c2*conj(b3)
	Matrix[19] = Matrix[19] + ( Be[5]*CM[2]*Be[1] + Be[6]*CM[0]*Be[2] )*Be[0]; //b4*c2*conj(b0) + b5*c0*conj(b1)
	Matrix[20] = Matrix[20] + ( Be[6]*CM[1]*Be[1] + Be[5]*CM[2]*Be[2] )*Be[0]; //b5*c1*conj(b0) + b4*c2*conj(b1)
	Matrix[21] = Matrix[21] + ( Be[7]*CM[0]*Be[7] + Be[8]*CM[2]*Be[8] )*Be[0]; //b6*c0*conj(b6) + b7*c2*conj(b7)
	Matrix[22] = Matrix[22] + ( Be[7]*CM[1]*Be[6] + Be[8]*CM[2]*Be[5] )*Be[0]; //b6*c1*conj(b5) + b7*c2*conj(b4)
	Matrix[23] = Matrix[23] + ( Be[7]*CM[0]*Be[5] + Be[8]*CM[2]*Be[6] )*Be[0]; //b6*c0*conj(b4) + b7*c2*conj(b5)
	Matrix[24] = Matrix[24] + ( Be[7]*CM[1]*Be[4] + Be[8]*CM[2]*Be[3] )*Be[0]; //b6*c1*conj(b3) + b7*c2*conj(b2)
	Matrix[25] = Matrix[25] + ( Be[7]*CM[0]*Be[3] + Be[8]*CM[2]*Be[4] )*Be[0]; //b6*c0*conj(b2) + b7*c2*conj(b3)
	Matrix[26] = Matrix[26] + ( Be[7]*CM[1]*Be[2] + Be[8]*CM[2]*Be[1] )*Be[0]; //b6*c1*conj(b1) + b7*c2*conj(b0)
	Matrix[27] = Matrix[27] + ( Be[7]*CM[0]*Be[1] + Be[8]*CM[2]*Be[2] )*Be[0]; //b6*c0*conj(b0) + b7*c2*conj(b1)
	Matrix[28] = Matrix[28] + ( Be[7]*CM[2]*Be[7] + Be[8]*CM[0]*Be[8] )*Be[0]; //b6*c2*conj(b6) + b7*c0*conj(b7)
	Matrix[29] = Matrix[29] + ( Be[8]*CM[1]*Be[7] + Be[7]*CM[2]*Be[8] )*Be[0]; //b7*c1*conj(b6) + b6*c2*conj(b7)
	Matrix[30] = Matrix[30] + ( Be[7]*CM[2]*Be[5] + Be[8]*CM[0]*Be[6] )*Be[0]; //b6*c2*conj(b4) + b7*c0*conj(b5)
	Matrix[31] = Matrix[31] + ( Be[8]*CM[1]*Be[5] + Be[7]*CM[2]*Be[6] )*Be[0]; //b7*c1*conj(b4) + b6*c2*conj(b5)
	Matrix[32] = Matrix[32] + ( Be[7]*CM[2]*Be[3] + Be[8]*CM[0]*Be[4] )*Be[0]; //b6*c2*conj(b2) + b7*c0*conj(b3)
	Matrix[33] = Matrix[33] + ( Be[8]*CM[1]*Be[3] + Be[7]*CM[2]*Be[4] )*Be[0]; //b7*c1*conj(b2) + b6*c2*conj(b3)
	Matrix[34] = Matrix[34] + ( Be[7]*CM[2]*Be[1] + Be[8]*CM[0]*Be[2] )*Be[0]; //b6*c2*conj(b0) + b7*c0*conj(b1)
	Matrix[35] = Matrix[35] + ( Be[8]*CM[1]*Be[1] + Be[7]*CM[2]*Be[2] )*Be[0]; //b7*c1*conj(b0) + b6*c2*conj(b1)


	

}



void CQ4::ElementStiffness(double* Matrix)
{
	GP = 0.57735;
	clear(Matrix, SizeOfStiffnessMatrix());
	double CM[3];
	//Constitute Matrix Element
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	CM[0] = material_->E / (1- material_->Nu * material_->Nu);   //1
	CM[1] = CM[0] * material_->Nu; // Nu
	CM[2] = (CM[0] - CM[1])/2;  // 1-Nu /2

	GaussPointLocalAssembly(CM,Matrix,-GP,-GP);
	GaussPointLocalAssembly(CM,Matrix,GP,-GP);
	GaussPointLocalAssembly(CM,Matrix,-GP,GP);
	GaussPointLocalAssembly(CM,Matrix,GP,GP);




}

void CQ4::ElementStress(double* stress, double* Displacement)
{
	double Be[9];
	double CM[3];
	GP = 0.57735;

	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);
	CM[0] = material_->E / (1- material_->Nu * material_->Nu);
	CM[1] = CM[0] * material_->Nu;
	CM[2] = (CM[0] - CM[1])/2;
	

	double strain_temp[3];
	memset(strain_temp,0,sizeof(strain_temp));
	ElementStrainMatrix(Be,-GP,-GP);

	for(int i = 1;i < 5;i++)
	{
		strain_temp[0] += Be[2*i-1]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0);
		strain_temp[1] += Be[2*i]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
		strain_temp[2] += Be[2*i]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0) + Be[2*i-1]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
	}


	stress[0] = CM[0]*strain_temp[0] + CM[1]*strain_temp[1];
	stress[1] = CM[1]*strain_temp[0] + CM[0]*strain_temp[1];
	stress[2] = CM[2]*strain_temp[2];

	memset(strain_temp,0,sizeof(strain_temp));
	ElementStrainMatrix(Be,GP,-GP);
	for(int i = 1;i < 5;i++)
	{
		strain_temp[0] += Be[2*i-1]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0);
		strain_temp[1] += Be[2*i]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
		strain_temp[2] += Be[2*i]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0) + Be[2*i-1]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
	}

	stress[3] = CM[0]*strain_temp[0] + CM[1]*strain_temp[1];
	stress[4] = CM[1]*strain_temp[0] + CM[0]*strain_temp[1];
	stress[5] = CM[2]*strain_temp[2];

	memset(strain_temp,0,sizeof(strain_temp));
	ElementStrainMatrix(Be,-GP,GP);
	for(int i = 1;i < 5;i++)
	{
		strain_temp[0] += Be[2*i-1]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0);
		strain_temp[1] += Be[2*i]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
		strain_temp[2] += Be[2*i]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0) + Be[2*i-1]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
	}

	stress[6] = CM[0]*strain_temp[0] + CM[1]*strain_temp[1];
	stress[7] = CM[1]*strain_temp[0] + CM[0]*strain_temp[1];
	stress[8] = CM[2]*strain_temp[2];

	memset(strain_temp,0,sizeof(strain_temp));
	ElementStrainMatrix(Be,GP,GP);
	for(int i = 1;i < 5;i++)
	{
		strain_temp[0] += Be[2*i-1]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0);
		strain_temp[1] += Be[2*i]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
		strain_temp[2] += Be[2*i]*Displacement[LocationMatrix_[2*i-2]-1]*(double)(!LocationMatrix_[2*i-2]==0) + Be[2*i-1]*Displacement[LocationMatrix_[2*i-1]-1]*(double)(!LocationMatrix_[2*i-1]==0);
	}

	stress[9] = CM[0]*strain_temp[0] + CM[1]*strain_temp[1];
	stress[10] = CM[1]*strain_temp[0] + CM[0]*strain_temp[1];
	stress[11] = CM[2]*strain_temp[2];
	

}

void CQ4::ElementCoordinates (double * coord)
{
	double eta;
	double psi;
	GP = 0.57735;
	eta = -GP; psi = -GP;
	coord[0]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[1]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = GP; psi = -GP;
	coord[2]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[3]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = -GP; psi = GP;
	coord[4]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[5]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;

	eta = GP; psi = GP;
	coord[6]=( nodes_[0]->XYZ[0]*(1-eta)*(1-psi) + nodes_[1]->XYZ[0]*(1-eta)*(1+psi) + nodes_[2]->XYZ[0]*(1+eta)*(1+psi) + nodes_[3]->XYZ[0]*(1+eta)*(1-psi) )/4;
	coord[7]=( nodes_[0]->XYZ[1]*(1-eta)*(1-psi) + nodes_[1]->XYZ[1]*(1-eta)*(1+psi) + nodes_[2]->XYZ[1]*(1+eta)*(1+psi) + nodes_[3]->XYZ[1]*(1+eta)*(1-psi) )/4;
}

	double CQ4::Gravity()
	{return 0;}