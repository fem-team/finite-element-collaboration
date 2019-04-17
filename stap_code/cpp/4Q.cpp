#include "4Q.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigen>
//constructor
using namespace std;
C4Q::C4Q()
{
    NEN=4;
    nodes= new CNode*[NEN];
	ND=12;
    LocationMatrix= new unsigned int[ND];
    ElementMaterial = nullptr;
}
//deconstructor
C4Q::~C4Q()
{
    delete [] nodes;
    delete [] LocationMatrix;
}

//	Read element data from stream Input
bool C4Q::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;
	// N to the element number
	Input >> N;
	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	// N turns to the number of node
    for(int i=0;i<NEN;i++){
        Input >> N;
        nodes[i]=&NodeList[N - 1];
    }
    
	//N turns to material property set number
    Input >> N;
    ElementMaterial = &(dynamic_cast<C4QMaterial*>(MaterialSets))[N - 1];
 
	return true;
}
//	Write element data to stream OutputFile
void C4Q::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele+1;
    for(int i=0;i<NEN;i++)
	{cout<< setw(9)<< nodes[i]->NodeNumber;}
    cout << setw(12) << ElementMaterial->nset << endl;

}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !

void C4Q::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
    
}
inline void C4Q::CalCttmat(){
    double E=(ElementMaterial)->E;
    double posi_rate=(dynamic_cast<C4QMaterial*>(ElementMaterial))->posi_rate;
    Cttmat(0,0)=E/(1-posi_rate*posi_rate);
    Cttmat(1,1)=E/(1-posi_rate*posi_rate);
    Cttmat(0,1)=posi_rate*Cttmat(0,0);
    Cttmat(1,0)=posi_rate*Cttmat(0,0);
    Cttmat(2,2)=(1-posi_rate)*Cttmat(0,0)/2;
}

void C4Q::CalBmatSingle(double* Bmat,double* JacobiDet,double* coor,double yita,double psi)
{
	MatrixXf GN(2,4)
    GN(0,0)=0.25*(yita-1);  
	GN(0,1)=-GN(0,0);       
	GN(0,2)=0.25*(yita+1);  
	GN(0,3)=-GN(0,2);
    GN(1,0)=0.25*(psi-1);   
	GN(1,1)=-0.25*(psi+1); 
	GN(1,2)=-GN(1,1);       
	GN(1,3)=-GN(1,0);
	
    MatrixXf Jacobi(2,2);
    Jacobi=GN*coor;
    *JacobiDet=Jacobi.determinant();
    
    MatrixXf Jacobi_inv(2,2)=Jacobi;
    Jacobi_inv=Jacobi_inv.inverse;
    
    MatrixXf Bmat(2,4);
    Bmat= Jacobi_inv*GN;
}

inline void C4Q::CalBmatTotal(){
    #define Bmat_size 8
    double item=1/sqrt(3);
    double yita[4]={-item,-item,item,item};
    double psi[4]={-item,item,-item,item};
    double coor[8];
    MatrixXf coor(4,2);
    coor(0,0)=nodes[0]->XYZ[0];  coor(0,1)=nodes[0]->XYZ[1];
    coor(1,0)=nodes[1]->XYZ[0];  coor(1,1)=nodes[1]->XYZ[1];
    coor(2,0)=nodes[2]->XYZ[0];  coor(2,1)=nodes[2]->XYZ[1];
    coor(3,0)=nodes[3]->XYZ[0];  coor(3,1)=nodes[3]->XYZ[1];

    for(int i=0;i<4;i++)
	{CalBmatSingle(BmatTotal+i*Bmat_size,JacobiDet+i,coor,yita[i],psi[i]);}

}
inline unsigned int C4Q::SizeOfStiffnessMatrix() { return 78; }

void C4Q::ElementStiffnessSingle(MatrixXf* Stiff,MatrixXf* Bmat,MatrixXf* Cttmat,double JacobiDet, double weight){
    CalCttmat();
    MatrixXf NBmat(3,8);
	NBmat= MatrixXf::Zero(3,8);
	NBmat(0,0)=Bmat(0,0);NBmat(0,2)=Bmat(0,1);NBmat(0,4)=Bmat(0,2);NBmat(0,6)=Bmat(0,3);
	NBmat(1,1)=Bmat(1,0);NBmat(1,3)=Bmat(1,1);NBmat(1,5)=Bmat(1,2);NBmat(1,7)=Bmat(1,3);
	NBmat(2,0)=Bmat(1,0);NBmat(2,1)=Bmat(0,0);NBmat(2,2)=Bmat(1,1);NBmat(2,3)=Bmat(0,1);
	NBmat(2,4)=Bmat(1,2);NBmat(2,5)=Bmat(0,2);NBmat(2,6)=Bmat(1,3);NBmat(2,7)=Bmat(0,3);
	MatrixXf NBmatT;
	NBmatT=NBmat.transpose(); 
    MatrixXf Stiff(8,8);
	Stiff=weight*JacobiDet*NBmatT*Cttmat*NBmat;
}

void C4Q::ElementStiffness(MatrixXf* Stiff){
    clear(Stiff, SizeOfStiffnessMatrix());
    CalCttmat();
    CalBmatTotal();
    ElementStiffnessSingle(Stiff,BmatTotal,Cttmat,JacobiDet[0],1);
    ElementStiffnessSingle(Stiff,BmatTotal+8,Cttmat,JacobiDet[1],1);
    ElementStiffnessSingle(Stiff,BmatTotal+16,Cttmat,JacobiDet[2],1);
    ElementStiffnessSingle(Stiff,BmatTotal+24,Cttmat,JacobiDet[3],1);
    
}

void C4Q::ElementGravity(double* bodyforce, double Gravity)
{
    clear(bodyforce,12);
}

void C4Q::ElementStress(double* stress, double* Displacement)
{
    double d[8];
    int item=0;
    for(int i=0;i<12;i++){
        if(LocationMatrix[i])d[item++]=Displacement[LocationMatrix[i]-1];
    }
	
    MatrixXf sigma(4,3);
	for (int i = 0; i < 4; i++)
	{
        double* BB=BmatTotal+i*Bmat_size;
        sigma(i,0)=BB[0]*d[0]+BB[1]*d[2]+BB[2]*d[4]+BB[3]*d[6];
        sigma(i,1)=BB[4]*d[1]+BB[5]*d[3]+BB[6]*d[5]+BB[7]*d[7];
        sigma(i,2)=BB[4]*d[0]+BB[0]*d[1]+BB[5]*d[2]+BB[1]*d[3]+BB[6]*d[4]+BB[2]*d[5]+BB[7]*d[6]+BB[3]*d[7];
        double item=sigma(i,0);
        sigma(i,0)=sigma(i,0)*Cttmat(0,0)+sigma(i,1)*Cttmat(0,1);
        sigma(i,1)=item*Cttmat(0,1)+sigma(i,1)*Cttmat(0,0);
        sigma(i,2)=sigma(i,2)*Cttmat(0,2);
	}
    double a1=1+sqrt(3)/2,a2=1-sqrt(3)/2,b=-0.5;
    MatrixXf sigma2(4,3);
    for (int j=0;j<3;j++){
        sigma2(0,j)=a1*sigma(0,j)+b*sigma(1,j)+a2*sigma(2,j)+b*sigma(3,j);
        sigma2(1,j)=b*sigma(0,j)+a1*sigma(1,j)+b*sigma(2,j)+a2*sigma(3,j);
        sigma2(2,j)=a2*sigma(0,j)+b*sigma(1,j)+a1*sigma(2,j)+b*sigma(3,j);
        sigma2(3,j)=b*sigma(0,j)+a2*sigma(1,j)+b*sigma(2,j)+a1*sigma(3,j);
    }

    for (int i=0;i<4;i++){
        stress[i]=sqrt(sigma2(i,0)*sigma2(i,0)+sigma2(i,1)*sigma2(i,1)+4*sigma2(i,2)*sigma2(i,2));
    }
