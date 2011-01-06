#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include"atoms/Atom.h"
#include"../mkl_wrapper/LapackWrapper.h"

using namespace std;
using namespace MolDS_base_atoms;

namespace MolDS_base{

class Molecule{
private:
   vector<Atom*>* atomVect;
   double* COMXyz;
   double* inertiaTensorOrigin;
   bool wasCalculatedCOMXyz;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   void CalcInertiaTensor(double** inertiaTensor);
   void FreeInertiaTensorMoments(double** inertiaTensor, double* inertiaMoments);
   void OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments);
   void OutputInertiaTensorOrigin();
   string messageTotalNumberAOs;
   string messageTotalNumberAtoms;
   string messageTotalNumberValenceElectrons;
   string messageConfiguration;
   string messageConfigurationTitleAU;
   string messageConfigurationTitleAng;
   string messageCOM;
   string messageCOMTitleAU;
   string messageCOMTitleAng;
   string messagePrincipalAxes;
   string messagePrincipalAxesTitleAU;
   string messagePrincipalAxesTitleAng;
   string messageInertiaTensorOrigin;
   string messageInertiaTensorOriginTitleAU;
   string messageInertiaTensorOriginTitleAng;
public:
   Molecule();
   ~Molecule();
   vector<Atom*>* GetAtomVect(); 
   double* GetCOMXyz();
   void CalcCOMXyz();
   int GetTotalNumberAOs();
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons();
   void CalcTotalNumberValenceElectrons();
   void OutputCOMXyz();
   void OutputTotalNumberAtomsAOsValenceelectrons();
   void OutputConfiguration();
   void SetInertiaTensorOrigin(double x, double y, double z);
   void CalcPrincipalAxes();
};

Molecule::Molecule(){
   this->atomVect = new vector<Atom*>;
   this->COMXyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->inertiaTensorOrigin = NULL;
   this->wasCalculatedCOMXyz = false;
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageConfiguration = "\tMolecular configration:\n";
   this->messageConfigurationTitleAU = "\t\t| i-th | atom type | x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageConfigurationTitleAng = "\t\t| i-th | atom type | x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageCOM = "\tCenter of Mass:\n";
   this->messageCOMTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageCOMTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messagePrincipalAxes = "\tPrincipal Axes:\n";
   this->messagePrincipalAxesTitleAU = "\t\t| inertia moments [a.u.] | x [a.u.] | y[a.u.] | z[a.u.] | (normalized)\n";
   this->messagePrincipalAxesTitleAng = "\t\t| inertia moments [g*angust**2/mol] | x [angst.] | y[angst.] | z[angst.] | (not normalized)\n";
   this->messageInertiaTensorOrigin = "\tInertia Tensor Origin:\n";
   this->messageInertiaTensorOriginTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageInertiaTensorOriginTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
}

Molecule::~Molecule(){
   if(this->atomVect != NULL){
      for(int i=0; i<this->atomVect->size(); i++){
         delete (*atomVect)[i];
         (*atomVect)[i] = NULL;
      }
      delete this->atomVect;
      this->atomVect = NULL;
      //cout << "atomVect deleted\n";
   }
   if(this->COMXyz != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(this->COMXyz);
      this->COMXyz = NULL;
      //cout << "COMXyz deleted\n";
   }
   if(this->inertiaTensorOrigin != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(this->inertiaTensorOrigin);
      this->inertiaTensorOrigin = NULL;
      //cout << "inertiaTensorOrigin deleted\n";
   }
}

vector<Atom*>* Molecule::GetAtomVect(){
   return this->atomVect;
}

double* Molecule::GetCOMXyz(){
   if(!this->wasCalculatedCOMXyz){
      this->CalcCOMXyz();
   }
   return this->COMXyz;
}

void Molecule::CalcCOMXyz(){
   if(!this->wasCalculatedCOMXyz){
      double totalAtomicMass;
      Atom* atom;
      double* atomicXyz;
      double atomicMass;
      for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         atomicXyz = atom->GetXyz();
         atomicMass = atom->GetAtomicMass();
         totalAtomicMass += atomicMass;
         for(int j=0; j<3; j++){
            this->COMXyz[j] += atomicXyz[j] * atomicMass;
         }
      }
      for(int i=0; i<3; i++){
         this->COMXyz[i]/=totalAtomicMass;
      }
   }
   this->wasCalculatedCOMXyz = true;
}

int Molecule::GetTotalNumberAOs(){
   return this->totalNumberAOs;
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAOs = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAOs);
      this->totalNumberAOs += (*this->atomVect)[i]->GetValence().size();
   }
}

int Molecule::GetTotalNumberValenceElectrons(){
   return this->totalNumberValenceElectrons;
}

void Molecule::CalcTotalNumberValenceElectrons(){
   this->totalNumberValenceElectrons = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      this->totalNumberValenceElectrons += (*this->atomVect)[i]->GetNumberValenceElectrons();
   }
}

void Molecule::OutputConfiguration(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageConfiguration;
   cout << this->messageConfigurationTitleAng;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0]/ang2AU, atom->GetXyz()[1]/ang2AU, atom->GetXyz()[2]/ang2AU);
   }
   cout << "\n";

   cout << this->messageConfigurationTitleAU;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0], atom->GetXyz()[1], atom->GetXyz()[2]);
   }
   cout << "\n";

}

void Molecule::OutputCOMXyz(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageCOM;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->COMXyz[0]/ang2AU,
                             this->COMXyz[1]/ang2AU,
                             this->COMXyz[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->COMXyz[0],
                             this->COMXyz[1],
                             this->COMXyz[2]);
   cout << "\n";

}

void Molecule::OutputTotalNumberAtomsAOsValenceelectrons(){
   cout << this->messageTotalNumberAtoms << this->atomVect->size() << "\n";
   cout << this->messageTotalNumberAOs << this->totalNumberAOs << "\n";
   cout << this->messageTotalNumberValenceElectrons << this->totalNumberValenceElectrons << "\n\n";
}

void Molecule::OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double gMolin2AU = Parameters::GetInstance()->GetGMolin2AU();

   cout << this->messagePrincipalAxes;
   cout << this->messagePrincipalAxesTitleAng;
   for(int i=0; i<3; i++){
      printf("\t\t%e\t%e\t%e\t%e\n",inertiaMoments[i]/gMolin2AU, 
                                    inertiaTensor[i][0]/ang2AU,
                                    inertiaTensor[i][1]/ang2AU,
                                    inertiaTensor[i][2]/ang2AU);
   }
   cout << "\n";

   cout << this->messagePrincipalAxesTitleAU;
   for(int i=0; i<3; i++){
      printf("\t\t%e\t%e\t%e\t%e\n",inertiaMoments[i], 
                                    inertiaTensor[i][0],
                                    inertiaTensor[i][1],
                                    inertiaTensor[i][2]);
   }
   cout << "\n";

}

void Molecule::OutputInertiaTensorOrigin(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();

   cout << this->messageInertiaTensorOrigin;
   cout << this->messageInertiaTensorOriginTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->inertiaTensorOrigin[0]/ang2AU,
                                 this->inertiaTensorOrigin[1]/ang2AU,
                                 this->inertiaTensorOrigin[2]/ang2AU);
   cout << "\n";

   cout << this->messageInertiaTensorOriginTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->inertiaTensorOrigin[0],
                                 this->inertiaTensorOrigin[1],
                                 this->inertiaTensorOrigin[2]);
   cout << "\n";

}

void Molecule::SetInertiaTensorOrigin(double x, double y, double z){
   if(this->inertiaTensorOrigin == NULL){
      this->inertiaTensorOrigin = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   }

   this->inertiaTensorOrigin[0] = x;
   this->inertiaTensorOrigin[1] = y;
   this->inertiaTensorOrigin[2] = z;

}

void Molecule::CalcPrincipalAxes(){

   if(this->inertiaTensorOrigin == NULL){
      this->SetInertiaTensorOrigin(this->COMXyz[0], this->COMXyz[1], this->COMXyz[2]);
   }

   double** inertiaTensor = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(3, 3);
   double*  inertiaMoments = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);

   try{
      this->CalcInertiaTensor(inertiaTensor);
      
      bool calcEigenVectors = true;
      MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(inertiaTensor,
                                                              inertiaMoments,
                                                              3,
                                                              calcEigenVectors);
      this->OutputPrincipalAxes(inertiaTensor, inertiaMoments);
      this->OutputInertiaTensorOrigin();
   }
   catch(MolDSException ex){
      this->FreeInertiaTensorMoments(inertiaTensor, inertiaMoments);
      throw ex;
   }

   this->FreeInertiaTensorMoments(inertiaTensor, inertiaMoments);
   
}

void Molecule::CalcInertiaTensor(double** inertiaTensor){

   Atom* atom;
   double x;
   double y;
   double z;
   double atomicMass;
   for(int a=0; a<this->atomVect->size(); a++){
      atom = (*this->atomVect)[a];
      atomicMass = atom->GetAtomicMass();
      x = atom->GetXyz()[0] - this->inertiaTensorOrigin[0];
      y = atom->GetXyz()[1] - this->inertiaTensorOrigin[1];
      z = atom->GetXyz()[2] - this->inertiaTensorOrigin[2];

      inertiaTensor[0][0] += atomicMass*(y*y + z*z);
      inertiaTensor[0][1] -= atomicMass*x*y;
      inertiaTensor[0][2] -= atomicMass*x*z;

      inertiaTensor[1][0] -= atomicMass*y*x;
      inertiaTensor[1][1] += atomicMass*(x*x + z*z);
      inertiaTensor[1][2] -= atomicMass*y*z;

      inertiaTensor[2][0] -= atomicMass*z*x;
      inertiaTensor[2][1] -= atomicMass*z*y;
      inertiaTensor[2][2] += atomicMass*(x*x + y*y);

   }
   
}

void Molecule::FreeInertiaTensorMoments(double** inertiaTensor, double* inertiaMoments){

   if(inertiaTensor != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(inertiaTensor, 3);
      inertiaTensor = NULL;
      //cout << "inertiaTensor deleted\n";
   }

   if(inertiaMoments != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(inertiaMoments);
      inertiaMoments = NULL;
      //cout << "inertiaMoments deleted\n";
   }

}



}
#endif





