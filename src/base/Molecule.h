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
public:
   Molecule();
   ~Molecule();
   vector<Atom*>* GetAtomVect(); 
   double* GetXyzCOM();
   void CalcXyzCOM();
   int GetTotalNumberAOs();
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons();
   void CalcTotalNumberValenceElectrons();
   double GetCoreRepulsionFirstDerivative(int indexAtomA, int indexAtomB, CartesianType axisA);
   void CalcTotalCoreRepulsionEnergy();
   double GetTotalCoreRepulsionEnergy();
   void OutputXyzCOM();
   void OutputTotalNumberAtomsAOsValenceelectrons();
   void OutputConfiguration();
   void OutputTotalCoreRepulsionEnergy();
   void CalcPrincipalAxes();
   void Rotate();
   void Translate();
   double GetDistanceAtoms(int atomAIndex, int atomBIndex);
   double GetDistanceAtoms(Atom* atomA, Atom* atomB);
private:
   vector<Atom*>* atomVect;
   double* xyzCOM;
   bool wasCalculatedXyzCOM;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   double totalCoreRepulsionEnergy;
   bool wasCalculatedTotalCoreRepulsionEnergy;
   void CalcInertiaTensor(double** inertiaTensor, double* inertiaTensorOrigin);
   void FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments);
   void Rotate(EularAngle eularAngle, double* rotatingOrigin, RotatedObjectType rotatedObj);
   void OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments);
   void OutputInertiaTensorOrigin(double* inertiaTensorOrigin);
   void OutputRotatingConditions(RotatingType rotatingType, double* rotatingOrigin, 
                                 double* rotatingAxis, double rotatingAngle, 
                                 EularAngle rotatingEularAngles);
   void OutputTranslatingConditions(double* translatingDifference);
   string messageTotalNumberAOs;
   string messageTotalNumberAtoms;
   string messageTotalNumberValenceElectrons;
   string messageConfiguration;
   string messageConfigurationTitleAU;
   string messageConfigurationTitleAng;
   string messageCoreRepulsion;
   string messageCoreRepulsionTitle;
   string messageCOM;
   string messageCOMTitleAU;
   string messageCOMTitleAng;
   string messageStartPrincipalAxes;
   string messageDonePrincipalAxes;
   string messagePrincipalAxes;
   string messagePrincipalAxesTitleAU;
   string messagePrincipalAxesTitleAng;
   string messageInertiaTensorOrigin;
   string messageInertiaTensorOriginTitleAU;
   string messageInertiaTensorOriginTitleAng;
   string messageStartRotate;
   string messageDoneRotate;
   string messageRotatingOrigin;
   string messageRotatingOriginTitleAU;
   string messageRotatingOriginTitleAng;
   string messageRotatingAxis;
   string messageRotatingAxisTitleAU;
   string messageRotatingAxisTitleAng;
   string messageRotatingAngle;
   string messageRotatingType;
   string messageRotatingEularAngles;
   string messageRotatingEularAnglesTitle;
   string messageStartTranslate;
   string messageDoneTranslate;
   string messageTranslatingDifference;
   string messageTranslatingDifferenceTitleAU;
   string messageTranslatingDifferenceTitleAng;
};

Molecule::Molecule(){
   this->atomVect = new vector<Atom*>;
   this->xyzCOM = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->wasCalculatedXyzCOM = false;
   this->wasCalculatedTotalCoreRepulsionEnergy = false;
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageConfiguration = "\tMolecular configration:\n";
   this->messageConfigurationTitleAU = "\t\t| i-th | atom type | x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageConfigurationTitleAng = "\t\t| i-th | atom type | x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageCoreRepulsion = "\tTotal core repulsion energy:\n";
   this->messageCoreRepulsionTitle = "\t\t| [a.u.] | [eV] |\n";
   this->messageCOM = "\tCenter of Mass:\n";
   this->messageCOMTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageCOMTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageStartPrincipalAxes = "**********  START: Principal Axes of Inertia  **********\n";
   this->messageDonePrincipalAxes =  "**********  DONE: Principal Axes of Inertia  ***********\n\n\n";
   this->messagePrincipalAxes = "\tPrincipal Axes:\n";
   this->messagePrincipalAxesTitleAU = "\t\t| inertia moments [a.u.] | x [a.u.] | y[a.u.] | z[a.u.] | (normalized)\n";
   this->messagePrincipalAxesTitleAng = "\t\t| inertia moments [g*angust**2/mol] | x [angst.] | y[angst.] | z[angst.] | (not normalized)\n";
   this->messageInertiaTensorOrigin = "\tInertia Tensor Origin:\n";
   this->messageInertiaTensorOriginTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageInertiaTensorOriginTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageStartRotate = "**********  START: Rotate molecule  **********\n";
   this->messageDoneRotate =  "**********  DONE: Rotate molecule  ***********\n\n\n";
   this->messageRotatingOrigin = "\tRotating Origin:\n";
   this->messageRotatingOriginTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageRotatingOriginTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAxis = "\tRotating Axis:\n";
   this->messageRotatingAxisTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageRotatingAxisTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAngle = "\tRotating Angle [degree]: ";
   this->messageRotatingType = "\tRotating Type: ";
   this->messageRotatingEularAngles = "\tRotating Eular Angles:\n";
   this->messageRotatingEularAnglesTitle = "\t\t| alpha[degree] | beta[degree] | gamma[degree] |\n";
   this->messageStartTranslate = "**********  START: Translate molecule  **********\n";
   this->messageDoneTranslate =  "**********  DONE: Translate molecule  ***********\n\n\n";
   this->messageTranslatingDifference = "\tTranslating Difference:\n";
   this->messageTranslatingDifferenceTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageTranslatingDifferenceTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
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
   if(this->xyzCOM != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->xyzCOM);
      //cout << "xyzCOM deleted\n";
   }
}

vector<Atom*>* Molecule::GetAtomVect(){
   return this->atomVect;
}

double* Molecule::GetXyzCOM(){
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   return this->xyzCOM;
}

void Molecule::CalcXyzCOM(){
   if(!this->wasCalculatedXyzCOM){
      double totalAtomicMass;
      Atom* atom;
      double* atomicXyz;
      double atomicMass;

      for(int j=0; j<3; j++){
         this->xyzCOM[j] = 0.0;
      }
      
      for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         atomicXyz = atom->GetXyz();
         atomicMass = atom->GetAtomicMass();
         totalAtomicMass += atomicMass;
         for(int j=0; j<3; j++){
            this->xyzCOM[j] += atomicXyz[j] * atomicMass;
         }
      }
      for(int i=0; i<3; i++){
         this->xyzCOM[i]/=totalAtomicMass;
      }
   }
   this->wasCalculatedXyzCOM = true;
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

// First derivative of the core repulsion related to the coordinate of atom A.
double Molecule::GetCoreRepulsionFirstDerivative(int indexAtomA, int indexAtomB, 
                                                 CartesianType axisA){
   double value=0.0;
   Atom* atomA = (*this->atomVect)[indexAtomA];
   Atom* atomB = (*this->atomVect)[indexAtomB];
   double distance = this->GetDistanceAtoms(indexAtomA, indexAtomB);
   value = atomA->GetCoreCharge()*atomB->GetCoreCharge();
   value *= (atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA])/distance;
   value *= -1.0/pow(distance,2.0);
   return value;
}

void Molecule::CalcTotalCoreRepulsionEnergy(){
   double energy = 0.0;
   double distance = 0.0;
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   for(int i=0; i<this->atomVect->size(); i++){
      atomA = (*this->atomVect)[i];
      for(int j=i+1; j<this->atomVect->size(); j++){
         atomB = (*this->atomVect)[j];
         distance = this->GetDistanceAtoms(i, j);
         energy += atomA->GetCoreCharge()*atomB->GetCoreCharge()/distance; 
      }
   }
   this->totalCoreRepulsionEnergy = energy;
   this->wasCalculatedTotalCoreRepulsionEnergy = true;
}

double Molecule::GetTotalCoreRepulsionEnergy(){
   if(!this->wasCalculatedTotalCoreRepulsionEnergy){
      this->CalcTotalCoreRepulsionEnergy();
   }
   return this->totalCoreRepulsionEnergy;
}

void Molecule::OutputTotalCoreRepulsionEnergy(){
   if(!this->wasCalculatedTotalCoreRepulsionEnergy){
      this->CalcTotalCoreRepulsionEnergy();
   }

   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   cout << this->messageCoreRepulsion;
   cout << this->messageCoreRepulsionTitle;
   printf("\t\t%e\t%e\n\n",this->totalCoreRepulsionEnergy, this->totalCoreRepulsionEnergy/eV2AU);

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

void Molecule::OutputXyzCOM(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageCOM;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOM[0]/ang2AU,
                             this->xyzCOM[1]/ang2AU,
                             this->xyzCOM[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOM[0],
                             this->xyzCOM[1],
                             this->xyzCOM[2]);
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

void Molecule::OutputInertiaTensorOrigin(double* inertiaTensorOrigin){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();

   cout << this->messageInertiaTensorOrigin;
   cout << this->messageInertiaTensorOriginTitleAng;
   printf("\t\t%e\t%e\t%e\n",inertiaTensorOrigin[0]/ang2AU,
                             inertiaTensorOrigin[1]/ang2AU,
                             inertiaTensorOrigin[2]/ang2AU);
   cout << "\n";

   cout << this->messageInertiaTensorOriginTitleAU;
   printf("\t\t%e\t%e\t%e\n",inertiaTensorOrigin[0],
                             inertiaTensorOrigin[1],
                             inertiaTensorOrigin[2]);
   cout << "\n";

}

void Molecule::CalcPrincipalAxes(){

   cout << this->messageStartPrincipalAxes;

   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double inertiaTensorOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetInertiaTensorOrigin() != NULL){
      inertiaTensorOrigin[0] = Parameters::GetInstance()->GetInertiaTensorOrigin()[0];
      inertiaTensorOrigin[1] = Parameters::GetInstance()->GetInertiaTensorOrigin()[1];
      inertiaTensorOrigin[2] = Parameters::GetInstance()->GetInertiaTensorOrigin()[2];
   }

   double** inertiaTensor = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(3, 3);
   double*  inertiaMoments = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);

   try{
      this->CalcInertiaTensor(inertiaTensor, inertiaTensorOrigin);
      
      bool calcEigenVectors = true;
      MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(inertiaTensor,
                                                              inertiaMoments,
                                                              3,
                                                              calcEigenVectors);
      this->OutputPrincipalAxes(inertiaTensor, inertiaMoments);
      this->OutputInertiaTensorOrigin(inertiaTensorOrigin);
   }
   catch(MolDSException ex){
      this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);
      throw ex;
   }
   this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);

   cout << this->messageDonePrincipalAxes;
   
}

void Molecule::CalcInertiaTensor(double** inertiaTensor, double* inertiaTensorOrigin){

   Atom* atom;
   double x;
   double y;
   double z;
   double atomicMass;
   for(int a=0; a<this->atomVect->size(); a++){
      atom = (*this->atomVect)[a];
      atomicMass = atom->GetAtomicMass();
      x = atom->GetXyz()[0] - inertiaTensorOrigin[0];
      y = atom->GetXyz()[1] - inertiaTensorOrigin[1];
      z = atom->GetXyz()[2] - inertiaTensorOrigin[2];

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

void Molecule::FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments){

   if(*inertiaTensor != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(inertiaTensor, 3);
      //cout << "inertiaTensor deleted\n";
   }

   if(*inertiaMoments != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(inertiaMoments);
      //cout << "inertiaMoments deleted\n";
   }

}

void Molecule::Rotate(){

   cout << this->messageStartRotate;

   // Default values are set if some conditions are not specified.
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double rotatingOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetRotatingOrigin() != NULL){
      rotatingOrigin[0] = Parameters::GetInstance()->GetRotatingOrigin()[0];
      rotatingOrigin[1] = Parameters::GetInstance()->GetRotatingOrigin()[1];
      rotatingOrigin[2] = Parameters::GetInstance()->GetRotatingOrigin()[2];
   }

   RotatingType rotatingType = Parameters::GetInstance()->GetRotatingType();
   double* rotatingAxis = Parameters::GetInstance()->GetRotatingAxis();
   EularAngle rotatingEularAngles = Parameters::GetInstance()->GetRotatingEularAngles();
   double rotatingAngle = Parameters::GetInstance()->GetRotatingAngle();

   this->OutputRotatingConditions(rotatingType, rotatingOrigin, 
                                  rotatingAxis, rotatingAngle, 
                                  rotatingEularAngles);

   // rotate
   if(rotatingType == Axis){
      EularAngle setZAxisEularAngles(rotatingAxis[0], rotatingAxis[1], rotatingAxis[2]);
      EularAngle angleAroundAxis;
      angleAroundAxis.SetAlpha(rotatingAngle);

      this->Rotate(setZAxisEularAngles, rotatingOrigin, Frame);
      this->Rotate(angleAroundAxis, rotatingOrigin, System);
      this->Rotate(setZAxisEularAngles, rotatingOrigin, System);
   }
   else if(rotatingType == Eular){
      this->Rotate(rotatingEularAngles, rotatingOrigin, System);
   }
   
   this->OutputConfiguration();
   cout << this->messageDoneRotate;
}

/***
 * rotatedObj == System: Molecule is rotated.
 * rotatedObj == Frame: De Cartesian is rotated.
 */
void Molecule::Rotate(EularAngle eularAngle, double* rotatingOrigin, RotatedObjectType rotatedObj){

   double rotatingMatrixAlpha[3][3];
   double rotatingMatrixBeta[3][3];
   double rotatingMatrixGamma[3][3];
   double inv = 1.0;
   if(rotatedObj == System){
      inv = -1.0;
   }

   CalcRotatingMatrix(rotatingMatrixAlpha, inv*eularAngle.GetAlpha(), ZAxis);
   CalcRotatingMatrix(rotatingMatrixBeta, inv*eularAngle.GetBeta(), YAxis);
   CalcRotatingMatrix(rotatingMatrixGamma, inv*eularAngle.GetGamma(), ZAxis);

   double temp1[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp1[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixGamma[k][j];
            }
            else if(rotatedObj == Frame){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixAlpha[k][j];
            }
         }
      }
   }

   double temp2[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp2[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp2[i][j] += rotatingMatrixAlpha[i][k] * temp1[k][j];
            }
            else if(rotatedObj == Frame){
               temp2[i][j] += rotatingMatrixGamma[i][k] * temp1[k][j];
            }
         }
      }
   }

   double rotatedXyz[3];
   Atom* atom;
   for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         for(int j=0; j<3; j++){
            rotatedXyz[j] = 0.0;
            for(int k=0; k<3; k++){
               rotatedXyz[j] += temp2[j][k] * (atom->GetXyz()[k] - rotatingOrigin[k]);
            }
         }
         for(int j=0; j<3; j++){
            atom->GetXyz()[j] = rotatedXyz[j] + rotatingOrigin[j];
         }
   }
}

void Molecule::OutputRotatingConditions(RotatingType rotatingType, double* rotatingOrigin, 
                                        double* rotatingAxis, double rotatingAngle, 
                                        EularAngle rotatingEularAngles){

   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double degree2Radian = Parameters::GetInstance()->GetDegree2Radian();

   // type
   printf("%s%s\n\n",this->messageRotatingType.c_str(),RotatingTypeStr(rotatingType));

   // rotating origin
   cout << this->messageRotatingOrigin;
   cout << this->messageRotatingOriginTitleAng;
   printf("\t\t%e\t%e\t%e\n\n",rotatingOrigin[0]/angst2AU,
                               rotatingOrigin[1]/angst2AU,
                               rotatingOrigin[2]/angst2AU);

   cout << this->messageRotatingOriginTitleAU;
   printf("\t\t%e\t%e\t%e\n\n",rotatingOrigin[0],
                               rotatingOrigin[1],
                               rotatingOrigin[2]);

   if(rotatingType == Axis){
      // rotating axis
      cout << this->messageRotatingAxis;
      cout << this->messageRotatingAxisTitleAng;
      printf("\t\t%e\t%e\t%e\n\n",rotatingAxis[0]/angst2AU,
                                  rotatingAxis[1]/angst2AU,
                                  rotatingAxis[2]/angst2AU);

      cout << this->messageRotatingAxisTitleAU;
      printf("\t\t%e\t%e\t%e\n\n",rotatingAxis[0],
                                  rotatingAxis[1],
                                  rotatingAxis[2]);

      // angle
      printf("%s%e\n\n",this->messageRotatingAngle.c_str(),rotatingAngle/degree2Radian);
   }
   else if (rotatingType == Eular){
      // Eular angles
      cout << this->messageRotatingEularAngles;
      cout << this->messageRotatingEularAnglesTitle;
      printf("\t\t%e\t%e\t%e\t\n\n",rotatingEularAngles.GetAlpha()/degree2Radian, 
                                    rotatingEularAngles.GetBeta()/degree2Radian,
                                    rotatingEularAngles.GetGamma()/degree2Radian);
   }

}


void Molecule::Translate(){

   cout << this->messageStartTranslate;

   double x = Parameters::GetInstance()->GetTranslatingDifference()[0];
   double y = Parameters::GetInstance()->GetTranslatingDifference()[1];
   double z = Parameters::GetInstance()->GetTranslatingDifference()[2];

   this->OutputTranslatingConditions(Parameters::GetInstance()->GetTranslatingDifference()); 

   Atom* atom;
   for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         atom->GetXyz()[0] += x;
         atom->GetXyz()[1] += y;
         atom->GetXyz()[2] += z;
   }
   
   this->wasCalculatedXyzCOM = false;
   this->CalcXyzCOM();

   this->OutputConfiguration();
   this->OutputXyzCOM();

   cout << this->messageDoneTranslate;
}

void Molecule::OutputTranslatingConditions(double* translatingDifference){

   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();

   cout << this->messageTranslatingDifference;
   cout << this->messageTranslatingDifferenceTitleAng;
   printf("\t\t%e\t%e\t%e\n\n",translatingDifference[0]/angst2AU,
                               translatingDifference[1]/angst2AU,
                               translatingDifference[2]/angst2AU);

   cout << this->messageTranslatingDifferenceTitleAU;
   printf("\t\t%e\t%e\t%e\n\n",translatingDifference[0],
                               translatingDifference[1],
                               translatingDifference[2]);

}

double Molecule::GetDistanceAtoms(int atomAIndex, int atomBIndex){
   Atom* atomA = (*this->atomVect)[atomAIndex];
   Atom* atomB = (*this->atomVect)[atomBIndex];
   return this->GetDistanceAtoms(atomA, atomB);
}

double Molecule::GetDistanceAtoms(Atom* atomA, Atom* atomB){

   double distance=0.0;
   distance = sqrt( pow(atomA->GetXyz()[0] - atomB->GetXyz()[0], 2.0)
                   +pow(atomA->GetXyz()[1] - atomB->GetXyz()[1], 2.0)
                   +pow(atomA->GetXyz()[2] - atomB->GetXyz()[2], 2.0) );
   return distance;

}





}
#endif





