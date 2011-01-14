#ifndef INCLUDED_PARAMETERS
#define INCLUDED_PARAMETERS

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>

using namespace std;
namespace MolDS_base{

// Parameters is singleton
class Parameters{
public:
   static Parameters* GetInstance();
   static void DeleteInstance();

   double GetThresholdSCF();
   void SetThresholdSCF(double thresholdSCF);
   int GetMaxIterationsSCF();
   void SetMaxIterationsSCF(int maxIterationsSCF);
   double GetEV2AU();
   double GetAngstrom2AU();
   double GetKayser2AU();
   double GetGMolin2AU();
   double GetDegree2Radian();
   double GetBondingAdjustParameterK();
   TheoryType GetCurrentTheory();
   void SetCurrentTheory(TheoryType theory);
   void SetTranslatingDifference(double x, double y, double z);
   double* GetTranslatingDifference();
   void SetInertiaTensorOrigin(double x, double y, double z);
   double* GetInertiaTensorOrigin();
   void SetRotatingOrigin(double x, double y, double z);
   double* GetRotatingOrigin();
   void SetRotatingType(RotatingType rotatingType);
   RotatingType GetRotatingType();
   void SetRotatingAxis(double x, double y, double z);
   double* GetRotatingAxis();
   void SetRotatingAngle(double rotatingAngle);
   double GetRotatingAngle();
   void SetRotatingEularAngles(double alpha, double beta, double gamma);
   EularAngle GetRotatingEularAngles();
   int GetActiveOccCIS();
   void SetActiveOccCIS(int activeOccCIS);
   int GetActiveVirCIS();
   void SetActiveVirCIS(int activeOccCIS);
   int GetNumberExcitedStatesCIS();
   void SetNumberExcitedStatesCIS(int nStates);
private:
   static Parameters* parameters;
   Parameters();
   Parameters(Parameters&);
   void operator = (Parameters&);
   ~Parameters();

   double thresholdSCF;
   int maxIterationsSCF;
   void SetDefaultValues();
   double eV2AU;
   double angstrom2AU;
   double kayser2AU;
   double gMolin2AU;
   double degree2Radian;
   double bondingAdjustParameterK; //see (3.79) in J. A. Pople book
   TheoryType currentTheory;
   double translatingDifference[3];
   double* inertiaTensorOrigin;
   double* rotatingOrigin;
   double rotatingAxis[3];
   double rotatingAngle;
   RotatingType rotatingType;
   EularAngle rotatingEularAngles;
   int activeOccCIS;
   int activeVirCIS;
   int numberExcitedStatesCIS;
};
Parameters* Parameters::parameters = NULL;

Parameters::Parameters(){
   this->SetDefaultValues();
}

Parameters::~Parameters(){
   if(this->inertiaTensorOrigin != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(this->inertiaTensorOrigin);
      this->inertiaTensorOrigin = NULL;
      //cout << "inertiaTensorOrigin deleted\n";
   }
   if(this->rotatingOrigin != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(this->rotatingOrigin);
      this->rotatingOrigin = NULL;
      //cout << "rotatingOrigin deleted\n";
   }

}

Parameters* Parameters::GetInstance(){
   if(parameters == NULL){
      parameters = new Parameters();
   }
   return parameters;
}

void Parameters::DeleteInstance(){
   if(parameters != NULL){
      delete parameters; 
   }
   parameters = NULL;
}

void Parameters::SetDefaultValues(){
   this->eV2AU = 0.03674903;
   this->angstrom2AU = 1.0/0.5291772;
   this->kayser2AU = 4.556336*pow(10.0,-6.0);
   this->thresholdSCF = pow(10.0, -8.0);
   this->maxIterationsSCF = 100;
   this->currentTheory = CNDO2;
   this->gMolin2AU = pow(10.0,5.0)/(6.0221415*9.1095);
   this->degree2Radian = M_PI / 180.0;
   this->translatingDifference[0] = 0.0;
   this->translatingDifference[1] = 0.0;
   this->translatingDifference[2] = 0.0;
   this->inertiaTensorOrigin = NULL;
   this->rotatingOrigin = NULL;
   this->rotatingAxis[0] = 0.0;
   this->rotatingAxis[1] = 0.0;
   this->rotatingAxis[2] = 1.0;
   this->rotatingType = Axis;
   this->rotatingEularAngles.SetAlpha(0.0);
   this->rotatingEularAngles.SetBeta(0.0);
   this->rotatingEularAngles.SetGamma(0.0);
   this->activeOccCIS = 10;
   this->activeVirCIS = 10;
   this->numberExcitedStatesCIS = 5;
}

double Parameters::GetThresholdSCF(){
   return this->thresholdSCF;
}

void Parameters::SetThresholdSCF(double thresholdSCF){
   this->thresholdSCF = thresholdSCF;
}

int Parameters::GetMaxIterationsSCF(){
   return this->maxIterationsSCF;
}

void Parameters::SetMaxIterationsSCF(int maxIterationsSCF){
   this->maxIterationsSCF = maxIterationsSCF;
}

double Parameters::GetEV2AU(){
   return this->eV2AU;
}

double Parameters::GetAngstrom2AU(){
   return this->angstrom2AU;
}

double Parameters::GetKayser2AU(){
   return this->kayser2AU;
}

double Parameters::GetGMolin2AU(){
   return this->gMolin2AU;
}

double Parameters::GetDegree2Radian(){
   return this->degree2Radian;
}

double Parameters::GetBondingAdjustParameterK(){
   return this->bondingAdjustParameterK;
}

TheoryType Parameters::GetCurrentTheory(){
   return this->currentTheory;
}

void Parameters::SetCurrentTheory(TheoryType theory){
   this->currentTheory = theory;
}

void Parameters::SetTranslatingDifference(double x, double y, double z){
   this->translatingDifference[0] = x;
   this->translatingDifference[1] = y;
   this->translatingDifference[2] = z;
}

double* Parameters::GetTranslatingDifference(){
   return this->translatingDifference;
}

void Parameters::SetInertiaTensorOrigin(double x, double y, double z){
   if(this->inertiaTensorOrigin == NULL){
      this->inertiaTensorOrigin = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   }

   this->inertiaTensorOrigin[0] = x;
   this->inertiaTensorOrigin[1] = y;
   this->inertiaTensorOrigin[2] = z;

}

double* Parameters::GetInertiaTensorOrigin(){
   return this->inertiaTensorOrigin;
}

void Parameters::SetRotatingOrigin(double x, double y, double z){
   if(this->rotatingOrigin == NULL){
      this->rotatingOrigin = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   }

   this->rotatingOrigin[0] = x;
   this->rotatingOrigin[1] = y;
   this->rotatingOrigin[2] = z;

}

double* Parameters::GetRotatingOrigin(){
   return this->rotatingOrigin;
}

void Parameters::SetRotatingType(RotatingType rotatingType){
   this->rotatingType = rotatingType;
}

RotatingType Parameters::GetRotatingType(){
   return this->rotatingType;
}

void Parameters::SetRotatingAxis(double x, double y, double z){
   this->rotatingAxis[0] = x;
   this->rotatingAxis[1] = y;
   this->rotatingAxis[2] = z;

}

double* Parameters::GetRotatingAxis(){
   return this->rotatingAxis;
}

void Parameters::SetRotatingAngle(double rotatingAngle){
   this->rotatingAngle = rotatingAngle;
}

double Parameters::GetRotatingAngle(){
   return this->rotatingAngle;
}

void Parameters::SetRotatingEularAngles(double alpha, double beta, double gamma){
   this->rotatingEularAngles.SetAlpha(alpha);
   this->rotatingEularAngles.SetBeta(beta);
   this->rotatingEularAngles.SetGamma(gamma);
}

EularAngle Parameters::GetRotatingEularAngles(){
   return this->rotatingEularAngles;
}

int Parameters::GetActiveOccCIS(){
   return this->activeOccCIS;
}
   
void Parameters::SetActiveOccCIS(int activeOccCIS){
   this->activeOccCIS = activeOccCIS;
}

int Parameters::GetActiveVirCIS(){
   return this->activeVirCIS;
}

void Parameters::SetActiveVirCIS(int activeVirCIS){
   this->activeVirCIS = activeVirCIS;
}

int Parameters::GetNumberExcitedStatesCIS(){
   return this->numberExcitedStatesCIS;
}

void Parameters::SetNumberExcitedStatesCIS(int nStates){
   this->numberExcitedStatesCIS = nStates;
}


}
#endif





