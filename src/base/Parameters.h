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
   double bondingAdjustParameterK; //see (3.79) in J. A. Pople book

public:
   static Parameters* GetInstance();
   static void DeleteInstance();

   double GetThresholdSCF();
   void SetThresholdSCF(double thresholdSCF);
   int GetMaxIterationsSCF();
   void SetMaxIterationsSCF(int maxIterationsSCF);
   double GetEV2AU();
   double GetAngstrom2AU();
   double GetBondingAdjustParameterK();

};
Parameters* Parameters::parameters = NULL;

Parameters::Parameters(){
   this->SetDefaultValues();
}

Parameters::~Parameters(){
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
   this->thresholdSCF = pow(10.0, -8.0);
   this->maxIterationsSCF = 100;
   this->bondingAdjustParameterK = 0.75;
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

double Parameters::GetBondingAdjustParameterK(){
   return this->bondingAdjustParameterK;
}

}
#endif





