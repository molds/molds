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
   double GetBondingAdjustParameterK();
   TheoryType GetCurrentTheory();
   void SetCurrentTheory(TheoryType theory);

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
   double bondingAdjustParameterK; //see (3.79) in J. A. Pople book
   TheoryType currentTheory;

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
   this->kayser2AU = 4.556336*pow(10.0,-6.0);
   this->thresholdSCF = pow(10.0, -8.0);
   this->maxIterationsSCF = 100;
   this->currentTheory = CNDO2;
   this->gMolin2AU = pow(10.0,5.0)/(6.0221415*9.1095);
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

double Parameters::GetBondingAdjustParameterK(){
   return this->bondingAdjustParameterK;
}

TheoryType Parameters::GetCurrentTheory(){
   return this->currentTheory;
}

void Parameters::SetCurrentTheory(TheoryType theory){
   this->currentTheory = theory;
}

}
#endif





