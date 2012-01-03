#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include"Uncopyable.h"
#include"Enums.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
using namespace std;
namespace MolDS_base{

Parameters* Parameters::parameters = NULL;
const double Parameters::Parameters::eV2AU = 0.03674903;
const double Parameters::kcalMolin2AU = 0.00159360175;
const double Parameters::angstrom2AU = 1.0/0.5291772;
const double Parameters::kayser2AU = 4.556336*pow(10.0,-6.0);
const double Parameters::fs2AU = 1.0/(2.418884326505*pow(10.0,-2.0));
const double Parameters::gMolin2AU = pow(10.0,5.0)/(6.0221415*9.1095);
const double Parameters::degree2Radian = M_PI / 180.0;

Parameters::Parameters(){
   this->SetDefaultValues();
   this->indecesMOPlot = new vector<int>;
}

Parameters::~Parameters(){
   if(this->inertiaTensorOrigin != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->inertiaTensorOrigin);
      //cout << "inertiaTensorOrigin deleted\n";
   }
   if(this->rotatingOrigin != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->rotatingOrigin);
      //cout << "rotatingOrigin deleted\n";
   }
   if(this->indecesMOPlot != NULL){
      delete this->indecesMOPlot;
      this->indecesMOPlot = NULL;
      //cout << "indecesMOPlot deleted\n";
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
   this->thresholdSCF = pow(10.0, -8.0);
   this->maxIterationsSCF = 100;
   this->dampingThreshSCF = 1.0;
   this->dampingWeightSCF = 0.8;
   this->diisNumErrorVectSCF = 5;
   this->diisStartErrorSCF = pow(10.0, -2.0);
   this->diisEndErrorSCF = pow(10.0, -8.0);
   this->fileNamePrefixMOPlot = "MO_";
   this->gridNumberMOPlot[XAxis] = 25;
   this->gridNumberMOPlot[YAxis] = 25;
   this->gridNumberMOPlot[ZAxis] = 25;
   this->frameLengthMOPlot[XAxis] = 20.0;
   this->frameLengthMOPlot[YAxis] = 20.0;
   this->frameLengthMOPlot[ZAxis] = 20.0;
   this->currentTheory = CNDO2;
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
   this->requiresCIS = false;
   this->isDavidsonCIS = true;
   this->maxIterationsCIS = 100;
   this->maxDimensionsCIS = 100;
   this->normToleranceCIS = pow(10.0, -6.0);
   this->electronicStateIndexMD = 0;
   this->totalStepsMD = 10;
   this->timeWidthMD = 0.1*this->fs2AU;
}

double Parameters::GetThresholdSCF() const{
   return this->thresholdSCF;
}

void Parameters::SetThresholdSCF(double thresholdSCF){
   this->thresholdSCF = thresholdSCF;
}

int Parameters::GetMaxIterationsSCF() const{
   return this->maxIterationsSCF;
}

void Parameters::SetMaxIterationsSCF(int maxIterationsSCF){
   this->maxIterationsSCF = maxIterationsSCF;
}

double Parameters::GetDampingThreshSCF() const{
   return this->dampingThreshSCF;
}

void Parameters::SetDampingThreshSCF(double dampingThreshSCF){
   this->dampingThreshSCF = dampingThreshSCF;
}

double Parameters::GetDampingWeightSCF() const{
   return this->dampingWeightSCF;
}

void Parameters::SetDampingWeightSCF(double dampingWeightSCF){
   this->dampingWeightSCF = dampingWeightSCF;
}

int Parameters::GetDiisNumErrorVectSCF() const{
   return this->diisNumErrorVectSCF;
}

void Parameters::SetDiisNumErrorVectSCF(int diisNumErrorVectSCF){
   this->diisNumErrorVectSCF = diisNumErrorVectSCF;
}

double Parameters::GetDiisStartErrorSCF() const{
   return this->diisStartErrorSCF;
}

void Parameters::SetDiisStartErrorSCF(double diisStartErrorSCF){
   this->diisStartErrorSCF = diisStartErrorSCF;
}

double Parameters::GetDiisEndErrorSCF() const{
   return this->diisEndErrorSCF;
}

void Parameters::SetDiisEndErrorSCF(double diisEndErrorSCF){
   this->diisEndErrorSCF = diisEndErrorSCF;
}

string Parameters::GetFileNamePrefixMOPlot() const{
   return this->fileNamePrefixMOPlot;
}

void Parameters::SetFileNamePrefixMOPlot(string fileNamePrefixMOPlot){
   this->fileNamePrefixMOPlot = fileNamePrefixMOPlot;
}

int* Parameters::GetGridNumberMOPlot() const{
   return (int*)this->gridNumberMOPlot;
}

void Parameters::SetGridNumberMOPlot(int Nx, int Ny, int Nz){
   this->gridNumberMOPlot[XAxis] = Nx;
   this->gridNumberMOPlot[YAxis] = Ny;
   this->gridNumberMOPlot[ZAxis] = Nz;
}

double* Parameters::GetFrameLengthMOPlot() const{
   return (double*)this->frameLengthMOPlot;
}

void Parameters::SetFrameLengthMOPlot(double lx, double ly, double lz){
   this->frameLengthMOPlot[XAxis] = lx;
   this->frameLengthMOPlot[YAxis] = ly;
   this->frameLengthMOPlot[ZAxis] = lz;
}

double Parameters::GetEV2AU() const{
   return this->eV2AU;
}

double Parameters::GetKcalMolin2AU() const{
   return this->kcalMolin2AU;
}

double Parameters::GetAngstrom2AU() const{
   return this->angstrom2AU;
}

double Parameters::GetKayser2AU() const{
   return this->kayser2AU;
}

double Parameters::GetGMolin2AU() const{
   return this->gMolin2AU;
}

double Parameters::GetDegree2Radian() const{
   return this->degree2Radian;
}

double Parameters::GetFs2AU() const{
   return this->fs2AU;
}

TheoryType Parameters::GetCurrentTheory() const{
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

double* Parameters::GetTranslatingDifference() const{
   return (double*)this->translatingDifference;
}

void Parameters::SetInertiaTensorOrigin(double x, double y, double z){
   if(this->inertiaTensorOrigin == NULL){
      this->inertiaTensorOrigin = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   }

   this->inertiaTensorOrigin[0] = x;
   this->inertiaTensorOrigin[1] = y;
   this->inertiaTensorOrigin[2] = z;

}

double* Parameters::GetInertiaTensorOrigin() const{
   return (double*)this->inertiaTensorOrigin;
}

void Parameters::SetRotatingOrigin(double x, double y, double z){
   if(this->rotatingOrigin == NULL){
      this->rotatingOrigin = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   }

   this->rotatingOrigin[0] = x;
   this->rotatingOrigin[1] = y;
   this->rotatingOrigin[2] = z;

}

double* Parameters::GetRotatingOrigin() const{
   return (double*)this->rotatingOrigin;
}

void Parameters::SetRotatingType(RotatingType rotatingType){
   this->rotatingType = rotatingType;
}

RotatingType Parameters::GetRotatingType() const{
   return this->rotatingType;
}

void Parameters::SetRotatingAxis(double x, double y, double z){
   this->rotatingAxis[0] = x;
   this->rotatingAxis[1] = y;
   this->rotatingAxis[2] = z;

}

double* Parameters::GetRotatingAxis() const{
   return (double*)this->rotatingAxis;
}

void Parameters::SetRotatingAngle(double rotatingAngle){
   this->rotatingAngle = rotatingAngle;
}

double Parameters::GetRotatingAngle() const{
   return this->rotatingAngle;
}

void Parameters::SetRotatingEularAngles(double alpha, double beta, double gamma){
   this->rotatingEularAngles.SetAlpha(alpha);
   this->rotatingEularAngles.SetBeta(beta);
   this->rotatingEularAngles.SetGamma(gamma);
}

EularAngle Parameters::GetRotatingEularAngles() const{
   return this->rotatingEularAngles;
}

int Parameters::GetActiveOccCIS() const{
   return this->activeOccCIS;
}
   
void Parameters::SetActiveOccCIS(int activeOccCIS){
   this->activeOccCIS = activeOccCIS;
}

int Parameters::GetActiveVirCIS() const{
   return this->activeVirCIS;
}

void Parameters::SetActiveVirCIS(int activeVirCIS){
   this->activeVirCIS = activeVirCIS;
}

int Parameters::GetNumberExcitedStatesCIS() const{
   return this->numberExcitedStatesCIS;
}

void Parameters::SetNumberExcitedStatesCIS(int nStates){
   this->numberExcitedStatesCIS = nStates;
}

bool Parameters::RequiresCIS() const{
   return this->requiresCIS;
}

void Parameters::SetRequiresCIS(bool requiresCIS){
   this->requiresCIS = requiresCIS;
}

vector<int>* Parameters::GetIndecesMOPlot() const{
   return this->indecesMOPlot;
}

void Parameters::AddIndexMOPlot(int moIndex){
   this->indecesMOPlot->push_back(moIndex);
}

bool Parameters::IsDavidsonCIS() const{
   return this->isDavidsonCIS;
}

void Parameters::SetIsDavidsonCIS(bool isDavidsonCIS){
   this->isDavidsonCIS = isDavidsonCIS;
}

int Parameters::GetMaxIterationsCIS() const{
   return this->maxIterationsCIS;
}

void Parameters::SetMaxIterationsCIS(int maxIterationsCIS){
   this->maxIterationsCIS = maxIterationsCIS;
}

int Parameters::GetMaxDimensionsCIS() const{
   return this->maxDimensionsCIS;
}

void Parameters::SetMaxDimensionsCIS(int maxDimensionsCIS){
   this->maxDimensionsCIS = maxDimensionsCIS;
}

double Parameters::GetNormToleranceCIS() const{
   return this->normToleranceCIS;
}

void Parameters::SetNormToleranceCIS(double normToleranceCIS){
   this->normToleranceCIS = normToleranceCIS;
}

bool Parameters::RequiresMD() const{
   return this->requiresMD;
}

void Parameters::SetRequiresMD(bool requiresMD){
   this->requiresMD = requiresMD;
}

int Parameters::GetElectronicStateIndexMD() const{
   return this->electronicStateIndexMD;
}

void Parameters::SetElectronicStateIndexMD(int electronicStateIndex){
   this->electronicStateIndexMD = electronicStateIndex;
}

int Parameters::GetTotalStepsMD() const{
   return this->totalStepsMD;
}

void Parameters::SetTotalStepsMD(int totalSteps){
   this->totalStepsMD = totalSteps;
}

double Parameters::GetTimeWidthMD() const{
   return this->timeWidthMD;
}

void Parameters::SetTimeWidthMD(double timeWidth){
   this->timeWidthMD = timeWidth;
}

}





