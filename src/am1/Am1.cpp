#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include"../base/Enums.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/atoms/Hatom.h"
#include"../base/atoms/Liatom.h"
#include"../base/atoms/Catom.h"
#include"../base/atoms/Natom.h"
#include"../base/atoms/Oatom.h"
#include"../base/atoms/Satom.h"
#include"../base/Molecule.h"
#include"../cndo/Cndo2.h"
#include"../zindo/ZindoS.h"
#include"../mndo/Mndo.h"
#include"Am1.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_am1{

/***
 *  Main Refferences for AM1 are [DZHS_1985, DY_1990]
 */
Am1::Am1() : MolDS_mndo::Mndo(){
   this->theory = AM1;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "Am1 created\n";
}

Am1::~Am1(){
   //cout << "Am1 deleted\n";
}

void Am1::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in am1::Am1::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in am1::Am1::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in am1::Am1::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in am1::Am1::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_am1::Am1::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_am1::Am1::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in am1::Am1::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in am1::Am1::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in am1::Am1::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: Atom A and B is same.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in am1::Am1::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in am1::Am1::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tAM1-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: AM1-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: AM1-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: AM1-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: AM1-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for AM1-CIS met convergence criterion(^^b\n\n\n";
}

void Am1::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double Am1::GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const{
   double energy = Mndo::GetDiatomCoreRepulsionEnergy(indexAtomA, indexAtomB);
   Atom* atomA = (*this->molecule->GetAtomVect())[indexAtomA];
   Atom* atomB = (*this->molecule->GetAtomVect())[indexAtomB];
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA = atomA->GetNddoAlpha(this->theory);
   double alphaB = atomB->GetNddoAlpha(this->theory);
   double temp = 0.0;
   for(int i=0; i<4; i++){
      double kA = atomA->GetNddoParameterK(this->theory, i);
      double lA = atomA->GetNddoParameterL(this->theory, i);
      double mA = atomA->GetNddoParameterM(this->theory, i);
      double kB = atomB->GetNddoParameterK(this->theory, i);
      double lB = atomB->GetNddoParameterL(this->theory, i);
      double mB = atomB->GetNddoParameterM(this->theory, i);
      temp += kA*exp(-lA*pow(distance-mA,2.0));
      temp += kB*exp(-lB*pow(distance-mB,2.0));
   }
   energy += atomA->GetCoreCharge()*atomB->GetCoreCharge()*temp/(distance/ang2AU);
   return energy;
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Am1::GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                  int atomBIndex, 
                                                  CartesianType axisA) const{
   double value = Mndo::GetDiatomCoreRepulsionFirstDerivative(atomAIndex,
                                                              atomBIndex,
                                                              axisA);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   double alphaA = atomA->GetNddoAlpha(this->theory);
   double alphaB = atomB->GetNddoAlpha(this->theory);
   double Rab = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double dRabDa = (atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA])/Rab;
   double temp1 = 0.0;
   double temp2 = 0.0;
   for(int i=0; i<4; i++){
      double kA = atomA->GetNddoParameterK(this->theory, i);
      double lA = atomA->GetNddoParameterL(this->theory, i);
      double mA = atomA->GetNddoParameterM(this->theory, i);
      double kB = atomB->GetNddoParameterK(this->theory, i);
      double lB = atomB->GetNddoParameterL(this->theory, i);
      double mB = atomB->GetNddoParameterM(this->theory, i);
      temp1 += kA*exp(-lA*pow(Rab-mA,2.0));
      temp1 += kB*exp(-lB*pow(Rab-mB,2.0));
      temp2 += -2.0*lA*(Rab-mA)*kA*exp(-lA*pow(Rab-mA,2.0));
      temp2 += -2.0*lB*(Rab-mB)*kB*exp(-lB*pow(Rab-mB,2.0));
   }
   value -= dRabDa
           *atomA->GetCoreCharge()
           *atomB->GetCoreCharge()
           *temp1
           /(pow(Rab,2.0)/ang2AU);
   value += dRabDa
           *atomA->GetCoreCharge()
           *atomB->GetCoreCharge()
           *temp2/(Rab/ang2AU);
   return value;
}

void Am1::CalcHFProperties(){
   MolDS_cndo::Cndo2::CalcHFProperties();
}

void Am1::OutputHFResults(double const* const* fockMatrix, 
                          double const* energiesMO, 
                          double const* atomicElectronPopulation, 
                          const Molecule& molecule) const{
   MolDS_cndo::Cndo2::OutputHFResults(fockMatrix, 
                                      energiesMO, 
                                      atomicElectronPopulation, 
                                      molecule);
}

}



