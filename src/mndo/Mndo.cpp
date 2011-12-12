#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include"../base/MolDSException.h"
#include"../mkl_wrapper/LapackWrapper.h"
#include"../base/Enums.h"
#include"../base/MallocerFreer.h"
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
#include"Mndo.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_mndo{

/***
 *  Main Refferences for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
Mndo::Mndo() : MolDS_zindo::ZindoS(){
   this->theory = MNDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   this->heatsFormation = 0.0;
   this->zMatrixForceElecStatesNum = 0;
   this->etaMatrixForceElecStatesNum = 0;
   this->zMatrixForce = NULL;
   //cout << "Mndo created\n";
}

Mndo::~Mndo(){
   if(this->twoElecTwoCore != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix6d(
                                    &this->twoElecTwoCore, 
                                    this->molecule->GetAtomVect()->size(),
                                    this->molecule->GetAtomVect()->size(),
                                    dxy,
                                    dxy,
                                    dxy);
      //cout << "twoElecTwoCore deleted\n";
   }
   if(this->zMatrixForce != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix3d(&this->zMatrixForce, 
                                                       this->zMatrixForceElecStatesNum,
                                                       this->molecule->GetTotalNumberAOs());
      //cout << "zMatrixForce deleted\n";
   }
   if(this->etaMatrixForce != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix3d(&this->etaMatrixForce, 
                                                       this->etaMatrixForceElecStatesNum,
                                                       this->molecule->GetTotalNumberAOs());
      //cout << "etaMatrixForce deleted\n";
   }
}

void Mndo::SetMolecule(Molecule* molecule){
   Cndo2::SetMolecule(molecule);
   this->twoElecTwoCore = MallocerFreer::GetInstance()->MallocDoubleMatrix6d(
                                                         molecule->GetAtomVect()->size(),
                                                         molecule->GetAtomVect()->size(),
                                                         dxy, dxy, dxy, dxy);
}
void Mndo::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in mndo::Mndo::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in mndo::Mndo::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in mndo::Mndo::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in mndo::Mndo::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_mndo::Mndo::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_mndo::Mndo::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in mndo::Mndo::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in mndo::Mndo::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageMultipoleA = "Multipole A is: ";
   this->errorMessageMultipoleB = "Multipole B is: ";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in mndo::Mndo::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix 
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreDiatomic: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreDiatomic: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreDiatomicFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreDiatomicFirstDerivatives: Atom A and B is same.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in mndo::Mndo::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in mndo::Mndo::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageCalcZMatrixForceEtaNull 
      = "Error in mndo::Mndo::CalcZMatrixForce: Nndo::etaMatrixForce is NULL. Call Mndo::CalcEtaMatrixForce before calling Mndo::CalcZMatrixForce.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tMNDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: MNDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: MNDO/S-SCF  **********\n\n\n";
   this->messageHeatsFormation = "\tHeats of formation:\n";
   this->messageHeatsFormationTitle = "\t\t| [a.u.] | [Kcal/mol] | \n";
   this->messageStartCIS = "**********  START: MNDO/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: MNDO/S-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for MNDO-CIS met convergence criterion(^^b\n\n\n";
}

void Mndo::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double Mndo::CalcDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB){
   Atom* atomA = (*this->molecule->GetAtomVect())[indexAtomA];
   Atom* atomB = (*this->molecule->GetAtomVect())[indexAtomB];
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA = atomA->GetNddoAlpha(this->theory);
   double alphaB = atomB->GetNddoAlpha(this->theory);
   double temp = 0.0;
   if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                    atomB->GetAtomType() == O)  ){
      temp = 1.0 + (distance/ang2AU)*exp(-alphaB*distance) + exp(-alphaA*distance);
   }
   else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                         atomA->GetAtomType() == O)  ){
      temp = 1.0 + (distance/ang2AU)*exp(-alphaA*distance) + exp(-alphaB*distance);
   }
   else{
      temp = 1.0 + exp(-alphaA*distance) + exp(-alphaB*distance);
   }
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   return  atomA->GetCoreCharge()*atomB->GetCoreCharge()*twoElecInt*temp; 
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Mndo::GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                   int atomBIndex, 
                                                   CartesianType axisA){
   double value =0.0;
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   double alphaA = atomA->GetNddoAlpha(this->theory);
   double alphaB = atomB->GetNddoAlpha(this->theory);
   double Rab = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double dRabDa = (atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA])/Rab;
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   double twoElecIntFirstDeri = this->GetNddoRepulsionIntegralFirstDerivative(
                                      atomA, s, s, atomB, s, s, axisA);
   double temp1 = 0.0;
   if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                    atomB->GetAtomType() == O)  ){
      temp1 = 1.0 + (Rab/ang2AU)*exp(-alphaB*Rab) + exp(-alphaA*Rab);
   }
   else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                         atomA->GetAtomType() == O)  ){
      temp1 = 1.0 + (Rab/ang2AU)*exp(-alphaA*Rab) + exp(-alphaB*Rab);
   }
   else{
      temp1 = 1.0 + exp(-alphaA*Rab) + exp(-alphaB*Rab);
   }

   double temp2 = 0.0;
   if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                    atomB->GetAtomType() == O)  ){
      temp2 = (1.0/ang2AU)*exp(-alphaB*Rab) 
             -alphaB*(Rab/ang2AU)*exp(-alphaB*Rab) 
             -alphaA*exp(-alphaA*Rab);
   }
   else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                         atomA->GetAtomType() == O)  ){
      temp2 = (1.0/ang2AU)*exp(-alphaA*Rab)
             -alphaA*(Rab/ang2AU)*exp(-alphaA*Rab) 
             -alphaB*exp(-alphaB*Rab);
   }
   else{
      temp2 = -alphaA*exp(-alphaA*Rab) 
              -alphaB*exp(-alphaB*Rab);
   }
   temp2 *= dRabDa;
   value = atomA->GetCoreCharge()*atomB->GetCoreCharge()
          *(twoElecIntFirstDeri*temp1 + twoElecInt*temp2); 
   return value;
}

void Mndo::CalcHeatsFormation(double* heatsFormation, Molecule* molecule){
   int groundState = 0;
   *heatsFormation = this->GetElectronicEnergy(groundState);
   for(int A=0; A<molecule->GetAtomVect()->size(); A++){
      Atom* atom = (*molecule->GetAtomVect())[A];
      *heatsFormation -= atom->GetMndoElecEnergyAtom();
      *heatsFormation += atom->GetMndoHeatsFormAtom();
   }
}

void Mndo::CalcHFProperties(){
   MolDS_cndo::Cndo2::CalcHFProperties();
   this->CalcHeatsFormation(&this->heatsFormation, this->molecule);
}

void Mndo::OutputHFResults(double** fockMatrix, 
                           double* energiesMO, 
                           double* atomicElectronPopulation, 
                           Molecule* molecule){
   MolDS_cndo::Cndo2::OutputHFResults(fockMatrix, 
                                      energiesMO, 
                                      atomicElectronPopulation, 
                                      molecule);
   // output heats of formation
   cout << this->messageHeatsFormation;
   cout << this->messageHeatsFormationTitle;
   printf("\t\t%e\t%e\n\n",this->heatsFormation,
                           this->heatsFormation/Parameters::GetInstance()->
                                                            GetKcalMolin2AU());
}

double Mndo::GetFockDiagElement(Atom* atomA, 
                                int atomAIndex, 
                                int mu, 
                                Molecule* molecule, 
                                double** gammaAB,
                                double** orbitalElectronPopulation, 
                                double* atomicElectronPopulation,
                                double****** twoElecTwoCore, 
                                bool isGuess){
   double value=0.0;
   int firstAOIndexA = atomA->GetFirstAOIndex();
   mu -= firstAOIndexA;
   value = atomA->GetCoreIntegral(atomA->GetValence()[mu], isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      OrbitalType orbitalMu = atomA->GetValence()[mu];
      for(int nu=0; nu<atomA->GetValence().size(); nu++){
         OrbitalType orbitalNu = atomA->GetValence()[nu];
         double coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);
         double exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
         temp += orbitalElectronPopulation[nu+firstAOIndexA]
                                          [nu+firstAOIndexA]
                *(coulomb - 0.5*exchange);
      }
      value += temp;

      temp = 0.0;
      for(int B=0; B<molecule->GetAtomVect()->size(); B++){
         if(B != atomAIndex){
            Atom* atomB = (*molecule->GetAtomVect())[B];
            int firstAOIndexB = atomB->GetFirstAOIndex();
            for(int lambda=0; lambda<atomB->GetValence().size(); lambda++){
               for(int sigma=0; sigma<atomB->GetValence().size(); sigma++){
                  temp += orbitalElectronPopulation[lambda+firstAOIndexB]
                                                   [sigma+firstAOIndexB]
                         *twoElecTwoCore[atomAIndex][B][mu][mu][lambda][sigma];
               }
            }
            temp += this->GetElectronCoreAttraction(atomAIndex, 
                                                    B, 
                                                    mu, 
                                                    mu, 
                                                    twoElecTwoCore);
         }
      }
      value += temp;
   }
   return value;
}

double Mndo::GetFockOffDiagElement(Atom* atomA, 
                                   Atom* atomB, 
                                   int atomAIndex, 
                                   int atomBIndex, 
                                   int mu, 
                                   int nu, 
                                   Molecule* molecule, 
                                   double** gammaAB, 
                                   double** overlap,
                                   double** orbitalElectronPopulation, 
                                   double****** twoElecTwoCore, 
                                   bool isGuess){
   double value = 0.0;
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   mu -= firstAOIndexA;
   nu -= firstAOIndexB;
   OrbitalType orbitalMu = atomA->GetValence()[mu];
   OrbitalType orbitalNu = atomB->GetValence()[nu];
   double bondParameter = 0.5*(atomA->GetBondingParameter(this->theory, orbitalMu) 
                              +atomB->GetBondingParameter(this->theory, orbitalNu)); 
   if(isGuess){
      value = bondParameter*overlap[mu+firstAOIndexA][nu+firstAOIndexB];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      double temp = 0.0;
      if(atomAIndex == atomBIndex){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         temp = (1.5*exchange - 0.5*coulomb)
               *orbitalElectronPopulation[mu+firstAOIndexA][nu+firstAOIndexB];
         for(int BB=0; BB<molecule->GetAtomVect()->size(); BB++){
            if(BB != atomAIndex){
               Atom* atomBB = (*molecule->GetAtomVect())[BB];
               int firstAOIndexBB = atomBB->GetFirstAOIndex();
               for(int lambda=0; lambda<atomBB->GetValence().size(); lambda++){
                  for(int sigma=0; sigma<atomBB->GetValence().size(); sigma++){
                     temp += orbitalElectronPopulation[lambda+firstAOIndexBB]
                                                      [sigma+firstAOIndexBB]
                            *twoElecTwoCore[atomAIndex][BB][mu][nu][lambda][sigma];
                  }
               }
               temp += this->GetElectronCoreAttraction(atomAIndex, 
                                                       BB, 
                                                       mu, 
                                                       nu, 
                                                       twoElecTwoCore);
            }
         }
      }
      else{
         temp = bondParameter*overlap[mu+firstAOIndexA][nu+firstAOIndexB];
         for(int sigma=0; sigma<atomA->GetValence().size(); sigma++){
            for(int lambda=0; lambda<atomB->GetValence().size(); lambda++){
               temp -= 0.5*orbitalElectronPopulation[lambda+firstAOIndexB]
                                                    [sigma+firstAOIndexA]
                      *twoElecTwoCore[atomAIndex][atomBIndex][mu][sigma][nu][lambda];
            }
         }
      }
      value += temp;
   }
   return value;
}

// NDDO Coulomb Interaction
double Mndo::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){
   double value=0.0;
   if( orbital1 == s && orbital2 == s){ 
      value = atom->GetNddoGss(this->theory);
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom->GetNddoGsp(this->theory);
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = this->GetCoulombInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetNddoGpp(this->theory);
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetNddoGpp2(this->theory);
   }   
   else{
      stringstream ss;
      ss << this->errorMessageCoulombInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   return value;
}

// NDDO Exchange Interaction
double Mndo::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){
   double value=0.0;
   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetNddoHsp(this->theory);
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = this->GetExchangeInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetNddoHpp(this->theory);
   }
   else{
      stringstream ss;
      ss << this->errorMessageExchangeInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   return value;
}

// electron in atom A (mu and nu) and core (atom B) attraction. 
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttraction(int atomAIndex, 
                                       int atomBIndex, 
                                       int mu, 
                                       int nu, 
                                       double****** twoElecTwoCore){
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   return -1.0*atomB->GetCoreCharge()*twoElecTwoCore[atomAIndex][atomBIndex][mu][nu][s][s];
}

// First derivative of electron in atom A (mu and nu) and core (atom B) attraction. 
// This derivative is related to the coordinate of atomA.
// Note that towtwoElecTwoCoreFirstDerivative is dioatomic one.
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttractionFirstDerivative(int atomAIndex, 
                                                      int atomBIndex, 
                                                      int mu, 
                                                      int nu, 
                                                      double***** twoElecTwoCoreFirstDerivative,
                                                      CartesianType axisA){
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   double value = -1.0*atomB->GetCoreCharge()
                  *twoElecTwoCoreFirstDerivative[mu][nu][s][s][axisA];
   return value;
}

void Mndo::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                              Atom* atomA, 
                                              Atom* atomB){
   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);
}

void Mndo::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                                          double** diatomicOverlapDeri, 
                                          Atom* atomA, 
                                          Atom* atomB){
   MolDS_cndo::Cndo2::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                        diatomicOverlapDeri,atomA, atomB);
}

// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double Mndo::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                         Molecule* molecule, 
                                         double** fockMatrix, 
                                         double** gammaAB){
   double value = 0.0;
   for(int A=0; A<molecule->GetAtomVect()->size(); A++){
      Atom* atomA = (*molecule->GetAtomVect())[A];
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValence().size();

      for(int B=A; B<molecule->GetAtomVect()->size(); B++){
         Atom* atomB = (*molecule->GetAtomVect())[B];
         int firstAOIndexB = atomB->GetFirstAOIndex();
         int numberAOsB = atomB->GetValence().size();

         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
                        OrbitalType orbitalSigma = atomB->GetValence()[sigma-firstAOIndexB];
                        gamma = this->twoElecTwoCore[A]
                                                    [B]
                                                    [mu-firstAOIndexA]
                                                    [nu-firstAOIndexA]
                                                    [lambda-firstAOIndexB]
                                                    [sigma-firstAOIndexB];

                        value += gamma*fockMatrix[moI][mu]
                                      *fockMatrix[moJ][nu]
                                      *fockMatrix[moK][lambda]
                                      *fockMatrix[moL][sigma];
                        value += gamma*fockMatrix[moI][lambda]
                                      *fockMatrix[moJ][sigma]
                                      *fockMatrix[moK][mu]
                                      *fockMatrix[moL][nu];
                        if(lambda != sigma){
                           value += gamma*fockMatrix[moI][mu]
                                         *fockMatrix[moJ][nu]
                                         *fockMatrix[moK][sigma]
                                         *fockMatrix[moL][lambda];
                           value += gamma*fockMatrix[moI][sigma]
                                         *fockMatrix[moJ][lambda]
                                         *fockMatrix[moK][mu]
                                         *fockMatrix[moL][nu];
                        }
                        if(mu != nu){
                           value += gamma*fockMatrix[moI][nu]
                                         *fockMatrix[moJ][mu]
                                         *fockMatrix[moK][lambda]
                                         *fockMatrix[moL][sigma];
                           value += gamma*fockMatrix[moI][lambda]
                                         *fockMatrix[moJ][sigma]
                                         *fockMatrix[moK][nu]
                                         *fockMatrix[moL][mu];
                        }
                        if(mu != nu && lambda != sigma){
                           value += gamma*fockMatrix[moI][nu]
                                         *fockMatrix[moJ][mu]
                                         *fockMatrix[moK][sigma]
                                         *fockMatrix[moL][lambda];
                           value += gamma*fockMatrix[moI][sigma]
                                         *fockMatrix[moJ][lambda]
                                         *fockMatrix[moK][nu]
                                         *fockMatrix[moL][mu];
                        }
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                           OrbitalType orbitalLambda = atomB->GetValence()[lambda-firstAOIndexB];
                           gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                           OrbitalType orbitalNu = atomA->GetValence()[nu-firstAOIndexA];
                           gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        }
                        else{
                           gamma = 0.0;
                        }
                        value += gamma*fockMatrix[moI][mu]
                                      *fockMatrix[moJ][nu]
                                      *fockMatrix[moK][lambda]
                                      *fockMatrix[moL][sigma];
                     }  
                  }
               }
            }
         }
      }
   }
   return value;
}

// right-upper part is only calculated by this method.
void Mndo::CalcCISMatrix(double** matrixCIS, int numberOcc, int numberVir){
   cout << this->messageStartCalcCISMatrix;
   double ompStartTime = omp_get_wtime();

   #pragma omp parallel for schedule(auto)
   for(int k=0; k<numberOcc*numberVir; k++){
      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->molecule->GetTotalNumberValenceElectrons()/2 - (k/numberVir) -1;
      int moA = this->molecule->GetTotalNumberValenceElectrons()/2 + (k%numberVir);

      for(int l=k; l<numberOcc*numberVir; l++){
         // single excitation from J-th (occupied)MO to B-th (virtual)MO
         int moJ = this->molecule->GetTotalNumberValenceElectrons()/2 - (l/numberVir) -1;
         int moB = this->molecule->GetTotalNumberValenceElectrons()/2 + (l%numberVir);
         double value=0.0;
          
         // Fast algorith, but this is not easy to read. Slow algorithm is alos written below.
         for(int A=0; A<molecule->GetAtomVect()->size(); A++){
            Atom* atomA = (*molecule->GetAtomVect())[A];
            int firstAOIndexA = atomA->GetFirstAOIndex();
            int numberAOsA = atomA->GetValence().size();

            for(int B=A; B<molecule->GetAtomVect()->size(); B++){
               Atom* atomB = (*molecule->GetAtomVect())[B];
               int firstAOIndexB = atomB->GetFirstAOIndex();
               int numberAOsB = atomB->GetValence().size();

               double gamma = 0.0;
               if(A!=B){
                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                        for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                           for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
                              OrbitalType orbitalSigma = atomB->GetValence()[sigma-firstAOIndexB];
                              gamma = this->twoElecTwoCore[A]
                                                          [B]
                                                          [mu-firstAOIndexA]
                                                          [nu-firstAOIndexA]
                                                          [lambda-firstAOIndexB]
                                                          [sigma-firstAOIndexB];

                              value += 2.0*gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][lambda]
                                                *fockMatrix[moB][sigma];
                              value += 2.0*gamma*fockMatrix[moA][lambda]
                                                *fockMatrix[moI][sigma]
                                                *fockMatrix[moJ][mu]
                                                *fockMatrix[moB][nu];
                              value -= gamma*fockMatrix[moA][mu]
                                            *fockMatrix[moB][nu]
                                            *fockMatrix[moI][lambda]
                                            *fockMatrix[moJ][sigma];
                              value -= gamma*fockMatrix[moA][lambda]
                                            *fockMatrix[moB][sigma]
                                            *fockMatrix[moI][mu]
                                            *fockMatrix[moJ][nu];
                              if(lambda != sigma){
                                 value += 2.0*gamma*fockMatrix[moA][mu]
                                                   *fockMatrix[moI][nu]
                                                   *fockMatrix[moJ][sigma]
                                                   *fockMatrix[moB][lambda];
                                 value += 2.0*gamma*fockMatrix[moA][sigma]
                                                   *fockMatrix[moI][lambda]
                                                   *fockMatrix[moJ][mu]
                                                   *fockMatrix[moB][nu];
                                 value -= gamma*fockMatrix[moA][mu]
                                               *fockMatrix[moB][nu]
                                               *fockMatrix[moI][sigma]
                                               *fockMatrix[moJ][lambda];
                                 value -= gamma*fockMatrix[moA][sigma]
                                               *fockMatrix[moB][lambda]
                                               *fockMatrix[moI][mu]
                                               *fockMatrix[moJ][nu];
                              }
                              if(mu != nu){
                                 value += 2.0*gamma*fockMatrix[moA][nu]
                                                   *fockMatrix[moI][mu]
                                                   *fockMatrix[moJ][lambda]
                                                   *fockMatrix[moB][sigma];
                                 value += 2.0*gamma*fockMatrix[moA][lambda]
                                                   *fockMatrix[moI][sigma]
                                                   *fockMatrix[moJ][nu]
                                                   *fockMatrix[moB][mu];
                                 value -= gamma*fockMatrix[moA][nu]
                                               *fockMatrix[moB][mu]
                                               *fockMatrix[moI][lambda]
                                               *fockMatrix[moJ][sigma];
                                 value -= gamma*fockMatrix[moA][lambda]
                                               *fockMatrix[moB][sigma]
                                               *fockMatrix[moI][nu]
                                               *fockMatrix[moJ][mu];
                              }
                              if(mu != nu && lambda != sigma){
                                 value += 2.0*gamma*fockMatrix[moA][nu]
                                                   *fockMatrix[moI][mu]
                                                   *fockMatrix[moJ][sigma]
                                                   *fockMatrix[moB][lambda];
                                 value += 2.0*gamma*fockMatrix[moA][sigma]
                                                   *fockMatrix[moI][lambda]
                                                   *fockMatrix[moJ][nu]
                                                   *fockMatrix[moB][mu];
                                 value -= gamma*fockMatrix[moA][nu]
                                               *fockMatrix[moB][mu]
                                               *fockMatrix[moI][sigma]
                                               *fockMatrix[moJ][lambda];
                                 value -= gamma*fockMatrix[moA][sigma]
                                               *fockMatrix[moB][lambda]
                                               *fockMatrix[moI][nu]
                                               *fockMatrix[moJ][mu];
                              }
                           }
                        }
                     }
                  }
               }
               else{
                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                        for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                           for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
                              if(mu==nu && lambda==sigma){
                                 OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                                 OrbitalType orbitalLambda = atomB->GetValence()[lambda-firstAOIndexB];
                                 gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                              }
                              else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                                 OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                                 OrbitalType orbitalNu = atomA->GetValence()[nu-firstAOIndexA];
                                 gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                              }
                              else{
                                 gamma = 0.0;
                              }
                              value += 2.0*gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][lambda]
                                                *fockMatrix[moB][sigma];
                              value -= gamma*fockMatrix[moA][mu]
                                            *fockMatrix[moB][nu]
                                            *fockMatrix[moI][lambda]
                                            *fockMatrix[moJ][sigma];
                           }  
                        }
                     }
                  }
               }
            }
         }
         // End of the fast algorith.
         
         /* 
         // Slow algorith, but this is easy to read. Fast altorithm is also written above.
         value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                                       this->molecule, this->fockMatrix, NULL)
                    -this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                       this->molecule, this->fockMatrix, NULL);
         // End of the slow algorith.
         */
         // Diagonal term
         if(k==l){
            value += this->energiesMO[moA] - this->energiesMO[moI];
         }
         matrixCIS[k][l] = value;
      }
   }
   double ompEndTime = omp_get_wtime();
   cout << this->messageOmpElapsedTimeCalcCISMarix;
   cout << ompEndTime - ompStartTime;
   cout << this->messageUnitSec << endl;
   cout << this->messageDoneCalcCISMatrix;
}

void Mndo::CheckZMatrixForce(vector<int> elecStates){
   // malloc or initialize Z matrix
   if(this->zMatrixForce == NULL){
      this->zMatrixForce = MallocerFreer::GetInstance()->
                           MallocDoubleMatrix3d(elecStates.size(),
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
      this->zMatrixForceElecStatesNum = elecStates.size();
   }
   else{
      MallocerFreer::GetInstance()->
      InitializeDoubleMatrix3d(this->zMatrixForce,
                               elecStates.size(),
                               this->molecule->GetTotalNumberAOs(), 
                               this->molecule->GetTotalNumberAOs());
   }
}

void Mndo::CheckEtaMatrixForce(vector<int> elecStates){
   // malloc or initialize eta matrix
   if(this->etaMatrixForce == NULL){
      this->etaMatrixForce = MallocerFreer::GetInstance()->
                             MallocDoubleMatrix3d(elecStates.size(),
                                                  this->molecule->GetTotalNumberAOs(), 
                                                  this->molecule->GetTotalNumberAOs());
      this->etaMatrixForceElecStatesNum = elecStates.size();
   }
   else{
      MallocerFreer::GetInstance()->
      InitializeDoubleMatrix3d(this->etaMatrixForce,
                               elecStates.size(),
                               this->molecule->GetTotalNumberAOs(), 
                               this->molecule->GetTotalNumberAOs());
   }
}

// see variable Q-vector in [PT_1996, PT_1997]
void Mndo::CalcActiveSetVariablesQ(vector<MoIndexPair>* nonRedundantQIndeces, 
                                   vector<MoIndexPair>* redundantQIndeces){
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   for(int moI=0; moI<numberOcc; moI++){
      bool isMoICIMO = numberOcc-numberActiveOcc<=moI ? true : false;
      for(int moJ=numberOcc; moJ<numberAOs; moJ++){
         bool isMoJCIMO = moJ<numberOcc+numberActiveVir ? true : false;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         nonRedundantQIndeces->push_back(moIndexPair);
      }
   }
   for(int moI=numberOcc-numberActiveOcc; moI<numberOcc; moI++){
      bool isMoICIMO = true;
      for(int moJ=moI; moJ<numberOcc; moJ++){
         bool isMoJCIMO = true;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         redundantQIndeces->push_back(moIndexPair);
      }
   }
   for(int moI=numberOcc; moI<numberOcc+numberActiveVir; moI++){
      bool isMoICIMO = true;
      for(int moJ=moI; moJ<numberOcc+numberActiveVir; moJ++){
         bool isMoJCIMO = true;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         redundantQIndeces->push_back(moIndexPair);
      }
   }
}

void Mndo::MallocTempMatrixForZMatrix(double** delta,
                                      double** q,
                                      double*** kNR,
                                      double*** kRDag,
                                      double** y,
                                      double*** transposedFockMatrix,
                                      double*** xiOcc,
                                      double*** xiVir,
                                      int sizeQNR,
                                      int sizeQR){
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberActiveMO = numberActiveOcc + numberActiveVir;
   int numberAOs = this->molecule->GetTotalNumberAOs();
   *delta = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(
                                          numberActiveMO);
   *q = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(
                                      sizeQNR+sizeQR);
   *kNR = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                        sizeQNR,
                                        sizeQNR);
   *kRDag = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                          sizeQNR,
                                          sizeQR);
   *y = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(
                                      sizeQNR);
   *transposedFockMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                                         numberAOs,
                                                         numberAOs);
   *xiOcc = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                          numberActiveOcc,
                                          numberAOs);
   *xiVir = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                          numberActiveVir,
                                          numberAOs);
}

void Mndo::FreeTempMatrixForZMatrix(double** delta,
                                    double** q,
                                    double*** kNR,
                                    double*** kRDag,
                                    double** y,
                                    double*** transposedFockMatrix,
                                    double*** xiOcc,
                                    double*** xiVir,
                                    int sizeQNR,
                                    int sizeQR){
   if(*delta != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(delta);
      //cout << "delta  deleted" << endl;
   }
   if(*q != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(q);
      //cout << "q  deleted" << endl;
   }
   if(*kNR != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(kNR, sizeQNR);
      //cout << "kNR  deleted" << endl;
   }
   if(*kRDag != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(kRDag, sizeQNR);
      //cout << "kRDag  deleted" << endl;
   }
   if(*y != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(y);
      //cout << "y  deleted" << endl;
   }
   if(*transposedFockMatrix != NULL){
      int numberAOs = this->molecule->GetTotalNumberAOs();
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(transposedFockMatrix,
                                                       numberAOs);
      //cout << "transposedFockMatrix  deleted" << endl;
   }
   if(*xiOcc != NULL){
      int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(xiOcc,
                                                       numberActiveOcc);
      //cout << "xiOcc deleted" << endl;
   }
   if(*xiVir != NULL){
      int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(xiVir,
                                                       numberActiveVir);
      //cout << "xiVir deleted" << endl;
   }
}

// \epsilon_{r}^{kl} in (1) in [PT_1997].
// k and l are index of CIS matrix.
double Mndo::GetCISCoefficientMOEnergy(int k, int l, int r, int numberActiveVir){
   double value=0.0;
   if(k==l){
      int moI = this->molecule->GetTotalNumberValenceElectrons()/2 - (k/numberActiveVir) -1;
      int moA = this->molecule->GetTotalNumberValenceElectrons()/2 + (k%numberActiveVir);
      if(r==moI){
         // r is index of occupied MO.
         value = -1.0;
      }
      else if(r==moA){
         // r is index of virtual MO.
         value = 1.0;
      }
   }
   return value;
}

// \f_{pqrs}^{lm} in (1) in [PT_1997].
// k and l are index of CIS matrix.
double Mndo::GetCISCoefficientTwoElecIntegral(int k, 
                                              int l, 
                                              int p, 
                                              int q, 
                                              int r, 
                                              int s,
                                              int numberActiveVir){
   double value=0.0;
   // single excitation from I-th (occupied)MO to A-th (virtual)MO
   int moI = this->molecule->GetTotalNumberValenceElectrons()/2 - (k/numberActiveVir) -1;
   int moA = this->molecule->GetTotalNumberValenceElectrons()/2 + (k%numberActiveVir);
   // single excitation from J-th (occupied)MO to B-th (virtual)MO
   int moJ = this->molecule->GetTotalNumberValenceElectrons()/2 - (l/numberActiveVir) -1;
   int moB = this->molecule->GetTotalNumberValenceElectrons()/2 + (l%numberActiveVir);
   if(p==moI && q==moA && r==moJ && s==moB ){
      value = 2.0;
   }
   else if(p==moI && q==moJ && r==moA && s==moB ){
      value = -1.0;
   }
   return value;
}

// see (40) in [PT_1996]
double Mndo::GetGammaNRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   if(moI==moK && moJ==moL){
      int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
      double nI = moI<numberOcc ? 2.0 : 0.0;
      double nJ = moJ<numberOcc ? 2.0 : 0.0;
      value = (this->energiesMO[moJ]-this->energiesMO[moI])/(nJ-nI);
   }
   return value;
}

// see (41) & (42) in [PT_1996]
double Mndo::GetGammaRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = moI==moJ ? 1.0 : this->energiesMO[moJ]-this->energiesMO[moI];
   }
   return value;
}

// see (43) in [PT_1996]
double Mndo::GetNNRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   if(moI==moK && moJ==moL){
      int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
      double nI = moI<numberOcc ? 2.0 : 0.0;
      double nJ = moJ<numberOcc ? 2.0 : 0.0;
      value = (nJ-nI);
   }
   return value;
}

// see (44) in [PT_1996]
double Mndo::GetNRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = 1.0;
   }
   return value;
}

// see (44) in [PT_1996]
double Mndo::GetKNRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;
   if(nI!=nJ && nK!=nL){
      value = 4.0*this->GetMolecularIntegralElement(moI, moJ, moK, moL, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL)
             -1.0*this->GetMolecularIntegralElement(moI, moK, moJ, moL, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL)
             -1.0*this->GetMolecularIntegralElement(moI, moL, moJ, moK, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL);
   }
   return 0.5*value;
}

// Dager of (45) in [PT_1996]. Note taht the (45) is real number.
double Mndo::GetKRDagerElement(int moI, int moJ, int moK, int moL){
   return this->GetKRElement(moK, moL, moI, moJ);
}

// see (45) in [PT_1996]
double Mndo::GetKRElement(int moI, int moJ, int moK, int moL){
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;
   if(nI==nJ && nK!=nL){
      value = 4.0*this->GetMolecularIntegralElement(moI, moJ, moK, moL, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL)
             -1.0*this->GetMolecularIntegralElement(moI, moK, moJ, moL, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL)
             -1.0*this->GetMolecularIntegralElement(moI, moL, moJ, moK, 
                                                    this->molecule, 
                                                    this->fockMatrix, NULL);
   }
   return 0.5*value;
}

// see (9) in [PT_1997]
void Mndo::CalcDeltaVector(double* delta, int exciteState){
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberActiveMO = numberActiveOcc + numberActiveVir;
   MallocerFreer::GetInstance()->InitializeDoubleMatrix1d(delta, numberActiveMO);
   #pragma omp parallel for schedule(auto)
   for(int r=0; r<numberActiveMO; r++){
      double value = 0.0;
      if(r<numberActiveOcc){
         // r is active occupied MO
         int rr=numberActiveOcc-(r+1);
         for(int a=0; a<numberActiveVir; a++){
            int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(rr,a);
            value -= pow(this->matrixCIS[exciteState][slaterDeterminantIndex],2.0);
         }
      }
      else{
         // r is active virtual MO
         int rr=r-numberActiveOcc;
         for(int i=0; i<numberActiveOcc; i++){
            int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,rr);
            value += pow(this->matrixCIS[exciteState][slaterDeterminantIndex],2.0);
         }
      }
      delta[r] = value;
   }
}

// see (18) in [PT_1977]
double Mndo::GetSmallQElement(int moI, 
                              int moP, 
                              double** xiOcc, 
                              double** xiVir, 
                              double** eta){
   double value = 0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   bool isMoPOcc = moP<numberOcc ? true : false;
   
   for(int A=0; A<molecule->GetAtomVect()->size(); A++){
      Atom* atomA = (*molecule->GetAtomVect())[A];
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValence().size();

      for(int B=A; B<molecule->GetAtomVect()->size(); B++){
         Atom* atomB = (*molecule->GetAtomVect())[B];
         int firstAOIndexB = atomB->GetFirstAOIndex();
         int numberAOsB = atomB->GetValence().size();

         if(A!=B){
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
                        double twoElecInt = 0.0;
                        twoElecInt = this->twoElecTwoCore[A]
                                                         [B]
                                                         [mu-firstAOIndexA]
                                                         [nu-firstAOIndexA]
                                                         [lambda-firstAOIndexB]
                                                         [sigma-firstAOIndexB];
                        double temp = 0.0;
                        if(isMoPOcc){
                           int p = numberOcc - (moP+1);
                           temp = 4.0*xiOcc[p][nu]*eta[lambda][sigma]
                                 -1.0*xiOcc[p][lambda]*eta[nu][sigma]
                                 -1.0*xiOcc[p][sigma]*eta[nu][lambda];
                           value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                           temp = 4.0*xiOcc[p][sigma]*eta[mu][nu]
                                 -1.0*xiOcc[p][mu]*eta[sigma][nu]
                                 -1.0*xiOcc[p][nu]*eta[sigma][mu];
                           value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                        }
                        else{
                           int p = moP - numberOcc;
                           temp = 4.0*xiVir[p][nu]*eta[lambda][sigma]
                                 -1.0*xiVir[p][lambda]*eta[sigma][nu]
                                 -1.0*xiVir[p][sigma]*eta[lambda][nu];
                           value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                           temp = 4.0*xiVir[p][sigma]*eta[mu][nu]
                                 -1.0*xiVir[p][mu]*eta[nu][sigma]
                                 -1.0*xiVir[p][nu]*eta[mu][sigma];
                           value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                        }
                         
                        if(lambda!=sigma){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][nu]*eta[sigma][lambda]
                                    -1.0*xiOcc[p][sigma]*eta[nu][lambda]
                                    -1.0*xiOcc[p][lambda]*eta[nu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                              temp = 4.0*xiOcc[p][lambda]*eta[mu][nu]
                                    -1.0*xiOcc[p][mu]*eta[lambda][nu]
                                    -1.0*xiOcc[p][nu]*eta[lambda][mu];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][nu]*eta[sigma][lambda]
                                    -1.0*xiVir[p][sigma]*eta[lambda][nu]
                                    -1.0*xiVir[p][lambda]*eta[sigma][nu];
                              value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                              temp = 4.0*xiVir[p][lambda]*eta[mu][nu]
                                    -1.0*xiVir[p][mu]*eta[nu][lambda]
                                    -1.0*xiVir[p][nu]*eta[mu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                        }
                        
                        if(mu!=nu){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][mu]*eta[lambda][sigma]
                                    -1.0*xiOcc[p][lambda]*eta[mu][sigma]
                                    -1.0*xiOcc[p][sigma]*eta[mu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiOcc[p][sigma]*eta[nu][mu]
                                    -1.0*xiOcc[p][nu]*eta[sigma][mu]
                                    -1.0*xiOcc[p][mu]*eta[sigma][nu];
                              value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][mu]*eta[lambda][sigma]
                                    -1.0*xiVir[p][lambda]*eta[sigma][mu]
                                    -1.0*xiVir[p][sigma]*eta[lambda][mu];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiVir[p][sigma]*eta[nu][mu]
                                    -1.0*xiVir[p][nu]*eta[mu][sigma]
                                    -1.0*xiVir[p][mu]*eta[nu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                           }
                        }

                        if(mu!=nu && lambda!=sigma){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][mu]*eta[sigma][lambda]
                                    -1.0*xiOcc[p][sigma]*eta[mu][lambda]
                                    -1.0*xiOcc[p][lambda]*eta[mu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiOcc[p][lambda]*eta[nu][mu]
                                    -1.0*xiOcc[p][nu]*eta[lambda][mu]
                                    -1.0*xiOcc[p][mu]*eta[lambda][nu];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][mu]*eta[sigma][lambda]
                                    -1.0*xiVir[p][sigma]*eta[lambda][mu]
                                    -1.0*xiVir[p][lambda]*eta[sigma][mu];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiVir[p][lambda]*eta[nu][mu]
                                    -1.0*xiVir[p][nu]*eta[mu][lambda]
                                    -1.0*xiVir[p][mu]*eta[nu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                        }
                        
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
                        double twoElecInt = 0.0;
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                           OrbitalType orbitalLambda = atomB->GetValence()[lambda-firstAOIndexB];
                           twoElecInt = this->GetCoulombInt(orbitalMu, 
                                                            orbitalLambda, 
                                                            atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                           OrbitalType orbitalNu = atomA->GetValence()[nu-firstAOIndexA];
                           twoElecInt = this->GetExchangeInt(orbitalMu, 
                                                             orbitalNu, 
                                                             atomA);
                        }
                        else{
                           twoElecInt = 0.0;
                        }

                        double temp = 0.0;
                        if(isMoPOcc){
                           int p = numberOcc - (moP+1);
                           temp = 4.0*xiOcc[p][nu]*eta[lambda][sigma]
                                 -1.0*xiOcc[p][lambda]*eta[nu][sigma]
                                 -1.0*xiOcc[p][sigma]*eta[nu][lambda];
                        }
                        else{
                           int p = moP - numberOcc;
                           temp = 4.0*xiVir[p][nu]*eta[lambda][sigma]
                                 -1.0*xiVir[p][lambda]*eta[sigma][nu]
                                 -1.0*xiVir[p][sigma]*eta[lambda][nu];
                        }
                        value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                     }  
                  }
               }
            }
         }
      }
   }
   return value;
}

// see (20) - (23) in [PT_1997]
void Mndo::CalcQVector(double* q, 
                       double* delta, 
                       double** xiOcc,
                       double** xiVir,
                       double** eta,
                       vector<MoIndexPair> nonRedundantQIndeces,
                       vector<MoIndexPair> redundantQIndeces){
   MallocerFreer::GetInstance()->InitializeDoubleMatrix1d(
                                 q,
                                 nonRedundantQIndeces.size()+redundantQIndeces.size());

   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      bool isMoICIMO = nonRedundantQIndeces[i].isMoICIMO;
      bool isMoJCIMO = nonRedundantQIndeces[i].isMoJCIMO;
      if(!isMoICIMO && isMoJCIMO){
         q[i] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta);
      }
      else if(isMoICIMO && !isMoJCIMO){
         q[i] = -1.0*this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
      }
      else if(isMoICIMO && isMoJCIMO){
         q[i] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta)
               -this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
      }
      else{
         q[i] = 0.0;
      }
   }
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<redundantQIndeces.size(); i++){
      int r = nonRedundantQIndeces.size() + i;
      int moI = redundantQIndeces[i].moI;
      int moJ = redundantQIndeces[i].moJ;
      if(moI == moJ){
         int rr = moI - (numberOcc-numberActiveOcc);
         q[r] = delta[rr];
      }
      else{
         q[r] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta)
               -this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
      }
   }
   /* 
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      printf("q[%d] = %e\n",i,q[i]);
   }
   for(int i=0; i<redundantQIndeces.size(); i++){
      int r = nonRedundantQIndeces.size() + i;
      printf("q[%d] = %e\n",r,q[r]);
   }
   */
}

// see (43) and (45) in [PT_1996].
// This method calculates "\Gamma_{NR} - K_{NR}".
// Note taht K_{NR} is not calculated.
void Mndo::CalcKNRMatrix(double** kNR, vector<MoIndexPair> nonRedundantQIndeces){
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      for(int j=i; j<nonRedundantQIndeces.size(); j++){
         int moK = nonRedundantQIndeces[j].moI;
         int moL = nonRedundantQIndeces[j].moJ;
         kNR[i][j] = this->GetGammaNRElement(moI, moJ, moK, moL)
                    -this->GetKNRElement(moI, moJ, moK, moL);
      }
   }
}

// see (44) and (46) in [PT_1996].
// This method calculates "K_{R}^{\dager} * \Gamma_{R}^{-1}".
// Note taht K_{R}^{\dager} is not calculated.
void Mndo::CalcKRDagerMatrix(double** kRDager, 
                             vector<MoIndexPair> nonRedundantQIndeces,
                             vector<MoIndexPair> redundantQIndeces){
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      for(int j=0; j<redundantQIndeces.size(); j++){
         int moK = redundantQIndeces[j].moI;
         int moL = redundantQIndeces[j].moJ;
         kRDager[i][j] = this->GetKRDagerElement(moI, moJ, moK, moL)
                        /this->GetGammaRElement(moK, moL, moK, moL);
      }
   }
}

// right hand side of (54) in [PT_1996]      
void Mndo::CalcAuxiliaryVector(double* y, 
                               double* q, 
                               double** kRDager, 
                               vector<MoIndexPair> nonRedundantQIndeces, 
                               vector<MoIndexPair> redundantQIndeces){
   MallocerFreer::GetInstance()->InitializeDoubleMatrix1d(
                                 y,
                                 nonRedundantQIndeces.size());
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      y[i] += q[i]/this->GetNNRElement(moI, moJ, moI, moJ);
      for(int j=0; j<redundantQIndeces.size(); j++){
         int k = nonRedundantQIndeces.size() + j; 
         y[i] += kRDager[i][j]*q[k];
      }
   }
}

void Mndo::TransposeFockMatrixMatrix(double** transposedFockMatrix){
   for(int i=0; i<this->molecule->GetTotalNumberAOs(); i++){
      for(int j=0; j<this->molecule->GetTotalNumberAOs(); j++){
         transposedFockMatrix[j][i] = this->fockMatrix[i][j];
      }
   }
}

// each element (mu, nu) of z matrix.
// see (57) in [PT_1996]
double Mndo::GetZMatrixForceElement(double* y,
                                    double* q,
                                    double** transposedFockMatrix,
                                    vector<MoIndexPair> nonRedundantQIndeces,
                                    vector<MoIndexPair> redundantQIndeces,
                                    int mu,
                                    int nu){
   double value=0.0;
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      value += y[i]
              *transposedFockMatrix[mu][moI]
              *transposedFockMatrix[nu][moJ];
   }
   for(int i=0; i<redundantQIndeces.size(); i++){
      int j = nonRedundantQIndeces.size() + i;
      int moI = redundantQIndeces[i].moI;
      int moJ = redundantQIndeces[i].moJ;
      value += (q[j]/this->GetGammaRElement(moI, moJ, moI, moJ))
              *transposedFockMatrix[mu][moI]
              *transposedFockMatrix[nu][moJ];
   }
   return value;
}

void Mndo::CalcXiMatrices(double** xiOcc, 
                          double** xiVir, 
                          int exciteState, 
                          double** transposedFockMatrix){
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   MallocerFreer::GetInstance()->InitializeDoubleMatrix2d(
                                 xiOcc, numberActiveOcc, numberAOs);
   MallocerFreer::GetInstance()->InitializeDoubleMatrix2d(
                                 xiVir, numberActiveVir, numberAOs);
   // xiOcc
   #pragma omp parallel for schedule(auto)
   for(int p=0; p<numberActiveOcc; p++){
      for(int mu=0; mu<numberAOs; mu++){
         for(int a=0; a<numberActiveVir; a++){
            int moA = numberOcc + a;
            int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(p,a);
            xiOcc[p][mu] += this->matrixCIS[exciteState][slaterDeterminantIndex]
                           *transposedFockMatrix[mu][moA];
         }
      }
   }
   // xiVir
   #pragma omp parallel for schedule(auto)
   for(int p=0; p<numberActiveVir; p++){
      for(int mu=0; mu<numberAOs; mu++){
         for(int i=0; i<numberActiveOcc; i++){
            int moI = numberOcc - (i+1);
            int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,p);
            xiVir[p][mu] += this->matrixCIS[exciteState][slaterDeterminantIndex]
                           *transposedFockMatrix[mu][moI];
         }
      }
   }
}

// see [PT_1996, PT_1997]
void Mndo::CalcZMatrixForce(vector<int> elecStates){
   if(this->etaMatrixForce == NULL){
      stringstream ss;
      ss << this->errorMessageCalcZMatrixForceEtaNull;
      throw MolDSException(ss.str());
   }
   this->CheckZMatrixForce(elecStates); 

   // creat MO-index-pair for Q variables. 
   vector<MoIndexPair> nonRedundantQIndeces;
   vector<MoIndexPair> redundantQIndeces;
   this->CalcActiveSetVariablesQ(&nonRedundantQIndeces, &redundantQIndeces);

   // malloc temporary arraies
   double* delta = NULL; // Delta matrix, see (9) in [PT_1997]
   double* q = NULL; //// Q-vector in (19) in [PT_1997]
   double** kNR = NULL; // K_{NR} matrix, see (45) in [PT_1996]
   double** kRDager = NULL; // Dagar of K_{R} matrix, see (46) in [PT_1996]
   double* y = NULL; // y-vector in (54) in [PT_1996]
   double** transposedFockMatrix = NULL; // transposed Fock matrix
   double** xiOcc = NULL;
   double** xiVir = NULL;
   this->MallocTempMatrixForZMatrix(&delta,
                                    &q,
                                    &kNR,
                                    &kRDager,
                                    &y,
                                    &transposedFockMatrix,
                                    &xiOcc,
                                    &xiVir,
                                    nonRedundantQIndeces.size(),
                                    redundantQIndeces.size());
   try{
      this->TransposeFockMatrixMatrix(transposedFockMatrix);
      this->CalcKNRMatrix(kNR, nonRedundantQIndeces);
      this->CalcKRDagerMatrix(kRDager, nonRedundantQIndeces,redundantQIndeces);
      int groundState=0;
      for(int n=0; n<elecStates.size(); n++){
         if(groundState < elecStates[n]){
            int exciteState = elecStates[n]-1;
            this->CalcDeltaVector(delta, exciteState);
            this->CalcXiMatrices(xiOcc, xiVir, exciteState, transposedFockMatrix);
            this->CalcQVector(q, 
                              delta, 
                              xiOcc, 
                              xiVir,
                              this->etaMatrixForce[n],
                              nonRedundantQIndeces, 
                              redundantQIndeces);
            this->CalcAuxiliaryVector(y, q, kRDager, nonRedundantQIndeces, redundantQIndeces);
            // solve (54) in [PT_1996]
            MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsysv(kNR, 
                                                                   y, 
                                                                   nonRedundantQIndeces.size());
            // calculate each element of Z matrix.
            #pragma omp parallel for schedule(auto)
            for(int mu=0; mu<this->molecule->GetTotalNumberAOs(); mu++){
               for(int nu=0; nu<this->molecule->GetTotalNumberAOs(); nu++){
                  this->zMatrixForce[n][mu][nu] = this->GetZMatrixForceElement(
                                                        y,
                                                        q,
                                                        transposedFockMatrix,
                                                        nonRedundantQIndeces,
                                                        redundantQIndeces,
                                                        mu,
                                                        nu);
               }
            }  

         }
      }
   }
   catch(MolDSException ex){
      this->FreeTempMatrixForZMatrix(&delta,
                                     &q,
                                     &kNR,
                                     &kRDager,
                                     &y,
                                     &transposedFockMatrix,
                                     &xiOcc,
                                     &xiVir,
                                     nonRedundantQIndeces.size(),
                                     redundantQIndeces.size());
      throw ex;
   }
   this->FreeTempMatrixForZMatrix(&delta,
                                  &q,
                                  &kNR,
                                  &kRDager,
                                  &y,
                                  &transposedFockMatrix,
                                  &xiOcc,
                                  &xiVir,
                                  nonRedundantQIndeces.size(),
                                  redundantQIndeces.size());
}

void Mndo::CalcEtaMatrixForce(vector<int> elecStates){
   this->CheckEtaMatrixForce(elecStates); 
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int groundState = 0;
   double** transposedFockMatrix = NULL; // transposed Fock matrix
   transposedFockMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                                        numberAOs,
                                                        numberAOs);
   try{
      this->TransposeFockMatrixMatrix(transposedFockMatrix);
      for(int n=0; n<elecStates.size(); n++){
         if(groundState < elecStates[n]){
            int exciteState = elecStates[n]-1;

            // calc each element
            #pragma omp parallel for schedule(auto)
            for(int mu=0; mu<numberAOs; mu++){
               for(int nu=0; nu<numberAOs; nu++){
                  for(int i=0; i<numberActiveOcc; i++){
                     int moI = numberOcc-(i+1);
                     for(int a=0; a<numberActiveVir; a++){
                        int moA = numberOcc+a;
                        int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,a);
                        this->etaMatrixForce[n][mu][nu] 
                                 += this->matrixCIS[exciteState][slaterDeterminantIndex]
                                   *transposedFockMatrix[mu][moI]
                                   *transposedFockMatrix[nu][moA];
                     }
                  }
               }
            }

         }
      }
   }
   catch(MolDSException ex){
      if(transposedFockMatrix != NULL){
         MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&transposedFockMatrix,numberAOs);
         //cout << "transposedFockMatrix  deleted" << endl;
      }
      throw ex;
   }
   if(transposedFockMatrix != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&transposedFockMatrix,numberAOs);
      //cout << "transposedFockMatrix  deleted" << endl;
   }
}

bool Mndo::RequiresExcitedStatesForce(vector<int> elecStates){
   bool requires = true;
   if(elecStates.size()==1 && elecStates[0]==0){
      requires = false;
   }
   return requires;
}

void Mndo::CalcForceHFElecCoreAttractionPart(double* force, 
                                             int atomAIndex, 
                                             int atomBIndex,
                                             double***** twoElecTwoCoreFirstDeriv){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->orbitalElectronPopulation[mu][nu]
                       *this->GetElectronCoreAttractionFirstDerivative(
                                   atomAIndex, 
                                   atomBIndex, 
                                   mu-firstAOIndexA, 
                                   nu-firstAOIndexA,
                                   twoElecTwoCoreFirstDeriv,
                                   (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceHFOverlapPart(double* force, 
                                  int atomAIndex, 
                                  int atomBIndex,
                                  double*** overlapDer){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   int numberAOsB = atomB->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
         double bondParameter = atomA->GetBondingParameter(
                                       this->theory, 
                                       atomA->GetValence()[mu-firstAOIndexA]) 
                               +atomB->GetBondingParameter(
                                       this->theory, 
                                       atomB->GetValence()[nu-firstAOIndexB]); 
         bondParameter*=0.5;
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->orbitalElectronPopulation[mu][nu]
                       *bondParameter
                       *overlapDer[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void Mndo::CalcForceHFTwoElecPart(double* force, 
                                  int atomAIndex, 
                                  int atomBIndex,
                                  double***** twoElecTwoCoreFirstDeriv){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   int numberAOsB = atomB->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= 0.5
                             *this->orbitalElectronPopulation[mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *twoElecTwoCoreFirstDeriv[mu-firstAOIndexA]
                                                      [nu-firstAOIndexA]
                                                      [lambda-firstAOIndexB]
                                                      [sigma-firstAOIndexB]
                                                      [(CartesianType)i];
                  force[i] += 0.25
                             *this->orbitalElectronPopulation[mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *twoElecTwoCoreFirstDeriv[mu-firstAOIndexA]
                                                      [nu-firstAOIndexA]
                                                      [lambda-firstAOIndexB]
                                                      [sigma-firstAOIndexB]
                                                      [(CartesianType)i];
               }
            }
         }
      }
   }
}

void Mndo::CalcForceExcitedStaticPart(double* force, 
                                      int elecStateIndex,
                                      int atomAIndex, 
                                      int atomBIndex,
                                      double***** twoElecTwoCoreFirstDeriv){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   int numberAOsB = atomB->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  double temp= 2.0*this->etaMatrixForce[elecStateIndex][mu][nu]
                                  *this->etaMatrixForce[elecStateIndex][lambda][sigma]
                              -1.0*this->etaMatrixForce[elecStateIndex][mu][lambda]
                                  *this->etaMatrixForce[elecStateIndex][nu][sigma];
                  force[i] += temp
                             *twoElecTwoCoreFirstDeriv[mu-firstAOIndexA]
                                                      [nu-firstAOIndexA]
                                                      [lambda-firstAOIndexB]
                                                      [sigma-firstAOIndexB]
                                                      [i];
               }
            }
         }
      }
   }
}

void Mndo::CalcForceExcitedElecCoreAttractionPart(double* force, 
                                                  int elecStateIndex,
                                                  int atomAIndex, 
                                                  int atomBIndex,
                                                  double***** twoElecTwoCoreFirstDeriv){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *this->GetElectronCoreAttractionFirstDerivative(
                                   atomAIndex, 
                                   atomBIndex, 
                                   mu-firstAOIndexA, 
                                   nu-firstAOIndexA,
                                   twoElecTwoCoreFirstDeriv,
                                   (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceExcitedOverlapPart(double* force, 
                                       int elecStateIndex,
                                       int atomAIndex, 
                                       int atomBIndex,
                                       double*** overlapDer){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   int numberAOsB = atomB->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
         double bondParameter = atomA->GetBondingParameter(
                                       this->theory, 
                                       atomA->GetValence()[mu-firstAOIndexA]) 
                               +atomB->GetBondingParameter(
                                       this->theory, 
                                       atomB->GetValence()[nu-firstAOIndexB]); 
         bondParameter *= 0.5;
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *bondParameter
                       *overlapDer[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void Mndo::CalcForceExcitedTwoElecPart(double* force, 
                                       int elecStateIndex,
                                       int atomAIndex, 
                                       int atomBIndex,
                                       double***** twoElecTwoCoreFirstDeriv){
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   int firstAOIndexA = atomA->GetFirstAOIndex();
   int firstAOIndexB = atomB->GetFirstAOIndex();
   int numberAOsA = atomA->GetValence().size();
   int numberAOsB = atomB->GetValence().size();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= this->zMatrixForce[elecStateIndex][mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *twoElecTwoCoreFirstDeriv[mu-firstAOIndexA]
                                                      [nu-firstAOIndexA]
                                                      [lambda-firstAOIndexB]
                                                      [sigma-firstAOIndexB]
                                                      [i];
                  force[i] += 0.50
                             *this->zMatrixForce[elecStateIndex][mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *twoElecTwoCoreFirstDeriv[mu-firstAOIndexA]
                                                      [nu-firstAOIndexA]
                                                      [lambda-firstAOIndexB]
                                                      [sigma-firstAOIndexB]
                                                      [(CartesianType)i];
               }
            }
         }
      }
   }
}

// electronicStateIndex is index of the electroinc eigen state.
// "electronicStateIndex = 0" means electronic ground state. 
void Mndo::CalcForce(vector<int> elecStates){
   this->CheckMatrixForce(elecStates);
   if(this->RequiresExcitedStatesForce(elecStates)){
      this->CalcEtaMatrixForce(elecStates);
      this->CalcZMatrixForce(elecStates);
   }
   #pragma omp parallel
   {
      double***** twoElecTwoCoreFirstDeriv = MallocerFreer::GetInstance()->MallocDoubleMatrix5d(
                                                                           dxy,
                                                                           dxy,
                                                                           dxy,
                                                                           dxy,
                                                                           CartesianType_end);
      double*** overlapDer = MallocerFreer::GetInstance()->MallocDoubleMatrix3d(
                                                           OrbitalType_end, 
                                                           OrbitalType_end, 
                                                           CartesianType_end);
      try{
         #pragma omp for schedule(auto)
         for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
            Atom* atomA = (*molecule->GetAtomVect())[a];
            int firstAOIndexA = atomA->GetFirstAOIndex();
            int numberAOsA = atomA->GetValence().size();
            for(int b=0; b<this->molecule->GetAtomVect()->size(); b++){
               if(a != b){
                  Atom* atomB = (*molecule->GetAtomVect())[b];
                  int firstAOIndexB = atomB->GetFirstAOIndex();
                  int numberAOsB = atomB->GetValence().size();

                  // calc. first derivative of overlap.
                  this->CalcDiatomicOverlapFirstDerivative(overlapDer, atomA, atomB);
                  // calc. first derivative of two elec two core interaction
                  this->CalcTwoElecTwoCoreDiatomicFirstDerivatives(twoElecTwoCoreFirstDeriv, 
                                                                   a, 
                                                                   b);

                  // core repulsion part
                  double coreRepulsion[CartesianType_end] = {0.0,0.0,0.0};
                  for(int i=0; i<CartesianType_end; i++){
                     coreRepulsion[i] += this->GetDiatomCoreRepulsionFirstDerivative(
                                               a, b, (CartesianType)i);
                  }  
                  // electron core attraction part (ground state)
                  double forceElecCoreAttPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceHFElecCoreAttractionPart(forceElecCoreAttPart,
                                                          a,
                                                          b,
                                                          twoElecTwoCoreFirstDeriv);
                  // overlap part (ground state)
                  double forceOverlapPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceHFOverlapPart(forceOverlapPart, 
                                               a,
                                               b,
                                               overlapDer);
                  // two electron part (ground state)
                  double forceTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceHFTwoElecPart(forceTwoElecPart,
                                               a,
                                               b,
                                               twoElecTwoCoreFirstDeriv);
                  // sum up contributions from each part (ground state)
                  #pragma omp critical
                  {
                     for(int n=0; n<elecStates.size(); n++){
                        for(int i=0; i<CartesianType_end; i++){
                           this->matrixForce[n][a][i] -= coreRepulsion[i];
                           this->matrixForce[n][a][i] += forceElecCoreAttPart[i];
                           this->matrixForce[n][a][i] += forceOverlapPart[i];
                           this->matrixForce[n][a][i] += forceTwoElecPart[i];
                           this->matrixForce[n][b][i] -= forceElecCoreAttPart[i];
                           this->matrixForce[n][b][i] -= forceOverlapPart[i];
                           this->matrixForce[n][b][i] -= forceTwoElecPart[i];
                        }
                     }
                  }
                  // excited state potential
                  for(int n=0; n<elecStates.size(); n++){
                     if(0<elecStates[n]){
                        // static part
                        double forceExcitedStaticPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedStaticPart(forceExcitedStaticPart,
                                                         n,
                                                         a,
                                                         b,
                                                         twoElecTwoCoreFirstDeriv);
                        // sum up contributions from static part (excited state)
                        for(int i=0; i<CartesianType_end; i++){
                           this->matrixForce[n][b][i] += forceExcitedStaticPart[i];
                           this->matrixForce[n][a][i] -= forceExcitedStaticPart[i];
                        }

                        // response part
                        // electron core attraction part (excited state)
                        double forceExcitedElecCoreAttPart[CartesianType_end]={0.0,0.0,0.0};
                        this->CalcForceExcitedElecCoreAttractionPart(
                                                   forceExcitedElecCoreAttPart,
                                                   n,
                                                   a,
                                                   b,
                                                   twoElecTwoCoreFirstDeriv);
                        // overlap part (excited state)
                        double forceExcitedOverlapPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedOverlapPart(forceExcitedOverlapPart, 
                                                          n,
                                                          a,
                                                          b,
                                                          overlapDer);
                        // two electron part (ground state)
                        double forceExcitedTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedTwoElecPart(forceExcitedTwoElecPart,
                                                          n,
                                                          a,
                                                          b,
                                                          twoElecTwoCoreFirstDeriv);
                        // sum up contributions from response part (excited state)
                        #pragma omp critical
                        {
                           for(int i=0; i<CartesianType_end; i++){
                              this->matrixForce[n][a][i] += forceExcitedElecCoreAttPart[i];
                              this->matrixForce[n][a][i] += forceExcitedOverlapPart[i];
                              this->matrixForce[n][a][i] += forceExcitedTwoElecPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedElecCoreAttPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedOverlapPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedTwoElecPart[i];
                           }
                        }

                     }
                  }
               }
            }
         }
      }
      catch(MolDSException ex){
         this->FreeCalcForceTempMatrices(&overlapDer, &twoElecTwoCoreFirstDeriv);
         throw ex;
      }
      this->FreeCalcForceTempMatrices(&overlapDer, &twoElecTwoCoreFirstDeriv);
   }
}

void Mndo::FreeCalcForceTempMatrices(double**** overlapDer, double****** twoElecTwoCoreFirstDeriv){
   if(*overlapDer != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix3d(overlapDer, 
                                                       OrbitalType_end,
                                                       OrbitalType_end);
      //cout << "overlapDer deleted\n";
   }
   if(*twoElecTwoCoreFirstDeriv != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix5d(twoElecTwoCoreFirstDeriv,
                                                       dxy,
                                                       dxy,
                                                       dxy,
                                                       dxy);
      //cout << "twoElecCoreFirstDeriv deleted\n";
   }
}

void Mndo::CalcTwoElecTwoCore(double****** twoElecTwoCore, Molecule* molecule){
   if(twoElecTwoCore == NULL){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix6d(twoElecTwoCore, 
                                                      molecule->GetAtomVect()->size(),
                                                      molecule->GetAtomVect()->size(),
                                                      dxy, dxy, dxy, dxy);
   } 

   #pragma omp parallel
   {
      double**** twoElecTwoCoreDiatomic = MallocerFreer::GetInstance()->MallocDoubleMatrix4d(
                                                                        dxy, dxy, dxy, dxy);
      try{
         // note that terms with condition a==b are not needed to calculate. 
         #pragma omp for schedule(auto)
         for(int a=0; a<molecule->GetAtomVect()->size(); a++){
            for(int b=a+1; b<molecule->GetAtomVect()->size(); b++){
               this->CalcTwoElecTwoCoreDiatomic(twoElecTwoCoreDiatomic, a, b);
               for(int mu=0; mu<dxy; mu++){
                  for(int nu=mu; nu<dxy; nu++){
                     for(int lambda=0; lambda<dxy; lambda++){
                        for(int sigma=lambda; sigma<dxy; sigma++){
                           double value = twoElecTwoCoreDiatomic[mu][nu][lambda][sigma];
                           twoElecTwoCore[a][b][mu][nu][lambda][sigma] = value;
                           twoElecTwoCore[a][b][mu][nu][sigma][lambda] = value;
                           twoElecTwoCore[a][b][nu][mu][lambda][sigma] = value;
                           twoElecTwoCore[a][b][nu][mu][sigma][lambda] = value;
                           twoElecTwoCore[b][a][lambda][sigma][mu][nu] = value;
                           twoElecTwoCore[b][a][lambda][sigma][nu][mu] = value;
                           twoElecTwoCore[b][a][sigma][lambda][mu][nu] = value;
                           twoElecTwoCore[b][a][sigma][lambda][nu][mu] = value;
                        }
                     }
                  }
               }
            }
         }
      }
      catch(MolDSException ex){
         MallocerFreer::GetInstance()->FreeDoubleMatrix4d(&twoElecTwoCoreDiatomic,
                                                          dxy, dxy, dxy);
         throw ex;
      }
      MallocerFreer::GetInstance()->FreeDoubleMatrix4d(&twoElecTwoCoreDiatomic,
                                                       dxy, dxy, dxy);
   }
}

// Calculation of two electrons two cores integral (mu, nu | lambda, sigma), 
// taht is, Eq. (9) in ref. [DT_1977-2].
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy] can be treatable.
void Mndo::CalcTwoElecTwoCoreDiatomic(double**** matrix, int atomAIndex, int atomBIndex){
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB->GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }
   else{
      atomA = (*this->molecule->GetAtomVect())[atomAIndex];
      atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix4d(matrix, dxy, dxy, dxy, dxy);
   } 

   // calclation in diatomic frame
   for(int mu=0; mu<atomA->GetValence().size(); mu++){
      for(int nu=0; nu<atomA->GetValence().size(); nu++){
         for(int lambda=0; lambda<atomB->GetValence().size(); lambda++){
            for(int sigma=0; sigma<atomB->GetValence().size(); sigma++){
               matrix[mu][nu][lambda][sigma] = this->GetNddoRepulsionIntegral(
                                               atomA, 
                                               atomA->GetValence()[mu],
                                               atomA->GetValence()[nu],
                                               atomB, 
                                               atomB->GetValence()[lambda],
                                               atomB->GetValence()[sigma]);
                     
            }
         }
      }
   }

   // rotate matirix into the space frame
   double** rotatingMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                            OrbitalType_end, OrbitalType_end);
   try{
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->RotateTwoElecTwoCoreDiatomicToSpaceFramegc(matrix, rotatingMatrix);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);

   /*
   printf("(mu, nu | lambda, sigma) matrix\n"); 
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               printf("mu=%d nu=%d lambda=%d sigma=%d $e\n",
                        mu,nu,lambda,sigma,matrix[mu][nu][lambda][sigma]);
            }
         }
      }
   }
   */
}

// Calculation of first derivatives of the two electrons two cores integral 
// (mu, nu | lambda, sigma), taht is, Eq. (9) in ref. [DT_1977-2].
// This derivative is related to the coordinates of atomA.
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy][CartesianType_end] can be treatable.
void Mndo::CalcTwoElecTwoCoreDiatomicFirstDerivatives(double***** matrix, 
                                                      int atomAIndex, 
                                                      int atomBIndex){
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB->GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }
   else{
      atomA = (*this->molecule->GetAtomVect())[atomAIndex];
      atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix5d(matrix, 
                                                             dxy, 
                                                             dxy, 
                                                             dxy, 
                                                             dxy, 
                                                             CartesianType_end);
   } 

   double** rotatingMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                             OrbitalType_end, OrbitalType_end);
   double*** rMatDeri = MallocerFreer::GetInstance()->MallocDoubleMatrix3d
                                   (OrbitalType_end, OrbitalType_end, CartesianType_end);
   double**** twoElecTwoCoreDiatomic = MallocerFreer::GetInstance()->MallocDoubleMatrix4d(
                                                                     dxy, dxy, dxy, dxy);
   try{
      // calclation in diatomic frame
      for(int mu=0; mu<atomA->GetValence().size(); mu++){
         for(int nu=0; nu<atomA->GetValence().size(); nu++){
            for(int lambda=0; lambda<atomB->GetValence().size(); lambda++){
               for(int sigma=0; sigma<atomB->GetValence().size(); sigma++){
                  for(int c=0; c<CartesianType_end; c++){
                     matrix[mu][nu][lambda][sigma][c] 
                        = this->GetNddoRepulsionIntegralFirstDerivative(
                                atomA, 
                                atomA->GetValence()[mu],
                                atomA->GetValence()[nu],
                                atomB, 
                                atomB->GetValence()[lambda],
                                atomB->GetValence()[sigma],
                                (CartesianType)c);
                  }  
                  twoElecTwoCoreDiatomic[mu][nu][lambda][sigma] 
                     = this->GetNddoRepulsionIntegral(
                             atomA, 
                             atomA->GetValence()[mu],
                             atomA->GetValence()[nu],
                             atomB, 
                             atomB->GetValence()[lambda],
                             atomB->GetValence()[sigma]);
               }
            }
         }
      }

      // rotate matirix into the space frame
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->CalcRotatingMatrixFirstDerivatives(rMatDeri, atomA, atomB);
      this->RotateTwoElecTwoCoreDiatomicFirstDerivativesToSpaceFramegc(matrix, 
                                                                       twoElecTwoCoreDiatomic,
                                                                       rotatingMatrix,
                                                                       rMatDeri);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);
      MallocerFreer::GetInstance()->FreeDoubleMatrix3d(&rMatDeri, 
                                                       OrbitalType_end,
                                                       OrbitalType_end);
      MallocerFreer::GetInstance()->FreeDoubleMatrix4d(&twoElecTwoCoreDiatomic,
                                                       dxy, dxy, dxy);
      throw ex;
   }
   MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);
   MallocerFreer::GetInstance()->FreeDoubleMatrix3d(&rMatDeri, 
                                                    OrbitalType_end,
                                                    OrbitalType_end);
   MallocerFreer::GetInstance()->FreeDoubleMatrix4d(&twoElecTwoCoreDiatomic,
                                                    dxy, dxy, dxy);
}

// Rotate 4-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateTwoElecTwoCoreDiatomicToSpaceFramegc(double**** matrix, double** rotatingMatrix){
   double oldMatrix[dxy][dxy][dxy][dxy];
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               oldMatrix[mu][nu][lambda][sigma] = matrix[mu][nu][lambda][sigma];
            }
         }
      }
   }
   
   // rotate
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = 0.0;
               for(int i=0; i<dxy; i++){
                  for(int j=0; j<dxy; j++){
                     for(int k=0; k<dxy; k++){
                        for(int l=0; l<dxy; l++){
                           matrix[mu][nu][lambda][sigma] += oldMatrix[i][j][k][l] 
                                                            *rotatingMatrix[mu][i] 
                                                            *rotatingMatrix[nu][j] 
                                                            *rotatingMatrix[lambda][k] 
                                                            *rotatingMatrix[sigma][l];
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

// Rotate 5-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateTwoElecTwoCoreDiatomicFirstDerivativesToSpaceFramegc(double***** matrix, 
                                                                      double**** twoElecTwoCoreDiatomic,
                                                                      double** rotatingMatrix,
                                                                      double*** rMatDeri){
   double oldMatrix[dxy][dxy][dxy][dxy][CartesianType_end];
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               for(int c=0; c<CartesianType_end; c++){
                  oldMatrix[mu][nu][lambda][sigma][c] = matrix[mu][nu][lambda][sigma][c];
               }
            }
         }
      }
   }
   
   // rotate
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               for(int c=0; c<CartesianType_end; c++){
                  matrix[mu][nu][lambda][sigma][c] = 0.0;
                  for(int i=0; i<dxy; i++){
                     for(int j=0; j<dxy; j++){
                        for(int k=0; k<dxy; k++){
                           for(int l=0; l<dxy; l++){
                              matrix[mu][nu][lambda][sigma][c] 
                                 += oldMatrix[i][j][k][l][c]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += twoElecTwoCoreDiatomic[i][j][k][l]
                                   *rMatDeri[mu][i][c]
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += twoElecTwoCoreDiatomic[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rMatDeri[nu][j][c]
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += twoElecTwoCoreDiatomic[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rMatDeri[lambda][k][c]
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += twoElecTwoCoreDiatomic[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rMatDeri[sigma][l][c];
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral(Atom* atomA, OrbitalType mu, OrbitalType nu,
                                      Atom* atomB, OrbitalType lambda, OrbitalType sigma){
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      value = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, mu, mu)
                  -this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, nu, nu));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegral;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// First derivative of NDDO repulsion integral.
// This derivation is related to the coordinate of atomA
// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegralFirstDerivative(
                                       Atom* atomA, OrbitalType mu, OrbitalType nu,
                                       Atom* atomB, OrbitalType lambda, OrbitalType sigma,
                                       CartesianType axisA){
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dRabDa = (atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA])/Rab;
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      value *= dRabDa;
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           mux, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muy, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      DA = atomA->GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB->GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA->GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB->GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegralFirstDerivative(
                         atomA, mu, mu, atomB, mu, mu, axisA)
                  -this->GetNddoRepulsionIntegralFirstDerivative(
                         atomA, mu, mu, atomB, nu, nu, axisA));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegralFirstDerivative;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteraction(MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab){
   double value = 0.0;
   double a = rhoA + rhoB;

   if(multipoleA == sQ && multipoleB == sQ){
      value = pow(pow(Rab,2.0) + pow(a,2.0), -0.5);
   }
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(Rab+DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleA, Qxx,
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/2.0 + pow(temp3,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(Rab,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(Rab-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(Rab+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(Rab+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             +pow(temp5,-0.5)/4.0 - pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0
             +pow(temp5,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(Rab-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 - pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/4.0 + pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, multipoleB, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/16.0 + pow(temp2,-0.5)/16.0 
             +pow(temp3,-0.5)/16.0 + pow(temp4,-0.5)/16.0
             -pow(temp5,-0.5)/8.0 - pow(temp6,-0.5)/8.0
             -pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0
             +pow(temp9,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(Rab-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 - pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/8.0 + pow(temp6,-0.5)/8.0
             +pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(Rab,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/2.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

// First derivative of semiempirical multipole-multipole interactions.
// This derivativ is related to the nuclear distance Rab.
// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                                                  MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab){
   double value = 0.0;
   double a = rhoA + rhoB;

   if(multipoleA == sQ && multipoleB == sQ){
      value = -1.0*Rab*pow(pow(Rab,2.0) + pow(a,2.0), -1.5);
   }
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(Rab+DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(a,2.0);
      value = (Rab+DB)*pow(temp1,-1.5)/2.0 
             -(Rab-DB)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      value = Rab*pow(temp1,-1.5)/2.0 
             -Rab*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleA, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      value = (Rab+2.0*DB)*pow(temp1,-1.5)/4.0 
             -(Rab)*pow(temp2,-1.5)/2.0 
             +(Rab-2.0*DB)*pow(temp3,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(Rab,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/2.0 
             -(Rab)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    mux, mux, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+DB,2.0) + pow(a,2.0);
      value = (Rab+DA-DB)*pow(temp1,-1.5)/4.0 
             -(Rab+DA+DB)*pow(temp2,-1.5)/4.0 
             -(Rab-DA-DB)*pow(temp3,-1.5)/4.0 
             +(Rab-DA+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(Rab-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-(Rab-DB)*pow(temp1,-1.5)/4.0 
             +(Rab-DB)*pow(temp2,-1.5)/4.0 
             +(Rab+DB)*pow(temp3,-1.5)/4.0 
             -(Rab+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    mux, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(Rab+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-(Rab+DA)*pow(temp1,-1.5)/4.0 
             +(Rab-DA)*pow(temp2,-1.5)/4.0 
             +(Rab+DA)*pow(temp3,-1.5)/4.0 
             -(Rab-DA)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    muz, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(Rab+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-(Rab+DA-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(Rab-DA-2.0*DB)*pow(temp2,-1.5)/8.0 
             -(Rab+DA+2.0*DB)*pow(temp3,-1.5)/8.0 
             +(Rab-DA+2.0*DB)*pow(temp4,-1.5)/8.0
             +(Rab+DA       )*pow(temp5,-1.5)/4.0 
             -(Rab-DA       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/8.0 
             +(Rab)*pow(temp2,-1.5)/8.0 
             -(Rab)*pow(temp3,-1.5)/4.0 
             -(Rab)*pow(temp4,-1.5)/4.0
             +(Rab)*pow(temp5,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/4.0 
             -(Rab)*pow(temp2,-1.5)/4.0 
             -(Rab)*pow(temp3,-1.5)/4.0 
             +(Rab)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(Rab-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(Rab+2.0*DB)*pow(temp2,-1.5)/8.0 
             -(Rab-2.0*DB)*pow(temp3,-1.5)/8.0 
             -(Rab+2.0*DB)*pow(temp4,-1.5)/8.0
             -(Rab       )*pow(temp5,-1.5)/4.0 
             +(Rab       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxx, multipoleB, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab+2.0*DA-2.0*DB)*pow(temp1,-1.5)/16.0 
             +(Rab+2.0*DA+2.0*DB)*pow(temp2,-1.5)/16.0 
             +(Rab-2.0*DA-2.0*DB)*pow(temp3,-1.5)/16.0 
             +(Rab-2.0*DA+2.0*DB)*pow(temp4,-1.5)/16.0
             -(Rab+2.0*DA)*pow(temp5,-1.5)/8.0 
             -(Rab-2.0*DA)*pow(temp6,-1.5)/8.0
             -(Rab+2.0*DB)*pow(temp7,-1.5)/8.0 
             -(Rab-2.0*DB)*pow(temp8,-1.5)/8.0
             +(Rab)*pow(temp9,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(Rab-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (Rab+DA-DB)*pow(temp1,-1.5)/8.0 
             -(Rab+DA-DB)*pow(temp2,-1.5)/8.0 
             -(Rab+DA+DB)*pow(temp3,-1.5)/8.0 
             +(Rab+DA+DB)*pow(temp4,-1.5)/8.0
             -(Rab-DA-DB)*pow(temp5,-1.5)/8.0 
             +(Rab-DA-DB)*pow(temp6,-1.5)/8.0
             +(Rab-DA+DB)*pow(temp7,-1.5)/8.0 
             -(Rab-DA+DB)*pow(temp8,-1.5)/8.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(Rab,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/4.0 
             +(Rab)*pow(temp2,-1.5)/4.0 
             -(Rab)*pow(temp3,-1.5)/2.0;
      value *= -1.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}
}

