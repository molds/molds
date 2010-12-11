#ifndef INCLUDED_ZINDOS
#define INCLUDED_ZINDOS

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
#include"../cndo/Cndo2.h"

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_zindo{

/***
 *  Refferences for 
 */
class ZindoS : public MolDS_cndo::Cndo2{
private:
   double GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom); // Apendix in [BZ_1979]
   double GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom); // Apendix in [BZ_1979]
protected:
   virtual void CalcGammaAB(double** gammaAB, Molecule* molecule);
   virtual void SetMessages();
   virtual double GetFockDiagElement(Atom* atomA, int atomAIndex, int firstAOIndexA, 
                                     int mu, Molecule* molecule, double** gammaAB,
                                     double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                     bool isGuess);
   virtual double GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                        int firstAOIndexA, int firstAOIndexB,
                                        int mu, int nu, Molecule* molecule, double** gammaAB, double** overelap,
                                        double** orbitalElectronPopulation, bool isGuess);
   virtual void SetEnableAtomTypes();
public:
   ZindoS();
   ~ZindoS();
};

ZindoS::ZindoS() : MolDS_cndo::Cndo2(){
   this->theory = ZINDOS;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "ZindoS created\n";
}

ZindoS::~ZindoS(){
   //cout << "ZindoS deleted\n";
}

void ZindoS::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in zindo::ZindoS::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in zindo::ZindoS::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in zindo::ZindoS::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in zindo::ZindoS::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_zindo::ZindoS::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_zindo::ZindoS::GetExchangeInt: Invalid orbitalType.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tZINDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: ZINDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: ZINDO/S-SCF  **********\n\n\n";
}

void ZindoS::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double ZindoS::GetFockDiagElement(Atom* atomA, int atomAIndex, int firstAOIndexA, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 bool isGuess){
   double value;
   value = atomA->GetCoreIntegral(atomA->GetValence()[mu-firstAOIndexA], 
                                     gammaAB[atomAIndex][atomAIndex], 
                                     isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
      for(int v=0; v<atomA->GetValence().size(); v++){
         OrbitalType orbitalLam = atomA->GetValence()[v];
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalLam, gammaAB[atomAIndex][atomAIndex], atomA);
         exchange = this->GetExchangeInt(orbitalMu, orbitalLam, gammaAB[atomAIndex][atomAIndex], atomA);
         lammda = firstAOIndexA + v;
         temp += orbitalElectronPopulation[lammda][lammda]*(coulomb - 0.5*exchange);
      }
      value += temp;
   
      temp = 0.0;
      for(int BB=0; BB<molecule->GetAtomVect()->size(); BB++){
         if(BB != atomAIndex){
            Atom* atomBB = (*molecule->GetAtomVect())[BB];
            temp += ( atomicElectronPopulation[BB] - atomBB->GetCoreCharge()  )
                     *gammaAB[atomAIndex][BB];
         }
      }
      value += temp;
   }

   return value;
}

double ZindoS::GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                    int firstAOIndexA, int firstAOIndexB,
                                    int mu, int nu, Molecule* molecule, double** gammaAB, double** overlap,
                                    double** orbitalElectronPopulation, bool isGuess){
   double value;
   double K = 1.0;
   if(m <= atomA->GetValenceShellType() || m <= atomB->GetValenceShellType()){
      K = 0.75;
   }
   double bondParameter = 0.5*K*(atomA->GetBondingParameter() + atomB->GetBondingParameter()); 

   if(isGuess){
      value = bondParameter*overlap[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(atomAIndex == atomBIndex){
         OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
         OrbitalType orbitalNu = atomA->GetValence()[nu-firstAOIndexA];
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, gammaAB[atomAIndex][atomAIndex], atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, gammaAB[atomAIndex][atomAIndex], atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlap[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]*gammaAB[atomAIndex][atomBIndex];
      }
   }

   return value;
}

void ZindoS::CalcGammaAB(double** gammaAB, Molecule* molecule){
   // Do nothing;
}

// Apendix in [BZ_1972]
// ZINDO Coulomb Interaction
double ZindoS::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom){

   double value=0.0;

   // ToDo: Coulomb interaction
   /*
   if( orbital1 == s && orbital2 == s){ 
      value = gamma;
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = gamma;
   }   
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz ) && orbital2 == s){ 
      value = gamma;
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = gamma + 4.0*atom->GetIndoF2()/25.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = gamma - 2.0*atom->GetIndoF2()/25.0;
   }   
   else{
      cout << this->errorMessageCoulombInt;
      cout << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      exit(EXIT_FAILURE);
   }   
   */
   return value;

}

// Apendix in [BZ_1972]
// ZINDO Exchange Interaction
double ZindoS::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom){

   double value=0.0;

   // ToDo: Exchange interaction
   /*
   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, gamma, atom);
   }   
   else if( (orbital1 == s) && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetIndoG1()/3.0;
   }   
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz) && orbital2 == s  ){  
      value = atom->GetIndoG1()/3.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = 3.0*atom->GetIndoF2()/25.0;
   }   
   else{
      cout << this->errorMessageExchangeInt;
      cout << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      exit(EXIT_FAILURE);
   }   
   */

   return value;
}

}
#endif



