#ifndef INCLUDED_INDO
#define INCLUDED_INDO

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
#include"../cndo/Cndo2.h"

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_indo{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class Indo : public MolDS_cndo::Cndo2{
private:
   double GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom); // Indo Coulomb Interaction, (3.87) - (3.91) in J. A. Pople book.
   double GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom); // Indo Exchange Interaction, (3.87) - (3.91) in J. A. Pople book.
protected:
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
   Indo();
   ~Indo();
};

Indo::Indo() : MolDS_cndo::Cndo2(){
   this->theory = INDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "Indo created\n";
}

Indo::~Indo(){
   //cout << "Indo deleted\n";
}

void Indo::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in indo::Indo::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in indo::Indo::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in indo::Indo::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in indo::Indo::CheckEnableAtomType: Not enable atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_indo::Indo::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_indo::Indo::GetExchangeInt: Invalid orbitalType.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tINDO-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: INDO-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: INDO-SCF  **********\n\n\n";
}

void Indo::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(Li);
   this->enableAtomTypes.push_back(Be);
   this->enableAtomTypes.push_back(B);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(F);
}

double Indo::GetFockDiagElement(Atom* atomA, int atomAIndex, int firstAOIndexA, int mu, 
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

double Indo::GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                    int firstAOIndexA, int firstAOIndexB,
                                    int mu, int nu, Molecule* molecule, double** gammaAB, double** overlap,
                                    double** orbitalElectronPopulation, bool isGuess){
   double value;
   double K = 1.0;  // = 1.0 or 0.75, see Eq. (3.79) in J. A. Pople book
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

// (3.87) - (3.91) in J. A. Pople book.
// Indo Coulomb Interaction
double Indo::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom){

   double value=0.0;
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

   return value;

}

// (3.87) - (3.91) in J. A. Pople book.
// Indo Exchange Interaction
double Indo::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, Atom* atom){

   double value=0.0;

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

   return value;

}


}
#endif



