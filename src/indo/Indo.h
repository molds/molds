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

class Indo : public MolDS_cndo::Cndo2{
public:
   Indo();
   ~Indo();
protected:
   void SetMessages();
   double GetFockDiagElement(Atom* atomA, int atomAIndex, int firstAOIndexA, 
                             int mu, Molecule* molecule, double** gammaAB,
                             double** orbitalElectronPopulation, double* atomicElectronPopulation,
                             bool isGuess);
   double GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                int firstAOIndexA, int firstAOIndexB,
                                int mu, int nu, Molecule* molecule, double** gammaAB, double** overelap,
                                double** orbitalElectronPopulation, bool isGuess);
};

Indo::Indo() : MolDS_cndo::Cndo2(){
   this->SetMessages();
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
   this->messageSCFMetConvergence = "\n\n\n\t\tINDO-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: INDO-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: INDO-SCF  **********\n\n\n";
}

double Indo::GetFockDiagElement(Atom* atomA, int atomAIndex, int firstAOIndexA, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 bool isGuess){
   double value;
   value = atomA->GetIndoCoreIntegral(atomA->GetValence()[mu-firstAOIndexA], 
                                     gammaAB[atomAIndex][atomAIndex], 
                                     isGuess);

   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
      for(int v=0; v<atomA->GetValence().size(); v++){
         OrbitalType orbitalLam = atomA->GetValence()[v];
         coulomb  = atomA->GetIndoCoulombInt(orbitalMu, orbitalLam, gammaAB[atomAIndex][atomAIndex]);
         exchange = atomA->GetIndoExchangeInt(orbitalMu, orbitalLam, gammaAB[atomAIndex][atomAIndex]);
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
         coulomb  = atomA->GetIndoCoulombInt(orbitalMu, orbitalNu, gammaAB[atomAIndex][atomAIndex]); 
         exchange = atomA->GetIndoExchangeInt(orbitalMu, orbitalNu, gammaAB[atomAIndex][atomAIndex]); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlap[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]*gammaAB[atomAIndex][atomBIndex];
      }
   }

   return value;
}



}
#endif



