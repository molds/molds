//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"../base/Enums.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"MC.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_mc{
MC::MC(){
   this->electronicStructure = NULL;
   this->SetMessages();
   //cout << "MC created \n";
}

MC::~MC(){
   //cout << "MC deleted\n";
}

void MC::SetTheory(ElectronicStructure* electronicStructure){
   // check enable electonic theory
   this->electronicStructure = electronicStructure;
}

void MC::DoMC(){
   cout << this->messageStartMC;

   int totalSteps = Parameters::GetInstance()->GetTotalStepsMC();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   double dr = Parameters::GetInstance()->GetStepWidthMC();
   double temperatur = Parameters::GetInstance()->GetTemperatureMC();
   bool requireGuess = false;
   //Molecule* molecule = this->electronicStructure->GetMolecule();

   // output initial conditions
   /*
   cout << this->messageinitialConditionMC;
   this->OutputEnergies();
   molecule->OutputConfiguration();
   molecule->OutputXyzCOM();
   molecule->OutputXyzCOC();
   molecule->OutputMomenta();
   */

   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepMC << s+1 << endl;
/*
      // create candidate
      Molecule candidate(*molecule);
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*candidate.GetAtomVect())[0];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dr;
         }
      }
      candidate.CalcXyzCOM();
      candidate.CalcXyzCOC();

      // calculate electronic structure of the candidate
      this->electronicStructure->DoSCF(requireGuess);
      if(elecState > 0){
         this->electronicStructure->DoCIS();
      }
      else if(elecState < 0){
         // ToDo: Error
      }

      delete molecule;
      molecule = &candidate;

      // update molecular basics
      molecule->CalcXyzCOM();
      molecule->CalcXyzCOC();

      // output results
      this->OutputEnergies();
      molecule->OutputConfiguration();
      molecule->OutputXyzCOM();
      molecule->OutputXyzCOC();
*/
      cout << this->messageEndStepMC << s+1 << endl;
   }

   cout << this->messageEndMC;
}

void MC::SetMessages(){
   this->messageStartMC = "**********  START: Monte Carlo  **********\n";
   this->messageEndMC = "**********  DONE: Monte Carlo  **********\n";
   this->messageinitialConditionMC = "\n\t========= Initial conditions \n";
   this->messageStartStepMC = "\n\t========== START: MC step ";
   this->messageEndStepMC =     "\t========== DONE: MC step ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageCoreRepulsionEnergy = "Core repulsion   ";
   this->messageElectronicEnergy = "Electronic\n\t\t(inc. core rep.)";
}

void MC::OutputEnergies(){
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   //Molecule* molecule = this->electronicStructure->GetMolecule();
   // output energies:
   cout << this->messageEnergies;
   cout << this->messageEnergiesTitle;
   printf("\t\t%s\t%e\t%e\n",this->messageCoreRepulsionEnergy.c_str(), 
                             this->electronicStructure->GetCoreRepulsionEnergy(),
                             this->electronicStructure->GetCoreRepulsionEnergy()
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageElectronicEnergy.c_str(), 
                             this->electronicStructure->GetElectronicEnergy(elecState),
                             this->electronicStructure->GetElectronicEnergy(elecState)
                             /Parameters::GetInstance()->GetEV2AU());
}

}



