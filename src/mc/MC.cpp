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
#include<boost/shared_ptr.hpp>
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"../base/Enums.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"../base/ElectronicStructureFactory.h"
#include"MC.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_mc{
MC::MC(){
   this->molecule = NULL;
   this->SetMessages();
   //cout << "MC created \n";
}

MC::~MC(){
   //cout << "MC deleted\n";
}

void MC::SetMolecule(Molecule* molecule){
   this->molecule = molecule;
}

void MC::DoMC(){
   cout << this->messageStartMC;

   // malloc electornic structure
   boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::GetInstance()->Create());
   electronicStructure->SetMolecule(this->molecule);

   int totalSteps = Parameters::GetInstance()->GetTotalStepsMC();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   double dr = Parameters::GetInstance()->GetStepWidthMC();
   double temperatur = Parameters::GetInstance()->GetTemperatureMC();
   bool requireGuess = false;

   // initial calculation
   electronicStructure->DoSCF();
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicStructure->DoCIS();
   }

   // output initial conditions
   cout << this->messageinitialConditionMC;
   this->OutputEnergies(electronicStructure);
   this->molecule->OutputConfiguration();
   this->molecule->OutputXyzCOM();
   this->molecule->OutputXyzCOC();
   this->molecule->OutputMomenta();

   for(int s=0; s<totalSteps; s++){
      requireGuess = (s==0) ? true : false;
      cout << this->messageStartStepMC << s+1 << endl;
/*
      // create candidate
      Molecule candidate(*this->molecule);
      for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*candidate.GetAtomVect())[0];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dr;
         }
      }
      candidate.CalcXyzCOM();
      candidate.CalcXyzCOC();
*/
      // calculate electronic structure of the candidate
      electronicStructure->DoSCF();
      if(Parameters::GetInstance()->RequiresCIS()){
         electronicStructure->DoCIS();
      }

      // output results
      this->OutputEnergies(electronicStructure);
      this->molecule->OutputConfiguration();
      this->molecule->OutputXyzCOM();
      this->molecule->OutputXyzCOC();

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

void MC::OutputEnergies(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure){
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   // output energies:
   cout << this->messageEnergies;
   cout << this->messageEnergiesTitle;
   printf("\t\t%s\t%e\t%e\n",this->messageCoreRepulsionEnergy.c_str(), 
                             electronicStructure->GetCoreRepulsionEnergy(),
                             electronicStructure->GetCoreRepulsionEnergy()
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageElectronicEnergy.c_str(), 
                             electronicStructure->GetElectronicEnergy(elecState),
                             electronicStructure->GetElectronicEnergy(elecState)
                             /Parameters::GetInstance()->GetEV2AU());
}

}



