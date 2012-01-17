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
#include<boost/random.hpp>
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

   int totalSteps = Parameters::GetInstance()->GetTotalStepsMC();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   double dr = Parameters::GetInstance()->GetStepWidthMC();
   double temperatur = Parameters::GetInstance()->GetTemperatureMC();
  
   // create real random generator
   unsigned long seed=100;
	boost::mt19937 realGenerator(seed);
	boost::uniform_real<> range(0, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > realRand( realGenerator, range );

   // create integer random generator
	boost::mt19937 intGenerator(seed);
   boost::uniform_smallint<> dst( 1, this->molecule->GetAtomVect()->size() );
	boost::variate_generator<boost::mt19937&, boost::uniform_smallint<> > intRand( intGenerator, dst );

   ElectronicStructure* currentES = NULL;
   ElectronicStructure* trialES = NULL;
   Molecule trialMolecule(*this->molecule);

   // malloc electornic structure
   //boost::shared_ptr<ElectronicStructure> electronicStructure1(ElectronicStructureFactory::GetInstance()->Create());
   //boost::shared_ptr<ElectronicStructure> electronicStructure2(ElectronicStructureFactory::GetInstance()->Create());

   // initial calculation
   /*
   currentES = electronicStructure1.get();
   currentES->SetMolecule(this->molecule);
   currentES->DoSCF();
   if(Parameters::GetInstance()->RequiresCIS()){
      currentES->DoCIS();
   }
   
   // output initial conditions
   cout << this->messageinitialConditionMC;
   this->OutputEnergies(currentES);
   this->molecule->OutputConfiguration();
   this->molecule->OutputXyzCOM();
   this->molecule->OutputXyzCOC();
   this->molecule->OutputMomenta();
   */

   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepMC << s+1 << endl;
      this->CreateTrialConfiguration(&trialMolecule, this->molecule, &realRand, &intRand);

/*
      // create candidate
      for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
         Atom* trialAtom = (*trialMolecule.GetAtomVect())[a];
         Atom* atom = (*this->molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            trialAtom->GetXyz()[i] = atom->GetXyz()[i] + dr;
         }
      }
      trialMolecule.CalcXyzCOM();
      trialMolecule.CalcXyzCOC();
      
      // calculate electronic structure of the candidate
      bool requireGuess = (s==0) ? true : false;
      trialES->SetMolecule(&trialMolecule);
      trialES->DoSCF(requireGuess);
      if(Parameters::GetInstance()->RequiresCIS()){
         trialES->DoCIS();
      }

      // which Electronic Structure is used?
      if(UsesTrial(currentES, trialES)){
         swap(currentES, trialES);
         // synchronous molecule
         this->SynchronousMolecularConfiguration(this->molecule, &trialMolecule);
      }
      
      // output results
      this->OutputEnergies(currentES);
      this->molecule->OutputConfiguration();
      this->molecule->OutputXyzCOM();
      this->molecule->OutputXyzCOC();
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

void MC::CreateTrialConfiguration(Molecule* trial,
                                  Molecule* current,
                                  boost::random::variate_generator<
                                     boost::random::mt19937&,
                                     boost::uniform_real<>
                                  > (*realRand),
                                  boost::random::variate_generator<
                                     boost::random::mt19937&,
                                     boost::uniform_smallint<>
                                  > (*intRand)){
   int changedAtomIndex = (*intRand)();
   double changedDistance = (*realRand)();
   cout << "atom: " << changedAtomIndex << "   dist: " << changedDistance << endl;
}

bool MC::UsesTrial(ElectronicStructure* currentES, ElectronicStructure* trialES) const {
   return true;
}

void MC::SynchronousMolecularConfiguration(Molecule* target, 
                                           Molecule* refference) const{
   for(int a=0; a<target->GetAtomVect()->size(); a++){
      Atom* targetAtom = (*target->GetAtomVect())[a];
      Atom* refferenceAtom = (*refference->GetAtomVect())[a];
      for(int i=0; i<CartesianType_end; i++){
         targetAtom->GetXyz()[i] = refferenceAtom->GetXyz()[i];
      }
   }
   target->CalcXyzCOM();
   target->CalcXyzCOC();
}

void MC::OutputEnergies(MolDS_base::ElectronicStructure* electronicStructure) const{
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


/*
#include <iostream>
using namespace std;

void CallRand(boost::random::variate_generator<
               boost::random::mt19937&, 
               boost::uniform_real<> 
              > (*randFunction)){
   cout << (*randFunction)() << endl;
}

int main()
{

	for( int i=0; i<10; ++i ){
		//cout << rand() << endl;
      CallRand(&rand);
	}
	cout << endl;

	return 0;
}
*/

