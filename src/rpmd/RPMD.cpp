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
#include"../mc/MC.h"
#include"RPMD.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_rpmd{
RPMD::RPMD(){
   this->SetMessages();
   this->SetEnableTheoryTypes();
   //cout << "RPMD created \n";
}

RPMD::~RPMD(){
   //cout << "RPMD deleted\n";
}

void RPMD::DoRPMD(const Molecule& refferenceMolecule){
   cout << this->messageStartRPMD;

   // check
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexRPMD();
   this->CheckEnableTheoryType(theory, elecState);

   double temperature = Parameters::GetInstance()->GetTemperatureRPMD();
   unsigned long seed = Parameters::GetInstance()->GetSeedRPMD();
   int totalSteps = Parameters::GetInstance()->GetTotalStepsRPMD();
   int numBeads = Parameters::GetInstance()->GetNumberBeadsRPMD();
   double dt = Parameters::GetInstance()->GetTimeWidthRPMD();
   double kB = Parameters::GetInstance()->GetBoltzmann();
   int numAtom = refferenceMolecule.GetAtomVect()->size();

   // create Beads
   vector<boost::shared_ptr<Molecule> > molecularBeads;
   vector<boost::shared_ptr<ElectronicStructure> > electronicStructureBeads;
   for(int b=0; b<numBeads; b++){
      // create molecular beads
      boost::shared_ptr<Molecule> molecule(new Molecule());
      *molecule = refferenceMolecule;
      molecularBeads.push_back(molecule);
      // create electronic structure beads
      boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::GetInstance()->Create());
      electronicStructure->SetMolecule(molecule.get());
      electronicStructureBeads.push_back(electronicStructure);
   }

   // initialize Beads
   for(int b=0; b<numBeads; b++){
      double stepWidth = 0.05;
      boost::shared_ptr<MolDS_mc::MC> mc(new MolDS_mc::MC());
      Molecule* molecule = molecularBeads[b].get();
      mc->SetMolecule(molecule);
      mc->DoMC(molecule->GetAtomVect()->size(), elecState, temperature, stepWidth, seed+b);
   }

   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepRPMD << s+1 << endl;

      // update momenta
      for(int b=0; b<numBeads; b++){
         int preB  = b==0 ? numBeads-1 : b-1;
         int postB = b==numBeads-1 ? 0 : b+1;
         double** electronicForceMatrix = electronicStructureBeads[b]->GetForce(elecState);;
         for(int a=0; a<numAtom; a++){
            Atom* atom = (*molecularBeads[b]->GetAtomVect())[a];
            Atom* preAtom = (*molecularBeads[preB]->GetAtomVect())[a];
            Atom* postAtom = (*molecularBeads[postB]->GetAtomVect())[a];
            double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
            for(int i=0; i<CartesianType_end; i++){
               double beadsForce = -1.0*coreMass*pow(kB*temperature*(double)numBeads,2.0)
                                  *(2.0*atom->GetXyz()[i] - preAtom->GetXyz()[i] - postAtom->GetXyz()[i]);
               double force = beadsForce + electronicForceMatrix[a][i];
               atom->GetPxyz()[i] += 0.5*dt*(force);
            }
         }
      }

   }
/*
   for(int i=0; i<numBeads; i++){
      Molecule* bead = molecularBeads[i];
      delete bead;
   }
*/

/*
   //Molecule trialMolecule(*this->molecule);
   // malloc electornic structure
   boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::GetInstance()->Create());
   electronicStructure->SetMolecule(this->molecule);

   int totalSteps = Parameters::GetInstance()->GetTotalStepsRPMD();
   double dt = Parameters::GetInstance()->GetTimeWidthRPMD();
   double time = 0.0;
   bool requireGuess = false;
   double** matrixForce = NULL;
   double initialEnergy;

   // initial calculation
   electronicStructure->DoSCF();
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicStructure->DoCIS();
   }
   matrixForce = electronicStructure->GetForce(elecState);

   // output initial conditions
   cout << this->messageinitialConditionRPMD;
   initialEnergy = this->OutputEnergies(electronicStructure);
   cout << endl;
   this->molecule->OutputConfiguration();
   this->molecule->OutputXyzCOC();
   this->molecule->OutputMomenta();

   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepRPMD << s+1 << endl;

      // update momenta
      for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*this->molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
         }
      }

      // update coordinates
      for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*this->molecule->GetAtomVect())[a];
         double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dt*atom->GetPxyz()[i]/coreMass;
         }
      }
      this->molecule->CalcXyzCOM();
      this->molecule->CalcXyzCOC();

      // update electronic structure
      electronicStructure->DoSCF(requireGuess);
      if(Parameters::GetInstance()->RequiresCIS()){
         electronicStructure->DoCIS();
      }

      // update force
      matrixForce = electronicStructure->GetForce(elecState);

      // update momenta
      for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*this->molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
         }
      }

      // output results
      this->OutputEnergies(electronicStructure, initialEnergy);
      this->molecule->OutputConfiguration();
      this->molecule->OutputXyzCOC();
      this->molecule->OutputMomenta();
      cout << this->messageTime << dt*((double)s+1)/Parameters::GetInstance()->GetFs2AU() << endl;
      cout << this->messageEndStepRPMD << s+1 << endl;
   }
*/
   cout << this->messageEndRPMD;
}

void RPMD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageElecState = "\tElectronic state = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in rpmd::RPMD::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartRPMD 
      = "**********  START: Ring Polymer Molecular dynamics  **********\n";
   this->messageEndRPMD 
      = "**********  DONE: Ring Polymer Molecular dynamics  **********\n";
   this->messageinitialConditionRPMD = "\n\t========= Initial conditions \n";
   this->messageStartStepRPMD = "\n\t========== START: RPMD step ";
   this->messageEndStepRPMD =     "\t========== DONE: RPMD step ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageCoreKineticEnergy =   "Core kinetic     ";
   this->messageCoreRepulsionEnergy = "Core repulsion   ";
   this->messageElectronicEnergy = "Electronic\n\t\t(inc. core rep.)";
   this->messageTotalEnergy =         "Total            ";
   this->messageErrorEnergy =         "Error            ";
   this->messageTime = "\tTime in [fs]: ";
}

/*
double RPMD::OutputEnergies(boost::shared_ptr<ElectronicStructure> electronicStructure){
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexRPMD();
   double coreKineticEnergy = 0.0;
   for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
      Atom* atom = (*this->molecule->GetAtomVect())[a];
      double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
      for(int i=0; i<CartesianType_end; i++){
      Atom* atom = (*this->molecule->GetAtomVect())[a];
         coreKineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
      }
   }  
   // output energies:
   cout << this->messageEnergies;
   cout << this->messageEnergiesTitle;
   printf("\t\t%s\t%e\t%e\n",this->messageCoreKineticEnergy.c_str(), 
                             coreKineticEnergy,
                             coreKineticEnergy/Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageCoreRepulsionEnergy.c_str(), 
                             electronicStructure->GetCoreRepulsionEnergy(),
                             electronicStructure->GetCoreRepulsionEnergy()
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageElectronicEnergy.c_str(), 
                             electronicStructure->GetElectronicEnergy(elecState),
                             electronicStructure->GetElectronicEnergy(elecState)
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageTotalEnergy.c_str(), 
                             (coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState)),
                             (coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState))
                             /Parameters::GetInstance()->GetEV2AU());

   return (coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState));
}

void RPMD::OutputEnergies(boost::shared_ptr<ElectronicStructure> electronicStructure, 
                        double initialEnergy){
   double energy = this->OutputEnergies(electronicStructure);
   printf("\t\t%s\t%e\t%e\n\n",this->messageErrorEnergy.c_str(), 
                             (initialEnergy - energy),
                             (initialEnergy - energy)
                             /Parameters::GetInstance()->GetEV2AU());
}
*/

void RPMD::SetEnableTheoryTypes(){
   // ground state
   this->enableGroundStateTheoryTypes.clear();
   this->enableGroundStateTheoryTypes.push_back(ZINDOS);
   this->enableGroundStateTheoryTypes.push_back(MNDO);
   this->enableGroundStateTheoryTypes.push_back(AM1);
   this->enableGroundStateTheoryTypes.push_back(PM3);
   this->enableGroundStateTheoryTypes.push_back(PM3PDDG);

   // excited state
   this->enableExcitedStateTheoryTypes.clear();
   this->enableExcitedStateTheoryTypes.push_back(MNDO);
   this->enableExcitedStateTheoryTypes.push_back(AM1);
   this->enableExcitedStateTheoryTypes.push_back(PM3);
   this->enableExcitedStateTheoryTypes.push_back(PM3PDDG);
}

void RPMD::CheckEnableTheoryType(TheoryType theoryType, int elecState){

   bool isEnable = false;
   int groundState = 0;
   if(elecState == groundState){
      for(int i=0; i<this->enableGroundStateTheoryTypes.size();i++){
         if(theoryType == this->enableGroundStateTheoryTypes[i]){
            isEnable = true;
            break;
         }
      }
   }
   else{
      for(int i=0; i<this->enableExcitedStateTheoryTypes.size();i++){
         if(theoryType == this->enableExcitedStateTheoryTypes[i]){
            isEnable = true;
            break;
         }
      }
   }
   if(!isEnable){
      stringstream ss;
      ss << this->errorMessageNotEnebleTheoryType;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      ss << this->errorMessageElecState << elecState << endl;
      throw MolDSException(ss.str());
   }
}

}



