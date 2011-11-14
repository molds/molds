#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include"../base/MolDSException.h"
#include"../base/Enums.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../cndo/Cndo2.h"
#include"../zindo/ZindoS.h"
#include"../mndo/Mndo.h"
#include"../am1/Am1.h"
#include"../pm3/Pm3.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_md{

/***
 *  Velocty Verlet is used here.
 */
class MD{
public:
   MD();
   ~MD();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
   void DoesMD();
private:
   string messageinitialConditionMD;
   string messageStartMD;
   string messageEndMD;
   string messageStartStepMD;
   string messageEndStepMD;
   string messageEnergies;
   string messageEnergiesTitle;
   string messageCoreKineticEnergy;
   string messageCoreRepulsionEnergy;
   string messageElectronicEnergy;
   string messageTotalEnergy;
   string messageErrorEnergy;
   string messageTime;
   string errorMessageNotEnebleTheoryType;
   string errorMessageTheoryType;
   MolDS_cndo::Cndo2* cndo;
   vector<TheoryType> enableTheoryTypes;
   void CheckEnableTheoryType(TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
   void OutputEnergies(double initialEnergy);
   double OutputEnergies();
};

MD::MD(){
   this->cndo = NULL;
   this->SetEnableTheoryTypes();
   this->SetMessages();
   //cout << "MD created \n";
}

MD::~MD(){
   //cout << "MD deleted\n";
}

void MD::SetTheory(MolDS_cndo::Cndo2* cndo){
   // check enable electonic theory
   this->CheckEnableTheoryType(cndo->GetTheoryType());
   this->cndo = cndo;
}

void MD::DoesMD(){
   cout << this->messageStartMD;

   int totalSteps = Parameters::GetInstance()->GetTotalStepsMD();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMD();
   double dt = Parameters::GetInstance()->GetTimeWidthMD();
   double time = 0.0;
   bool requireGuess = false;
   Molecule* molecule = this->cndo->GetMolecule();
   double** matrixForce = this->cndo->GetForce(elecState);
   double initialEnergy;

   // output initial conditions
   cout << this->messageinitialConditionMD;
   initialEnergy = this->OutputEnergies();
   cout << endl;
   molecule->OutputConfiguration();
   molecule->OutputXyzCOM();
   molecule->OutputXyzCOC();
   molecule->OutputMomenta();

   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepMD << s+1 << endl;

      // update momenta
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
         }
      }

      // update coordinates
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*molecule->GetAtomVect())[a];
         double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dt*atom->GetPxyz()[i]/coreMass;
         }
      }

      // update molecular basics
      molecule->CalcXyzCOM();
      molecule->CalcXyzCOC();

      // update electronic structure
      this->cndo->DoesSCF(requireGuess);
      if(elecState > 0){
         this->cndo->DoesCIS();
      }
      else if(elecState < 0){
         // ToDo: Error
      }

      // update force
      matrixForce = this->cndo->GetForce(elecState);

      // update momenta
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
         }
      }

      // output results
      this->OutputEnergies(initialEnergy);
      molecule->OutputConfiguration();
      molecule->OutputXyzCOM();
      molecule->OutputXyzCOC();
      molecule->OutputMomenta();
      cout << this->messageTime << dt*((double)s+1)/Parameters::GetInstance()->GetFs2AU() << endl;
      cout << this->messageEndStepMD << s+1 << endl;
   }

   cout << this->messageEndMD;
}

void MD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in md::MD::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartMD = "**********  START: Molecular dynamics  **********\n";
   this->messageEndMD = "**********  DONE: Molecular dynamics  **********\n";
   this->messageinitialConditionMD = "\n\t========= Initial conditions \n";
   this->messageStartStepMD = "\n\t========== START: MD step ";
   this->messageEndStepMD =     "\t========== DONE: MD step ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageCoreKineticEnergy =   "Core kinetic     ";
   this->messageCoreRepulsionEnergy = "Core repulsion   ";
   this->messageElectronicEnergy = "Electronic\n\t\t(inc. core rep.)";
   this->messageTotalEnergy =         "Total            ";
   this->messageErrorEnergy =         "Error            ";
   this->messageTime = "\tTime in [fs]: ";
}

double MD::OutputEnergies(){
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMD();
   Molecule* molecule = this->cndo->GetMolecule();
   double coreKineticEnergy = 0.0;
   for(int a=0; a<molecule->GetAtomVect()->size(); a++){
      Atom* atom = (*molecule->GetAtomVect())[a];
      double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
      for(int i=0; i<CartesianType_end; i++){
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
                             this->cndo->GetCoreRepulsionEnergy(),
                             this->cndo->GetCoreRepulsionEnergy()
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageElectronicEnergy.c_str(), 
                             this->cndo->GetElectronicEnergy(elecState),
                             this->cndo->GetElectronicEnergy(elecState)
                             /Parameters::GetInstance()->GetEV2AU());
   printf("\t\t%s\t%e\t%e\n",this->messageTotalEnergy.c_str(), 
                             (coreKineticEnergy + this->cndo->GetElectronicEnergy(elecState)),
                             (coreKineticEnergy + this->cndo->GetElectronicEnergy(elecState))
                             /Parameters::GetInstance()->GetEV2AU());

   return (coreKineticEnergy + this->cndo->GetElectronicEnergy(elecState));
}

void MD::OutputEnergies(double initialEnergy){
   double energy = this->OutputEnergies();
   printf("\t\t%s\t%e\t%e\n\n",this->messageErrorEnergy.c_str(), 
                             (initialEnergy - energy),
                             (initialEnergy - energy)
                             /Parameters::GetInstance()->GetEV2AU());
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
   this->enableTheoryTypes.push_back(MNDO);
   this->enableTheoryTypes.push_back(AM1);
   this->enableTheoryTypes.push_back(PM3);
}

void MD::CheckEnableTheoryType(TheoryType theoryType){

   bool isEnable = false;
   for(int i=0; i<this->enableTheoryTypes.size();i++){
      if(theoryType == this->enableTheoryTypes[i]){
         isEnable = true;
         break;
      }
   }
   if(!isEnable){
      stringstream ss;
      ss << this->errorMessageNotEnebleTheoryType;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      throw MolDSException(ss.str());
   }
}

}



