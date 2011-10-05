#ifndef INCLUDED_MD
#define INCLUDED_MD

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_md{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class MD{
public:
   MD();
   ~MD();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
   void DoesMD();
private:
   string messageStartMD;
   string messageEndMD;
   string messageStartStepMD;
   string messageEndStepMD;
   string messageZindoSMD;
   MolDS_cndo::Cndo2* cndo;
   vector<TheoryType> enableTheoryTypes;
   string errorMessageNotEnebleTheoryType;
   string errorMessageTheoryType;
   void CheckEnableTheoryType(TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
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
   double** matrixForce;
   bool requireGuess = false;
   Molecule* molecule = this->cndo->GetMolecule();

   matrixForce = this->cndo->GetForce(elecState);

   // calc initial kinetic energy 
   double kineticEnergy = 0.0;
   for(int a=0; a<molecule->GetAtomVect()->size(); a++){
      Atom* atom = (*molecule->GetAtomVect())[a];
      double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
      for(int i=0; i<CartesianType_end; i++){
         kineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
      }
   }
   cout << "initial kinetic Energy = " << kineticEnergy << endl; 
   cout << "initial electronic energy = " << cndo->GetElectronicEnergy() << endl;
   cout << "initial total energy = " << kineticEnergy + cndo->GetElectronicEnergy() << endl;
   double iniEne =kineticEnergy + cndo->GetElectronicEnergy(); 
   for(int s=0; s<totalSteps; s++){
      cout << this->messageStartStepMD << s << endl;

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
      molecule->CalcTotalCoreRepulsionEnergy();

      molecule->OutputConfiguration();
      molecule->OutputXyzCOM();
      molecule->OutputXyzCOC();
      molecule->OutputTotalCoreRepulsionEnergy();

      // calc electronic structure and force
      this->cndo->DoesSCF(requireGuess);
      if(elecState > 0){
         this->cndo->DoesCIS();
      }
      else if(elecState < 0){
         // ToDo: Error
      }
      matrixForce = this->cndo->GetForce(elecState);

      // update momenta
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*molecule->GetAtomVect())[a];
         for(int i=0; i<CartesianType_end; i++){
            atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
         }
      }

      kineticEnergy = 0.0;
      for(int a=0; a<molecule->GetAtomVect()->size(); a++){
         Atom* atom = (*molecule->GetAtomVect())[a];
         double coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
         for(int i=0; i<CartesianType_end; i++){
            kineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
         }
      }  
      cout << "kinetic Energy = " << kineticEnergy << endl; 
      cout << "electronic energy = " << cndo->GetElectronicEnergy() << endl;
      cout << "total energy = " << kineticEnergy + cndo->GetElectronicEnergy() << endl;
      cout << "gosa energy [au] = " << kineticEnergy + cndo->GetElectronicEnergy() - iniEne << endl;
      //printf("total energy = %.10lf\n",kineticEnergy + cndo->GetElectronicEnergy());
      time = dt*((double)s+1)/Parameters::GetInstance()->GetFs2AU();
      cout << "Time: " << time << "[fs.]" << endl << endl << endl;
      cout << this->messageEndStepMD << s << endl;
   }


   cout << this->messageEndMD;
}

void MD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in md::MD::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartMD = "**********  START: Molecular dynamics  **********\n";
   this->messageEndMD = "**********  DONE: Molecular dynamics  **********\n";
   this->messageStartStepMD = "\n\t========== START: MD step ";
   this->messageEndStepMD =     "\t========== DONE: MD step ";
   this->messageZindoSMD = "\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!! A L A R T !!!!!!!!!!!!!!!!!!!!!!!!!\n\tNote that this MD algorythm can not work correctly with ZINDO/S. In this MD algorythm, the overlap matrix between AO can not calculated correctry. The reason is the using of the GTO expansion techniques for the first derivative of the overlap integrals. If you are one of the developpers, see ZndoS::CalcDiatomicOverlapInDiatomicFrame and comments in there. \n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
}

void MD::CheckEnableTheoryType(TheoryType theoryType){

   if(theoryType == ZINDOS){
      cout << this->messageZindoSMD;
   }

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
#endif



