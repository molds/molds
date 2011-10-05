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
   cout << "dt = " << dt << endl;
   cout << "initial kinetic Energy = " << kineticEnergy << endl; 
   cout << "initial electronic energy = " << cndo->GetElectronicEnergy() << endl;
   cout << "initial total energy = " << kineticEnergy + cndo->GetElectronicEnergy() << endl;
   double iniEne =kineticEnergy + cndo->GetElectronicEnergy(); 
   for(int s=0; s<totalSteps; s++){
      cout << "Starting MD step: " << s << endl;

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
      molecule->CalcXyzCOM();
      molecule->CalcXyzCOC();
      molecule->CalcTotalCoreRepulsionEnergy();

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
      cout << "Time: " << time << endl << endl << endl;
   }



}

void MD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in md::MD::CheckEnableTheoryType: Non available theory is set.\n";
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
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
#endif



