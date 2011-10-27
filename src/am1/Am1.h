#ifndef INCLUDED_AM1
#define INCLUDED_AM1

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_am1{

/***
 *  Main Refferences for AM1 are [DZHS_1985, DY_1990]
 */
class Am1 : public MolDS_mndo::Mndo{
public:
   Am1();
   virtual ~Am1();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcCoreRepulsionEnergy();
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        CartesianType axisA);
   virtual void CalcHFProperties();
   virtual void OutputHFResults(double** fockMatrix, double* energiesMO, 
                                double* atomicElectronPopulation, Molecule* molecule);
private:
   string errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
   string errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles;
   string errorMessageGetNddoRepulsionIntegral;
   string errorMessageGetNddoRepulsionIntegralFirstDerivative;
   string errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix;
   string errorMessageCalcTwoElecTwoCoreNullMatrix;
   string errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms;
   string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix;
   string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms;
};

Am1::Am1() : MolDS_mndo::Mndo(){
   this->theory = AM1;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "Am1 created\n";
}

Am1::~Am1(){
   if(this->twoElecTwoCore != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix6d(
                                    &this->twoElecTwoCore, 
                                    this->molecule->GetAtomVect()->size(),
                                    this->molecule->GetAtomVect()->size(),
                                    dxy,
                                    dxy,
                                    dxy);
      //cout << "twoElecTwoCore deleted\n";
   }
   //cout << "Am1 deleted\n";
}

void Am1::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in am1::Am1::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in am1::Am1::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in am1::Am1::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in am1::Am1::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_am1::Am1::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_am1::Am1::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in am1::Am1::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in am1::Am1::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in am1::Am1::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: Atom A and B is same.\n"; 
   this->messageSCFMetConvergence = "\n\n\n\t\tAM1/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: AM1/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: AM1/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: AM1/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: AM1/S-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for AM1-CIS met convergence criterion(^^b\n\n\n";
}

void Am1::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

void Am1::CalcCoreRepulsionEnergy(){
   double energy = 0.0;
   double distance = 0.0;
   double twoElecInt = 0.0;
   double alphaA = 0.0;
   double alphaB = 0.0;
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   double temp = 0.0;
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   for(int i=0; i<this->molecule->GetAtomVect()->size(); i++){
      atomA = (*this->molecule->GetAtomVect())[i];
      alphaA = atomA->GetNddoAlpha(this->theory);
      for(int j=i+1; j<this->molecule->GetAtomVect()->size(); j++){
         atomB = (*this->molecule->GetAtomVect())[j];
         alphaB = atomB->GetNddoAlpha(this->theory);
         distance = this->molecule->GetDistanceAtoms(i, j);
         if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                          atomB->GetAtomType() == O)  ){
            temp = 1.0 + (distance/ang2AU)*exp(-alphaB*distance) + exp(-alphaA*distance);
         }
         else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                               atomA->GetAtomType() == O)  ){
            temp = 1.0 + (distance/ang2AU)*exp(-alphaA*distance) + exp(-alphaB*distance);
         }
         else{
            temp = 1.0 + exp(-alphaA*distance) + exp(-alphaB*distance);
         }
         twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
         energy += atomA->GetCoreCharge()*atomB->GetCoreCharge()*twoElecInt*temp; 
         temp = 0.0;
         for(int i=0; i<4; i++){
            double kA = atomA->GetNddoParameterK(this->theory, i);
            double lA = atomA->GetNddoParameterL(this->theory, i);
            double mA = atomA->GetNddoParameterM(this->theory, i);
            double kB = atomB->GetNddoParameterK(this->theory, i);
            double lB = atomB->GetNddoParameterL(this->theory, i);
            double mB = atomB->GetNddoParameterM(this->theory, i);
            temp += kA*exp(-lA*pow(distance-mA,2.0));
            temp += kB*exp(-lB*pow(distance-mB,2.0));
         }
         energy += atomA->GetCoreCharge()*atomB->GetCoreCharge()*temp/(distance/ang2AU);
      }
   }
   this->coreRepulsionEnergy = energy;
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Am1::GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                   int atomBIndex, 
                                                   CartesianType axisA){
   double value =0.0;
   /*
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   Atom* atomA = (*this->molecule->GetAtomVect())[atomAIndex];
   Atom* atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   double alphaA = atomA->GetNddoAlpha(this->theory);
   double alphaB = atomB->GetNddoAlpha(this->theory);
   double Rab = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double dRabDa = (atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA])/Rab;
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   double twoElecIntFirstDeri = this->GetNddoRepulsionIntegralFirstDerivative(
                                      atomA, s, s, atomB, s, s, axisA);
   double temp1 = 0.0;
   if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                    atomB->GetAtomType() == O)  ){
      temp1 = 1.0 + (Rab/ang2AU)*exp(-alphaB*Rab) + exp(-alphaA*Rab);
   }
   else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                         atomA->GetAtomType() == O)  ){
      temp1 = 1.0 + (Rab/ang2AU)*exp(-alphaA*Rab) + exp(-alphaB*Rab);
   }
   else{
      temp1 = 1.0 + exp(-alphaA*Rab) + exp(-alphaB*Rab);
   }

   double temp2 = 0.0;
   if(atomA->GetAtomType() == H && (atomB->GetAtomType() == N || 
                                    atomB->GetAtomType() == O)  ){
      temp2 = (1.0/ang2AU)*exp(-alphaB*Rab) 
             -alphaB*(Rab/ang2AU)*exp(-alphaB*Rab) 
             -alphaA*exp(-alphaA*Rab);
   }
   else if(atomB->GetAtomType() == H && (atomA->GetAtomType() == N || 
                                         atomA->GetAtomType() == O)  ){
      temp2 = (1.0/ang2AU)*exp(-alphaA*Rab)
             -alphaA*(Rab/ang2AU)*exp(-alphaA*Rab) 
             -alphaB*exp(-alphaB*Rab);
   }
   else{
      temp2 = -alphaA*exp(-alphaA*Rab) 
              -alphaB*exp(-alphaB*Rab);
   }
   temp2 *= dRabDa;
   value = atomA->GetCoreCharge()*atomB->GetCoreCharge()
          *(twoElecIntFirstDeri*temp1 + twoElecInt*temp2); 
   */
   return value;
}

void Am1::CalcHFProperties(){
   MolDS_cndo::Cndo2::CalcHFProperties();
}

void Am1::OutputHFResults(double** fockMatrix, double* energiesMO, 
                          double* atomicElectronPopulation, Molecule* molecule){
   MolDS_cndo::Cndo2::OutputHFResults(fockMatrix, 
                                      energiesMO, 
                                      atomicElectronPopulation, 
                                      molecule);
}

}
#endif



