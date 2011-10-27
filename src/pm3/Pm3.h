#ifndef INCLUDED_PM3
#define INCLUDED_PM3

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_pm3{

/***
 *  Main Refferences for PM3 are [S_1989, S_1989-2, S_1991, S_2004, S_2007]
 */
class Pm3 : public MolDS_am1::Am1{
public:
   Pm3();
   virtual ~Pm3();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
private:
};

Pm3::Pm3() : MolDS_am1::Am1(){
   this->theory = PM3;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "Pm3 created\n";
}

Pm3::~Pm3(){
   //cout << "Pm3 deleted\n";
}

void Pm3::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in pm3::Pm3::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in pm3::Pm3::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in pm3::Pm3::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in pm3::Pm3::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_pm3::Pm3::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_pm3::Pm3::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in pm3::Pm3::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in pm3::Pm3::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in pm3:: Pm3::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in pm3:: Pm3::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in pm3::Pm3::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in pm3::Pm3::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in pm3::Pm3::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix 
      = "Error in pm3::Pm3::CalcTwoElecTwoCoreDiatomic: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms
      = "Error in pm3::Pm3::CalcTwoElecTwoCoreDiatomic: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix
      = "Error in pm3::Pm3::CalcTwoElecTwoCoreDiatomicFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms
      = "Error in pm3::Pm3::CalcTwoElecTwoCoreDiatomicFirstDerivatives: Atom A and B is same.\n"; 
   this->messageSCFMetConvergence = "\n\n\n\t\tPM3/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: PM3/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: PM3/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: PM3/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: PM3/S-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for PM3-CIS met convergence criterion(^^b\n\n\n";
}

void Pm3::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

}
#endif



