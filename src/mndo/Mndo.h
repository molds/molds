#ifndef INCLUDED_MNDO
#define INCLUDED_MNDO

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_mndo{

/***
 *  Main Refferences for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
class Mndo : public MolDS_zindo::ZindoS{
public:
   Mndo();
   ~Mndo();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcCoreRepulsionEnergy();
   virtual double GetFockDiagElement(Atom* atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     Molecule* molecule, 
                                     double** gammaAB,
                                     double** orbitalElectronPopulation, 
                                     double* atomicElectronPopulation,
                                     double****** twoElecTwoCore,
                                     bool isGuess);
   virtual double GetFockOffDiagElement(Atom* atomA, 
                                        Atom* atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
                                        int mu, int nu, 
                                        Molecule* molecule, 
                                        double** gammaAB, 
                                        double** overelap,
                                        double** orbitalElectronPopulation, 
                                        double****** twoElecTwoCore,
                                        bool isGuess);
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                   Atom* atomA, 
                                                   Atom* atomB);
   virtual void CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                                                double** diatomicOverlapDeri, 
                                                Atom* atomA, Atom* atomB);
   virtual double GetCoulombInt(OrbitalType orbital1, 
                                OrbitalType orbital2, 
                                Atom* atom); 
   virtual double GetExchangeInt(OrbitalType orbital1, 
                                 OrbitalType orbital2, 
                                 Atom* atom); 
   virtual double GetElectronCoreAttraction(Atom* atomA, Atom* atomB, 
                                            OrbitalType mu, OrbitalType nu,
                                            double**** twoElecTwoCoreMatrixTwoAtoms);
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              Molecule* molecule, 
                                              double** fockMatrix, 
                                              double** gammaAB);
   virtual void CalcCISMatrix(double** matrixCIS, int numberOcc, int numberVir);
   virtual void CalcForce(int electronicStateIndex);
private:
   string errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
   string errorMessageMultipoleA;
   string errorMessageMultipoleB;
   string errorMessageGetNddoRepulsionIntegral;
   string errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsNullMatrix;
   string errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsSameAtoms;
   void CalcTwoElecTwoCoreMatrixTwoAtoms(double**** matrix, int atomAIndex, int atomBIndex);
   void RotateTwoElecTwoCoreMatrixToSpaceFrame(double**** matrix, double** rotatingMatrix);
   double GetNddoRepulsionIntegral(Atom* atomA, OrbitalType mu, OrbitalType nu,
                                   Atom* atomB, OrbitalType lambda, OrbitalType sigma);
   double GetSemiEmpiricalMultipoleInteraction(MultipoleType multipoleA,
                                               MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab);
};

Mndo::Mndo() : MolDS_zindo::ZindoS(){
   this->theory = MNDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "Mndo created\n";
}

Mndo::~Mndo(){
   //cout << "Mndo deleted\n";
}

void Mndo::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in mndo::Mndo::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in mndo::Mndo::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in mndo::Mndo::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in mndo::Mndo::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_mndo::Mndo::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_mndo::Mndo::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in mndo::Mndo::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in mndo::Mndo::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageMultipoleA = "Multipole A is: ";
   this->errorMessageMultipoleB = "Multipole B is: ";
   this->errorMessageGetNddoRepulsionIntegral = "Error in mndo::Mndo::GetNddoRepulsionIntegral: Bad orbital is set.\n";

   this->errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsNullMatrix 
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreMatrixTwoAtoms: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsSameAtoms
      = "Error in mndo::Mndo::CalcTwoElecTwoCoreMatrixTwoAtoms: The matrix is NULL.\n"; 
   this->messageSCFMetConvergence = "\n\n\n\t\tMNDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: MNDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: MNDO/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: MNDO/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: MNDO/S-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for MNDO-CIS met convergence criterion(^^b\n\n\n";
}

void Mndo::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

void Mndo::CalcCoreRepulsionEnergy(){
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
      alphaA = atomA->GetMndoAlpha();
      for(int j=i+1; j<this->molecule->GetAtomVect()->size(); j++){
         atomB = (*this->molecule->GetAtomVect())[j];
         alphaB = atomB->GetMndoAlpha();
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
      }
   }
   this->coreRepulsionEnergy = energy;
}

double Mndo::GetFockDiagElement(Atom* atomA, int atomAIndex, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 double****** twoElecTwoCore, bool isGuess){
   double value=0.0;
   /*
   int firstAOIndexA = atomA->GetFirstAOIndex();
   value = atomA->GetCoreIntegral(atomA->GetValence()[mu-firstAOIndexA], 
                                  isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      int totalNumberAOs = this->molecule->GetTotalNumberAOs();
      double orbitalElectronPopulationDiagPart[totalNumberAOs];

      for(int i=0; i<totalNumberAOs; i++){
         orbitalElectronPopulationDiagPart[i] = orbitalElectronPopulation[i][i];
      }

      OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
      OrbitalType orbitalLam;
      int atomANumberValence = atomA->GetValence().size();
      for(int v=0; v<atomANumberValence; v++){
         orbitalLam = atomA->GetValence()[v];
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalLam, atomA);
         exchange = this->GetExchangeInt(orbitalMu, orbitalLam, atomA);
         lammda = v + firstAOIndexA;
         temp += orbitalElectronPopulationDiagPart[lammda]*(coulomb - 0.5*exchange);
      }
      value += temp;
   
      temp = 0.0;
      int totalNumberAtoms = molecule->GetAtomVect()->size();
      for(int B=0; B<totalNumberAtoms; B++){
         if(B != atomAIndex){
            Atom* atomB = (*molecule->GetAtomVect())[B];
            OrbitalType orbitalSigma;
            int sigma;
            int atomBNumberValence = atomB->GetValence().size();
            for(int i=0; i<atomBNumberValence; i++){
               sigma = i + atomB->GetFirstAOIndex();
               orbitalSigma = atomB->GetValence()[i];
               temp += orbitalElectronPopulationDiagPart[sigma]
                      *this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalSigma);
            }
            temp -= atomB->GetCoreCharge() 
                   *this->GetNishimotoMatagaTwoEleInt(atomA, s, atomB, s);
         }
      }
      value += temp;
   }
   */
   return value;
}

double Mndo::GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                    int mu, int nu, Molecule* molecule, double** gammaAB, double** overlap,
                                    double** orbitalElectronPopulation, 
                                    double****** twoElecTwoCore, bool isGuess){
   
   double value = 0.0;
   /*
   OrbitalType orbitalMu = atomA->GetValence()[mu-atomA->GetFirstAOIndex()];
   OrbitalType orbitalNu = atomB->GetValence()[nu-atomB->GetFirstAOIndex()];

   double bondParameter = 0.5*(atomA->GetBondingParameter(this->theory, orbitalMu) 
                              +atomB->GetBondingParameter(this->theory, orbitalNu)); 

   if(isGuess){
      value = bondParameter*overlap[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(atomAIndex == atomBIndex){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlap[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]
                  *this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
      }
   }
   */
   return value;
}

// MNDO Coulomb Interaction
double Mndo::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;
    
   if( orbital1 == s && orbital2 == s){ 
      value = atom->GetMndoGss();
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom->GetMndoGsp();
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = this->GetCoulombInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetMndoGpp();
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetMndoGpp2();
   }   
   
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( orbital1 == s && ( orbital2 == dxy || 
                               orbital2 == dyz || 
                               orbital2 == dzz || 
                               orbital2 == dzx || 
                               orbital2 == dxxyy )){ 
      value = atom->GetZindoF0sdLower();
   }   
   else if( orbital2 == s && ( orbital1 == dxy || 
                               orbital1 == dyz || 
                               orbital1 == dzz || 
                               orbital1 == dzx || 
                               orbital1 == dxxyy )){ 
      value = atom->GetZindoF0sdLower();
   }
   else if( orbital1 == dzz && (orbital2 == px || orbital2==py) ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*2.0;
   }
   else if( orbital2 == dzz && (orbital1 == px || orbital1==py) ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dzz && orbital2 == pz) ||
            (orbital2 == dzz && orbital1 == pz) ){
      value = atom->GetZindoF0sdLower()
             +atom->GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == orbital2) && ( orbital1 == dxy || 
                                        orbital1 == dyz || 
                                        orbital1 == dzz || 
                                        orbital1 == dzx || 
                                        orbital1 == dxxyy )){ 
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*36.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == px) ||
            (orbital2 == dxxyy && orbital1 == px) || 
            (orbital1 == dxxyy && orbital2 == py) ||
            (orbital2 == dxxyy && orbital1 == py) || 
            (orbital1 == dxy && orbital2 == px) ||
            (orbital2 == dxy && orbital1 == px) || 
            (orbital1 == dxy && orbital2 == py) ||
            (orbital2 == dxy && orbital1 == py) ||
            (orbital1 == dzx && orbital2 == px) ||
            (orbital2 == dzx && orbital1 == px) || 
            (orbital1 == dzx && orbital2 == pz) ||
            (orbital2 == dzx && orbital1 == pz) || 
            (orbital1 == dyz && orbital2 == py) ||
            (orbital2 == dyz && orbital1 == py) || 
            (orbital1 == dyz && orbital2 == pz) ||
            (orbital2 == dyz && orbital1 == pz) ){
      value = atom->GetZindoF0sdLower()
             +atom->GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == pz) ||
            (orbital2 == dxxyy && orbital1 == pz) || 
            (orbital1 == dxy && orbital2 == pz) ||
            (orbital2 == dxy && orbital1 == pz) ||
            (orbital1 == dzx && orbital2 == py) ||
            (orbital2 == dzx && orbital1 == py) ||
            (orbital1 == dyz && orbital2 == px) ||
            (orbital2 == dyz && orbital1 == px)  ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzz) ||
            (orbital2 == dxxyy && orbital1 == dzz) || 
            (orbital1 == dxy && orbital2 == dzz) ||
            (orbital2 == dxy && orbital1 == dzz)  ){
      value = atom->GetZindoF0ddLower()
             -atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*6.0;
   }
   else if( (orbital1 == dxy && orbital2 == dxxyy) ||
            (orbital2 == dxy && orbital1 == dxxyy)  ){
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*4.0
             -atom->GetZindoF4ddLower()*34.0;
   }
   else if( (orbital1 == dzx && orbital2 == dzz) ||
            (orbital2 == dzx && orbital1 == dzz) || 
            (orbital1 == dyz && orbital2 == dzz) ||
            (orbital2 == dyz && orbital1 == dzz)  ){
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*2.0
             -atom->GetZindoF4ddLower()*24.0;
   }
   else if( (orbital1 == dzx && orbital2 == dxxyy) ||
            (orbital2 == dzx && orbital1 == dxxyy) || 
            (orbital1 == dzx && orbital2 == dxy) || 
            (orbital2 == dzx && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dxxyy) || 
            (orbital2 == dyz && orbital1 == dxxyy) || 
            (orbital1 == dyz && orbital2 == dxy) || 
            (orbital2 == dyz && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dzx) || 
            (orbital2 == dyz && orbital1 == dzx) ){
      value = atom->GetZindoF0ddLower()
             -atom->GetZindoF2ddLower()*2.0
             -atom->GetZindoF4ddLower()*4.0;
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageCoulombInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   
   return value;

}

// MNDO Exchange Interaction
double Mndo::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;
   
   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetMndoHsp();
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = this->GetExchangeInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetMndoHpp();
   }
   
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( (orbital1 == s) && (orbital2 == dxy || 
                                orbital2 == dyz || 
                                orbital2 == dzz || 
                                orbital2 == dzx || 
                                orbital2 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital2 == s) && (orbital1 == dxy || 
                                orbital1 == dyz || 
                                orbital1 == dzz || 
                                orbital1 == dzx || 
                                orbital1 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital1 == px && orbital2 == dzz) ||
            (orbital2 == px && orbital1 == dzz) ||
            (orbital1 == py && orbital2 == dzz) ||
            (orbital2 == py && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()
             +atom->GetZindoG3pdLower()*18.0;
   }
   else if( (orbital1 == px && orbital2 == dxxyy) ||
            (orbital2 == px && orbital1 == dxxyy) ||
            (orbital1 == px && orbital2 == dxy) ||
            (orbital2 == px && orbital1 == dxy) ||
            (orbital1 == px && orbital2 == dzx) ||
            (orbital2 == px && orbital1 == dzx) ||
            (orbital1 == py && orbital2 == dxxyy) ||
            (orbital2 == py && orbital1 == dxxyy) ||
            (orbital1 == py && orbital2 == dxy) ||
            (orbital2 == py && orbital1 == dxy) ||
            (orbital1 == py && orbital2 == dyz) ||
            (orbital2 == py && orbital1 == dyz) ||
            (orbital1 == pz && orbital2 == dzx) ||
            (orbital2 == pz && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dyz) ||
            (orbital2 == pz && orbital1 == dyz) ){
      value = atom->GetZindoG1pdLower()*3.0
             +atom->GetZindoG3pdLower()*24.0;
   }
   else if( (orbital1 == px && orbital2 == dyz) ||
            (orbital2 == px && orbital1 == dyz) ||
            (orbital1 == py && orbital2 == dzx) ||
            (orbital2 == py && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dxxyy) ||
            (orbital2 == pz && orbital1 == dxxyy) ||
            (orbital1 == pz && orbital2 == dxy) ||
            (orbital2 == pz && orbital1 == dxy) ){
      value = atom->GetZindoG3pdLower()*15.0;
   }
   else if( (orbital1 == pz && orbital2 == dzz) ||
            (orbital2 == pz && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()*4.0
             +atom->GetZindoG3pdLower()*27.0;
   }
   else if( (orbital1 == dzz && orbital2 == dxxyy) ||
            (orbital2 == dzz && orbital1 == dxxyy) ||
            (orbital1 == dzz && orbital2 == dxy) ||
            (orbital2 == dzz && orbital1 == dxy) ){
      value = atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*15.0;
   }
   else if( (orbital1 == dzz && orbital2 == dzx) ||
            (orbital2 == dzz && orbital1 == dzx) ||
            (orbital1 == dzz && orbital2 == dyz) ||
            (orbital2 == dzz && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()
             +atom->GetZindoF4ddLower()*30.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dxy) ||
            (orbital2 == dxxyy && orbital1 == dxy) ){
      value = atom->GetZindoF4ddLower()*35.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzx) ||
            (orbital2 == dxxyy && orbital1 == dzx) ||
            (orbital1 == dxxyy && orbital2 == dyz) ||
            (orbital2 == dxxyy && orbital1 == dyz) ||
            (orbital1 == dxy && orbital2 == dzx) ||
            (orbital2 == dxy && orbital1 == dzx) ||
            (orbital1 == dxy && orbital2 == dyz) ||
            (orbital2 == dxy && orbital1 == dyz) ||
            (orbital1 == dzx && orbital2 == dyz) ||
            (orbital2 == dzx && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()*3.0
             +atom->GetZindoF4ddLower()*20.0;
   }
   */
   
   else{
      stringstream ss;
      ss << this->errorMessageExchangeInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   
   return value;
}

// electron in atom A (mu and nu) and core (atom B) attraction. 
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttraction(Atom* atomA, Atom* atomB, 
                                          OrbitalType mu, OrbitalType nu,
                                          double**** twoElecTwoCoreMatrixTwoAtoms){
   return -1.0*atomB->GetCoreCharge()*twoElecTwoCoreMatrixTwoAtoms[mu][nu][s][s];
}

void Mndo::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, Atom* atomA, Atom* atomB){

   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         printf("diatomicOverlap[%d][%d]=%lf\n",i,j,diatomicOverlap[i][j]);
      }
   }
   */

}

void Mndo::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                                                double** diatomicOverlapDeri, 
                                                Atom* atomA, Atom* atomB){

   MolDS_cndo::Cndo2::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                        diatomicOverlapDeri,atomA, atomB);

}
// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double Mndo::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                          Molecule* molecule, double** fockMatrix, double** gammaAB){
   double value = 0.0;
   return value;
}

void Mndo::CalcCISMatrix(double** matrixCIS, int numberOcc, int numberVir){
}

// electronicStateIndex is index of the electroinc eigen state.
// "electronicStateIndex = 0" means electronic ground state. 
void Mndo::CalcForce(int electronicStateIndex){
}

// Calculation of two electrons two cores integral (mu, nu | lambda, sigma), 
// taht is, Eq. (9) in ref. [DT_1977-2].
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy] can be treatable.
void Mndo::CalcTwoElecTwoCoreMatrixTwoAtoms(double**** matrix, int atomAIndex, int atomBIndex){

   Atom* atomA = NULL;
   Atom* atomB = NULL;
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB->GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }
   else{
      atomA = (*this->molecule->GetAtomVect())[atomAIndex];
      atomB = (*this->molecule->GetAtomVect())[atomBIndex];
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreMatrixTwoAtomsNullMatrix;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB->GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix4d(matrix, dxy, dxy, dxy, dxy);
   } 

   // calclation in diatomic frame
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = this->GetNddoRepulsionIntegral(
                                               atomA, 
                                               atomA->GetValence()[mu],
                                               atomA->GetValence()[nu],
                                               atomB, 
                                               atomB->GetValence()[lambda],
                                               atomB->GetValence()[sigma]);
                     
            }
         }
      }
   }

   // rotate matirix into the space frame
   double** rotatingMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(
                                             OrbitalType_end, OrbitalType_end);
   try{
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->RotateTwoElecTwoCoreMatrixToSpaceFrame(matrix, rotatingMatrix);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&rotatingMatrix, OrbitalType_end);

   /*
   printf("(mu, nu | lambda, sigma) matrix\n"); 
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               printf("mu=%d nu=%d lambda=%d sigma=%d $e\n",
                        mu,nu,lambda,sigma,matrix[mu][nu][lambda][sigma]);
            }
         }
      }
   }
   */

}

// Rotate 4-d matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateTwoElecTwoCoreMatrixToSpaceFrame(double**** matrix, double** rotatingMatrix){

   double oldMatrix[dxy][dxy][dxy][dxy];
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               oldMatrix[mu][nu][lambda][sigma] = matrix[mu][nu][lambda][sigma];
            }
         }
      }
   }
   
   // rotate
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = 0.0;

               for(int i=0; i<dxy; i++){
                  for(int j=0; j<dxy; j++){
                     for(int k=0; k<dxy; k++){
                        for(int l=0; l<dxy; l++){
                           matrix[mu][nu][lambda][sigma] += oldMatrix[i][j][k][l] 
                                                            *rotatingMatrix[mu][i] 
                                                            *rotatingMatrix[nu][j] 
                                                            *rotatingMatrix[lambda][k] 
                                                            *rotatingMatrix[sigma][l];
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral(Atom* atomA, OrbitalType mu, OrbitalType nu,
                                      Atom* atomB, OrbitalType lambda, OrbitalType sigma){
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      value = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(0);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(0);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(0);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(1);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(1);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(1);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      DA = atomA->GetMndoDerivedParameterD(2);
      DB = atomB->GetMndoDerivedParameterD(2);
      rhoA = atomA->GetMndoDerivedParameterRho(2);
      rhoB = atomB->GetMndoDerivedParameterRho(2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, mu, mu)
                  -this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, nu, nu));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      // ToDo: error log.
      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegral;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB->GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteraction(MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab){
   double value = 0.0;
   double a = rhoA + rhoB;

   if(multipoleA == sQ && multipoleB == sQ){
      value = pow(pow(Rab,2.0) + pow(a,2.0), -0.5);
   }
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(Rab+DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleA, Qxx,
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/2.0 + pow(temp3,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(Rab,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(Rab-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(Rab+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(Rab+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             +pow(temp5,-0.5)/4.0 - pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0
             +pow(temp5,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(Rab-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 - pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/4.0 + pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, multipoleB, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/16.0 + pow(temp2,-0.5)/16.0 
             +pow(temp3,-0.5)/16.0 + pow(temp4,-0.5)/16.0
             -pow(temp5,-0.5)/8.0 - pow(temp6,-0.5)/8.0
             -pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0;
             +pow(temp9,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(Rab-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 - pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/8.0 + pow(temp6,-0.5)/8.0
             +pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(Rab,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/2.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

}
#endif



