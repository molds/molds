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
   virtual double GetFockDiagElement(Atom* atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     Molecule* molecule, 
                                     double** gammaAB,
                                     double** orbitalElectronPopulation, 
                                     double* atomicElectronPopulation,
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
   double GetNddoRepulsionIntegral(Atom* atomA, OrbitalType Mu, OrbitalType Nu,
                                   Atom* atomB, OrbitalType Lambda, OrbitalType Sigma);
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

double Mndo::GetFockDiagElement(Atom* atomA, int atomAIndex, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 bool isGuess){
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
                                    double** orbitalElectronPopulation, bool isGuess){
   
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

// See Apendix in [DT_1977]
// Orbital Mu and Nu belong atom A, 
// orbital Lambda and Sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral(Atom* atomA, OrbitalType Mu, OrbitalType Nu,
                                      Atom* atomB, OrbitalType Lambda, OrbitalType Sigma){
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(Mu == s && Nu == s && Lambda == s && Sigma == s){
   }
   // (29) in [DT_1977]
   else if(Mu == s && Nu == s && Lambda == px && Sigma == px){
   }
   else if(Mu == s && Nu == s && Lambda == py && Sigma == py){
   }
   // (30) in [DT_1977]
   else if(Mu == s && Nu == s && Lambda == pz && Sigma == pz){
   }
   // (31) in [DT_1977]
   else if(Mu == px && Nu == px && Lambda == s && Sigma == s){
   }
   else if(Mu == py && Nu == py && Lambda == s && Sigma == s){
   }
   // (32) in [DT_1977]
   else if(Mu == pz && Nu == pz && Lambda == s && Sigma == s){
   }
   // (33) in [DT_1977]
   else if(Mu == px && Nu == px && Lambda == px && Sigma == px){
   }
   else if(Mu == py && Nu == py && Lambda == py && Sigma == py){
   }
   // (34) in [DT_1977]
   else if(Mu == px && Nu == px && Lambda == py && Sigma == py){
   }
   else if(Mu == py && Nu == py && Lambda == px && Sigma == px){
   }
   // (35) in [DT_1977]
   else if(Mu == px && Nu == px && Lambda == pz && Sigma == pz){
   }
   else if(Mu == py && Nu == py && Lambda == pz && Sigma == pz){
   }
   // (36) in [DT_1977]
   else if(Mu == pz && Nu == pz && Lambda == px && Sigma == px){
   }
   else if(Mu == pz && Nu == pz && Lambda == py && Sigma == py){
   }
   // (37) in [DT_1977]
   else if(Mu == pz && Nu == pz && Lambda == pz && Sigma == pz){
   }
   // (38) in [DT_1977]
   else if(Mu == s && Nu == pz && Lambda == s && Sigma == s){
   }
   else if(Mu == pz && Nu == s && Lambda == s && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   // (39) in [DT_1977]
   else if(Mu == s && Nu == pz && Lambda == px && Sigma == px){
   }
   else if(Mu == pz && Nu == s && Lambda == px && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == pz && Lambda == py && Sigma == py){
   }
   else if(Mu == pz && Nu == s && Lambda == py && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   // (40) in [DT_1977]
   else if(Mu == s && Nu == pz && Lambda == pz && Sigma == pz){
   }
   else if(Mu == pz && Nu == s && Lambda == pz && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   // (41) in [DT_1977]
   else if(Mu == s && Nu == s && Lambda == s && Sigma == pz){
   }
   else if(Mu == s && Nu == s && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   // (42) in [DT_1977]
   else if(Mu == px && Nu == px && Lambda == s && Sigma == pz){
   }
   else if(Mu == px && Nu == px && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == py && Lambda == s && Sigma == pz){
   }
   else if(Mu == py && Nu == py && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   // (43) in [DT_1977]
   else if(Mu == pz && Nu == pz && Lambda == s && Sigma == pz){
   }
   else if(Mu == pz && Nu == pz && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   // (44) in [DT_1977]
   else if(Mu == s && Nu == px && Lambda == s && Sigma == px){
   }
   else if(Mu == px && Nu == s && Lambda == s && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == px && Lambda == px && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == px && Nu == s && Lambda == px && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   else if(Mu == s && Nu == py && Lambda == s && Sigma == py){
   }
   else if(Mu == py && Nu == s && Lambda == s && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == py && Lambda == py && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == s && Lambda == py && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // (45) in [DT_1977]
   else if(Mu == s && Nu == pz && Lambda == s && Sigma == pz){
   }
   else if(Mu == pz && Nu == s && Lambda == s && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == pz && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == pz && Nu == s && Lambda == pz && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // (46) in [DT_1977]
   else if(Mu == s && Nu == px && Lambda == px && Sigma == pz){
   }
   else if(Mu == px && Nu == s && Lambda == px && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == px && Lambda == pz && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == px && Nu == s && Lambda == pz && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   else if(Mu == s && Nu == py && Lambda == py && Sigma == pz){
   }
   else if(Mu == py && Nu == s && Lambda == py && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == s && Nu == py && Lambda == pz && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == s && Lambda == pz && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // (47) in [DT_1977]
   else if(Mu == px && Nu == pz && Lambda == s && Sigma == px){
   }
   else if(Mu == pz && Nu == px && Lambda == s && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == px && Nu == pz && Lambda == px && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == pz && Nu == px && Lambda == px && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == pz && Lambda == s && Sigma == py){
   }
   else if(Mu == pz && Nu == py && Lambda == s && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == py && Nu == pz && Lambda == py && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == pz && Nu == py && Lambda == py && Sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // (48) in [DT_1977]
   else if(Mu == px && Nu == pz && Lambda == px && Sigma == pz){
   }
   else if(Mu == pz && Nu == px && Lambda == px && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == px && Nu == pz && Lambda == pz && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == pz && Nu == px && Lambda == pz && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == pz && Lambda == py && Sigma == pz){
   }
   else if(Mu == pz && Nu == py && Lambda == py && Sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == py && Nu == pz && Lambda == pz && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == pz && Nu == py && Lambda == pz && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(Mu == px && Nu == py && Lambda == px && Sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral(atomA, Mu, Mu, atomB, Mu, Mu)
                  -this->GetNddoRepulsionIntegral(atomA, Mu, Mu, atomB, Nu, Nu));
   }
   else if(Mu == py && Nu == px && Lambda == px && Sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Lambda, Sigma);
   }
   else if(Mu == px && Nu == py && Lambda == py && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Mu, Nu, atomB, Sigma, Lambda);
   }
   else if(Mu == py && Nu == px && Lambda == py && Sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, Nu, Mu, atomB, Sigma, Lambda);
   }
   // d-orbitals
   else if(Mu == dxy || Mu == dyz || Mu == dzz || Mu == dzx || Mu == dxxyy ||
           Nu == dxy || Nu == dyz || Nu == dzz || Nu == dzx || Nu == dxxyy ||
           Lambda == dxy || Lambda == dyz || Lambda == dzz || Lambda  == dzx || Lambda == dxxyy ||
           Sigma == dxy || Sigma == dyz || Sigma == dzz || Sigma  == dzx || Sigma == dxxyy){

      // ToDo: error log.
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

}
#endif



