#ifndef INCLUDED_ZINDOS
#define INCLUDED_ZINDOS

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_zindo{

/***
 *  Main Refference for Zindo is [RZ_1973]
 */
class ZindoS : public MolDS_cndo::Cndo2{
public:
   ZindoS();
   ~ZindoS();
   virtual void DoesCIS();
protected:
   virtual void CalcGammaAB(double** gammaAB, Molecule* molecule);
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
   double** matrixCIS;
   double* excitedEnergies;
   int matrixCISdimension;
   double GetCoulombInt(OrbitalType orbital1, 
                        OrbitalType orbital2, 
                        Atom* atom); // Apendix in [BZ_1979]
   double GetExchangeInt(OrbitalType orbital1, 
                         OrbitalType orbital2, 
                         Atom* atom); // Apendix in [BZ_1979]
   double GetNishimotoMatagaTwoEleInt(Atom* atomA, OrbitalType orbitalA, 
                                      Atom* atomB, OrbitalType orbitalB); // ref. [MN_1957] and (5a) in [AEZ_1986]
   double GetNishimotoMatagaTwoEleIntFirstDerivative(Atom* atomA, OrbitalType orbitalA, 
                                                     Atom* atomB, OrbitalType orbitalB,
                                                     CartesianType axisA);// ref. [MN_1957] and (5a) in [AEZ_1986]
   double nishimotoMatagaParamA;
   double nishimotoMatagaParamB;
   double overlapCorrectionSigma;
   double overlapCorrectionPi;
   void DoesCISDirect();
   void DoesCISDavidson();
   void CalcRitzVector(double* ritzVector, 
                       double** expansionVectors, 
                       double** interactionMatrix, 
                       int interactionMatrixDimension, 
                       int ritzVectorIndex);
   void CalcResidualVectorAndNorm(double* residualVector, 
                                  double* norm, 
                                  double* ritzVector, 
                                  double* interactionEigenEnergies, 
                                  int residualVectorIndex);
   void SortSingleExcitationSlaterDeterminants(vector<MoEnergy>* moEnergies);
   void UpdateExpansionVectors(double** expansionVectors, 
                               double* interactionEigenEnergies, 
                               double* residualVector,
                               int interactionMatrixDimension, 
                               int* notConvergedStates, 
                               int residualVectorIndex);
   void CalcInteractionMatrix(double** interactionMatrix, 
                              double** expansionVectors, 
                              int interactionMatrixDimension);
   void FreeDavidsonCISTemporaryMtrices(double*** expansionVectors, 
                                        double** residualVector, 
                                        double** ritzVector);
   void FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix, 
                                            double interactionMatrixDimension, 
                                            double** interactionEigenEnergies);
   string errorMessageNishimotoMataga;
   string errorMessageDavidsonNotConverged;
   string errorMessageDavidsonMaxIter;
   string errorMessageDavidsonMaxDim;
   string errorMessageCalcCISMatrix;
   string messageStartCalcCISMatrix;
   string messageOmpElapsedTimeCalcCISMarix;
   string messageOmpElapsedTimeCIS;
   string messageDoneCalcCISMatrix;
   string messageStartCIS;
   string messageDoneCIS;
   string messageStartDirectCIS;
   string messageDoneDirectCIS;
   string messageStartDavidsonCIS;
   string messageDoneDavidsonCIS;
   string messageNumIterCIS;
   string messageResidualNorm;
   string messageDavidsonConverge;
   string messageDavidsonReachCISMatrix;
   string messageDavidsonGoToDirect;
   string messageExcitedStatesEnergies;
   string messageExcitedStatesEnergiesTitle;
};

ZindoS::ZindoS() : MolDS_cndo::Cndo2(){
   this->theory = ZINDOS;
   this->SetMessages();
   this->SetEnableAtomTypes();
   this->matrixCIS = NULL;
   this->excitedEnergies = NULL;
   this->matrixCISdimension = 0;
   this->nishimotoMatagaParamA = 1.2;
   this->nishimotoMatagaParamB = 2.4;
   this->matrixForce = NULL;
   this->overlapCorrectionSigma = 1.267;
   this->overlapCorrectionPi = 0.585;
   //cout << "ZindoS created\n";
}

ZindoS::~ZindoS(){
   if(this->matrixCIS != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&this->matrixCIS, this->matrixCISdimension);
      //cout << "matrixCIS deleted\n";
   }
   if(this->excitedEnergies != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->excitedEnergies);
      //cout << "exceitedEnergies deleted\n";
   }
   if(this->matrixForce != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(&this->matrixForce, this->molecule->GetAtomVect()->size());
      //cout << "matrixForce deleted\n";
   }
   //cout << "ZindoS deleted\n";
}

void ZindoS::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in zindo::ZindoS::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in zindo::ZindoS::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in zindo::ZindoS::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in zindo::ZindoS::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_zindo::ZindoS::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_zindo::ZindoS::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageNishimotoMataga = "Error in base_zindo::ZindoS::GetNishimotoMatagaTwoEleInt: Invalid orbitalType.\n";
   this->errorMessageMolecularIntegralElement
      = "Error in zindo::ZindoS::GetMolecularIntegralElement: Non available orbital is contained.\n";
   this->errorMessageCalcCISMatrix
      = "Error in zindo::ZindoS::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in zindo::ZindoS::DoesCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageDavidsonMaxIter = "Davidson roop reaches max_iter=";
   this->errorMessageDavidsonMaxDim = "Dimension of the expansion vectors reaches max_dim=";
   this->messageSCFMetConvergence = "\n\n\n\t\tZINDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: ZINDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: ZINDO/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: ZINDO/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: ZINDO/S-CIS  **********\n\n\n";
   this->messageOmpElapsedTimeCalcCISMarix = "\tElapsed time(omp) for the calc. of the CIS matrix = ";
   this->messageOmpElapsedTimeCIS = "\tElapsed time(omp) for the CIS = ";
   this->messageStartCalcCISMatrix = "----------- START: Calculation of the CIS matrix -----------\n";
   this->messageDoneCalcCISMatrix  = "----------- DONE: Calculation of the CIS matrix -----------\n\n";
   this->messageStartDirectCIS = "\t======  START: Direct-CIS  =====\n\n";
   this->messageDoneDirectCIS =  "\t======  DONE: Direct-CIS  =====\n\n\n";
   this->messageStartDavidsonCIS = "\t======  START: Davidson-CIS  =====\n";
   this->messageDoneDavidsonCIS =  "\t======  DONE: Davidson-CIS  =====\n\n\n";
   this->messageNumIterCIS = "\tDavidson iter=";
   this->messageResidualNorm = "-th excited: norm of the residual = ";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for ZINDO/S-CIS met convergence criterion(^^b\n\n\n";
   this->messageDavidsonReachCISMatrix = "\n\t\tDimension of the expansion vectors reaches to the dimension of the CIS-matrix.\n";
   this->messageDavidsonGoToDirect = "\t\tHence, we go to the Direct-CIS.\n\n";
   this->messageExcitedStatesEnergies = "\tExcitation energies:\n";
   this->messageExcitedStatesEnergiesTitle = "\t\t| i-th | e[a.u.] | e[eV] | \n";
}

void ZindoS::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double ZindoS::GetFockDiagElement(Atom* atomA, int atomAIndex, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 bool isGuess){
   double value=0.0;
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

   return value;
}

double ZindoS::GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                    int mu, int nu, Molecule* molecule, double** gammaAB, double** overlap,
                                    double** orbitalElectronPopulation, bool isGuess){
   double value = 0.0;
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
   
   return value;
}

void ZindoS::CalcGammaAB(double** gammaAB, Molecule* molecule){
   // Do nothing;
}

// Apendix in [BZ_1972]
// ZINDO Coulomb Interaction
double ZindoS::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;
   
   if( orbital1 == s && orbital2 == s){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetZindoF0ssLower()
             +atom->GetZindoF2ppLower()*4.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoF0ssLower()
             -atom->GetZindoF2ppLower()*2.0;
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

// Apendix in [BZ_1972]
// ZINDO Exchange Interaction
double ZindoS::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;

   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoG1spLower();
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = atom->GetZindoG1spLower();
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoF2ppLower()*3.0;
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

// ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt(Atom* atomA, OrbitalType orbitalA, 
                                           Atom* atomB, OrbitalType orbitalB){
   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   double gammaAA;
   if(orbitalA == s || 
      orbitalA == px ||
      orbitalA == py ||
      orbitalA == pz ){
      gammaAA = atomA->GetZindoF0ss();
   }
   /*
   else if(orbitalA == dxy ||
           orbitalA == dyz ||
           orbitalA == dzz ||
           orbitalA == dzx ||
           orbitalA == dxxyy ){
      gammaAA = atomA->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalA) << "\n";
      throw MolDSException(ss.str());
   }   

   double gammaBB;
   if(orbitalB == s || 
      orbitalB == px ||
      orbitalB == py ||
      orbitalB == pz ){
      gammaBB = atomB->GetZindoF0ss();
   }
   /*
   else if(orbitalB == dxy ||
           orbitalB == dyz ||
           orbitalB == dzz ||
           orbitalB == dzx ||
           orbitalB == dxxyy ){
      gammaBB = atomB->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomB->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalB) << "\n";
      throw MolDSException(ss.str());
   }  

   return this->nishimotoMatagaParamA/( r+this->nishimotoMatagaParamB/(gammaAA+gammaBB) );

}

// First derivative of Nishimoto-Mataga related to the coordinate of atom A.
// For Nishimoto-Mataga, See ZindoS::GetNishimotoMatagaTwoEleInt 
// or ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleIntFirstDerivative(Atom* atomA, OrbitalType orbitalA, 
                                                          Atom* atomB, OrbitalType orbitalB,
                                                          CartesianType axisA){
   double gammaAA;
   if(orbitalA == s || 
      orbitalA == px ||
      orbitalA == py ||
      orbitalA == pz ){
      gammaAA = atomA->GetZindoF0ss();
   }
   /*
   else if(orbitalA == dxy ||
           orbitalA == dyz ||
           orbitalA == dzz ||
           orbitalA == dzx ||
           orbitalA == dxxyy ){
      gammaAA = atomA->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalA) << "\n";
      throw MolDSException(ss.str());
   }   

   double gammaBB;
   if(orbitalB == s || 
      orbitalB == px ||
      orbitalB == py ||
      orbitalB == pz ){
      gammaBB = atomB->GetZindoF0ss();
   }
   /*
   else if(orbitalB == dxy ||
           orbitalB == dyz ||
           orbitalB == dzz ||
           orbitalB == dzx ||
           orbitalB == dxxyy ){
      gammaBB = atomB->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomB->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalB) << "\n";
      throw MolDSException(ss.str());
   }  

   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dCartesian = atomA->GetXyz()[axisA] - atomB->GetXyz()[axisA];
   double value = -1.0*dCartesian/r;
   value *= this->nishimotoMatagaParamA;
   value *= pow( r+this->nishimotoMatagaParamB/(gammaAA+gammaBB) ,-2.0);
   return value;

}

void ZindoS::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, Atom* atomA, Atom* atomB){

   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);

   // see (4f) in [AEZ_1986]
   diatomicOverlap[pz][pz] *= this->overlapCorrectionSigma;
   diatomicOverlap[py][py] *= this->overlapCorrectionPi;
   diatomicOverlap[px][px] *= this->overlapCorrectionPi;
   
   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         printf("diatomicOverlap[%d][%d]=%lf\n",i,j,diatomicOverlap[i][j]);
      }
   }
   */

}

void ZindoS::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                                                double** diatomicOverlapDeri, 
                                                Atom* atomA, Atom* atomB){

   MolDS_cndo::Cndo2::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                        diatomicOverlapDeri,atomA, atomB);

   // see (4f) in [AEZ_1986] like as overlap integlral
   diatomicOverlapDeri[pz][pz] *= this->overlapCorrectionSigma;
   diatomicOverlapDeri[py][py] *= this->overlapCorrectionPi;
   diatomicOverlapDeri[px][px] *= this->overlapCorrectionPi;

}
// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double ZindoS::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                          Molecule* molecule, double** fockMatrix, double** gammaAB){
   double value = 0.0;
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   int firstAOIndexA;
   int firstAOIndexB;
   int numberAOsA;
   int numberAOsB;
   double gamma;
   double exchange;
   double coulomb;
   OrbitalType orbitalMu;
   OrbitalType orbitalNu;

   for(int A=0; A<molecule->GetAtomVect()->size(); A++){
      atomA = (*molecule->GetAtomVect())[A];
      firstAOIndexA = atomA->GetFirstAOIndex();
      numberAOsA = atomA->GetValence().size();

      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         orbitalMu = atomA->GetValence()[mu-firstAOIndexA];

         // CNDO term
         for(int B=A; B<molecule->GetAtomVect()->size(); B++){
            atomB = (*molecule->GetAtomVect())[B];
            firstAOIndexB = atomB->GetFirstAOIndex();
            numberAOsB = atomB->GetValence().size();

            for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
               orbitalNu = atomB->GetValence()[nu-firstAOIndexB];

               if(A<B){
                  gamma = this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
                  value += gamma*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];
                  value += gamma*fockMatrix[moI][nu]*fockMatrix[moJ][nu]*fockMatrix[moK][mu]*fockMatrix[moL][mu];
               }
               else{
                  gamma = atomA->GetZindoF0ss();
                  value += gamma*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];
               }  

            }
         }

         // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
         for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
            orbitalNu = atomA->GetValence()[nu-firstAOIndexA];

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
               value += exchange*fockMatrix[moI][mu]*fockMatrix[moJ][nu]*fockMatrix[moK][nu]*fockMatrix[moL][mu];
               value += exchange*fockMatrix[moI][mu]*fockMatrix[moJ][nu]*fockMatrix[moK][mu]*fockMatrix[moL][nu];
            }

            coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

            if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                  gamma = atomA->GetZindoF0ss();
            }
            else{
               stringstream ss;
               ss << this->errorMessageMolecularIntegralElement;
               ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
               throw MolDSException(ss.str());
            }   

            value += (coulomb-gamma)*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];

         }

      }
   }

   return value;
}

void ZindoS::DoesCIS(){
   cout << this->messageStartCIS;
   double ompStartTime = omp_get_wtime();

   int numberOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberVir = Parameters::GetInstance()->GetActiveVirCIS();

   // malloc or initialize  CIS matrix
   if(this->matrixCIS == NULL){
      this->matrixCISdimension = numberOcc*numberVir;
      this->matrixCIS = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(this->matrixCISdimension, 
                                                                           this->matrixCISdimension);
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix2d(this->matrixCIS, 
                                                             this->matrixCISdimension, 
                                                             this->matrixCISdimension);
   }

   // malloc or initialize CIS eigen vector
   if(this->excitedEnergies == NULL){
      this->excitedEnergies = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(this->matrixCISdimension);
   }
   else{
      MallocerFreer::GetInstance()->InitializeDoubleMatrix1d(this->excitedEnergies, this->matrixCISdimension);
   }

   // calculate CIS matrix
   this->CalcCISMatrix(matrixCIS, numberOcc, numberVir);

   // calculate excited energies
   if(Parameters::GetInstance()->IsDavidsonCIS()){
      this->DoesCISDavidson();
   }
   else{
      this->DoesCISDirect();
   }

   // output eigen energies
   cout << this->messageExcitedStatesEnergies;
   cout << this->messageExcitedStatesEnergiesTitle;
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      printf("\t\t %d\t%e\t%e\n",k+1, this->excitedEnergies[k], this->excitedEnergies[k]/eV2AU);
   }
   cout << endl;
   double ompEndTime = omp_get_wtime();
   cout << this->messageOmpElapsedTimeCIS;
   cout << ompEndTime - ompStartTime;
   cout << this->messageUnitSec << endl;
   cout << this->messageDoneCIS;
}

void ZindoS::SortSingleExcitationSlaterDeterminants(vector<MoEnergy>* moEnergies){

   int numberOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberVir = Parameters::GetInstance()->GetActiveVirCIS();
   for(int k=0; k<numberOcc*numberVir; k++){
      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->molecule->GetTotalNumberValenceElectrons()/2 - (k/numberVir) -1;
      int moA = this->molecule->GetTotalNumberValenceElectrons()/2 + (k%numberVir);

      MoEnergy moEnergy = {this->energiesMO[moA] - this->energiesMO[moI], moI, moA, k};
      moEnergies->push_back(moEnergy);
   }
   sort(moEnergies->begin(), moEnergies->end(), LessMoEnergy());
}

// This method is used for Davidson
void ZindoS::CalcRitzVector(double* ritzVector, double** expansionVectors, double** interactionMatrix, 
                            int interactionMatrixDimension, int ritzVectorIndex){

   for(int j=0; j<this->matrixCISdimension; j++){
      ritzVector[j] = 0.0;
      for(int k=0; k<interactionMatrixDimension; k++){
         ritzVector[j] += expansionVectors[j][k]*interactionMatrix[ritzVectorIndex][k];
      }
   }
}

// This method is used for Davidson
void ZindoS::CalcResidualVectorAndNorm(double* residualVector, double* norm, double* ritzVector, 
                                       double* interactionEigenEnergies, int residualVectorIndex){

   double sqNorm = 0.0;
   for(int j=0; j<this->matrixCISdimension; j++){
      residualVector[j] = interactionEigenEnergies[residualVectorIndex] * ritzVector[j];
      for(int k=0; k<this->matrixCISdimension; k++){
         double value = j<=k ? this->matrixCIS[j][k] : this->matrixCIS[k][j];
         residualVector[j] -= value*ritzVector[k];
      }
      sqNorm += pow(residualVector[j],2.0);
   }
   *norm = sqrt(sqNorm);
}

// This method is used for Davidson
void ZindoS::UpdateExpansionVectors(double** expansionVectors, double* interactionEigenEnergies, double* residualVector,
                                    int interactionMatrixDimension, int* notConvergedStates, int residualVectorIndex){

   double newExpansionVector[this->matrixCISdimension];

   // calculate new expansion vector from residual vector
   for(int j=0; j<this->matrixCISdimension; j++){
      double temp = interactionEigenEnergies[residualVectorIndex]-this->matrixCIS[j][j];
      if(temp == 0.0){
         // prevent dividing by 0.
         temp = pow(10,-100);
      }
      newExpansionVector[j]=pow(temp, -1.0)*residualVector[j];
   }

   // orthonormalize old expansion vectors and new expansion vector
   for(int k=0; k<interactionMatrixDimension+*notConvergedStates; k++){
      double overlap=0.0;
      for(int j=0; j<this->matrixCISdimension; j++){
         overlap += expansionVectors[j][k] * newExpansionVector[j];
      }
      for(int j=0; j<this->matrixCISdimension; j++){
         newExpansionVector[j] -= overlap*expansionVectors[j][k];
      }
   }

   // add new expansion vector to old expansion vectors
   double sqNormNewExpVect = 0.0;
   for(int j=0; j<this->matrixCISdimension; j++){
      sqNormNewExpVect += pow(newExpansionVector[j],2.0);
   }
   double normNewExpVect = sqrt(sqNormNewExpVect);
   for(int j=0; j<this->matrixCISdimension; j++){
      expansionVectors[j][interactionMatrixDimension+*notConvergedStates] 
            = pow(normNewExpVect,-1.0)*newExpansionVector[j];
   }
   *notConvergedStates += 1;
}

// This method is used for Davidson
void ZindoS::CalcInteractionMatrix(double** interactionMatrix, double** expansionVectors, int interactionMatrixDimension){

   #pragma omp parallel for schedule(auto)
   for(int k=0; k<interactionMatrixDimension*interactionMatrixDimension; k++){
      int i = k/interactionMatrixDimension;
      int j = k%interactionMatrixDimension;
      if(i<=j){
         for(int a=0; a<this->matrixCISdimension; a++){
            for(int b=0; b<this->matrixCISdimension; b++){
               double value = a<=b ? this->matrixCIS[a][b] : this->matrixCIS[b][a];
               interactionMatrix[i][j] += expansionVectors[a][i] 
                                          *value
                                          *expansionVectors[b][j];
            }
         }
      }
   }
}

void ZindoS::DoesCISDavidson(){
   cout << this->messageStartDavidsonCIS;

   int numberOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberExcitedStates = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   int maxIter = Parameters::GetInstance()->GetMaxIterationsCIS();
   int maxDim  = Parameters::GetInstance()->GetMaxDimensionsCIS();
   double normTol = Parameters::GetInstance()->GetNormToleranceCIS();
   bool convergeExcitedStates[numberExcitedStates];
   int interactionMatrixDimension;
   bool reachMaxDim;
   bool allConverged;
   int notConvergedStates;
   bool goToDirectCIS;
   double** expansionVectors = NULL;
   double*  ritzVector = NULL;
   double*  residualVector = NULL;

   // malloc
   expansionVectors = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(this->matrixCISdimension, maxDim);
   ritzVector = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(this->matrixCISdimension);
   residualVector = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(this->matrixCISdimension);

   try{
      // sort single excitation slater determinants
      vector<MoEnergy> moEnergies;
      this->SortSingleExcitationSlaterDeterminants(&moEnergies);

      // set initial expansion vectors and initial conveged vectors
      for(int k=0; k<numberExcitedStates; k++){
         expansionVectors[moEnergies[k].slaterIndex][k] = 1.0;
         convergeExcitedStates[k] = false;
      }

      interactionMatrixDimension = 0;
      reachMaxDim = false;
      goToDirectCIS = false;
      // Davidson roop
      for(int k=0; k<maxIter; k++){
         cout << messageNumIterCIS << k << endl;

         // calculate dimension of the interaction matrix (= number of the expansion vectors).
         for(int i=0; i<numberExcitedStates; i++){
            if(!convergeExcitedStates[i]){
               interactionMatrixDimension += 1;
            }
         }

         // malloc interaction matrix and etc.
         double** interactionMatrix = NULL;
         double* interactionEigenEnergies = NULL;
         interactionMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d
                              (interactionMatrixDimension, interactionMatrixDimension);
         interactionEigenEnergies = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(interactionMatrixDimension);

         try{

            // calculate interaction matrix
            this->CalcInteractionMatrix(interactionMatrix, expansionVectors, interactionMatrixDimension);

            // diagonalize interaction matrix
            bool calcEigenVectors = true;
            MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(interactionMatrix,
                                                                    interactionEigenEnergies, 
                                                                    interactionMatrixDimension, 
                                                                    calcEigenVectors);

            // check convergence of all excited states
            notConvergedStates=0;
            allConverged = true;
            for(int i=0; i<numberExcitedStates; i++){
       
               // calculate i-th ritz vector
               this->CalcRitzVector(ritzVector, expansionVectors, interactionMatrix, interactionMatrixDimension, i);

               // calculate i-th residual vector and the norm of the residual vector
               double norm = 0.0;
               this->CalcResidualVectorAndNorm(residualVector, &norm, ritzVector, interactionEigenEnergies, i);

               // output norm of residual vector
               cout << "\t  " << i+1 << this->messageResidualNorm << norm << endl;
               if(i == numberExcitedStates-1){
                  cout << endl;
               }
        
               // check tolerance for the norm of the residual vector.
               if(norm < normTol){
                  convergeExcitedStates[i] = true;
               }
               else{
                  convergeExcitedStates[i] = false;
                  allConverged = false;
                  if(interactionMatrixDimension+notConvergedStates == maxDim && maxDim !=this->matrixCISdimension){
                     reachMaxDim = true;
                     break;
                  }
                  else if(interactionMatrixDimension+notConvergedStates == this->matrixCISdimension){
                     goToDirectCIS = true;
                     break;
                  }
            
                  // update expansion vectors
                  this->UpdateExpansionVectors(expansionVectors, interactionEigenEnergies, residualVector,
                                                interactionMatrixDimension, &notConvergedStates, i);

               }
            } 

            if(allConverged){
               // copy to cis eigen vector and value
               for(int i=0; i<numberExcitedStates; i++){
                  this->excitedEnergies[i] = interactionEigenEnergies[i];
                  this->CalcRitzVector(ritzVector, expansionVectors, interactionMatrix, interactionMatrixDimension, i);
                  for(int j=0; j<this->matrixCISdimension; j++){
                     this->matrixCIS[i][j] = ritzVector[j];
                  }
               }
            }
         }
         catch(MolDSException ex){
            this->FreeDavidsonRoopCISTemporaryMtrices(&interactionMatrix, 
                                                      interactionMatrixDimension, 
                                                      &interactionEigenEnergies);
            throw ex;
         }
         this->FreeDavidsonRoopCISTemporaryMtrices(&interactionMatrix, 
                                                   interactionMatrixDimension, 
                                                   &interactionEigenEnergies);

         // stop the Davidson roop
         if(allConverged){
            cout << this->messageDavidsonConverge;
            break;
         }
         else if(!allConverged && goToDirectCIS){
            cout << this->messageDavidsonReachCISMatrix;
            cout << this->messageDavidsonGoToDirect;
            break;
         }
         else if(!allConverged && reachMaxDim){
            stringstream ss;
            ss << endl;
            ss << this->errorMessageDavidsonNotConverged;
            ss << this->errorMessageDavidsonMaxDim << maxDim << endl;
            throw MolDSException(ss.str());
         }
         else if(!allConverged && k==maxIter-1){
            stringstream ss;
            ss << this->errorMessageDavidsonNotConverged;
            ss << this->errorMessageDavidsonMaxIter << maxIter << endl;
            throw MolDSException(ss.str());
         }

      }// end Davidson roop
   }
   catch(MolDSException ex){
      this->FreeDavidsonCISTemporaryMtrices(&expansionVectors, &residualVector, &ritzVector);
      throw ex;
   }
   this->FreeDavidsonCISTemporaryMtrices(&expansionVectors, &residualVector, &ritzVector);

   cout << this->messageDoneDavidsonCIS;
   // change algorithm from Davidso to direct
   if(goToDirectCIS){
      this->DoesCISDirect();
   }

}

void ZindoS::FreeDavidsonCISTemporaryMtrices(double*** expansionVectors, double** residualVector, double** ritzVector){

   if(*expansionVectors != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(expansionVectors, this->matrixCISdimension);
      //cout << "expansion vectors deleted\n";
   }
   if(*residualVector != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(residualVector);
      //cout << "residual vector deleted\n";
   }
   if(*ritzVector != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(ritzVector);
      //cout << "ritz vector deleted\n";
   }

}

void ZindoS::FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix, 
                                             double interactionMatrixDimension, 
                                             double** interactionEigenEnergies){
   if(*interactionMatrix != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(interactionMatrix, interactionMatrixDimension);
      //cout << "interactionMatrix deleted\n";
   }
   if(*interactionEigenEnergies != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(interactionEigenEnergies);
      //cout << "interactionEigenEnergies deleted\n";
   }

}

void ZindoS::DoesCISDirect(){

   cout << this->messageStartDirectCIS;
   bool calcEigenVectors = true;
   MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(this->matrixCIS,
                                                           this->excitedEnergies, 
                                                           this->matrixCISdimension, 
                                                           calcEigenVectors);
   cout << this->messageDoneDirectCIS;
}

void ZindoS::CalcCISMatrix(double** matrixCIS, int numberOcc, int numberVir){
   cout << this->messageStartCalcCISMatrix;
   double ompStartTime = omp_get_wtime();

   #pragma omp parallel for schedule(auto)
   for(int k=0; k<numberOcc*numberVir; k++){
      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->molecule->GetTotalNumberValenceElectrons()/2 - (k/numberVir) -1;
      int moA = this->molecule->GetTotalNumberValenceElectrons()/2 + (k%numberVir);

      for(int l=k; l<numberOcc*numberVir; l++){
         // single excitation from J-th (occupied)MO to B-th (virtual)MO
         int moJ = this->molecule->GetTotalNumberValenceElectrons()/2 - (l/numberVir) -1;
         int moB = this->molecule->GetTotalNumberValenceElectrons()/2 + (l%numberVir);
         double value=0.0;

         // Fast algorith, but this is not easy to read. Slow algorithm is alos written below.
         Atom* atomA = NULL;
         Atom* atomB = NULL;
         int firstAOIndexA;
         int firstAOIndexB;
         int numberAOsA;
         int numberAOsB;
         double gamma;
         double exchange;
         double coulomb;
         OrbitalType orbitalMu;
         OrbitalType orbitalNu;
         // Off diagonal term (right upper)
         if(k<l){
            for(int A=0; A<molecule->GetAtomVect()->size(); A++){
               atomA = (*molecule->GetAtomVect())[A];
               firstAOIndexA = atomA->GetFirstAOIndex();
               numberAOsA = atomA->GetValence().size();

               for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                  orbitalMu = atomA->GetValence()[mu-firstAOIndexA];

                  // CNDO term
                  for(int B=A; B<molecule->GetAtomVect()->size(); B++){
                     atomB = (*molecule->GetAtomVect())[B];
                     firstAOIndexB = atomB->GetFirstAOIndex();
                     numberAOsB = atomB->GetValence().size();

                     for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                        orbitalNu = atomB->GetValence()[nu-firstAOIndexB];

                        if(A<B){
                           gamma = this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
                           value += 2.0*gamma*fockMatrix[moA][mu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moJ][nu]
                                             *fockMatrix[moB][nu];
                           value -=     gamma*fockMatrix[moA][mu]
                                             *fockMatrix[moB][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][nu];
                           value += 2.0*gamma*fockMatrix[moA][nu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][mu]
                                             *fockMatrix[moB][mu];
                           value -=     gamma*fockMatrix[moA][nu]
                                             *fockMatrix[moB][nu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moJ][mu];
                        }
                        else{
                           gamma = atomA->GetZindoF0ss();
                           value += 2.0*gamma*fockMatrix[moA][mu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moJ][nu]
                                             *fockMatrix[moB][nu];
                           value -=     gamma*fockMatrix[moA][mu]
                                             *fockMatrix[moB][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][nu];
                        }  
                     }
                  }

                  // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
                  for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                     orbitalNu = atomA->GetValence()[nu-firstAOIndexA];

                     if(mu!=nu){
                        exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        value += 2.0*exchange*fockMatrix[moA][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][nu]
                                             *fockMatrix[moB][mu];
                        value += 2.0*exchange*fockMatrix[moA][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][mu]
                                             *fockMatrix[moB][nu];
                        value -=     exchange*fockMatrix[moA][mu]
                                             *fockMatrix[moB][nu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moJ][mu];
                        value -=     exchange*fockMatrix[moA][mu]
                                             *fockMatrix[moB][nu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moJ][nu];
                     }

                     coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

                     if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                         (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                           gamma = atomA->GetZindoF0ss();
                     }
                     else{
                        stringstream ss;
                        ss << this->errorMessageCalcCISMatrix;
                        ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
                        ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
                        ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
                        throw MolDSException(ss.str());
                     }   

                     value += 2.0*(coulomb-gamma)*fockMatrix[moA][mu]
                                                 *fockMatrix[moI][mu]
                                                 *fockMatrix[moJ][nu]
                                                 *fockMatrix[moB][nu];
                     value -=     (coulomb-gamma)*fockMatrix[moA][mu]
                                                 *fockMatrix[moB][mu]
                                                 *fockMatrix[moI][nu]
                                                 *fockMatrix[moJ][nu];
                  }
               }
            }
         }
         // Diagonal term
         else if(k==l){
            value = this->energiesMO[moA] - this->energiesMO[moI];
            for(int A=0; A<molecule->GetAtomVect()->size(); A++){
               atomA = (*molecule->GetAtomVect())[A];
               firstAOIndexA = atomA->GetFirstAOIndex();
               numberAOsA = atomA->GetValence().size();

               for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                  orbitalMu = atomA->GetValence()[mu-firstAOIndexA];

                  // CNDO term
                  for(int B=A; B<molecule->GetAtomVect()->size(); B++){
                     atomB = (*molecule->GetAtomVect())[B];
                     firstAOIndexB = atomB->GetFirstAOIndex();
                     numberAOsB = atomB->GetValence().size();

                     for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                        orbitalNu = atomB->GetValence()[nu-firstAOIndexB];

                        if(A<B){
                           gamma = this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
                           value += 2.0*gamma*fockMatrix[moI][mu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moI][nu];
                           value -=     gamma*fockMatrix[moI][mu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][nu];
                           value += 2.0*gamma*fockMatrix[moI][nu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moI][mu];
                           value -=     gamma*fockMatrix[moI][nu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moA][mu];
                        }
                        else{
                           gamma = atomA->GetZindoF0ss();
                           value += 2.0*gamma*fockMatrix[moI][mu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moI][nu];
                           value -=     gamma*fockMatrix[moI][mu]
                                             *fockMatrix[moI][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][nu];
                        }  
                     }
                  }

                  // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
                  for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                     orbitalNu = atomA->GetValence()[nu-firstAOIndexA];

                     if(mu!=nu){
                        exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        value += 2.0*exchange*fockMatrix[moI][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moI][mu];
                        value += 2.0*exchange*fockMatrix[moI][mu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moI][nu];
                        value -=     exchange*fockMatrix[moI][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moA][nu]
                                             *fockMatrix[moA][mu];
                        value -=     exchange*fockMatrix[moI][mu]
                                             *fockMatrix[moI][nu]
                                             *fockMatrix[moA][mu]
                                             *fockMatrix[moA][nu];
                     }

                     coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

                     if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                         (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                           gamma = atomA->GetZindoF0ss();
                     }
                     else{
                        stringstream ss;
                        ss << this->errorMessageCalcCISMatrix;
                        ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
                        ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
                        ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
                        throw MolDSException(ss.str());
                     }   

                     value += 2.0*(coulomb-gamma)*fockMatrix[moI][mu]
                                                 *fockMatrix[moA][mu]
                                                 *fockMatrix[moA][nu]
                                                 *fockMatrix[moI][nu];
                     value -=     (coulomb-gamma)*fockMatrix[moI][mu]
                                                 *fockMatrix[moI][mu]
                                                 *fockMatrix[moA][nu]
                                                 *fockMatrix[moA][nu];
                  }
               }
            }
         }
         // End of the fast algorith.

         /*
         // Slow algorith, but this is easy to read. Fast altorithm is also written above.
         // diagonal term
         if(k==l){
            value = this->energiesMO[moA] - this->energiesMO[moI] 
                     +2.0*this->GetMolecularIntegralElement(moI, moA, moA, moI, 
                                          this->molecule, this->fockMatrix, NULL)
                     -    this->GetMolecularIntegralElement(moI, moI, moA, moA, 
                                          this->molecule, this->fockMatrix, NULL);

         }
         // Off diagonal term (right upper)
         else if(k<l){
            value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                          this->molecule, this->fockMatrix, NULL)
                     -    this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                          this->molecule, this->fockMatrix, NULL);
         }
         // End of the slow algorith.
         */

         matrixCIS[k][l] = value;
      }
   }
   double ompEndTime = omp_get_wtime();
   cout << this->messageOmpElapsedTimeCalcCISMarix;
   cout << ompEndTime - ompStartTime;
   cout << this->messageUnitSec << endl;
   cout << this->messageDoneCalcCISMatrix;
}

// electronicStateIndex is index of the electroinc eigen state.
// "electronicStateIndex = 0" means electronic ground state. 
void ZindoS::CalcForce(int electronicStateIndex){

   // malloc or initialize Force matrix
   if(this->matrixForce == NULL){
      this->matrixForce = MallocerFreer::GetInstance()->
                          MallocDoubleMatrix2d(this->molecule->GetAtomVect()->size(), 
                                               CartesianType_end);
   }
   else{
      MallocerFreer::GetInstance()->
      InitializeDoubleMatrix2d(this->matrixForce,
                               this->molecule->GetAtomVect()->size(),
                               CartesianType_end);
   }

   #pragma omp parallel for schedule(auto)
   for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
      Atom* atomA = (*molecule->GetAtomVect())[a];
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValence().size();
      double coreRepulsion[CartesianType_end] = {0.0,0.0,0.0};
      double electronicForce1[CartesianType_end] = {0.0,0.0,0.0};
      double electronicForce2[CartesianType_end] = {0.0,0.0,0.0};
      double electronicForce3[CartesianType_end] = {0.0,0.0,0.0};
      double*** overlapDer = MallocerFreer::GetInstance()->MallocDoubleMatrix3d
                                      (OrbitalType_end, OrbitalType_end, CartesianType_end);
      for(int b=0; b<this->molecule->GetAtomVect()->size(); b++){
         if(a != b){
            Atom* atomB = (*molecule->GetAtomVect())[b];
            int firstAOIndexB = atomB->GetFirstAOIndex();
            int numberAOsB = atomB->GetValence().size();

            // calc. first derivative of overlap.
            this->CalcDiatomicOverlapFirstDerivative(overlapDer, atomA, atomB);

            for(int i=0; i<CartesianType_end; i++){
               coreRepulsion[i] += this->GetDiatomCoreRepulsionFirstDerivative
                                                (a, b, (CartesianType)i);
               electronicForce1[i] += ( atomA->GetCoreCharge()*atomicElectronPopulation[b]
                                       +atomB->GetCoreCharge()*atomicElectronPopulation[a])
                                       *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                          (atomA, s, atomB, s, (CartesianType)i);
            }

            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
               for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                  OrbitalType orbitalNu = atomB->GetValence()[nu-firstAOIndexB];
                  double bondParameter = 0.5*(atomA->GetBondingParameter(
                                                      this->theory, orbitalMu) 
                                             +atomB->GetBondingParameter(
                                                      this->theory, orbitalNu)); 
                  for(int i=0; i<CartesianType_end; i++){
                     electronicForce2[i] += 2.0*this->orbitalElectronPopulation[mu][nu]
                                           *bondParameter
                                           *overlapDer[mu-firstAOIndexA][nu-firstAOIndexB][i];
                     electronicForce3[i] += (this->orbitalElectronPopulation[mu][mu]
                                            *this->orbitalElectronPopulation[nu][nu]
                                            -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                            *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                                (atomA, orbitalMu, atomB, orbitalNu,
                                                (CartesianType)i);
                     }
                  }
               }

            for(int i=0; i<CartesianType_end; i++){
               this->matrixForce[a][i] = -1.0*(coreRepulsion[i]
                                          - electronicForce1[i] 
                                          + electronicForce2[i]
                                          + electronicForce3[i]);
            }

         }
      }
      if(overlapDer != NULL){
         MallocerFreer::GetInstance()->FreeDoubleMatrix3d(&overlapDer, 
                                                                OrbitalType_end,
                                                                OrbitalType_end);
      }
   }
   
  
   /*
   // Calculate force. First derivative of overlap integral is
   // calculated with GTO expansion technique.
   //#pragma omp parallel for schedule(auto)
   for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
      Atom* atomA = (*molecule->GetAtomVect())[a];
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValence().size();
      for(int i=0; i<CartesianType_end; i++){

         double coreRepulsion = 0.0;
         double electronicForce1 = 0.0;
         double electronicForce2 = 0.0;
         double electronicForce3 = 0.0;
         for(int b=0; b<this->molecule->GetAtomVect()->size(); b++){
            if(a != b){
               Atom* atomB = (*molecule->GetAtomVect())[b];
               int firstAOIndexB = atomB->GetFirstAOIndex();
               int numberAOsB = atomB->GetValence().size();

               // Calculation of core repusion force
               coreRepulsion += this->GetDiatomCoreRepulsionFirstDerivative
                                                (a, b, (CartesianType)i);

               // Calculate force arise from electronic part.
               electronicForce1 += ( atomA->GetCoreCharge()*atomicElectronPopulation[b]
                                    +atomB->GetCoreCharge()*atomicElectronPopulation[a])
                                    *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                          (atomA, s, atomB, s, (CartesianType)i);

               for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                  OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
                  for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                     OrbitalType orbitalNu = atomB->GetValence()[nu-firstAOIndexB];

                     double bondParameter = 0.5*(atomA->GetBondingParameter(this->theory, orbitalMu) 
                                                +atomB->GetBondingParameter(this->theory, orbitalNu)); 

                     electronicForce2 += 2.0*this->orbitalElectronPopulation[mu][nu]
                                         *bondParameter
                                         *this->GetOverlapElementFirstDerivativeByGTOExpansion
                                                (atomA, mu-firstAOIndexA, 
                                                 atomB, nu-firstAOIndexB,
                                                 STO6G, (CartesianType)i);

                     electronicForce3 += (this->orbitalElectronPopulation[mu][mu]
                                         *this->orbitalElectronPopulation[nu][nu]
                                         -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                         *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                                (atomA, orbitalMu, atomB, orbitalNu,
                                                (CartesianType)i);
                  }
               }
            }
         }

         this->matrixForce[a][i] = -1.0*(coreRepulsion 
                                          - electronicForce1 
                                          + electronicForce2
                                          + electronicForce3);
      }
   }
   */

   /*  
   // checking of calculated force
   cout << "chek the force\n";
   double checkSumForce[3] = {0.0, 0.0, 0.0};
   for(int a=0; a<this->molecule->GetAtomVect()->size(); a++){
      for(int i=0; i<CartesianType_end; i++){
         cout << this->matrixForce[a][i] << " ";
         checkSumForce[i] += this->matrixForce[a][i];
      }
      cout << endl;
   }
   cout << endl << endl;
   for(int i=0; i<CartesianType_end; i++){
      cout << "force: " << i << " "  << checkSumForce[i] << endl;
   }
   */
}

}
#endif



