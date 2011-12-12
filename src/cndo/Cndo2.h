#ifndef INCLUDED_CNDO
#define INCLUDED_CNDO
namespace MolDS_cndo{

/***
 *  Refferences for Cndo2 are [PB_1970], [PSS_1965], and [PS_1965].
 */
class Cndo2{
public:
   Cndo2();
   virtual ~Cndo2();
   MolDS_base::TheoryType GetTheoryType();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   MolDS_base::Molecule* GetMolecule();
   void DoesSCF();
   void DoesSCF(bool requiresGuess);
   virtual void DoesCIS();
   double** GetForce(int elecState);
   double*** GetForce(std::vector<int> elecStates);
   double GetElectronicEnergy(int elecState);
   double GetCoreRepulsionEnergy();
protected:
   std::string errorMessageAtomA;
   std::string errorMessageAtomB;
   std::string errorMessageAtomType;
   std::string errorMessageOrbitalType;
   std::string errorMessageSCFNotConverged;
   std::string errorMessageMoleculeNotSet;
   std::string errorMessageOddTotalValenceElectrions;
   std::string errorMessageNotEnebleAtomType;
   std::string errorMessageCoulombInt;
   std::string errorMessageExchangeInt;
   std::string errorMessageMolecularIntegralElement;
   std::string errorMessageGetGaussianOverlapOrbitalD;
   std::string errorMessageGetGaussianOverlapFirstDerivativeOrbitalD;
   std::string errorMessageCISNotImplemented;
   std::string errorMessageCalcForceNotImplemented;
   std::string errorMessageGetElectronicEnergyNULLCISEnergy;
   std::string errorMessageGetElectronicEnergyEnergyNotCalculated;
   std::string errorMessageGetElectronicEnergyNumberCISStates;
   std::string errorMessageGetElectronicEnergySetElecState;
   std::string messageSCFMetConvergence;
   std::string messageStartSCF;
   std::string messageDoneSCF;
   std::string messageOmpElapsedTimeSCF;
   std::string messageUnitSec; 
   std::vector<MolDS_base::AtomType> enableAtomTypes;
   double coreRepulsionEnergy;
   virtual void CalcHFProperties();
   double GetBondingAdjustParameterK(MolDS_base::ShellType shellA, 
                                     MolDS_base::ShellType shellB);
   virtual double CalcDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB);
   virtual double GetDiatomCoreRepulsionFirstDerivative(int indexAtomA, 
                                                        int indexAtomB, 
                                                        MolDS_base::CartesianType axisA);
   double** orbitalElectronPopulation; //P_{\mu\nu} of (2.50) in J. A. Pople book.
   double*   atomicElectronPopulation; //P_{AB} of (3.21) in J. A. Pople book.
   double GetReducedOverlap(int na, int la, int m, 
                            int nb, int lb, double alpha, double beta);
   double GetReducedOverlap(int na, int nb, double alpha, double beta);
   double GetReducedOverlapFirstDerivativeAlpha(int na, int la, int m, 
                                                int nb, int lb, double alpha, double beta);
   double GetReducedOverlapFirstDerivativeBeta(int na, int la, int m, 
                                               int nb, int lb, double alpha, double beta);
   double GetOverlapElementFirstDerivativeByGTOExpansion(MolDS_base_atoms::Atom* atomA, 
                                                         int valenceIndexA, 
                                                         MolDS_base_atoms::Atom* atomB, 
                                                         int valenceIndexB,
                                                         MolDS_base::STOnGType stonG, 
                                                         MolDS_base::CartesianType axisA); // See [DY_1977].
   void CalcRotatingMatrix(double** rotatingMatrix, 
                           MolDS_base_atoms::Atom* atomA, 
                           MolDS_base_atoms::Atom* atomB);
   virtual void CalcGammaAB(double** gammaAB, MolDS_base::Molecule* molecule);
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetFockDiagElement(MolDS_base_atoms::Atom* atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     MolDS_base::Molecule* molecule, 
                                     double** gammaAB,
                                     double** orbitalElectronPopulation, 
                                     double* atomicElectronPopulation,
                                     double****** twoElecTwoCore, 
                                     bool isGuess);
   virtual double GetFockOffDiagElement(MolDS_base_atoms::Atom* atomA, 
                                        MolDS_base_atoms::Atom* atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
                                        int mu, 
                                        int nu, 
                                        MolDS_base::Molecule* molecule, 
                                        double** gammaAB, 
                                        double** overlap,
                                        double** orbitalElectronPopulation, 
                                        double****** twoElecTwoCore, 
                                        bool isGuess);
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                   MolDS_base_atoms::Atom* atomA, 
                                                   MolDS_base_atoms::Atom* atomB);
   virtual void CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                                                double** diatomicOverlapDeri, 
                                                MolDS_base_atoms::Atom* atomA, 
                                                MolDS_base_atoms::Atom* atomB);
   void CalcDiatomicOverlapFirstDerivative(double*** overlapFirstDeri, 
                                           MolDS_base_atoms::Atom* atomA, 
                                           MolDS_base_atoms::Atom* atomB);
   void FreeDiatomicOverlapDeriTemps(double*** diatomicOverlap, 
                                     double*** rotatingMatrix,
                                     double*** diaOverlapDeriR,
                                     double**** rMatDeri);
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              MolDS_base::Molecule* molecule, 
                                              double** fockMatrix, 
                                              double** gammaAB);
   virtual void CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                                   MolDS_base::Molecule* molecule);
   virtual void CalcForce(std::vector<int> elecStates);
   virtual void OutputHFResults(double** fockMatrix, 
                                double* energiesMO, 
                                double* atomicElectronPopulation, 
                                MolDS_base::Molecule* molecule);
   MolDS_base::TheoryType theory;
   MolDS_base::Molecule* molecule;
   double** fockMatrix;
   double* energiesMO;
   double*** matrixForce;
   double****** twoElecTwoCore;
   double** matrixCIS;
   double* excitedEnergies;
   int matrixCISdimension;
   void CalcRotatingMatrixFirstDerivatives(double*** rMatFirstDeri, 
                                           MolDS_base_atoms::Atom* atomA,
                                           MolDS_base_atoms::Atom* atomB);
   struct MoEnergyGap{
      double energyGap;
      int occIndex;
      int virIndex;
      int slaterIndex;
   };
   struct LessMoEnergyGap { 
      bool operator()(const MoEnergyGap& rLeft, const MoEnergyGap& rRight) 
      const { return rLeft.energyGap < rRight.energyGap; } 
   };
private:
   std::string messageEnergiesMOs;
   std::string messageEnergiesMOsTitle;
   std::string messageMullikenAtoms;
   std::string messageMullikenAtomsTitle;
   std::string messageElecEnergy;
   std::string messageElecEnergyTitle;
   std::string messageOcc;
   std::string messageUnOcc;
   std::string messageCoreRepulsionTitle;
   std::string messageCoreRepulsion;
   double elecHFEnergy;
   double** gammaAB;
   double** overlap;
   double bondingAdjustParameterK[2]; //see (3.79) in J. A. Pople book

   // use Y[na][nb][la][lb][m][i][j] 
   // as Y_{ij\lammda} in (B.20) in Pople book for give na, nb, la, lb, m, i, and j.
   static const double Y[MolDS_base::ShellType_end+1]
                        [MolDS_base::ShellType_end+1]
                        [MolDS_base::ShellType_end]
                        [MolDS_base::ShellType_end]
                        [MolDS_base::ShellType_end]
                        [2*MolDS_base::ShellType_end+1]
                        [2*MolDS_base::ShellType_end+1];
   // use Z[na][nb][k] as Z_{k} in (B.30) in Pople book for give na, nb, and k. 
   static const double Z[2*MolDS_base::ShellType_end]
                        [2*MolDS_base::ShellType_end]
                        [4*MolDS_base::ShellType_end-1];
   void CalcCoreRepulsionEnergy();
   bool SatisfyConvergenceCriterion(double** oldOrbitalElectronPopulation, 
                                    double** orbitalElectronPopulation,
                                    int numberAOs, 
                                    double* rmsDensity, 
                                    int times);
   void UpdateOldOrbitalElectronPopulation(double** oldOrbitalElectronPopulation, 
                                           double** orbitalElectronPopulation,
                                           int numberAOs);
   void CalcOrbitalElectronPopulation(double** orbitalElectronPopulation, 
                                      MolDS_base::Molecule* molecule, 
                                      double** fockMatrix);
   void CalcAtomicElectronPopulation(double* atomicElectronPopulation,
                                     double** orbitalElectronPopulation, 
                                     MolDS_base::Molecule* molecule);
   void CalcOverlap(double** overlap, MolDS_base::Molecule* molecule);
   void CalcOverlapByGTOExpansion(double** overlap, 
                                  MolDS_base::Molecule* molecule, 
                                  MolDS_base::STOnGType stonG); //See [DY_1977]
   double GetOverlapElementByGTOExpansion(MolDS_base_atoms::Atom* atomA, 
                                          int valenceIndexA, 
                                          MolDS_base_atoms::Atom* atomB, 
                                          int valenceIndexB,
                                          MolDS_base::STOnGType stonG); // see [DY_1977]
   double GetGaussianOverlap(MolDS_base::AtomType atomTypeA, 
                             MolDS_base::OrbitalType valenceOrbitalA, 
                             double gaussianExponentA, 
                             MolDS_base::AtomType atomTypeB, 
                             MolDS_base::OrbitalType valenceOrbitalB, 
                             double gaussianExponentB, 
                             double dx, 
                             double dy, 
                             double dz, 
                             double Rab); // see [DY_1977]
   double GetGaussianOverlapSaSb(double gaussianExponentA, 
                                 double gaussianExponentB, 
                                 double Rab); // see [DY_1977]
   double GetGaussianOverlapFirstDerivative(MolDS_base::AtomType atomTypeA, 
                                            MolDS_base::OrbitalType valenceOrbitalA, 
                                            double gaussianExponentA, 
                                            MolDS_base::AtomType atomTypeB, 
                                            MolDS_base::OrbitalType valenceOrbitalB, 
                                            double gaussianExponentB, 
                                            double dx, 
                                            double dy, 
                                            double dz, 
                                            double Rab, 
                                            MolDS_base::CartesianType axisA);// see [DY_1977]
   void CalcFockMatrix(double** fockMatrix, 
                       MolDS_base::Molecule* molecule, 
                       double** overlap, 
                       double** gammaAB,
                       double** orbitalElectronPopulation, 
                       double* atomicElectronPopulation,
                       double****** twoElecTwoCore,
                       bool isGuess);
   void RotateDiatmicOverlapToSpaceFrame(double** diatomicOverlap, 
                                         double** rotatingMatrix);
   void SetOverlapElement(double** overlap, 
                          double** diatomicOverlap, 
                          MolDS_base_atoms::Atom* atomA, 
                          MolDS_base_atoms::Atom* atomB);
   double GetAuxiliaryA(int k, double rho);
   double GetAuxiliaryB(int k, double rho);
   double GetAuxiliaryD(int la, int lb, int m);
   double GetAuxiliaryAFirstDerivative(int k, double rho);
   double GetAuxiliaryBFirstDerivative(int k, double rho);
   void DoesDamp(double rmsDensity, 
                 double** orbitalElectronPopulation, 
                 double** oldOrbitalElectronPopulation, 
                 MolDS_base::Molecule* molecule);
   void DoesDIIS(double** orbitalElectronPopulation,
                 double** oldOrbitalElectronPopulation,
                 double*** diisStoredDensityMatrix,
                 double*** diisStoredErrorVect,
                 double** diisErrorProducts,
                 double* diisErrorCoefficients,
                 int diisNumErrorVect,
                 MolDS_base::Molecule* molecule, 
                 int step);
   void CheckEnableAtomType(MolDS_base::Molecule* molecule);
   void CheckNumberValenceElectrons(MolDS_base::Molecule* molecule);
   void FreeDiatomicOverlapAndRotatingMatrix(double*** diatomicOverlap, 
                                             double*** rotatingMatrix);
   void CalcElecEnergy(double* elecHFEnergy, 
                       MolDS_base::Molecule* molecule, 
                       double* energiesMO, 
                       double** fockMatrix, 
                       double** gammaAB, 
                       double coreRepulsionEnergy);
   void FreeElecEnergyMatrices(double*** fMatrix, 
                               double*** hMatrix, 
                               double*** dammyOrbitalElectronPopulation, 
                               double**  dammyAtomicElectronPopulation );
   void FreeSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                 double**** diisStoredDensityMatrix,
                                 double**** diisStoredErrorVect,
                                 double*** diisErrorProducts,
                                 double** diisErrorCoefficients);
};


}
#endif



