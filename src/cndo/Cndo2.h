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
   MolDS_base::TheoryType GetTheoryType() const;
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   MolDS_base::Molecule* GetMolecule();
   void DoesSCF();
   void DoesSCF(bool requiresGuess);
   virtual void DoesCIS();
   double** GetForce(int elecState);
   double*** GetForce(std::vector<int> elecStates);
   double GetElectronicEnergy(int elecState) const;
   double GetCoreRepulsionEnergy() const;
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
   MolDS_base::TheoryType theory;
   MolDS_base::Molecule* molecule;
   double coreRepulsionEnergy;
   double** fockMatrix;
   double* energiesMO;
   double*** matrixForce;
   double****** twoElecTwoCore;
   double** orbitalElectronPopulation; //P_{\mu\nu} of (2.50) in J. A. Pople book.
   double*   atomicElectronPopulation; //P_{AB} of (3.21) in J. A. Pople book.
   double** matrixCIS;
   double* excitedEnergies;
   int matrixCISdimension;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcHFProperties();
   double GetBondingAdjustParameterK(MolDS_base::ShellType shellA, 
                                     MolDS_base::ShellType shellB) const;
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int indexAtomA, 
                                                        int indexAtomB, 
                                                        MolDS_base::CartesianType axisA) const;
   double GetReducedOverlap(int na, int la, int m, 
                            int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlap(int na, int nb, double alpha, double beta) const;
   double GetReducedOverlapFirstDerivativeAlpha(int na, int la, int m, 
                                                int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapFirstDerivativeBeta(int na, int la, int m, 
                                               int nb, int lb, double alpha, double beta) const;
   double GetOverlapElementFirstDerivativeByGTOExpansion(const MolDS_base_atoms::Atom& atomA, 
                                                         int valenceIndexA, 
                                                         const MolDS_base_atoms::Atom& atomB, 
                                                         int valenceIndexB,
                                                         MolDS_base::STOnGType stonG, 
                                                         MolDS_base::CartesianType axisA) const; // See [DY_1977].
   void CalcRotatingMatrix(double** rotatingMatrix, 
                           const MolDS_base_atoms::Atom& atomA, 
                           const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcGammaAB(double** gammaAB, const MolDS_base::Molecule& molecule) const;
   virtual double GetFockDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     const MolDS_base::Molecule& molecule, 
                                     double const* const* gammaAB,
                                     double const* const* orbitalElectronPopulation, 
                                     double const* atomicElectronPopulation,
                                     double const* const* const* const* const* const* twoElecTwoCore, 
                                     bool isGuess) const;
   virtual double GetFockOffDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                        const MolDS_base_atoms::Atom& atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
                                        int mu, 
                                        int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overlap,
                                        double const* const* orbitalElectronPopulation, 
                                        double const* const* const* const* const* const* twoElecTwoCore, 
                                        bool isGuess) const;
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                   const MolDS_base_atoms::Atom& atomA, 
                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                                  const MolDS_base_atoms::Atom& atomA, 
                                                                  const MolDS_base_atoms::Atom& atomB) const;
   void CalcDiatomicOverlapFirstDerivative(double*** overlapFirstDeri, 
                                           const MolDS_base_atoms::Atom& atomA, 
                                           const MolDS_base_atoms::Atom& atomB) const;
   void FreeDiatomicOverlapDeriTemps(double*** diatomicOverlap, 
                                     double*** rotatingMatrix,
                                     double*** diaOverlapDeriR,
                                     double**** rMatDeri) const;
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                                   const MolDS_base::Molecule& molecule) const;
   virtual void CalcForce(std::vector<int> elecStates);
   virtual void OutputHFResults(double const* const* fockMatrix, 
                                double const* energiesMO, 
                                double const* atomicElectronPopulation, 
                                const MolDS_base::Molecule& molecule) const;
   void CalcRotatingMatrixFirstDerivatives(double*** rMatFirstDeri, 
                                           const MolDS_base_atoms::Atom& atomA,
                                           const MolDS_base_atoms::Atom& atomB) const;
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
   bool SatisfyConvergenceCriterion(double const* const* oldOrbitalElectronPopulation, 
                                    double const* const* orbitalElectronPopulation,
                                    int numberAOs, 
                                    double* rmsDensity, 
                                    int times) const;
   void UpdateOldOrbitalElectronPopulation(double** oldOrbitalElectronPopulation, 
                                           double const* const* orbitalElectronPopulation,
                                           int numberAOs) const;
   void CalcOrbitalElectronPopulation(double** orbitalElectronPopulation, 
                                      const MolDS_base::Molecule& molecule, 
                                      double const* const* fockMatrix) const;
   void CalcAtomicElectronPopulation(double* atomicElectronPopulation,
                                     double const* const* orbitalElectronPopulation, 
                                     const MolDS_base::Molecule& molecule) const;
   void CalcOverlap(double** overlap, const MolDS_base::Molecule& molecule) const;
   void CalcOverlapByGTOExpansion(double** overlap, 
                                  const MolDS_base::Molecule& molecule, 
                                  MolDS_base::STOnGType stonG) const; //See [DY_1977]
   double GetOverlapElementByGTOExpansion(const MolDS_base_atoms::Atom& atomA, 
                                          int valenceIndexA, 
                                          const MolDS_base_atoms::Atom& atomB, 
                                          int valenceIndexB,
                                          MolDS_base::STOnGType stonG) const; // see [DY_1977]
   double GetGaussianOverlap(MolDS_base::AtomType atomTypeA, 
                             MolDS_base::OrbitalType valenceOrbitalA, 
                             double gaussianExponentA, 
                             MolDS_base::AtomType atomTypeB, 
                             MolDS_base::OrbitalType valenceOrbitalB, 
                             double gaussianExponentB, 
                             double dx, 
                             double dy, 
                             double dz, 
                             double Rab) const; // see [DY_1977]
   double GetGaussianOverlapSaSb(double gaussianExponentA, 
                                 double gaussianExponentB, 
                                 double Rab) const; // see [DY_1977]
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
                                            MolDS_base::CartesianType axisA) const;// see [DY_1977]
   void CalcFockMatrix(double** fockMatrix, 
                       const MolDS_base::Molecule& molecule, 
                       double const* const* overlap, 
                       double const* const* gammaAB,
                       double const* const* orbitalElectronPopulation, 
                       double const* atomicElectronPopulation,
                       double const* const* const* const* const* const* twoElecTwoCore,
                       bool isGuess) const;
   void RotateDiatmicOverlapToSpaceFrame(double** diatomicOverlap, 
                                         double const* const* rotatingMatrix) const;
   void SetOverlapElement(double** overlap, 
                          double const* const* diatomicOverlap, 
                          const MolDS_base_atoms::Atom& atomA, 
                          const MolDS_base_atoms::Atom& atomB) const;
   double GetAuxiliaryA(int k, double rho) const;
   double GetAuxiliaryB(int k, double rho) const;
   double GetAuxiliaryD(int la, int lb, int m) const;
   double GetAuxiliaryAFirstDerivative(int k, double rho) const;
   double GetAuxiliaryBFirstDerivative(int k, double rho) const;
   void DoesDamp(double rmsDensity, 
                 double** orbitalElectronPopulation, 
                 double const* const* oldOrbitalElectronPopulation, 
                 const MolDS_base::Molecule& molecule) const;
   void DoesDIIS(double** orbitalElectronPopulation,
                 double const* const* oldOrbitalElectronPopulation,
                 double*** diisStoredDensityMatrix,
                 double*** diisStoredErrorVect,
                 double** diisErrorProducts,
                 double* diisErrorCoefficients,
                 int diisNumErrorVect,
                 const MolDS_base::Molecule& molecule, 
                 int step) const;
   void CheckEnableAtomType(const MolDS_base::Molecule& molecule) const;
   void CheckNumberValenceElectrons(const MolDS_base::Molecule& molecule) const;
   void FreeDiatomicOverlapAndRotatingMatrix(double*** diatomicOverlap, 
                                             double*** rotatingMatrix) const;
   void CalcElecHFEnergy(double* elecHFEnergy, 
                         const MolDS_base::Molecule& molecule, 
                         double const* energiesMO, 
                         double const* const* fockMatrix, 
                         double const* const* gammaAB, 
                         double coreRepulsionEnergy) const;
   void FreeElecEnergyMatrices(double*** fMatrix, 
                               double*** hMatrix, 
                               double*** dammyOrbitalElectronPopulation, 
                               double**  dammyAtomicElectronPopulation ) const;
   void FreeSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                 double**** diisStoredDensityMatrix,
                                 double**** diisStoredErrorVect,
                                 double*** diisErrorProducts,
                                 double** diisErrorCoefficients) const;
};


}
#endif



