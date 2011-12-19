#ifndef INCLUDED_MNDO
#define INCLUDED_MNDO
namespace MolDS_mndo{

/***
 *  Main Refferences for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
class Mndo : public MolDS_zindo::ZindoS{
public:
   Mndo();
   virtual ~Mndo();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
protected:
   std::string errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
   std::string errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles;
   std::string errorMessageGetNddoRepulsionIntegral;
   std::string errorMessageGetNddoRepulsionIntegralFirstDerivative;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms;
   std::string errorMessageCalcZMatrixForceEtaNull;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA) const;
   virtual void CalcHFProperties();
   virtual void OutputHFResults(double const* const* fockMatrix, 
                                double const* energiesMO, 
                                double const* atomicElectronPopulation, 
                                const MolDS_base::Molecule& molecule) const;
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
                                        int mu, int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overelap,
                                        double const* const* orbitalElectronPopulation, 
                                        double const* const* const* const* const* const* twoElecTwoCore,
                                        bool isGuess) const;
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                   const MolDS_base_atoms::Atom& atomA, 
                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                                  const MolDS_base_atoms::Atom& atomA, 
                                                                  const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                                MolDS_base::OrbitalType orbital2, 
                                const MolDS_base_atoms::Atom& atom) const; 
   virtual double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                                 MolDS_base::OrbitalType orbital2, 
                                 const MolDS_base_atoms::Atom& atom) const; 
   virtual void CalcTwoElecTwoCore(double****** twoElecTwoCore, MolDS_base::Molecule* molecule);
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcCISMatrix(double** matrixCIS, int numberActiveOcc, int numberActiveVir) const;
   virtual void CalcForce(std::vector<int> elecStates);
   double GetNddoRepulsionIntegral(const MolDS_base_atoms::Atom& atomA, 
                                   MolDS_base::OrbitalType mu, 
                                   MolDS_base::OrbitalType nu,
                                   const MolDS_base_atoms::Atom& atomB, 
                                   MolDS_base::OrbitalType lambda, 
                                   MolDS_base::OrbitalType sigma) const;
   double GetNddoRepulsionIntegralFirstDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                  MolDS_base::OrbitalType mu, 
                                                  MolDS_base::OrbitalType nu,
                                                  const MolDS_base_atoms::Atom& atomB, 
                                                  MolDS_base::OrbitalType lambda, 
                                                  MolDS_base::OrbitalType sigma,
                                                  MolDS_base::CartesianType axisA) const;
private:
   std::string errorMessageMultipoleA;
   std::string errorMessageMultipoleB;
   std::string messageHeatsFormation;
   std::string messageHeatsFormationTitle;
   struct MoIndexPair{int moI; int moJ; bool isMoICIMO; bool isMoJCIMO;};
   double*** zMatrixForce;
   double*** etaMatrixForce;
   int zMatrixForceElecStatesNum;
   int etaMatrixForceElecStatesNum;
   double heatsFormation;
   double GetGammaNRElement(int moI, int moJ, int moK, int moL) const;
   double GetGammaRElement(int moI, int moJ, int moK, int moL) const;
   double GetNNRElement(int moI, int moJ, int moK, int moL) const;
   double GetNRElement(int moI, int moJ, int moK, int moL) const;
   double GetKNRElement(int moI, int moJ, int moK, int moL) const;
   double GetKRElement(int moI, int moJ, int moK, int moL) const;
   double GetKRDagerElement(int moI, int moJ, int moK, int moL) const;
   void MallocTempMatrixForZMatrix(double** delta,
                                   double** q,
                                   double*** kNR, 
                                   double*** kRDag,
                                   double** y,
                                   double*** transposedFockMatrix,
                                   double*** xiOcc,
                                   double*** xiVir,
                                   int sizeQNR,
                                   int sizeQR) const;
   void FreeTempMatrixForZMatrix(double** delta,
                                 double** q,
                                 double*** kNR, 
                                 double*** kRDag,
                                 double** y,
                                 double*** transposedFockMatrix,
                                 double*** xiOcc,
                                 double*** xiVir,
                                 int sizeQNR,
                                 int sizeQR) const;
   void CalcDeltaVector(double* delta, int exciteState) const;
   double GetSmallQElement(int moI, 
                           int moP, 
                           double const* const* xiOcc, 
                           double const* const* xiVir,
                           double const* const* eta) const;
   void CalcQVector(double* q, 
                    double const* delta, 
                    double const* const* xiOcc,
                    double const* const* xiVir,
                    double const* const* eta,
                    std::vector<MoIndexPair> nonRedundantQIndeces,
                    std::vector<MoIndexPair> redundantQIndeces) const;
   void TransposeFockMatrixMatrix(double** transposedFockMatrix) const;
   void CalcKNRMatrix(double** kNR, 
                      std::vector<MoIndexPair> nonRedundantQIndeces) const;
   void CalcKRDagerMatrix(double** kRDager, 
                          std::vector<MoIndexPair> nonRedundantQIndeces,
                          std::vector<MoIndexPair> redundantQIndeces) const;
   void CalcAuxiliaryVector(double* y,
                            double const* q,
                            double const* const* kRDager,
                            std::vector<MoIndexPair> nonRedundantQIndeces,
                            std::vector<MoIndexPair> redundantQIndeces) const;
   void CalcXiMatrices(double** xiOcc, 
                       double** xiVir, 
                       int exciteState,
                       double const* const* transposedFockMatrix) const;
   double GetZMatrixForceElement(double const* y,
                                 double const* q,
                                 double const* const* transposedFockMatrix,
                                 std::vector<MoIndexPair> nonRedundantQIndeces,
                                 std::vector<MoIndexPair> redundantQIndeces,
                                 int mu, 
                                 int nu) const;
   void CheckZMatrixForce(std::vector<int> elecStates);
   void CheckEtaMatrixForce(std::vector<int> elecStates);
   void CalcZMatrixForce(std::vector<int> elecStates);
   void CalcEtaMatrixForce(std::vector<int> elecStates);
   bool RequiresExcitedStatesForce(std::vector<int> elecStates) const;
   double GetCISCoefficientMOEnergy(int k, 
                                    int l, 
                                    int r, 
                                    int numberActiveVir) const;
   double GetCISCoefficientTwoElecIntegral(int k, 
                                           int l, 
                                           int p, 
                                           int q, 
                                           int r, 
                                           int s, 
                                           int numberActiveVir) const;
   void CalcActiveSetVariablesQ(std::vector<MoIndexPair>* nonRedundantQIndeces, 
                                std::vector<MoIndexPair>* redundantQIndeces) const;
   void CalcHeatsFormation(double* heatsFormation, 
                           const MolDS_base::Molecule& molecule) const;
   double GetElectronCoreAttraction(int atomAIndex, 
                                    int atomBIndex, 
                                    int mu, 
                                    int nu, 
                                    double const* const* const* const* const* const* twoElecTwoCore) const;
   double GetElectronCoreAttractionFirstDerivative(int atomAIndex, 
                                                   int atomBIndex, 
                                                   int mu, 
                                                   int nu, 
                                                   double const* const* const* const* const* twoElecTwoCoreFirstDerivative,
                                                   MolDS_base::CartesianType axisA) const;
   void CalcTwoElecTwoCoreDiatomic(double**** matrix, int atomAIndex, int atomBIndex) const;
   void CalcTwoElecTwoCoreDiatomicFirstDerivatives(double***** matrix, 
                                                   int atomAIndex, 
                                                   int atomBIndex) const;
   void RotateTwoElecTwoCoreDiatomicToSpaceFramegc(double**** matrix, 
                                                   double const* const* rotatingMatrix) const;
   void RotateTwoElecTwoCoreDiatomicFirstDerivativesToSpaceFramegc(
        double***** matrix, 
        double const* const* const* const* twoElecTwoCoreDiatomic,
        double const* const* rotatingMatrix,
        double const* const* const* rMatDeri) const;
   double GetSemiEmpiricalMultipoleInteraction(MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab) const;
   double GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                                               MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab) const;
   void FreeCalcForceTempMatrices(double**** overlapDer, 
                                  double****** twoElecTwoCoreFirstDeriv) const;
   void CalcForceHFElecCoreAttractionPart(double* force, 
                                          int atomAIndex,
                                          int atomBIndex,
                                          double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;
   void CalcForceHFOverlapPart(double* force, 
                               int atomAIndex,
                               int atomBIndex,
                               double*** overlapDer);
   void CalcForceHFTwoElecPart(double* force, 
                               int atomAIndex,
                               int atomBIndex,
                               double***** twoElecTwoCoreFirstDeriv);
   void CalcForceExcitedStaticPart(double* force, 
                                   int elecStateIndex,
                                   int atomAIndex,
                                   int atomBIndex,
                                   double***** twoElecTwoCoreFirstDeriv);
   void CalcForceExcitedElecCoreAttractionPart(double* force, 
                                               int elecStateIndex,
                                               int atomAIndex,
                                               int atomBIndex,
                                               double***** twoElecTwoCoreFirstDeriv);
   void CalcForceExcitedOverlapPart(double* force, 
                                    int elecStateIndex,
                                    int atomAIndex,
                                    int atomBIndex,
                                    double*** overlapDer);
   void CalcForceExcitedTwoElecPart(double* force, 
                                    int elecStateIndex,
                                    int atomAIndex,
                                    int atomBIndex,
                                    double***** twoElecTwoCoreFirstDeriv);

};

}
#endif



