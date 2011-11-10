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
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double CalcDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB);
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA);
   virtual void CalcHFProperties();
   virtual void OutputHFResults(double** fockMatrix, 
                                double* energiesMO, 
                                double* atomicElectronPopulation, 
                                MolDS_base::Molecule* molecule);
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
                                        int mu, int nu, 
                                        MolDS_base::Molecule* molecule, 
                                        double** gammaAB, 
                                        double** overelap,
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
   virtual double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                                MolDS_base::OrbitalType orbital2, 
                                MolDS_base_atoms::Atom* atom); 
   virtual double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                                 MolDS_base::OrbitalType orbital2, 
                                 MolDS_base_atoms::Atom* atom); 
   virtual void CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                                   MolDS_base::Molecule* molecule);
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              MolDS_base::Molecule* molecule, 
                                              double** fockMatrix, 
                                              double** gammaAB);
   virtual void CalcCISMatrix(double** matrixCIS, int numberOcc, int numberVir);
   virtual void CalcForce(std::vector<int> elecStates);
   double GetNddoRepulsionIntegral(MolDS_base_atoms::Atom* atomA, 
                                   MolDS_base::OrbitalType mu, 
                                   MolDS_base::OrbitalType nu,
                                   MolDS_base_atoms::Atom* atomB, 
                                   MolDS_base::OrbitalType lambda, 
                                   MolDS_base::OrbitalType sigma);
   double GetNddoRepulsionIntegralFirstDerivative(MolDS_base_atoms::Atom* atomA, 
                                                  MolDS_base::OrbitalType mu, 
                                                  MolDS_base::OrbitalType nu,
                                                  MolDS_base_atoms::Atom* atomB, 
                                                  MolDS_base::OrbitalType lambda, 
                                                  MolDS_base::OrbitalType sigma,
                                                  MolDS_base::CartesianType axisA);
private:
   std::string errorMessageMultipoleA;
   std::string errorMessageMultipoleB;
   std::string messageHeatsFormation;
   std::string messageHeatsFormationTitle;
   double heatsFormation;
   void CalcHeatsFormation(double* heatsFormation, MolDS_base::Molecule* molecule);
   double GetElectronCoreAttraction(int atomAIndex, 
                                    int atomBIndex, 
                                    int mu, 
                                    int nu, 
                                    double****** twoElecTwoCore);
   double GetElectronCoreAttractionFirstDerivative(
                                    int atomAIndex, 
                                    int atomBIndex, 
                                    int mu, 
                                    int nu, 
                                    double***** twoElecTwoCoreFirstDerivative,
                                    MolDS_base::CartesianType axisA);
   void CalcTwoElecTwoCoreDiatomic(double**** matrix, int atomAIndex, int atomBIndex);
   void CalcTwoElecTwoCoreDiatomicFirstDerivatives(double***** matrix, 
                                                   int atomAIndex, 
                                                   int atomBIndex);
   void RotateTwoElecTwoCoreDiatomicToSpaceFramegc(double**** matrix, 
                                                   double** rotatingMatrix);
   void RotateTwoElecTwoCoreDiatomicFirstDerivativesToSpaceFramegc(
                                             double***** matrix, 
                                             double**** twoElecTwoCoreDiatomic,
                                             double** rotatingMatrix,
                                             double*** rMatDeri);
   double GetSemiEmpiricalMultipoleInteraction(MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab);
   double GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                                               MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab);
   void FreeCalcForceTempMatrices(double**** overlapDer, 
                                  double****** twoElecTwoCoreFirstDeriv);
};

}
#endif



