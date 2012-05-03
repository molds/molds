//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#ifndef INCLUDED_CNDO
#define INCLUDED_CNDO
namespace MolDS_cndo{

/***
 *  Refferences for Cndo2 are [PB_1970], [PSS_1965], and [PS_1965].
 */
class Cndo2 : public MolDS_base::ElectronicStructure{
public:
   Cndo2();
   virtual ~Cndo2();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   void DoSCF();
   void DoSCF(bool requiresGuess);
   virtual void OutputSCFResults() const;
   virtual void DoCIS();
   virtual void OutputCISResults() const;
   double** GetForce(int elecState);
   double*** GetForce(const std::vector<int>& elecStates);
   double GetElectronicEnergy(int elecState) const;
   double GetCoreRepulsionEnergy() const;
   double GetVdWCorrectionEnergy() const;
   MolDS_base::TheoryType GetTheoryType() const;
protected:
   std::string errorMessageAtomA;
   std::string errorMessageAtomB;
   std::string errorMessageAtomType;
   std::string errorMessageOrbitalType;
   std::string errorMessageCartesianType;
   std::string errorMessageSCFNotConverged;
   std::string errorMessageMoleculeNotSet;
   std::string errorMessageOddTotalValenceElectrions;
   std::string errorMessageNotEnebleAtomType;
   std::string errorMessageCoulombInt;
   std::string errorMessageExchangeInt;
   std::string errorMessageMolecularIntegralElement;
   std::string errorMessageGetGaussianCartesianMatrixBadOrbital;
   std::string errorMessageGetGaussianOverlapBadOrbital;
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
   double vdWCorrectionEnergy;
   double** fockMatrix;
   double* energiesMO;
   double*** matrixForce;
   double****** twoElecTwoCore;
   double** orbitalElectronPopulation; //P_{\mu\nu} of (2.50) in J. A. Pople book.
   double*   atomicElectronPopulation; //P_{AB} of (3.21) in J. A. Pople book.
   double** matrixCIS;
   double* excitedEnergies;
   double* freeExcitonEnergiesCIS;
   int matrixCISdimension;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcSCFProperties();
   virtual void CalcCISProperties();
   double GetBondingAdjustParameterK(MolDS_base::ShellType shellA, 
                                     MolDS_base::ShellType shellB) const;
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int indexAtomA, 
                                                        int indexAtomB, 
                                                        MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomVdWCorrectionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomVdWCorrectionFirstDerivative(int indexAtomA, 
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
   virtual void CalcForce(const std::vector<int>& elecStates);
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
   struct CISEigenVectorCoefficient{
      double coefficient;
      int occIndex;
      int virIndex;
      int slaterIndex;
   };
   struct MoreCISEigenVectorCoefficient { 
      bool operator()(const CISEigenVectorCoefficient& rLeft, const CISEigenVectorCoefficient& rRight) 
      const { return fabs(rLeft.coefficient) > fabs(rRight.coefficient); } 
   };
private:
   std::string errorMessageCalDiaOverlapDiaFrameNullMatrix;
   std::string errorMessageCalcRotatingMatrixNullRotMatrix;
   std::string errorMessageRotDiaOverlapToSpaceFrameNullDiaMatrix;
   std::string errorMessageRotDiaOverlapToSpaceFrameNullRotMatrix;
   std::string errorMessageSetOverlapElementNullDiaMatrix;
   std::string messageIterSCF;
   std::string messageDensityRMS;
   std::string messageEnergiesMOs;
   std::string messageEnergiesMOsTitle;
   std::string messageMullikenAtoms;
   std::string messageMullikenAtomsTitle;
   std::string messageElecEnergy;
   std::string messageElecEnergyVdW;
   std::string messageElecEnergyTitle;
   std::string messageOcc;
   std::string messageUnOcc;
   std::string messageCoreRepulsionTitle;
   std::string messageCoreRepulsion;
   std::string messageVdWCorrectionTitle;
   std::string messageVdWCorrection;
   std::string messageElectronicDipoleTitleAU;
   std::string messageElectronicDipoleTitleDebye;
   std::string messageElectronicDipole;
   std::string messageCoreDipoleTitleAU;
   std::string messageCoreDipoleTitleDebye;
   std::string messageCoreDipole;
   std::string messageTotalDipoleTitleAU;
   std::string messageTotalDipoleTitleDebye;
   std::string messageTotalDipole;
   double elecSCFEnergy;
   double** gammaAB;
   double** overlap; // overlap integral between AOs
   double*** cartesianMatrix; // cartesian matrix represented by AOs
   double** electronicDipoleEigenStates; // electronic dipole moment of eqch eigen states.
   double* coreDipole; // electronic dipole moment of eqch eigen states.
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
   void OutputMOEnergies() const;
   void OutputSCFEnergies() const;
   void OutputSCFDipole() const;
   void OutputSCFMulliken() const;
   void CalcCoreRepulsionEnergy();
   void CalcVdWCorrectionEnergy();
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
   void CalcCoreDipole(double* coreDipole,
                       const MolDS_base::Molecule& molecule) const;
   void CalcElectronicDipoleGroundState(double** electronicDipoleEigenStates,
                                        double const* const* const* cartesianMatrix,
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* orbitalElectronPopulation) const;
   void CalcCartesianMatrixByGTOExpansion(double*** cartesianMatrix,
                                          const MolDS_base::Molecule& molecule, 
                                          MolDS_base::STOnGType stonG) const; 
   double GetCartesianMatrixElementByGTOExpansion(const MolDS_base_atoms::Atom& atomA, 
                                                  int valenceIndexA, 
                                                  const MolDS_base_atoms::Atom& atomB, 
                                                  int valenceIndexB,
                                                  MolDS_base::CartesianType axis,
                                                  MolDS_base::STOnGType stonG) const;
   double GetGaussianCartesianMatrix(MolDS_base::AtomType atomTypeA, 
                                     MolDS_base::OrbitalType valenceOrbitalA, 
                                     double gaussianExponentA, 
                                     double const* xyzA,
                                     MolDS_base::AtomType atomTypeB, 
                                     MolDS_base::OrbitalType valenceOrbitalB, 
                                     double gaussianExponentB, 
                                     double const* xyzB,
                                     double Rab,
                                     MolDS_base::CartesianType axis) const;
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
   void DoDamp(double rmsDensity, 
               double** orbitalElectronPopulation, 
               double const* const* oldOrbitalElectronPopulation, 
               const MolDS_base::Molecule& molecule) const;
   void DoDIIS(double** orbitalElectronPopulation,
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
   void CalcElecSCFEnergy(double* elecSCFEnergy, 
                         const MolDS_base::Molecule& molecule, 
                         double const* energiesMO, 
                         double const* const* fockMatrix, 
                         double const* const* gammaAB, 
                         double coreRepulsionEnergy,
                         double vdWCorrectionEnergy) const;
   void FreeElecEnergyMatrices(double*** fMatrix, 
                               double*** hMatrix, 
                               double*** dammyOrbitalElectronPopulation, 
                               double**  dammyAtomicElectronPopulation ) const;
   void FreeSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                 double**** diisStoredDensityMatrix,
                                 double**** diisStoredErrorVect,
                                 double*** diisErrorProducts,
                                 double** diisErrorCoefficients) const;
   void MallocSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                   double**** diisStoredDensityMatrix,
                                   double**** diisStoredErrorVect,
                                   double*** diisErrorProducts,
                                   double** diisErrorCoefficients);
};


}
#endif



