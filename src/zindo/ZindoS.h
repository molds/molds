#ifndef INCLUDED_ZINDOS
#define INCLUDED_ZINDOS
namespace MolDS_zindo{

/***
 *  Main Refference for Zindo is [RZ_1973]
 */
class ZindoS : public MolDS_cndo::Cndo2{
public:
   ZindoS();
   virtual ~ZindoS();
   virtual void DoesCIS();
protected:
   std::string errorMessageDavidsonNotConverged;
   std::string errorMessageCalcCISMatrix;
   std::string messageStartCIS;
   std::string messageDoneCIS;
   std::string messageDavidsonConverge;
   std::string messageStartCalcCISMatrix;
   std::string messageOmpElapsedTimeCalcCISMarix;
   std::string messageOmpElapsedTimeCIS;
   std::string messageDoneCalcCISMatrix;
   virtual void CalcGammaAB(double** gammaAB, const MolDS_base::Molecule& molecule) const;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
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
                                const MolDS_base_atoms::Atom& atom) const; // Apendix in [BZ_1979]
   virtual double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                                 MolDS_base::OrbitalType orbital2, 
                                 const MolDS_base_atoms::Atom& atom) const; // Apendix in [BZ_1979]
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcCISMatrix(double** matrixCIS, 
                              int numberActiveOcc, 
                              int numberActiveVir) const;
   virtual void CalcForce(std::vector<int> elecStates);
   int GetSlaterDeterminantIndex(int activeOccIndex, int activeVirIndex) const;
   void CheckMatrixForce(std::vector<int> elecStates);
private:
   std::string errorMessageCalcForceNotGroundState;
   std::string errorMessageElecState;
   std::string errorMessageNishimotoMataga;
   std::string errorMessageDavidsonMaxIter;
   std::string errorMessageDavidsonMaxDim;
   std::string messageStartDirectCIS;
   std::string messageDoneDirectCIS;
   std::string messageStartDavidsonCIS;
   std::string messageDoneDavidsonCIS;
   std::string messageNumIterCIS;
   std::string messageResidualNorm;
   std::string messageDavidsonReachCISMatrix;
   std::string messageDavidsonGoToDirect;
   std::string messageExcitedStatesEnergies;
   std::string messageExcitedStatesEnergiesTitle;
   int matrixForceElecStatesNum;
   double nishimotoMatagaParamA;
   double nishimotoMatagaParamB;
   double overlapCorrectionSigma;
   double overlapCorrectionPi;
   double GetNishimotoMatagaTwoEleInt(const MolDS_base_atoms::Atom& atomA, 
                                      MolDS_base::OrbitalType orbitalA, 
                                      const MolDS_base_atoms::Atom& atomB, 
                                      MolDS_base::OrbitalType orbitalB) const; // ref. [MN_1957] and (5a) in [AEZ_1986]
   double GetNishimotoMatagaTwoEleIntFirstDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                     MolDS_base::OrbitalType orbitalA, 
                                                     const MolDS_base_atoms::Atom& atomB, 
                                                     MolDS_base::OrbitalType orbitalB,
                                                     MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
   void DoesCISDirect();
   void DoesCISDavidson();
   void CalcRitzVector(double* ritzVector, 
                       double const* const* expansionVectors, 
                       double const* const* interactionMatrix, 
                       int interactionMatrixDimension, 
                       int ritzVectorIndex) const;
   void CalcResidualVectorAndNorm(double* residualVector, 
                                  double* norm, 
                                  double const* ritzVector, 
                                  double const* interactionEigenEnergies, 
                                  int residualVectorIndex) const;
   void SortSingleExcitationSlaterDeterminants(std::vector<MoEnergyGap>* moEnergyGaps) const;
   void UpdateExpansionVectors(double** expansionVectors, 
                               int* notConvergedStates, 
                               double const* interactionEigenEnergies, 
                               double const* residualVector,
                               int interactionMatrixDimension, 
                               int residualVectorIndex) const;
   void CalcInteractionMatrix(double** interactionMatrix, 
                              double const* const* expansionVectors, 
                              int interactionMatrixDimension) const;
   void FreeDavidsonCISTemporaryMtrices(double*** expansionVectors, 
                                        double** residualVector, 
                                        double** ritzVector) const;
   void FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix, 
                                            double interactionMatrixDimension, 
                                            double** interactionEigenEnergies) const;
};

}
#endif



