#ifndef INCLUDED_PARAMETERS
#define INCLUDED_PARAMETERS

namespace MolDS_base{

// Parameters is singleton
class Parameters{
public:
   static Parameters* GetInstance();
   static void DeleteInstance();

   double GetThresholdSCF();
   void SetThresholdSCF(double thresholdSCF);
   int GetMaxIterationsSCF();
   void SetMaxIterationsSCF(int maxIterationsSCF);
   double GetDampingThreshSCF();
   void SetDampingThreshSCF(double dampingThreshSCF);
   double GetDampingWeightSCF();
   void SetDampingWeightSCF(double dampingWeightSCF);
   int GetDiisNumErrorVectSCF();
   void SetDiisNumErrorVectSCF(int diisNumErrorVectSCF);
   double GetDiisStartErrorSCF();
   void SetDiisStartErrorSCF(double diisStartErrorSCF);
   double GetDiisEndErrorSCF();
   void SetDiisEndErrorSCF(double diisEndErrorSCF);
   double GetEV2AU();
   double GetKcalMolin2AU();
   double GetAngstrom2AU();
   double GetKayser2AU();
   double GetGMolin2AU();
   double GetDegree2Radian();
   double GetFs2AU();
   double GetBondingAdjustParameterK();
   TheoryType GetCurrentTheory();
   void SetCurrentTheory(TheoryType theory);
   void SetTranslatingDifference(double x, double y, double z);
   double* GetTranslatingDifference();
   void SetInertiaTensorOrigin(double x, double y, double z);
   double* GetInertiaTensorOrigin();
   void SetRotatingOrigin(double x, double y, double z);
   double* GetRotatingOrigin();
   void SetRotatingType(RotatingType rotatingType);
   RotatingType GetRotatingType();
   void SetRotatingAxis(double x, double y, double z);
   double* GetRotatingAxis();
   void SetRotatingAngle(double rotatingAngle);
   double GetRotatingAngle();
   void SetRotatingEularAngles(double alpha, double beta, double gamma);
   EularAngle GetRotatingEularAngles();
   int GetActiveOccCIS();
   void SetActiveOccCIS(int activeOccCIS);
   int GetActiveVirCIS();
   void SetActiveVirCIS(int activeOccCIS);
   int GetNumberExcitedStatesCIS();
   void SetNumberExcitedStatesCIS(int nStates);
   bool RequiresCIS();
   void SetRequiresCIS(bool requiresCIS);
   bool IsDavidsonCIS();
   void SetIsDavidsonCIS(bool isDavidsonCIS);
   int GetMaxIterationsCIS();
   void SetMaxIterationsCIS(int maxIterationsCIS);
   int GetMaxDimensionsCIS();
   void SetMaxDimensionsCIS(int maxDimensionsCIS);
   double GetNormToleranceCIS();
   void SetNormToleranceCIS(double normToleranceCIS);
   bool RequiresMD();
   void SetRequiresMD(bool requiresMD);
   int GetElectronicStateIndexMD();
   void SetElectronicStateIndexMD(int electronicStateIndex);
   int GetTotalStepsMD();
   void SetTotalStepsMD(int totalSteps);
   double GetTimeWidthMD();
   void SetTimeWidthMD(double timeWidth);
private:
   static Parameters* parameters;
   Parameters();
   Parameters(Parameters&);
   void operator = (Parameters&);
   ~Parameters();

   double thresholdSCF;
   int maxIterationsSCF;
   double dampingThreshSCF;
   double dampingWeightSCF;
   int diisNumErrorVectSCF;
   double diisStartErrorSCF;
   double diisEndErrorSCF;
   void SetDefaultValues();
   double eV2AU;
   double kcalMolin2AU;
   double angstrom2AU;
   double kayser2AU;
   double gMolin2AU;
   double degree2Radian;
   double fs2AU;
   double bondingAdjustParameterK; //see (3.79) in J. A. Pople book
   TheoryType currentTheory;
   double translatingDifference[3];
   double* inertiaTensorOrigin;
   double* rotatingOrigin;
   double rotatingAxis[3];
   double rotatingAngle;
   RotatingType rotatingType;
   EularAngle rotatingEularAngles;
   int activeOccCIS;
   int activeVirCIS;
   int numberExcitedStatesCIS;
   int maxIterationsCIS;
   int maxDimensionsCIS;
   double normToleranceCIS;
   bool requiresCIS;
   bool isDavidsonCIS;
   bool requiresMD;
   int electronicStateIndexMD;
   int totalStepsMD;
   double timeWidthMD;
};

}
#endif





