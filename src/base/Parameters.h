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
#ifndef INCLUDED_PARAMETERS
#define INCLUDED_PARAMETERS

namespace MolDS_base{

// Parameters is singleton
class Parameters: private Uncopyable{
public:
   static Parameters* GetInstance();
   static void DeleteInstance();

   SimulationType GetCurrentSimulation() const;
   void SetCurrentSimulation(SimulationType simulation);
   TheoryType GetCurrentTheory() const;
   void SetCurrentTheory(TheoryType theory);
   // Pysical constants
   double GetEV2AU() const;
   double GetKcalMolin2AU() const;
   double GetAngstrom2AU() const;
   double GetKayser2AU() const;
   double GetGMolin2AU() const;
   double GetDegree2Radian() const;
   double GetFs2AU() const;
   double GetBoltzmann() const;
   // SCF
   double GetThresholdSCF() const;
   void SetThresholdSCF(double thresholdSCF);
   int GetMaxIterationsSCF() const;
   void SetMaxIterationsSCF(int maxIterationsSCF);
   double GetDampingThreshSCF() const;
   void SetDampingThreshSCF(double dampingThreshSCF);
   double GetDampingWeightSCF() const;
   void SetDampingWeightSCF(double dampingWeightSCF);
   int GetDiisNumErrorVectSCF() const;
   void SetDiisNumErrorVectSCF(int diisNumErrorVectSCF);
   double GetDiisStartErrorSCF() const;
   void SetDiisStartErrorSCF(double diisStartErrorSCF);
   double GetDiisEndErrorSCF() const;
   void SetDiisEndErrorSCF(double diisEndErrorSCF);
   // MOPlot
   std::string GetFileNamePrefixMOPlot() const;
   void SetFileNamePrefixMOPlot(std::string fileNamePrefixMOPlot);
   int* GetGridNumberMOPlot() const;
   void SetGridNumberMOPlot(int Nx, int Ny, int Nz);
   double* GetFrameLengthMOPlot() const;
   void SetFrameLengthMOPlot(double lx, double ly, double lz);
   std::vector<int>* GetIndecesMOPlot() const;
   void AddIndexMOPlot(int moIndex);
   bool RequiresMOPlot() const;
   // Translation
   void SetTranslatingDifference(double x, double y, double z);
   double* GetTranslatingDifference() const;
   // Principal axes
   void SetInertiaTensorOrigin(double x, double y, double z);
   double* GetInertiaTensorOrigin() const;
   // Rotation
   void SetRotatingOrigin(double x, double y, double z);
   double* GetRotatingOrigin() const;
   void SetRotatingType(RotatingType rotatingType);
   RotatingType GetRotatingType() const;
   void SetRotatingAxis(double x, double y, double z);
   double* GetRotatingAxis() const;
   void SetRotatingAngle(double rotatingAngle);
   double GetRotatingAngle() const;
   void SetRotatingEularAngles(double alpha, double beta, double gamma);
   EularAngle GetRotatingEularAngles() const;
   // CIS
   int GetActiveOccCIS() const;
   void SetActiveOccCIS(int activeOccCIS);
   int GetActiveVirCIS() const;
   void SetActiveVirCIS(int activeOccCIS);
   int GetNumberExcitedStatesCIS() const;
   void SetNumberExcitedStatesCIS(int nStates);
   bool RequiresCIS() const;
   void SetRequiresCIS(bool requiresCIS);
   bool IsDavidsonCIS() const;
   void SetIsDavidsonCIS(bool isDavidsonCIS);
   int GetMaxIterationsCIS() const;
   void SetMaxIterationsCIS(int maxIterationsCIS);
   int GetMaxDimensionsCIS() const;
   void SetMaxDimensionsCIS(int maxDimensionsCIS);
   double GetNormToleranceCIS() const;
   void SetNormToleranceCIS(double normToleranceCIS);
   // Memory
   double GetLimitHeapMemory() const;
   void SetLimitHeapMemory(double limitHeap);
   // MD
   int GetElectronicStateIndexMD() const;
   void SetElectronicStateIndexMD(int electronicStateIndex);
   int GetTotalStepsMD() const;
   void SetTotalStepsMD(int totalSteps);
   double GetTimeWidthMD() const;
   void SetTimeWidthMD(double timeWidth);
   // MC
   int GetElectronicStateIndexMC() const;
   void SetElectronicStateIndexMC(int electronicStateIndex);
   int GetTotalStepsMC() const;
   void SetTotalStepsMC(int totalSteps);
   double GetTemperatureMC() const;
   void SetTemperatureMC(double temperature);
   double GetStepWidthMC() const;
   void SetStepWidthMC(double stepWidth);
   unsigned long GetSeedMC() const;
   void SetSeedMC(unsigned long seed);
   // RPMD
   int GetElectronicStateIndexRPMD() const;
   void SetElectronicStateIndexRPMD(int electronicStateIndex);
   int GetNumberElectronicStatesRPMD() const;
   void SetNumberElectronicStatesRPMD(int NumberElectronicStates);
   int GetTotalStepsRPMD() const;
   void SetTotalStepsRPMD(int totalSteps);
   double GetTemperatureRPMD() const;
   void SetTemperatureRPMD(double temperature);
   double GetTimeWidthRPMD() const;
   void SetTimeWidthRPMD(double stepWidth);
   int GetNumberBeadsRPMD() const;
   void SetNumberBeadsRPMD(int numberBeads);
   unsigned long GetSeedRPMD() const;
   void SetSeedRPMD(unsigned long seed);
   // Opt (steepest descent)
   int GetLineReturnTimesSteepestDescent();
   void SetLineReturnTimesSteepestDescent(int lineReturnTimes);
   int GetStepsSteepestDescent() const;
   void SetStepsSteepestDescent(int steps);
   int GetElectronicStateIndexSteepestDescent() const;
   void SetElectronicStateIndexSteepestDescent(int electronicStateIndex);
   double GetMaxGradientSteepestDescent() const;
   void SetMaxGradientSteepestDescent(double maxGradient);
   double GetRmsGradientSteepestDescent() const;
   void SetRmsGradientSteepestDescent(double rmsGradient);
   double GetTimeWidthSteepestDescent() const;
   void SetTimeWidthSteepestDescent(double timeWidth);
private:
   static Parameters* parameters;
   Parameters();
   ~Parameters();
   std::string errorMessageGetIndecesMOPlotNull;
   SimulationType currentSimulation;
   TheoryType currentTheory;
   // Physical constants
   static const double eV2AU;
   static const double kcalMolin2AU;
   static const double angstrom2AU;
   static const double kayser2AU;
   static const double gMolin2AU;
   static const double degree2Radian;
   static const double fs2AU;
   static const double boltzmann;
   // SCF
   double thresholdSCF;
   int maxIterationsSCF;
   double dampingThreshSCF;
   double dampingWeightSCF;
   int diisNumErrorVectSCF;
   double diisStartErrorSCF;
   double diisEndErrorSCF;
   // MOPlot
   std::string fileNamePrefixMOPlot;
   int gridNumberMOPlot[CartesianType_end];
   double frameLengthMOPlot[CartesianType_end];
   std::vector<int>* indecesMOPlot;
   // Translation
   double translatingDifference[3];
   // Principal axes
   double* inertiaTensorOrigin;
   // Rotation
   double* rotatingOrigin;
   double rotatingAxis[3];
   double rotatingAngle;
   RotatingType rotatingType;
   EularAngle rotatingEularAngles;
   // CIS
   int activeOccCIS;
   int activeVirCIS;
   int numberExcitedStatesCIS;
   int maxIterationsCIS;
   int maxDimensionsCIS;
   double normToleranceCIS;
   bool requiresCIS;
   bool isDavidsonCIS;
   // Memory
   double limitHeapMemory;
   // MD
   int electronicStateIndexMD;
   int totalStepsMD;
   double timeWidthMD;
   // MC
   int electronicStateIndexMC;
   int totalStepsMC;
   double temperatureMC;
   double stepWidthMC;
   unsigned long seedMC;
   // RPMD
   int electronicStateIndexRPMD;
   int numberElectronicStatesRPMD;
   int totalStepsRPMD;
   double temperatureRPMD;
   double timeWidthRPMD;
   int numberBeadsRPMD;
   unsigned long seedRPMD;
   // Opt (steepest descent)
   int lineReturnTimesSteepestDescent;
   int stepsSteepestDescent;
   int electronicStateIndexSteepestDescent;
   double maxGradientSteepestDescent;
   double rmsGradientSteepestDescent;
   double timeWidthSteepestDescent;
   // Other
   void SetDefaultValues();
   void SetMessages();
};

}
#endif





