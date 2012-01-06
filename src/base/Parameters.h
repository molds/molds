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
   std::string GetFileNamePrefixMOPlot() const;
   void SetFileNamePrefixMOPlot(std::string fileNamePrefixMOPlot);
   int* GetGridNumberMOPlot() const;
   void SetGridNumberMOPlot(int Nx, int Ny, int Nz);
   double* GetFrameLengthMOPlot() const;
   void SetFrameLengthMOPlot(double lx, double ly, double lz);
   std::vector<int>* GetIndecesMOPlot() const;
   void AddIndexMOPlot(int moIndex);
   double GetEV2AU() const;
   double GetKcalMolin2AU() const;
   double GetAngstrom2AU() const;
   double GetKayser2AU() const;
   double GetGMolin2AU() const;
   double GetDegree2Radian() const;
   double GetFs2AU() const;
   TheoryType GetCurrentTheory() const;
   void SetCurrentTheory(TheoryType theory);
   void SetTranslatingDifference(double x, double y, double z);
   double* GetTranslatingDifference() const;
   void SetInertiaTensorOrigin(double x, double y, double z);
   double* GetInertiaTensorOrigin() const;
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
   double GetLimitHeapMemory() const;
   void SetLimitHeapMemory(double limitHeap);
   bool RequiresMD() const;
   void SetRequiresMD(bool requiresMD);
   bool RequiresMOPlot() const;
   int GetElectronicStateIndexMD() const;
   void SetElectronicStateIndexMD(int electronicStateIndex);
   int GetTotalStepsMD() const;
   void SetTotalStepsMD(int totalSteps);
   double GetTimeWidthMD() const;
   void SetTimeWidthMD(double timeWidth);
private:
   static Parameters* parameters;
   Parameters();
   ~Parameters();
   std::string errorMessageGetIndecesMOPlotNull;
   void SetDefaultValues();
   void SetMessages();
   static const double eV2AU;
   static const double kcalMolin2AU;
   static const double angstrom2AU;
   static const double kayser2AU;
   static const double gMolin2AU;
   static const double degree2Radian;
   static const double fs2AU;
   double thresholdSCF;
   int maxIterationsSCF;
   double dampingThreshSCF;
   double dampingWeightSCF;
   int diisNumErrorVectSCF;
   double diisStartErrorSCF;
   double diisEndErrorSCF;
   std::string fileNamePrefixMOPlot;
   int gridNumberMOPlot[CartesianType_end];
   double frameLengthMOPlot[CartesianType_end];
   std::vector<int>* indecesMOPlot;
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
   double limitHeapMemory;
   bool requiresMD;
   int electronicStateIndexMD;
   int totalStepsMD;
   double timeWidthMD;
};

}
#endif





