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
#ifndef INCLUDED_OPTIMIZER
#define INCLUDED_OPTIMIZER
namespace MolDS_optimization{

class Optimizer : public MolDS_base::PrintController{
public:
   Optimizer();
   virtual ~Optimizer();
   void Optimize(MolDS_base::Molecule& molecule);
protected:
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageGeometyrOptimizationNotConverged;
   std::string messageLineSearchSteps;
   virtual void SetMessages();
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double** matrixForce, double dt) const;
   void UpdateElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                  MolDS_base::Molecule& molecule,
                                  bool requireGuess, 
                                  bool printsLogs) const;
   bool SatisfiesConvergenceCriterion(double** matrixForce, 
                                      const MolDS_base::Molecule& molecule,
                                      double oldEnergy,
                                      double currentEnergy,
                                      double maxGradientThreshold,
                                      double rmsGradientThreshold) const;
   void OutputMoleculeElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                          MolDS_base::Molecule& molecule,
                                          bool printsLogs) const;
private:
   std::string errorMessageTheoryType;
   std::string errorMessageTotalSteps;
   std::string messageGeometyrOptimizationMetConvergence;
   std::string messageStartGeometryOptimization;
   std::string messageEndGeometryOptimization;
   std::string messageReducedTimeWidth;
   std::string messageOptimizationLog;
   std::string messageEnergyDifference;
   std::string messageMaxGradient;
   std::string messageRmsGradient;
   std::string messageAu;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void SetEnableTheoryTypes();
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType) const;
   void ClearMolecularMomenta(MolDS_base::Molecule& molecule) const;
   virtual void LineSearch(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                   MolDS_base::Molecule& molecule,
                   double* lineSearchedEnergy,
                   bool* obainesOptimizedStructure) const = 0;
};

}
#endif


