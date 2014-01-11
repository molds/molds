//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
protected:
   class OptimizerState{
   protected:
      double currentEnergy;
      double initialEnergy;
      double const* const* matrixForce;
      std::string errorMessageFailedToDowncastState;
      virtual void SetMessages();
   public:
      OptimizerState():
         currentEnergy(0.0), initialEnergy(0.0), matrixForce(NULL){this->SetMessages();}
      virtual ~OptimizerState(){}
      double& GetCurrentEnergyRef(){return this->currentEnergy;}
      double GetCurrentEnergy(){return this->currentEnergy;}
      double GetInitialEnergy(){return this->initialEnergy;}
      double const* const*  GetMatrixForce(){return this->matrixForce;}
      double const* const** GetMatrixForcePtr(){return &this->matrixForce;}
      void SetCurrentEnergy(double currentEnergy){this->currentEnergy = currentEnergy;}
      void SetInitialEnergy(double initialEnergy){this->initialEnergy = initialEnergy;}
      void SetMatrixForce(double const* const* matrixForce){this->matrixForce = matrixForce;}
      template<class State>
      State& CastRef(){
         try{
            return dynamic_cast<State&>(*this);
         }
         catch(std::bad_cast& ex){
            throw MolDS_base::MolDSException(this->errorMessageFailedToDowncastState);
         }
      }
   };
public:
   Optimizer();
   virtual ~Optimizer();
   void Optimize(MolDS_base::Molecule& molecule);
protected:
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageGeometyrOptimizationNotConverged;
   std::string messageLineSearchSteps;
   virtual void SetMessages();
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double const* const* matrixForce, double dt) const;
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double const* const* matrixForce) const;
   void UpdateElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                  MolDS_base::Molecule& molecule,
                                  bool requireGuess, 
                                  bool printsLogs) const;
   bool SatisfiesConvergenceCriterion(double const* const* matrixForce, 
                                      const MolDS_base::Molecule& molecule,
                                      double oldEnergy,
                                      double currentEnergy,
                                      double maxGradientThreshold,
                                      double rmsGradientThreshold) const;
   void OutputMoleculeElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                          MolDS_base::Molecule& molecule,
                                          bool printsLogs) const;
   void LineSearch(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                   MolDS_base::Molecule& molecule,
                   double &lineSearchCurrentEnergy,
                   double const* const* matrixForce,
                   int elecState,
                   double dt) const;
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
   virtual void SearchMinimum(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                              MolDS_base::Molecule& molecule,
                              double* lineSearchedEnergy,
                              bool* obainesOptimizedStructure) const = 0;
   virtual void InitializeState(OptimizerState &state, const MolDS_base::Molecule& molecule) const = 0;
};

}
#endif



