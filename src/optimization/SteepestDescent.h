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
#ifndef INCLUDED_STEEPEST_DESCENT
#define INCLUDED_STEEPEST_DESCENT
namespace MolDS_optimization{

class SteepestDescent : public MolDS_optimization::Optimizer{
public:
   SteepestDescent();
   ~SteepestDescent();
protected:
   void SetMessages();
private:
   std::string messageStartSteepestDescentStep;
   const std::string& OptimizationStepMessage() const{
      return this->messageStartSteepestDescentStep;
   }
   void SearchMinimum(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                      MolDS_base::Molecule& molecule,
                      double* lineSearchedEnergy,
                      bool* obainesOptimizedStructure) const;
   void InitializeState(OptimizerState&, const MolDS_base::Molecule&) const{}
   virtual void PrepareState(OptimizerState& state,
                             const MolDS_base::Molecule& molecule,
                             const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                             const int elecState) const{};
   void UpdateSearchDirection(OptimizerState& state,
                              boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                              const MolDS_base::Molecule& molecule,
                              int elecState) const;
};

}
#endif



