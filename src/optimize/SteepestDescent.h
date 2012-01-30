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
#ifndef INCLUDED_STEEPEST_DESCENT
#define INCLUDED_STEEPEST_DESCENT
namespace MolDS_optimize{

/***
 *  Line search algorythm is used in the early stage of the optimization.
 *  Then, steepest descent algorythm is used to get final optimized struture.
 */
class SteepestDescent : public MolDS_base::PrintController{
public:
   SteepestDescent();
   ~SteepestDescent();
   void Optimize(MolDS_base::Molecule& molecule);
private:
   std::string messageStartSteepestDescent;
   std::string messageEndSteepestDescent;
   std::string messageStartStepSteepestDescent;
   std::string messageEndStepSteepestDescent;
   std::string messageStartLineSearch;
   std::string messageEndLineSearch;
   std::string messageStartLineReturnTimes;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreRepulsionEnergy;
   std::string messageElectronicEnergy;
   std::string messageTotalEnergy;
   std::string messageDifferentEnergy;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void SetMessages();
   void SetEnableTheoryTypes();
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType) const;
   void ClearMolecularMomenta(MolDS_base::Molecule& molecule) const;
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double** matrixForce, double dt) const;
   void UpdateElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                  bool requireGuess, 
                                  bool printsLogs) const;
   void LineSearch(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, MolDS_base::Molecule& molecule) const;
   void SteepestDescentSearch(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, MolDS_base::Molecule& molecule) const;
};

}
#endif



