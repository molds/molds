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
#ifndef INCLUDED_RPMD
#define INCLUDED_RPMD
namespace MolDS_rpmd{

/***
 *  Velocty Verlet is used here.
 */
class RPMD{
public:
   RPMD();
   ~RPMD();
   void DoRPMD(const MolDS_base::Molecule&);
private:
   std::string messageinitialConditionRPMD;
   std::string messageStartRPMD;
   std::string messageEndRPMD;
   std::string messageStartStepRPMD;
   std::string messageEndStepRPMD;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreKineticEnergy;
   std::string messageCoreRepulsionEnergy;
   std::string messageElectronicEnergy;
   std::string messageTotalEnergy;
   std::string messageErrorEnergy;
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   std::string errorMessageElecState;
   std::vector<MolDS_base::TheoryType> enableGroundStateTheoryTypes;
   std::vector<MolDS_base::TheoryType> enableExcitedStateTheoryTypes;
   void SetMessages();
   void SetEnableTheoryTypes();
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType, int elecState);
   //void OutputEnergies(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, double initialEnergy);
   //double OutputEnergies(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure);
};

}
#endif


