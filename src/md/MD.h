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
#ifndef INCLUDED_MD
#define INCLUDED_MD
namespace MolDS_md{

/***
 *  Velocty Verlet is used here.
 */
class MD{
public:
   MD();
   ~MD();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
   void DoMD();
private:
   std::string messageinitialConditionMD;
   std::string messageStartMD;
   std::string messageEndMD;
   std::string messageStartStepMD;
   std::string messageEndStepMD;
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
   MolDS_cndo::Cndo2* cndo;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
   void OutputEnergies(double initialEnergy);
   double OutputEnergies();
};

}
#endif



