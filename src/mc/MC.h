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
#ifndef INCLUDED_MC
#define INCLUDED_MC
namespace MolDS_mc{

/***
 *  Velocty Verlet is used here.
 */
class MC{
public:
   MC();
   ~MC();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
   void DoMC();
private:
   std::string messageinitialConditionMC;
   std::string messageStartMC;
   std::string messageEndMC;
   std::string messageStartStepMC;
   std::string messageEndStepMC;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreRepulsionEnergy;
   std::string messageElectronicEnergy;
   std::string messageTotalEnergy;
   std::string errorMessageNotEnebleExcitedTheoryType;
   std::string errorMessageTheoryType;
   MolDS_cndo::Cndo2* cndo;
   void SetMessages();
   void OutputEnergies();
};

}
#endif



