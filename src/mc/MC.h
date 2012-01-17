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
   void SetMolecule(MolDS_base::Molecule* molecule);
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
   MolDS_base::Molecule* molecule;
   void SetMessages();
   void CreateTrialConfiguration(MolDS_base::Molecule* trial,
                                 MolDS_base::Molecule* current,
                                 boost::random::variate_generator<
                                    boost::random::mt19937&,
                                    boost::uniform_real<>
                                 > (*realRand),
                                 boost::random::variate_generator<
                                    boost::random::mt19937&,
                                    boost::uniform_smallint<>
                                 > (*intRand));
   void SynchronousMolecularConfiguration(MolDS_base::Molecule* target, 
                                          MolDS_base::Molecule* refference) const;
   bool UsesTrial(MolDS_base::ElectronicStructure* currentES, 
                  MolDS_base::ElectronicStructure* trialES) const;
   void OutputEnergies(MolDS_base::ElectronicStructure* electronicStructure) const;
};

}
#endif



