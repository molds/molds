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
#ifndef INCLUDED_HOLE_DENSITY_LOGGER
#define INCLUDED_HOLE_DENSITY_LOGGER
namespace MolDS_base{

class HoleDensityLogger: public PrintController{
public:
   HoleDensityLogger(const MolDS_base::Molecule& molecule, 
                     double const* const* fockMatrix, 
                     double const* const* cisMatrix, 
                     MolDS_base::TheoryType theory);
   void DrawDensity(int elecStateIndex);
   void DrawDensity(std::vector<int> elecStateIndeces);
private:
   std::string errorMessageCISMatrixNULL;
   std::string errorMessageFockMatrixNULL;
   std::string messageCubeHeaderComment1;
   std::string messageCubeHeaderComment2;
   std::string messageStartHoleDensityPlot;
   std::string messageEndHoleDensityPlot;
   std::string messageSkippedElecStateIndex;
   std::string stringCubeExtension;
   HoleDensityLogger();
   MolDS_base::Molecule const* molecule;
   double const* const* fockMatrix;
   double const* const* cisMatrix;
   MolDS_base::TheoryType theory;
   void MatricesNullCheck() const;
   void OutputHeaderToFile(std::ofstream& ofs, 
                           double const* origin, 
                           int const* gridNumber,
                           double dx, 
                           double dy,
                           double dz) const;
   void OutputMoleculeToFile(std::ofstream& ofs, const MolDS_base::Molecule& molecule)const ;
   double GetMOValue(int moIndex, const MolDS_base::Molecule& molecule, double x, double y, double z) const;
   void SetMessage();
};

}
#endif
