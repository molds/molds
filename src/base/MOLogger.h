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
#ifndef INCLUDED_MOLOGGER
#define INCLUDED_MOLOGGER
namespace MolDS_base{

class MOLogger: public PrintController{
public:
   MOLogger(const MolDS_base::Molecule& molecule, double const* const* fockMatrix, MolDS_base::TheoryType theory);
   void DrawMO(int moIndex);
   void DrawMO(std::vector<int> moIndeces);
private:
   std::string stringCubeExtension;
   std::string messageCubeHeaderComment1;
   std::string messageCubeHeaderComment2;
   std::string messageStartMOPlot;
   std::string messageEndMOPlot;
   std::string messageSkippedMOIndex;
   MOLogger();
   MolDS_base::Molecule const* molecule;
   double const* const* fockMatrix;
   MolDS_base::TheoryType theory;
   void SetMessage();
};

}
#endif
