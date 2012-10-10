//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   //
// Copyright (C) 2012-2012 Katsuhiko Nishimra                             //
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
#ifndef INCLUDED_GDIIS
#define INCLUDED_GDIIS
namespace MolDS_optimization{

class GDIIS : public MolDS_base::PrintController{
public:
   GDIIS(int sizeErrorVector);
	 ~GDIIS();
   bool DoGDIIS(double* vectorError,
                double* vectorPosition,
                double const* vectorRefStep);
   bool DoGDIIS(double* vectorError,
                MolDS_base::Molecule& molecule,
                double const* vectorRefStep);

   void DiscardPrevious();
   void DiscardOldest();
   typedef std::list<double*> list;
   typedef std::list<double*>::iterator iterator;
private:
   const int sizeErrorVector;
   const int maxnumErrors;
   double** matrixGDIIS;
   list listErrors;
   list listPositions;
   double MinCosine();
   std::string messageTakingGDIISStep;
   std::string messageSingularGDIISMatrix;
   boost::format formatTooSmallLagrangeMultiplier;
   boost::format formatTooLargeGDIISStep;
   boost::format formatWrongDirection;
};

}
#endif
