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

class GDIISException : public MolDS_base::MolDSException {
   public:
   GDIISException(const std::string& cause)
      : MolDS_base::MolDSException(boost::format("GDIISException: %s") % cause){
   }

   GDIISException(const boost::format& cause)
      : MolDS_base::MolDSException(boost::format("GDIISException: %s") % cause.str()){
   }
};

class GDIIS : public MolDS_base::PrintController{
public:
   GDIIS(int sizeErrorVector);
	 ~GDIIS();
   void DoGDIIS(double* vectorError,
                double* vectorPosition,
                double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException);
   void DoGDIIS(double* vectorError,
                MolDS_base::Molecule& molecule,
                double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException);

   void DiscardPrevious();
   void DiscardOldest();
   typedef std::list<double*> list;
   typedef std::list<double*>::iterator iterator;
private:
   void Update(double const* vectorError,
               double const* vectorPosition);
   void CalcGDIIS(double*       vectorError,
                  double*       vectorPosition,
                  double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException);
   void RecalcGDIIS(double*       vectorError,
                    double*       vectorPosition,
                    double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException);
   void PopOldest(double** vectorError, double** vectorPosition);
   void PushOldest(double* vectorError, double* vectorPosition);
   const int sizeErrorVector;
   const int maxnumErrors;
   double** matrixGDIIS;
   list listErrors;
   list listPositions;
   double MinCosine();
   std::string messageTakingGDIISStep;
   std::string messageSingularGDIISMatrix;
   std::string messageOnlyOneErrorVector;
   std::string messageRecalcGDIISStep;
   boost::format formatTooSmallLagrangeMultiplier;
   boost::format formatTooLargeGDIISStep;
   boost::format formatWrongDirection;
};

}
#endif
