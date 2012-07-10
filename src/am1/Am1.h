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
#ifndef INCLUDED_AM1
#define INCLUDED_AM1
namespace MolDS_am1{

/***
 *  Main Refferences for AM1 are [DZHS_1985, DY_1990]
 */
class Am1 : public MolDS_mndo::Mndo{
public:
   Am1();
   virtual ~Am1();
   virtual void OutputSCFResults() const;
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcSCFProperties();
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomCoreRepulsionSecondDerivative(int atomAIndex,
                                                         int atomBIndex, 
                                                         MolDS_base::CartesianType axisA1,
                                                         MolDS_base::CartesianType axisA2) const;
private:
   double GetAdditionalDiatomCoreRepulsionTerm(double k, double l, double m, double distance) const;
   double GetAdditionalDiatomCoreRepulsionTermFirstDerivative(double k, double l, double m, double distance) const;
   double GetAdditionalDiatomCoreRepulsionTermSecondDerivative(double k, double l, double m, double distance) const;
};

}
#endif



