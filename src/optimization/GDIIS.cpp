//************************************************************************//
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

#include<string>
#include<iostream>
#include<math.h>
#include<vector>
#include<stdexcept>
#include<list>
#include<boost/format.hpp>
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"../wrappers/Lapack.h"
#include"../base/Enums.h"
#include"../base/MallocerFreer.h"
#include"../base/EularAngle.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"GDIIS.h"

#define MAXNUMERRORS (5)

using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_optimization;

GDIIS::GDIIS(int sizeErrorVector):
   sizeErrorVector(sizeErrorVector),
   maxnumErrors(MAXNUMERRORS),
   matrixGDIIS(NULL),
   listErrors(),
   listPositions(),
   messageTakingGDIISStep("Taking GDIIS step.\n"),
   messageSingularGDIISMatrix("Error while solving GDIIS equation. Discarding current data.\n"),
   formatTooSmallLagrangeMultiplier("GDIIS: Lagrange Multiplier is too small. (%e)\n"),
   formatTooLargeGDIISStep("GDIIS: GDIIS step is too large. (gdiis:%e, reference:%e)\n"),
   formatWrongDirection("GDIIS: GDIIS step direction is too far from reference step. (cosine: %+f)\n")
{
   MallocerFreer::GetInstance()->Malloc(&this->matrixGDIIS,     maxnumErrors+1, maxnumErrors+1);
}

GDIIS::~GDIIS(){
   MallocerFreer::GetInstance()->Free(&this->matrixGDIIS,     maxnumErrors+1, maxnumErrors+1);
   for(GDIIS::iterator it = listErrors.begin(); it != listErrors.end(); it++){
      MallocerFreer::GetInstance()->Free(&*it, sizeErrorVector);
   }
   for(GDIIS::iterator it = listPositions.begin(); it != listPositions.end(); it++){
      MallocerFreer::GetInstance()->Free(&*it, sizeErrorVector);
   }
}

bool GDIIS::DoGDIIS(double* vectorError,
                    double* vectorPosition,
                    double const* vectorRefStep){
   this->Update(vectorError, vectorPosition);
   return this->CalcGDIIS(vectorError, vectorPosition, vectorRefStep);
}

void GDIIS::Update(double const* vectorError,
                   double const* vectorPosition){
   // Prepare GDIIS parameters
   double *tmp = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc(&tmp, sizeErrorVector);
      for(int i=0;i<sizeErrorVector;i++){
         tmp[i]=vectorError[i];
      }
      listErrors.push_back(tmp);
      tmp = NULL;
   }catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&tmp, sizeErrorVector);
      throw ex;
   }

   try{
      MallocerFreer::GetInstance()->Malloc(&tmp, sizeErrorVector);
      for(int i=0;i<sizeErrorVector;i++){
         tmp[i]=vectorPosition[i];
      }
      listPositions.push_back(tmp);
      tmp = NULL;
   }catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&tmp, sizeErrorVector);
      throw ex;
   }

   // Discard the oldest data if the size of the list exceeds the maximum
   if(listErrors.size() > maxnumErrors){
      this->DiscardOldest();
   }
}

bool GDIIS::CalcGDIIS(double* vectorError,
                      double* vectorPosition,
                      double const* vectorRefStep){
   // Prepare GDIIS matrix
   GDIIS::iterator it=listErrors.begin();
   for(int i=0; it!=listErrors.end();i++,it++){
      GDIIS::iterator it2=it;
      for(int j=i;it2!=listErrors.end();j++,it2++){
         matrixGDIIS[i][j] = 0.0;
         for(int k=0;k<this->sizeErrorVector;k++){
            matrixGDIIS[i][j] += (*it)[k] * (*it2)[k];
         }
         matrixGDIIS[j][i] = matrixGDIIS[i][j];
      }
   }
   const int numErrors = this->listErrors.size();
   for(int i=0;i<numErrors;i++){
      matrixGDIIS[i][numErrors] = matrixGDIIS[numErrors][i] = 1.0;
   }
   matrixGDIIS[numErrors][numErrors] = 0;

   // If only one error vector is given, following routine is meaningless.
   if(numErrors <= 1){
     return false;
   }

   double*  vectorCoefs = NULL;
   try{
      // Solve DIIS equation
      MallocerFreer::GetInstance()->Malloc(&vectorCoefs, numErrors+1);
      vectorCoefs[numErrors]=1.0;
      try{
         MolDS_wrappers::Lapack::GetInstance()->Dsysv(matrixGDIIS,
                                                      vectorCoefs,
                                                      numErrors+1);
      }
      catch(MolDSException ex){
         // Assume all errors to be due to singular GDIIS matrix.
         // Remove the newest data to eliminate singularity.
         this->DiscardPrevious();
         this->OutputLog(messageSingularGDIISMatrix);
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         return false;
      }

      // If Lagrange multiplier is too small;
      if(-vectorCoefs[numErrors] < 1e-8){
         this->OutputLog((formatTooSmallLagrangeMultiplier % -vectorCoefs[numErrors]).str());
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         // Recalculate GDIIS step without the oldest data.
         this->DiscardOldest();
         return CalcGDIIS(vectorError,vectorPosition,vectorRefStep);
      }

      // Interpolate error vectors and positions
      for(int i=0;i<sizeErrorVector;i++){
         vectorError[i] = vectorPosition[i] = 0;
         GDIIS::iterator it=listErrors.begin(),it2=listPositions.begin();
         for(int j=0; it!=listErrors.end();it++,it2++,j++){
            vectorError[i]    += vectorCoefs[j] * (*it)[i];
            vectorPosition[i] += vectorCoefs[j] * (*it2)[i];
         }
      }

      // Calculate cosine of the angle between GDIIS step and vectorRefStep
      // and lengths of the vectors.
      double innerprod = 0, normSquaregdiis = 0, normSquareref = 0;
      for(int i=0;i<this->sizeErrorVector;i++){
         double diff = vectorPosition[i] - listPositions.back()[i];
         innerprod  += diff * vectorRefStep[i];
         normSquaregdiis += diff * diff;
         normSquareref   += vectorRefStep[i] * vectorRefStep[i];
      }
      double cosine = innerprod/sqrt(normSquaregdiis*normSquareref);

      // If length of the GDIIS step is larger than reference step * 10
      if(normSquaregdiis >= normSquareref * 100){
         // Rollback vectorPosition and vectorError to original value
         for(int i=0; i<this->sizeErrorVector; i++){
            vectorError[i]    = listErrors.back()[i];
            vectorPosition[i] = listPositions.back()[i];
         }
         this->OutputLog((formatTooLargeGDIISStep % sqrt(normSquaregdiis) % sqrt(normSquareref)).str());
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         // and recalculate GDIIS step without the oldest data
         this->DiscardOldest();
         return CalcGDIIS(vectorError,vectorPosition,vectorRefStep);
      }

      // If the calculated cos(theta) on Eq. 8 of [FS_2002] is below the minimum tolerant value
      if(cosine < this->MinCosine()){
         // Rollback vectorPosition and vectorError to original value
         for(int i=0; i<this->sizeErrorVector; i++){
            vectorError[i]    = listErrors.back()[i];
            vectorPosition[i] = listPositions.back()[i];
         }
         this->OutputLog((formatWrongDirection % cosine).str());
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         // and recalculate GDIIS step without the oldest data
         this->DiscardOldest();
         return CalcGDIIS(vectorError,vectorPosition,vectorRefStep);
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
   this->OutputLog(messageTakingGDIISStep);
   return true;
}

bool GDIIS::DoGDIIS(double *vectorError, Molecule& molecule, double const* vectorRefStep){
   double** matrixPosition = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            matrixPosition[i][j] = atom->GetXyz()[j];
         }
      }
      if(this->DoGDIIS(vectorError,&matrixPosition[0][0],vectorRefStep)){
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            Atom* atom = molecule.GetAtom(i);
            for(int j=0;j<CartesianType_end;j++){
               atom->GetXyz()[j] = matrixPosition[i][j];
            }
         }
      }
      else{
         MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
         return false;
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
   return true;
}

double GDIIS::MinCosine(){
   static const double inf = std::numeric_limits<double>::infinity();
   // Taken from [FS_2002], p 12.
   static const double mincos[] = {inf, inf, 0.97, 0.84, 0.71, 0.67, 0.62, 0.56, 0.49, 0.41};
   static const int nummincos = sizeof(mincos)/sizeof(mincos[0]);
   int numErrors = listErrors.size();
   return numErrors >= nummincos ? 0.0 : mincos[numErrors];
}

void GDIIS::DiscardPrevious(){
   if(listErrors.size()==0) return;
   double* tmp = listErrors.back();
   MallocerFreer::GetInstance()->Free(&tmp, this->sizeErrorVector);
   listErrors.pop_back();
   tmp = listPositions.back();
   MallocerFreer::GetInstance()->Free(&tmp, this->sizeErrorVector);
   listPositions.pop_back();
}


void GDIIS::DiscardOldest(){
   if(listErrors.size()==0) return;
   double *tmp = listErrors.front();
   MallocerFreer::GetInstance()->Free(&tmp, sizeErrorVector);
   listErrors.pop_front();
   tmp = listPositions.front();
   MallocerFreer::GetInstance()->Free(&tmp, sizeErrorVector);
   listPositions.pop_front();
}
