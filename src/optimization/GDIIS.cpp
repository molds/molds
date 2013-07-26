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
#include"../base/Uncopyable.h"
#include"../mpi/MpiProcess.h"
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
   messageOnlyOneErrorVector("There is only one error vector.\n"),
   messageRecalcGDIISStep("Recalculate GDIIS step without the oldest error vector.\n"),
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

void GDIIS::DoGDIIS(double* vectorError,
                    double* vectorPosition,
                    double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException){
   this->Update(vectorError, vectorPosition);
   this->CalcGDIIS(vectorError, vectorPosition, vectorRefStep);
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

void GDIIS::CalcGDIIS(double* vectorError,
                      double* vectorPosition,
                      double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException){
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
      this->OutputLog(messageOnlyOneErrorVector);
      throw GDIISException(this->messageOnlyOneErrorVector);
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
         this->OutputLog(this->messageSingularGDIISMatrix);
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         throw GDIISException(this->messageSingularGDIISMatrix);
      }

      // If Lagrange multiplier is too small;
      if(-vectorCoefs[numErrors] < 1e-8){
         this->OutputLog((this->formatTooSmallLagrangeMultiplier % -vectorCoefs[numErrors]).str());
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         // Recalculate GDIIS step without the oldest data.
         return this->RecalcGDIIS(vectorError,vectorPosition,vectorRefStep);
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
         this->OutputLog((this->formatTooLargeGDIISStep % sqrt(normSquaregdiis) % sqrt(normSquareref)).str());
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         // and recalculate GDIIS step without the oldest data
         return this->RecalcGDIIS(vectorError,vectorPosition,vectorRefStep);
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
         return this->RecalcGDIIS(vectorError,vectorPosition,vectorRefStep);
      }
   }
   catch(GDIISException ex){
      MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
      throw ex;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
   this->OutputLog(messageTakingGDIISStep);
}

void GDIIS::RecalcGDIIS(double* vectorError,
                        double* vectorPosition,
                        double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException){
   double *vectorErrorOldest    = NULL;
   double *vectorPositionOldest = NULL;

   this->PopOldest(&vectorErrorOldest, &vectorPositionOldest);

   this->OutputLog(messageRecalcGDIISStep);

   try{
      this->CalcGDIIS(vectorError, vectorPosition, vectorRefStep);
   }
   catch(GDIISException ex){
      this->PushOldest(vectorErrorOldest, vectorPositionOldest);
      throw ex;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&vectorErrorOldest,   this->sizeErrorVector);
      MallocerFreer::GetInstance()->Free(&vectorPositionOldest, this->sizeErrorVector);
      throw ex;
   }
}

void GDIIS::DoGDIIS(double *vectorError, Molecule& molecule, double const* vectorRefStep) throw(GDIISException, MolDS_base::MolDSException){
   double** matrixPosition = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            matrixPosition[i][j] = atom->GetXyz()[j];
         }
      }
      try{
         this->DoGDIIS(vectorError,&matrixPosition[0][0],vectorRefStep);
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            Atom* atom = molecule.GetAtom(i);
            for(int j=0;j<CartesianType_end;j++){
               atom->GetXyz()[j] = matrixPosition[i][j];
            }
         }
      }
      catch(GDIISException ex){
         MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
         throw ex;
      }
   }
   catch(GDIISException ex){
      MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
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
   if(this->listErrors.size()==0) return;
   double* tmp = listErrors.back();
   MallocerFreer::GetInstance()->Free(&tmp, this->sizeErrorVector);
   this->listErrors.pop_back();
   tmp = listPositions.back();
   MallocerFreer::GetInstance()->Free(&tmp, this->sizeErrorVector);
   this->listPositions.pop_back();
}

void GDIIS::PopOldest(double** vectorError, double** vectorPosition){
   if(this->listErrors.size()==0) return;
   *vectorError    = listErrors.front();
   *vectorPosition = listPositions.front();
   this->listErrors.pop_front();
   this->listPositions.pop_front();
}

void GDIIS::PushOldest(double* vectorError, double* vectorPosition){
   this->listErrors.push_front(vectorError);
   this->listPositions.push_front(vectorPosition);
}

void GDIIS::DiscardOldest(){
   if(this->listErrors.size()==0) return;
   double *vectorErrorOldest    = NULL;
   double *vectorPositionOldest = NULL;
   this->PopOldest(&vectorErrorOldest, &vectorPositionOldest);
   MallocerFreer::GetInstance()->Free(&vectorErrorOldest,    sizeErrorVector);
   MallocerFreer::GetInstance()->Free(&vectorPositionOldest, sizeErrorVector);
}
