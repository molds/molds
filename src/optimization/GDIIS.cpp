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
   numErrors(0),
   nextError(0),
   matrixGDIIS(NULL),
   matrixErrors(NULL),
   matrixPositions(NULL)
{
   MallocerFreer::GetInstance()->Malloc(&this->matrixGDIIS,     maxnumErrors+1, maxnumErrors+1);
   MallocerFreer::GetInstance()->Malloc(&this->matrixErrors,    maxnumErrors,   sizeErrorVector);
   MallocerFreer::GetInstance()->Malloc(&this->matrixPositions, maxnumErrors,   sizeErrorVector);
   double *tmp = this->matrixErrors[0];
}

GDIIS::~GDIIS(){
   MallocerFreer::GetInstance()->Free(&this->matrixGDIIS,     maxnumErrors+1, maxnumErrors+1);
   MallocerFreer::GetInstance()->Free(&this->matrixErrors,    maxnumErrors,   sizeErrorVector);
   MallocerFreer::GetInstance()->Free(&this->matrixPositions, maxnumErrors,   sizeErrorVector);
}

void GDIIS::DoGDIIS(double* vectorError,
                    double* vectorPosition,
                    double const* vectorRefStep){
   // Prepare GDIIS parameters
   const int current = this->nextError++;
   this->nextError %= this->maxnumErrors;
   if(this->numErrors < this->maxnumErrors){
      this->numErrors++;
   }
   for(int i=0;i<this->sizeErrorVector;i++){
      this->matrixErrors[current][i]    = vectorError[i];
      this->matrixPositions[current][i] = vectorPosition[i];
   }

   // Prepare GDIIS matrix
   for(int i=0;i<this->numErrors;i++){
      matrixGDIIS[i][current] = 0.0;
      for(int j=0;j<this->sizeErrorVector;j++){
         matrixGDIIS[i][current] += vectorError[j] * matrixErrors[i][j];
      }
      matrixGDIIS[current][i] = matrixGDIIS[i][current];
   }
   for(int i=0;i<this->numErrors;i++){
      matrixGDIIS[i][numErrors] = matrixGDIIS[numErrors][i] = 1.0;
   }
   matrixGDIIS[numErrors][numErrors] = 0;

   double*  vectorCoefs = NULL;
   try{
      // Solve DIIS equation
      MallocerFreer::GetInstance()->Malloc(&vectorCoefs, this->numErrors+1);
      for(int i=0;i<this->numErrors;i++){
         vectorCoefs[i]=0.0;
      }
      vectorCoefs[numErrors]=1.0;
      MolDS_wrappers::Lapack::GetInstance()->Dsysv(matrixGDIIS,
                                                   vectorCoefs,
                                                   numErrors+1);

      // If only one error vector is given, following routine is meaningless.
      if(numErrors <= 1){
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         return;
      }

      // If Lagrange multiplier is too small, don't take GDIIS step.
      if(-vectorCoefs[numErrors] < 1e-8){
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         return;
      }

      // Interpolate error vectors and positions
      for(int i=0;i<sizeErrorVector;i++){
         vectorError[i] = vectorPosition[i] = 0;
         for(int j=0;j<numErrors;j++){
            vectorError[i]    += vectorCoefs[j] * this->matrixErrors[j][i];
            vectorPosition[i] += vectorCoefs[j] * this->matrixPositions[j][i];
         }
      }

      // Calculate cosine of the angle between GDIIS step and vectorRefStep
      // and lengths of the vectors.
      double innerprod = 0, norm2gdiis = 0, norm2ref = 0;
      for(int i=0;i<this->sizeErrorVector;i++){
         double diff = vectorPosition[i] - matrixPositions[current][i];
         innerprod  += diff * vectorRefStep[i];
         norm2gdiis += diff * diff;
         norm2ref   += vectorRefStep[i] * vectorRefStep[i];
      }
      double cosine = innerprod/sqrt(norm2gdiis*norm2ref);

      // If length of the GDIIS step is larger than reference step * 10
      if(norm2gdiis >= norm2ref * 100){
         for(int i=0; i<this->sizeErrorVector; i++){
            vectorError[i]    = matrixErrors[current][i];
            vectorPosition[i] = matrixPositions[current][i];
         }
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         return;
      }

      // If the calculated cosine value is below the minimum tolerant value
      if(cosine < this->MinCosine()){
         // Rollback vectorPosition and vectorError to original value
         for(int i=0; i<this->sizeErrorVector; i++){
            vectorError[i]    = matrixErrors[current][i];
            vectorPosition[i] = matrixPositions[current][i];
         }
         MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
         return;
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&vectorCoefs, numErrors+1);
}

void GDIIS::DoGDIIS(double *vectorError, Molecule& molecule, double const* vectorRefStep){
   double** matrixPosition = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            matrixPosition[i][j] = atom->GetXyz()[j];
         }
      }
      this->DoGDIIS(vectorError,&matrixPosition[0][0],vectorRefStep);
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            atom->GetXyz()[j] = matrixPosition[i][j];
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixPosition, molecule.GetNumberAtoms(), CartesianType_end);
}

double GDIIS::MinCosine(){
   static const double mincos[] = {1.0/0.0, 1.0/0.0, 0.97, 0.84, 0.71, 0.67, 0.62, 0.56, 0.49, 0.41};
   return this->numErrors >= sizeof(mincos)/sizeof(mincos[0]) ? 0.0 : mincos[this->numErrors];
}
