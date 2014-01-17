//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/RealSphericalHarmonicsIndex.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"Optimizer.h"
#include"ConjugateGradient.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_optimization{

ConjugateGradient::ConjugateGradientState::ConjugateGradientState(Molecule& molecule,
                                                                  boost::shared_ptr<ElectronicStructure>& electronicStructure):
   OptimizerState(molecule, electronicStructure),
   oldMatrixForce(NULL),
   matrixSearchDirection(NULL),
   numAtoms(molecule.GetAtomVect().size()){
   MallocerFreer::GetInstance()->Malloc<double>(&this->oldMatrixForce       , this->numAtoms, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(&this->matrixSearchDirection, this->numAtoms, CartesianType_end);
}

ConjugateGradient::ConjugateGradientState::~ConjugateGradientState(){
   MallocerFreer::GetInstance()->Free<double>(&this->oldMatrixForce       , this->numAtoms, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->matrixSearchDirection, this->numAtoms, CartesianType_end);
}

ConjugateGradient::ConjugateGradient(){
   this->SetMessages();
   //this->OutputLog("ConjugateGradient created\n");
}

ConjugateGradient::~ConjugateGradient(){
   //this->OutputLog("ConjugateGradient deleted\n");
}

void ConjugateGradient::SetMessages(){
   Optimizer::SetMessages();
   this->errorMessageNotEnebleTheoryType  
      = "Error in optimization::ConjugateGradient::CheckEnableTheoryType: Non available theory is set.\n";
   this->errorMessageGeometyrOptimizationNotConverged 
      = "Error in optimization::ConjugateGradient::Optimize: Optimization did not met convergence criterion.\n";
   this->messageStartConjugateGradientStep = "\n==========  START: Conjugate gradient step ";
}

void ConjugateGradient::SearchMinimum(boost::shared_ptr<ElectronicStructure> electronicStructure,
                                      Molecule& molecule,
                                      double* lineSearchedEnergy,
                                      bool* obtainesOptimizedStructure) const{
   ConjugateGradientState state(molecule,electronicStructure);

   // initial calculation
   bool requireGuess = true;
   this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, this->CanOutputLogs());
   state.SetCurrentEnergy(electronicStructure->GetElectronicEnergy(state.GetElecState()));

   requireGuess = false;
   state.SetMatrixForce(electronicStructure->GetForce(state.GetElecState()));

   this->InitializeState(state, molecule);

   // conugate gradient loop
   for(int s=0; s<state.GetTotalSteps(); s++){
      this->OutputOptimizationStepMessage(s);

      this->PrepareState(state, molecule, electronicStructure, state.GetElecState());

      this->CalcNextStepGeometry(molecule, state, electronicStructure, state.GetElecState(), state.GetDeltaT());

      this->UpdateSearchDirection(state, electronicStructure, molecule, state.GetElecState());

      // check convergence
      if(this->SatisfiesConvergenceCriterion(state.GetMatrixForce(),
                                             molecule,
                                             state.GetInitialEnergy(),
                                             state.GetCurrentEnergy(),
                                             state.GetMaxGradientThreshold(),
                                             state.GetRmsGradientThreshold())){
         *obtainesOptimizedStructure = true;
         break;
      }
   }

   *lineSearchedEnergy = state.GetCurrentEnergy();
}

void ConjugateGradient::InitializeState(OptimizerState &stateOrig, const Molecule& molecule) const{
   ConjugateGradientState& state = stateOrig.CastRef<ConjugateGradientState>();
   for(int a=0;a<molecule.GetAtomVect().size();a++){
      for(int i=0; i<CartesianType_end; i++){
         state.GetMatrixSearchDirection()[a][i] = state.GetMatrixForce()[a][i];
      }
   }
}

void ConjugateGradient::CalcNextStepGeometry(Molecule &molecule,
                                             OptimizerState& stateOrig,
                                             boost::shared_ptr<ElectronicStructure> electronicStructure,
                                             const int elecState,
                                             const double dt) const{
   ConjugateGradientState& state = stateOrig.CastRef<ConjugateGradientState>();

   state.SetInitialEnergy(state.GetCurrentEnergy());

   this->LineSearch(electronicStructure, molecule, state.GetCurrentEnergyRef(), state.GetMatrixSearchDirection(), elecState, dt);
}

void ConjugateGradient::UpdateSearchDirection(OptimizerState& stateOrig,
                                              boost::shared_ptr<ElectronicStructure> electronicStructure, 
                                              const MolDS_base::Molecule& molecule,
                                              int elecState) const{
   ConjugateGradientState& state = stateOrig.CastRef<ConjugateGradientState>();
   for(int a=0;a<molecule.GetAtomVect().size();a++){
      for(int i=0; i<CartesianType_end; i++){
         state.GetOldMatrixForce()[a][i] = state.GetMatrixForce()[a][i];
      }
   }
   state.SetMatrixForce(electronicStructure->GetForce(elecState));
   double beta=0.0;
   double temp=0.0;
   for(int a=0;a<molecule.GetAtomVect().size();a++){
      for(int i=0; i<CartesianType_end; i++){
         temp += pow(state.GetOldMatrixForce()[a][i],2.0);
         beta += (state.GetMatrixForce()[a][i] - state.GetOldMatrixForce()[a][i])*state.GetMatrixForce()[a][i];
      }
   }
   beta /= temp;
   for(int a=0;a<molecule.GetAtomVect().size();a++){
      for(int i=0; i<CartesianType_end; i++){
         state.GetMatrixSearchDirection()[a][i] *= beta;
         state.GetMatrixSearchDirection()[a][i] += state.GetMatrixForce()[a][i];
      }
   }
}


}



