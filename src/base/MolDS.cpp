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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<time.h>
#include<list>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/shared_ptr.hpp>
#include<boost/random.hpp>
#include"mkl.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"../mkl_wrapper/LapackWrapper.h"
#include"Utilities.h"
#include"Enums.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"AtomFactory.h"
#include"Molecule.h"
#include"InputParser.h"
#include"GTOExpansionSTO.h"
#include"ElectronicStructure.h"
#include"ElectronicStructureFactory.h"
#include"../md/MD.h"
#include"../mc/MC.h"
#include"../rpmd/RPMD.h"
#include"../optimize/SteepestDescent.h"
#include"MolDS.h"
using namespace std;
namespace MolDS_base{
void MolDS::Run(){
   // Welcome Messages
   this->OutputLog(Utilities::GetWelcomeMessage());
   
   //timer set
   time_t startTime;
   time(&startTime);
   clock_t startTick = clock();
   double ompStartTime = omp_get_wtime();

   Molecule* molecule = NULL;
   bool runingNormally = true;
   try{
      // declare 
      MallocerFreer::GetInstance();
      AtomFactory::GetInstance();
      InputParser::GetInstance();
      molecule = new Molecule();
      Parameters::GetInstance();
      MolDS_mkl_wrapper::LapackWrapper::GetInstance();
      GTOExpansionSTO::GetInstance();
      // Parse input
      InputParser::GetInstance()->Parse(molecule);
   }
   catch(MolDSException ex){
      this->OutputLog((boost::format("%s\n") % ex.what()).str());
      runingNormally = false;
   }
   
   // once electronic structure calculation
   if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == Once){
      ElectronicStructure* electronicStructure = NULL;
      try{
         electronicStructure = ElectronicStructureFactory::GetInstance()->Create();
         electronicStructure->SetMolecule(molecule);
         electronicStructure->DoSCF();
         if(Parameters::GetInstance()->RequiresCIS()){
            electronicStructure->DoCIS();
         }
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }
      if(electronicStructure != NULL){
         delete electronicStructure;
      }
   }

   // MD
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == MD){
      MolDS_md::MD* md = NULL;
      try{
         md = new MolDS_md::MD();
         md->SetMolecule(molecule);
         md->DoMD();
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }
      if(md != NULL){
         delete md;
      }
   }

   // MC
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == MC){
      MolDS_mc::MC* mc = NULL;
      try{
         mc = new MolDS_mc::MC();
         mc->SetMolecule(molecule);
         mc->DoMC();
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }
      if(mc != NULL){
         delete mc;
      }
   }

   // RPMD
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == RPMD){
      MolDS_rpmd::RPMD* rpmd = NULL;
      try{
         rpmd = new MolDS_rpmd::RPMD();
         rpmd->DoRPMD(*molecule);
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }
      if(rpmd != NULL){
         delete rpmd;
      }
   }

   // Optimize (Steepest Descent)
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == Optimize){
      MolDS_optimize::SteepestDescent* steepestDescent = NULL;
      try{
         steepestDescent = new MolDS_optimize::SteepestDescent();
         steepestDescent->Optimize(*molecule);
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }
      if(steepestDescent != NULL){
         delete steepestDescent;
      }
   }

   // Diagonalize Inertia Tensor
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == PrincipalAxes ){
      try{
         molecule->CalcPrincipalAxes();
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }

   }

   // Translate molecule
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == Translate){
      try{
         molecule->Translate();
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }

   }

   // Rotate molecule
   else if(runingNormally && Parameters::GetInstance()->GetCurrentSimulation() == Rotate){
      try{
         molecule->Rotate();
      }
      catch(MolDSException ex){
         this->OutputLog((boost::format("%s\n") % ex.what()).str());
         runingNormally = false;
      }

   }

   //Free 
   GTOExpansionSTO::DeleteInstance();
   MolDS_mkl_wrapper::LapackWrapper::DeleteInstance(); 
   Parameters::DeleteInstance();
   delete molecule;
   InputParser::DeleteInstance();
   AtomFactory::DeleteInstance();
   MallocerFreer::DeleteInstance();

   // Farewell Messages
   this->OutputLog(Utilities::GetFarewellMessage(startTime, startTick, ompStartTime, runingNormally));
}
}

