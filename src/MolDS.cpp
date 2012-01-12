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
#include"mkl.h"
#include"base/MolDSException.h"
#include"base/Uncopyable.h"
#include"mkl_wrapper/LapackWrapper.h"
#include"base/Utilities.h"
#include"base/Enums.h"
#include"base/MallocerFreer.h"
#include"base/EularAngle.h"
#include"base/Parameters.h"
#include"base/atoms/Atom.h"
#include"base/AtomFactory.h"
#include"base/Molecule.h"
#include"base/InputParser.h"
#include"base/GTOExpansionSTO.h"
#include"cndo/Cndo2.h"
#include"indo/Indo.h"
#include"zindo/ZindoS.h"
#include"mndo/Mndo.h"
#include"am1/Am1.h"
#include"pm3/Pm3.h"
#include"pm3/Pm3Pddg.h"
#include"md/MD.h"
#include"mc/MC.h"
using namespace std;
using namespace MolDS_base;

int main(){

   try{
      // Welcome Messages
      OutputWelcomeMessage();
      
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
         cout << ex.what() << endl;
         runingNormally = false;
      }
      
      // electronic structure calculation
      if(runingNormally && (Parameters::GetInstance()->GetCurrentTheory() == CNDO2 ||
                            Parameters::GetInstance()->GetCurrentTheory() == INDO || 
                            Parameters::GetInstance()->GetCurrentTheory() == ZINDOS ||
                            Parameters::GetInstance()->GetCurrentTheory() == MNDO ||
                            Parameters::GetInstance()->GetCurrentTheory() == AM1 ||
                            Parameters::GetInstance()->GetCurrentTheory() == PM3 ||
                            Parameters::GetInstance()->GetCurrentTheory() == PM3PDDG)){
         MolDS_cndo::Cndo2* electronicStructure = NULL;
         if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2 ){
            electronicStructure = new MolDS_cndo::Cndo2();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == INDO ){
            electronicStructure = new MolDS_indo::Indo();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS ){
            electronicStructure = new MolDS_zindo::ZindoS();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == MNDO ){
            electronicStructure = new MolDS_mndo::Mndo();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == AM1 ){
            electronicStructure = new MolDS_am1::Am1();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == PM3 ){
            electronicStructure = new MolDS_pm3::Pm3();
         }
         else if(Parameters::GetInstance()->GetCurrentTheory() == PM3PDDG ){
            electronicStructure = new MolDS_pm3::Pm3Pddg();
         }
         else{
         }
         try{
            electronicStructure->SetMolecule(molecule);
            electronicStructure->DoSCF();
            if(Parameters::GetInstance()->RequiresCIS()){
               electronicStructure->DoCIS();
            }
            if(Parameters::GetInstance()->RequiresMD()){
               MolDS_md::MD* md = new MolDS_md::MD();
               md->SetTheory(electronicStructure);
               md->DoMD();
               delete md;
            }
            if(Parameters::GetInstance()->RequiresMC()){
               MolDS_mc::MC* mc = new MolDS_mc::MC();
               mc->SetTheory(electronicStructure);
               mc->DoMC();
               delete mc;
            }
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }
         delete electronicStructure;
      }

      // Diagonalize Inertia Tensor
      else if(Parameters::GetInstance()->GetCurrentTheory() == PrincipalAxes && runingNormally){
         try{
            molecule->CalcPrincipalAxes();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }

      }

      // Translate molecule
      else if(Parameters::GetInstance()->GetCurrentTheory() == Translate && runingNormally){
         try{
            molecule->Translate();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }

      }

      // Rotate molecule
      else if(Parameters::GetInstance()->GetCurrentTheory() == Rotate && runingNormally){
         try{
            molecule->Rotate();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
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
      OutputFarewellMessage(startTime, startTick, ompStartTime, runingNormally);
   }
   catch(exception ex){
      cout << ex.what();
   }

   return 0;
}

