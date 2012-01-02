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
#include"mkl_wrapper/LapackWrapper.h"
#include"base/Utilities.h"
#include"base/Enums.h"
#include"base/MallocerFreer.h"
#include"base/EularAngle.h"
#include"base/Parameters.h"
#include"base/atoms/Atom.h"
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
using namespace std;
using namespace MolDS_base;

int main(){

   try{
      bool runingNormally = true;
      // Welcome Messages
      time_t startTime;
      struct tm *ltm;
      char s[50];
      clock_t startTick = clock();
      double ompStartTime = omp_get_wtime();
      time(&startTime);
      ltm = localtime(&startTime);
      fmttm(s, ltm);
      cout << "\n\n     >>>>>  Welcome to the MolDS world at " << s << "  <<<<<\n\n\n";

      // declare
      InputParser::GetInstance();
      Molecule* molecule = new Molecule();
      MallocerFreer::GetInstance();
      Parameters::GetInstance();
      MolDS_mkl_wrapper::LapackWrapper::GetInstance();
      GTOExpansionSTO::GetInstance();

      // Parse input
      try{
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
            electronicStructure->DoesSCF();
            if(Parameters::GetInstance()->RequiresCIS()){
               electronicStructure->DoesCIS();
            }
            if(Parameters::GetInstance()->RequiresMD()){
               MolDS_md::MD* md = new MolDS_md::MD();
               md->SetTheory(electronicStructure);
               md->DoesMD();
               delete md;
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
      MallocerFreer::DeleteInstance();
      delete molecule;
      InputParser::DeleteInstance();

      // Farewell Messages
      time_t endTime;
      time(&endTime);
      clock_t endTick = clock();
      double consumedTime = (double)(endTick - startTick)/(double)CLOCKS_PER_SEC;
      double ompEndTime = omp_get_wtime();
      if(runingNormally){
         cout << "\n\n     >>>>>  The MolDS finished normally!  <<<<<\n";
      }
      else{
         cout << "\n\n     >>>>>  The MolDS finished abnormally..............  <<<<<\n";
      }
      cout << "     >>>>>  CPU time: " << consumedTime << "[s].  <<<<<\n";
      cout << "     >>>>>  Elapsed time: " << endTime - startTime << "[s].  <<<<<\n";
      cout << "     >>>>>  Elapsed time(OMP): " << ompEndTime - ompStartTime << "[s].  <<<<<\n";
      cout << "     >>>>>  See you.  <<<<<\n\n\n";

   }
   catch(exception ex){
      cout << ex.what();
   }

   return 0;
}

