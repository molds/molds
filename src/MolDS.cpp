#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<time.h>
#include<list>
#include<vector>
#include<omp.h>
#include"base/Utilities.h"
#include"base/MolDSException.h"
#include"base/Enums.h"
#undef INCLUDED_ENUMS
#define RENUMSTR_BODY 1  
#include"base/Enums.h"
#include"base/MathUtilities.h"
#include"base/EularAngle.h"
#include"base/Molecule.h"
#include"base/atoms/Atom.h"
#include"base/atoms/Hatom.h"
#include"base/atoms/Catom.h"
#include"base/atoms/Natom.h"
#include"base/atoms/Oatom.h"
#include"base/atoms/Liatom.h"
#include"base/atoms/Satom.h"
#include"base/MallocerFreer.h"
#include"base/InputParser.h"
#include"base/Parameters.h"
#include"base/GTOExpansionSTO.h"
#include"cndo/Cndo2.h"
#include"indo/Indo.h"
#include"zindo/ZindoS.h"
#include"mkl_wrapper/LapackWrapper.h"



using namespace std;
using namespace MolDS_base;
using namespace MolDS_mkl_wrapper;


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
      LapackWrapper::GetInstance();
      GTOExpansionSTO::GetInstance();


      // Parse input
      try{
         InputParser::GetInstance()->Parse(molecule);
      }
      catch(MolDSException ex){
         cout << ex.what() << endl;
         runingNormally = false;
      }
      
      // CNDO/2
      if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2 && runingNormally){
         MolDS_cndo::Cndo2* cndo2 = new MolDS_cndo::Cndo2();
         try{
            cndo2->SetMolecule(molecule);
            cndo2->DoesSCF();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }
         delete cndo2;
      }

      // INDO
      else if(Parameters::GetInstance()->GetCurrentTheory() == INDO && runingNormally){
         MolDS_indo::Indo* indo = new MolDS_indo::Indo();
         try{
            indo->SetMolecule(molecule);
            indo->DoesSCF();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }
         delete indo;
      }

      // ZINDO/S
      else if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS && runingNormally){
         MolDS_zindo::ZindoS* zindoS = new MolDS_zindo::ZindoS();
         try{
            zindoS->SetMolecule(molecule);
            zindoS->DoesSCF();
            zindoS->DoesCIS();
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }
         delete zindoS;
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
      LapackWrapper::DeleteInstance(); 
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
   catch(MolDSException ex){
      cout << ex.what();
   }

   return 0;
}













