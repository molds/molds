#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<time.h>
#include<list>
#include<vector>
#include"base/MolDSException.h"
#include"base/Enums.h"
#undef INCLUDED_ENUMS
#define RENUMSTR_BODY 1  
#include"base/Enums.h"
#include"base/Molecule.h"
#include"base/atoms/Atom.h"
#include"base/atoms/Hatom.h"
#include"base/atoms/Catom.h"
#include"base/atoms/Liatom.h"
#include"base/atoms/Satom.h"
#include"base/MallocerFreer.h"
#include"base/InputParser.h"
#include"base/Utilities.h"
#include"base/EularAngle.h"
#include"base/Parameters.h"
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


      // Parse input
      InputParser::GetInstance()->Parse(molecule);

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
         }
         catch(MolDSException ex){
            cout << ex.what() << endl;
            runingNormally = false;
         }
         delete zindoS;
      }


      //Free 
      LapackWrapper::DeleteInstance(); 
      Parameters::DeleteInstance();
      MallocerFreer::DeleteInstance();
      delete molecule;
      InputParser::DeleteInstance();


      // Farewell Messages
      clock_t endTick = clock();
      double consumedTime = (double)(endTick - startTick)/(double)CLOCKS_PER_SEC;
      if(runingNormally){
         cout << "\n\n     >>>>>  The MolDS finished normally!  <<<<<\n";
      }
      else{
         cout << "\n\n     >>>>>  The MolDS finished abnormally..............  <<<<<\n";
      }
      cout <<     "     >>>>>  Consumed time (CPU time): " << consumedTime << "[s].  <<<<<\n";
      cout <<     "     >>>>>  See you.  <<<<<\n\n\n";

   }
   catch(MolDSException ex){
      cout << ex.what();
   }

   return 0;
}













