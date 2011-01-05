#ifndef INCLUDED_INPUT_PARSER
#define INCLUDED_INPUT_PARSER

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>
#include"Molecule.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"atoms/Hatom.h"
#include"atoms/Liatom.h"
#include"atoms/Catom.h"

using namespace std;
using namespace MolDS_base_atoms;

namespace MolDS_base{

// InputParser is singleton
class InputParser{
private:
   static InputParser* inputParser;
   InputParser();
   InputParser(InputParser&);
   void operator = (InputParser&);
   ~InputParser();
   string messageStartParseInput;
   string messageDoneParseInput;
   string messageTotalNumberAOs;
   string messageTotalNumberAtoms;
   string messageTotalNumberValenceElectrons;
   string messageConfigration;
   string messageConfigrationTitleAU;
   string messageConfigrationTitleAng;
   string messageCOM;
   string messageCOMTitleAU;
   string messageCOMTitleAng;
   string messageInputTerms;
   string messageScfConditions;
   string messageScfMaxIterations;
   string messageScfRmsDensity;
   string stringSpace;
   string stringCommentOut;
   string stringGeometry;
   string stringGeometryEnd;
   string stringTheory;
   string stringTheoryEnd;
   string stringTheoryCNDO2;
   string stringTheoryINDO;
   string stringTheoryZINDOS;
   string stringTheoryNONE;
   string stringScf;
   string stringScfEnd;
   string stringScfMaxIter;
   string stringScfRmsDensity;
   void CalcMolecularBasics(Molecule* molecule);
   void OutputMolecularBasics(Molecule* molecule);
   void OutputScfConditions();
   void OutputInputTerms(vector<string> inputTerms);
   bool IsCommentOut(string str);
   vector<string> GetInputTerms();

public:
   static InputParser* GetInstance();
   static void DeleteInstance();
   void Parse(Molecule* molecule);

};
InputParser* InputParser::inputParser = NULL;

InputParser::InputParser(){
   this->messageStartParseInput = "**********  START: Parse input  **********\n";
   this->messageDoneParseInput =  "**********  DONE: Parse input  ***********\n\n\n";
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageConfigration = "\tMolecular configration:\n";
   this->messageConfigrationTitleAU = "\t\t| i-th | atom type | x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageConfigrationTitleAng = "\t\t| i-th | atom type | x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageCOM = "\tCenter of Mass:\n";
   this->messageCOMTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageCOMTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageInputTerms = "Input terms:\n";
   this->messageScfConditions = "\tSCF conditions:\n";
   this->messageScfMaxIterations = "\t\tMax iterations: ";
   this->messageScfRmsDensity = "\t\tRMS density: ";
   this->stringSpace = " ";
   this->stringCommentOut = "//";
   this->stringGeometry =    "geometry";
   this->stringGeometryEnd = "geometry_end";
   this->stringTheory = "theory";
   this->stringTheoryEnd = "theory_end";
   this->stringScf = "scf";
   this->stringScfEnd = "scf_end";
   this->stringScfMaxIter = "max_iter";
   this->stringScfRmsDensity = "rms_density";
   this->stringTheoryCNDO2 = "cndo/2";
   this->stringTheoryINDO = "indo";
   this->stringTheoryZINDOS = "zindo/s";
   this->stringTheoryNONE = "none";
}

InputParser::~InputParser(){
}

InputParser* InputParser::GetInstance(){
   if(inputParser == NULL){
      inputParser = new InputParser();
   }
   return inputParser;
}

void InputParser::DeleteInstance(){
   if(inputParser != NULL){
      delete inputParser; 
   }
   inputParser = NULL;
}

vector<string> InputParser::GetInputTerms(){

   string str;
   string inputTerm;
   vector<string> inputTerms;
   bool isPreCharSpace = false;
   while(getline(cin, str)){

      // check comment out 
      if(this->IsCommentOut(str)){
         continue;
      }

      // get input terms
      for(int i=0; i<str.length(); i++){
         if(str.data()[i] != stringSpace.data()[0]){
            // change to lower case.
            inputTerm += tolower(str.data()[i]);
            isPreCharSpace = false;
         }
         else{
            if(!isPreCharSpace){
               inputTerms.push_back(inputTerm);
               inputTerm = "";
               isPreCharSpace = true;
            }
         }
      }
      if(inputTerm.length()>0){
         inputTerms.push_back(inputTerm);
         inputTerm = "";
      }
      isPreCharSpace = true;
   }

   return inputTerms;

}

void InputParser::Parse(Molecule* molecule){

   cout << messageStartParseInput;

   // read input
   vector<string> inputTerms = this->GetInputTerms();

   // parse input
   for(int i=0; i<inputTerms.size();i++){

      // molecular geometry
      if(inputTerms[i].compare(this->stringGeometry) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringGeometryEnd) != 0){
            double x = atof(inputTerms[j+1].c_str());
            double y = atof(inputTerms[j+2].c_str());
            double z = atof(inputTerms[j+3].c_str());
            if(inputTerms[j] == "h"){
               molecule->GetAtomVect()->push_back(new Hatom(x, y, z));
            }
            else if(inputTerms[j] == "li"){
               molecule->GetAtomVect()->push_back(new Liatom(x, y, z));
            }
            else if(inputTerms[j] == "c"){
               molecule->GetAtomVect()->push_back(new Catom(x, y, z));
            }
            else if(inputTerms[j] == "s"){
               molecule->GetAtomVect()->push_back(new Satom(x, y, z));
            }
            j += 4;
         }
         i = j;
      }

      // scf condition
      if(inputTerms[i].compare(this->stringScf) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringScfEnd) != 0){
            // max iterations
            if(inputTerms[j].compare(this->stringScfMaxIter) == 0){
               Parameters::GetInstance()->SetMaxIterationsSCF(atoi(inputTerms[j+1].c_str()));
               j++;
            }
            // RMS density 
            if(inputTerms[j].compare(this->stringScfRmsDensity) == 0){
               Parameters::GetInstance()->SetThresholdSCF(atof(inputTerms[j+1].c_str()));
               j++;
            }
            j++;   
         }
         i = j;
      }
      
      // theory
      if(inputTerms[i].compare(this->stringTheory) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringTheoryEnd) != 0){

            // CNDO/2
            if(inputTerms[j].compare(this->stringTheoryCNDO2) == 0){
               Parameters::GetInstance()->SetCurrentTheory(CNDO2);
            }

            // INDO
            else if(inputTerms[j].compare(this->stringTheoryINDO) == 0){
               Parameters::GetInstance()->SetCurrentTheory(INDO);
            }

            // ZINDO/S
            else if(inputTerms[j].compare(this->stringTheoryZINDOS) == 0){
               Parameters::GetInstance()->SetCurrentTheory(ZINDOS);
            }

            // NONE 
            else if(inputTerms[j].compare(this->stringTheoryNONE) == 0){
               Parameters::GetInstance()->SetCurrentTheory(NONE);
            }
            j++;
         }
         i = j;
      }

   }

   this->CalcMolecularBasics(molecule);
   this->OutputMolecularBasics(molecule);
   this->OutputScfConditions();
   this->OutputInputTerms(inputTerms);
   cout << messageDoneParseInput;

}

void InputParser::CalcMolecularBasics(Molecule* molecule){

   molecule->CalcTotalNumberAOs();
   molecule->CalcTotalNumberValenceElectrons();
   molecule->CalcCOMXyz();

}

void InputParser::OutputMolecularBasics(Molecule* molecule){

   // Total number of Atoms
   cout << this->messageTotalNumberAtoms << molecule->GetAtomVect()->size() << "\n";
   cout << this->messageTotalNumberAOs << molecule->GetTotalNumberAOs() << "\n";
   cout << this->messageTotalNumberValenceElectrons << molecule->GetTotalNumberValenceElectrons() << "\n\n";

   // configuration
   cout << this->messageConfigration;
   cout << this->messageConfigrationTitleAng;
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   for(int a=0; a<molecule->GetAtomVect()->size(); a++){
      Atom* atom = (*molecule->GetAtomVect())[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0]/ang2AU, atom->GetXyz()[1]/ang2AU, atom->GetXyz()[2]/ang2AU);
   }
   cout << "\n";

   cout << this->messageConfigrationTitleAU;
   for(int a=0; a<molecule->GetAtomVect()->size(); a++){
      Atom* atom = (*molecule->GetAtomVect())[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0], atom->GetXyz()[1], atom->GetXyz()[2]);
   }
   cout << "\n";

   // Center of Mass
   cout << this->messageCOM;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",molecule->GetCOMXyz()[0]/ang2AU, 
                             molecule->GetCOMXyz()[1]/ang2AU, 
                             molecule->GetCOMXyz()[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",molecule->GetCOMXyz()[0], 
                             molecule->GetCOMXyz()[1], 
                             molecule->GetCOMXyz()[2]);
   cout << "\n";

}

void InputParser::OutputScfConditions(){

   cout << this->messageScfConditions;
   printf("%s%d\n",this->messageScfMaxIterations.c_str(),Parameters::GetInstance()->GetMaxIterationsSCF());
   printf("%s%e\n",this->messageScfRmsDensity.c_str(),Parameters::GetInstance()->GetThresholdSCF());
   cout << "\n";


}
void InputParser::OutputInputTerms(vector<string> inputTerms){
   
   // output input terms
   cout << this->messageInputTerms;
   for(int i=0; i<inputTerms.size();i++){
      cout << inputTerms[i] << " | ";
      if(i%10 == 9){
         cout << "\n";
      }
   }
   cout << endl << endl;

}

/****
 *
 *  # or // are treated as comment out 
 *
 ****/
bool InputParser::IsCommentOut(string tempStr){

   string str = TrimString(tempStr);

   string commentPrefix1 = "#";
   string prefix1;
   if(str.length()>=1){
      prefix1 += str.data()[0];
   }

   string commentPrefix2 = "//";
   string prefix2;
   if(str.length()>=2){
      prefix2 += str.data()[0];
      prefix2 += str.data()[1];
   }

   return 0==prefix1.compare(commentPrefix1) || 0==prefix2.compare(commentPrefix2) ;

}

}
#endif





