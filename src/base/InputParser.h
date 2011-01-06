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
public:
   static InputParser* GetInstance();
   static void DeleteInstance();
   void Parse(Molecule* molecule);
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
   string messageInputTerms;
   string messageScfConditions;
   string messageScfMaxIterations;
   string messageScfRmsDensity;
   string stringSpace;
   string stringCommentOut;
   string stringTheory;
   string stringTheoryEnd;
   string stringTheoryCNDO2;
   string stringTheoryINDO;
   string stringTheoryZINDOS;
   string stringTheoryPrincipalAxes;
   string stringTheoryTranslate;
   string stringTheoryRotate;
   string stringTheoryNONE;
   string stringGeometry;
   string stringGeometryEnd;
   string stringScf;
   string stringScfEnd;
   string stringScfMaxIter;
   string stringScfRmsDensity;
   string stringInertiaTensor;
   string stringInertiaTensorEnd;
   string stringInertiaTensorOrigin;
   string stringRotate;
   string stringRotateEnd;
   string stringRotatingOrigin;
   string stringRotatingAxis;
   string stringRotatingAngle;
   string stringRotatingAngles;
   string stringRotatingType;
   string stringRotatingTypeAxis;
   string stringRotatingTypeEularAngle;
   string stringTranslate;
   string stringTranslateEnd;
   string stringTranslatingDifference;
   void CalcMolecularBasics(Molecule* molecule);
   void OutputMolecularBasics(Molecule* molecule);
   void OutputScfConditions();
   void OutputInputTerms(vector<string> inputTerms);
   bool IsCommentOut(string str);
   vector<string> GetInputTerms();
};
InputParser* InputParser::inputParser = NULL;

InputParser::InputParser(){
   this->messageStartParseInput = "**********  START: Parse input  **********\n";
   this->messageDoneParseInput =  "**********  DONE: Parse input  ***********\n\n\n";
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageInputTerms = "Input terms:\n";
   this->messageScfConditions = "\tSCF conditions:\n";
   this->messageScfMaxIterations = "\t\tMax iterations: ";
   this->messageScfRmsDensity = "\t\tRMS density: ";
   this->stringSpace = " ";
   this->stringCommentOut = "//";
   this->stringTheoryCNDO2 = "cndo/2";
   this->stringTheoryINDO = "indo";
   this->stringTheoryZINDOS = "zindo/s";
   this->stringTheoryPrincipalAxes = "principal_axes";
   this->stringTheoryTranslate = "translate";
   this->stringTheoryRotate = "rotate";
   this->stringTheoryNONE = "none";
   this->stringGeometry =    "geometry";
   this->stringGeometryEnd = "geometry_end";
   this->stringTheory = "theory";
   this->stringTheoryEnd = "theory_end";
   this->stringScf = "scf";
   this->stringScfEnd = "scf_end";
   this->stringScfMaxIter = "max_iter";
   this->stringScfRmsDensity = "rms_density";
   this->stringInertiaTensor = "inertia";
   this->stringInertiaTensorEnd = "inertia_end";
   this->stringInertiaTensorOrigin = "origin";
   this->stringRotate = "rotate";
   this->stringRotateEnd = "rotate_end";
   this->stringRotatingOrigin = "origin";
   this->stringRotatingAxis = "axis";
   this->stringRotatingAngle = "angle";
   this->stringRotatingAngles = "angles";
   this->stringRotatingType = "type";
   this->stringRotatingTypeAxis = "axis";
   this->stringRotatingTypeEularAngle = "eular_angle";
   this->stringTranslate = "translate";
   this->stringTranslateEnd = "translate_end";
   this->stringTranslatingDifference = "difference";
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
            double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
            double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
            double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
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
      
      // inertia tensor condition
      if(inputTerms[i].compare(this->stringInertiaTensor) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringInertiaTensorEnd) != 0){
            // origin
            if(inputTerms[j].compare(this->stringInertiaTensorOrigin) == 0){
               double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               molecule->SetInertiaTensorOrigin(x, y, z);
               j+=3;
            }
            j++;   
         }
         i = j;
      }
      
      // translating condition
      if(inputTerms[i].compare(this->stringTranslate) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringTranslateEnd) != 0){
            // origin
            if(inputTerms[j].compare(this->stringTranslatingDifference) == 0){
               double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               molecule->SetTranslatingDifference(x, y, z);
               j+=3;
            }
            j++;   
         }
         i = j;
      }
      
      // rotating condition
      if(inputTerms[i].compare(this->stringRotate) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringRotateEnd) != 0){
            // origin
            if(inputTerms[j].compare(this->stringRotatingOrigin) == 0){
               double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               molecule->SetRotatingOrigin(x, y, z);
               j+=3;
            }
            // axis
            else if(inputTerms[j].compare(this->stringRotatingAxis) == 0){
               double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               molecule->SetRotatingAxis(x, y, z);
               j+=3;
            }
            // angle 
            else if(inputTerms[j].compare(this->stringRotatingAngle) == 0){
               double angle = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               molecule->SetRotatingAngle(angle);
               j++;
            }
            // angles (EularAngle)
            else if(inputTerms[j].compare(this->stringRotatingAngles) == 0){
               double alpha = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               double beta  = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               double gamma = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               molecule->SetRotatingEularAngles(alpha, beta, gamma);
               j += 3;
            }
            // type
            else if(inputTerms[j].compare(this->stringRotatingType) == 0){
               if(inputTerms[j+1].compare(this->stringRotatingTypeAxis) == 0){
                  molecule->SetRotatingType(Axis);
               }
               else if(inputTerms[j+1].compare(this->stringRotatingTypeEularAngle) == 0){
                  molecule->SetRotatingType(Eular);
               }
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

            // Princepal axes
            else if(inputTerms[j].compare(this->stringTheoryPrincipalAxes) == 0){
               Parameters::GetInstance()->SetCurrentTheory(PrincipalAxes);
            }

            // Translate
            else if(inputTerms[j].compare(this->stringTheoryTranslate) == 0){
               Parameters::GetInstance()->SetCurrentTheory(Translate);
            }

            // Rotate
            else if(inputTerms[j].compare(this->stringTheoryRotate) == 0){
               Parameters::GetInstance()->SetCurrentTheory(Rotate);
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

   molecule->OutputTotalNumberAtomsAOsValenceelectrons();
   molecule->OutputConfiguration();
   molecule->OutputCOMXyz();
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





