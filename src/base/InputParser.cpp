#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include"Utilities.h"
#include"Enums.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"atoms/Hatom.h"
#include"atoms/Liatom.h"
#include"atoms/Catom.h"
#include"atoms/Natom.h"
#include"atoms/Oatom.h"
#include"atoms/Satom.h"
#include"Molecule.h"
#include"InputParser.h"
using namespace std;
using namespace MolDS_base_atoms;
namespace MolDS_base{

InputParser* InputParser::inputParser = NULL;

InputParser::InputParser(){
   this->SetMessages();
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

void InputParser::SetMessages(){
   this->messageStartParseInput = "**********  START: Parse input  **********\n";
   this->messageDoneParseInput =  "**********  DONE: Parse input  ***********\n\n\n";
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageInputTerms = "Input terms:\n";
   this->messageScfConditions = "\tSCF conditions:\n";
   this->messageScfMaxIterations = "\t\tMax iterations: ";
   this->messageScfRmsDensity = "\t\tRMS density: ";
   this->messageScfDampingThresh = "\t\tDamping threshold: ";
   this->messageScfDampingWeight = "\t\tDamping weight: ";
   this->messageScfDiisNumErrorVect = "\t\tDIIS number of error vectors: ";
   this->messageScfDiisStartError = "\t\tDIIS starting error: ";
   this->messageScfDiisEndError = "\t\tDIIS ending error: ";
   this->messageCisConditions = "\tCIS conditions:\n";
   this->messageCisNumberActiveOcc = "\t\tNumber of active Occ.: ";
   this->messageCisNumberActiveVir = "\t\tNumber of active Vir.: ";
   this->messageCisNumberExcitedStates = "\t\tNumber of excited states: ";
   this->messageCisDavidson = "\t\tCIS-Davidson: ";
   this->messageCisNormTolerance = "\t\tNorm tolerance for the residual of the Davidson: ";
   this->messageCisMaxIterations = "\t\tMax iterations for the Davidson: ";
   this->messageCisMaxDimensions = "\t\tMax dimensions for the Davidson: ";
   this->messageMdConditions = "\tMD conditions:\n";
   this->messageMdTotalSteps = "\t\tTotal steps: ";
   this->messageMdElecState = "\t\tElectronic eigen state: ";
   this->messageMdTimeWidth = "\t\tTime width(dt): ";
   this->messageMOPlotConditions = "\tMO plot conditions:\n";
   this->messageMOPlotIndex = "\t\tMO index: ";
   this->messageMOPlotGridNumber = "\t\tNumber of grid(x, y, z): ";
   this->messageMOPlotFrameLength = "\t\tFrame length[angst.](x, y, z): ";
   this->messageMOPlotFilePrefix = "\t\tFile name prefix: ";
   this->messageFs = "[fs]";
   this->stringYES = "yes";
   this->stringNO = "no";
   this->stringSpace = " ";
   this->stringCommentOut = "//";
   this->stringTheoryCNDO2 = "cndo/2";
   this->stringTheoryINDO = "indo";
   this->stringTheoryZINDOS = "zindo/s";
   this->stringTheoryMNDO = "mndo";
   this->stringTheoryAM1 = "am1";
   this->stringTheoryPM3 = "pm3";
   this->stringTheoryPM3PDDG = "pm3/pddg";
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
   this->stringScfDampingThresh = "damping_thresh";
   this->stringScfDampingWeight = "damping_weight";
   this->stringScfDiisNumErrorVect = "diis_num_error_vect";
   this->stringScfDiisStartError = "diis_start_error";
   this->stringScfDiisEndError = "diis_end_error";
   this->stringMO = "mo";
   this->stringMOPlot = "moplot";
   this->stringMOPlotEnd = "moplot_end";
   this->stringMOPlotGridNumber = "grid_number";
   this->stringMOPlotFrameLength = "frame_length";
   this->stringMOPlotFilePrefix = "file_prefix";
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
   this->stringCIS = "cis";
   this->stringCISEnd = "cis_end";
   this->stringCISActiveOcc = "active_occ";
   this->stringCISActiveVir = "active_vir";
   this->stringCISNStates = "nstates";
   this->stringCISDavidson = "davidson";
   this->stringCISMaxIter = "max_iter";
   this->stringCISMaxDimensions = "max_dim";
   this->stringCISNormTolerance = "norm_tol";
   this->stringMD = "md";
   this->stringMDEnd = "md_end";
   this->stringMDTotalSteps = "total_steps";
   this->stringMDElecState = "electronic_state";
   this->stringMDTimeWidth = "dt";
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
            else if(inputTerms[j] == "n"){
               molecule->GetAtomVect()->push_back(new Natom(x, y, z));
            }
            else if(inputTerms[j] == "o"){
               molecule->GetAtomVect()->push_back(new Oatom(x, y, z));
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
            // Damping Threshold 
            if(inputTerms[j].compare(this->stringScfDampingThresh) == 0){
               Parameters::GetInstance()->SetDampingThreshSCF(atof(inputTerms[j+1].c_str()));
               j++;
            }
            // Damping Weight
            if(inputTerms[j].compare(this->stringScfDampingWeight) == 0){
               Parameters::GetInstance()->SetDampingWeightSCF(atof(inputTerms[j+1].c_str()));
               j++;
            }
            // DIIS number of stored error vectors
            if(inputTerms[j].compare(this->stringScfDiisNumErrorVect) == 0){
               Parameters::GetInstance()->SetDiisNumErrorVectSCF(atoi(inputTerms[j+1].c_str()));
               j++;
            }
            // DIIS starting error
            if(inputTerms[j].compare(this->stringScfDiisStartError) == 0){
               Parameters::GetInstance()->SetDiisStartErrorSCF(atof(inputTerms[j+1].c_str()));
               j++;
            }
            // DIIS ending error
            if(inputTerms[j].compare(this->stringScfDiisEndError) == 0){
               Parameters::GetInstance()->SetDiisEndErrorSCF(atof(inputTerms[j+1].c_str()));
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
               Parameters::GetInstance()->SetInertiaTensorOrigin(x, y, z);
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
               Parameters::GetInstance()->SetTranslatingDifference(x, y, z);
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
               Parameters::GetInstance()->SetRotatingOrigin(x, y, z);
               j+=3;
            }
            // axis
            else if(inputTerms[j].compare(this->stringRotatingAxis) == 0){
               double x = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double y = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               double z = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
               Parameters::GetInstance()->SetRotatingAxis(x, y, z);
               j+=3;
            }
            // angle 
            else if(inputTerms[j].compare(this->stringRotatingAngle) == 0){
               double angle = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               Parameters::GetInstance()->SetRotatingAngle(angle);
               j++;
            }
            // angles (EularAngle)
            else if(inputTerms[j].compare(this->stringRotatingAngles) == 0){
               double alpha = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               double beta  = atof(inputTerms[j+2].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               double gamma = atof(inputTerms[j+3].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
               Parameters::GetInstance()->SetRotatingEularAngles(alpha, beta, gamma);
               j += 3;
            }
            // type
            else if(inputTerms[j].compare(this->stringRotatingType) == 0){
               if(inputTerms[j+1].compare(this->stringRotatingTypeAxis) == 0){
                  Parameters::GetInstance()->SetRotatingType(Axis);
               }
               else if(inputTerms[j+1].compare(this->stringRotatingTypeEularAngle) == 0){
                  Parameters::GetInstance()->SetRotatingType(Eular);
               }
               j++;
            }
            j++;   
         }
         i = j;
      }
      
   this->stringMOPlotFilePrefix = "file_prefix";
      // mo plot condition
      if(inputTerms[i].compare(this->stringMOPlot) == 0){
         int j=i+1;
         while(inputTerms[j].compare(this->stringMOPlotEnd) != 0){
            // Frame length
            if(inputTerms[j].compare(this->stringMOPlotFrameLength) == 0){
               double lx = atof(inputTerms[j+1].c_str()) 
                          *Parameters::GetInstance()->GetAngstrom2AU();
               double ly = atof(inputTerms[j+2].c_str()) 
                          *Parameters::GetInstance()->GetAngstrom2AU();
               double lz = atof(inputTerms[j+3].c_str()) 
                          *Parameters::GetInstance()->GetAngstrom2AU();
               Parameters::GetInstance()->SetFrameLengthMOPlot(lx, ly, lz);
               j += 3;
            }
            // Grid number
            if(inputTerms[j].compare(this->stringMOPlotGridNumber) == 0){
               int nx = atof(inputTerms[j+1].c_str());
               int ny = atof(inputTerms[j+2].c_str());
               int nz = atof(inputTerms[j+3].c_str());
               Parameters::GetInstance()->SetGridNumberMOPlot(nx, ny, nz);
               j += 3;
            }
            // mo index
            if(inputTerms[j].compare(this->stringMO) == 0){
               int moIndex = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->AddIndexMOPlot(moIndex);
               j++;
            }
            // file prefix
            if(inputTerms[j].compare(this->stringMOPlotFilePrefix) == 0){
               string filePrefix(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetFileNamePrefixMOPlot(filePrefix);
               j++;
            }
            j++;   
         }
         i = j;
      }
      // cis condition
      if(inputTerms[i].compare(this->stringCIS) == 0){
         Parameters::GetInstance()->SetRequiresCIS(true);
         int j=i+1;
         while(inputTerms[j].compare(this->stringCISEnd) != 0){
            // number of active occupied orbitals
            if(inputTerms[j].compare(this->stringCISActiveOcc) == 0){
               int activeOccCIS = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetActiveOccCIS(activeOccCIS);
               j++;
            }
            // number of active virtual orbitals
            if(inputTerms[j].compare(this->stringCISActiveVir) == 0){
               int activeVirCIS = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetActiveVirCIS(activeVirCIS);
               j++;
            }
            // number of excited states
            if(inputTerms[j].compare(this->stringCISNStates) == 0){
               int nStates = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetNumberExcitedStatesCIS(nStates);
               j++;
            }
            // Davidson is used or not
            if(inputTerms[j].compare(this->stringCISDavidson) == 0){
               if(inputTerms[j+1].compare(this->stringYES) == 0){
                  Parameters::GetInstance()->SetIsDavidsonCIS(true);
               }
               else{
                  Parameters::GetInstance()->SetIsDavidsonCIS(false);
               }
               j++;
            }
            // max iterations for the Davidson roop
            if(inputTerms[j].compare(this->stringCISMaxIter) == 0){
               int maxIter = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetMaxIterationsCIS(maxIter);
               j++;
            }
            // max dimensions for the Davidson expansion
            if(inputTerms[j].compare(this->stringCISMaxDimensions) == 0){
               int maxDim = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetMaxDimensionsCIS(maxDim);
               j++;
            }
            // nolm tolerance for the norm of the resiudal vectors of the Davidson.
            if(inputTerms[j].compare(this->stringCISNormTolerance) == 0){
               double normTol = atof(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetNormToleranceCIS(normTol);
               j++;
            }
            j++;   
         }
         i = j;
      }

      // MD condition
      if(inputTerms[i].compare(this->stringMD) == 0){
         Parameters::GetInstance()->SetRequiresMD(true);
         int j=i+1;
         while(inputTerms[j].compare(this->stringMDEnd) != 0){
            // number of total steps 
            if(inputTerms[j].compare(this->stringMDTotalSteps) == 0){
               int totalSteps = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetTotalStepsMD(totalSteps);
               j++;
            }
            // index of electronic eigen state on whichi MD runs. 
            if(inputTerms[j].compare(this->stringMDElecState) == 0){
               int elecStateIndex = atoi(inputTerms[j+1].c_str());
               Parameters::GetInstance()->SetElectronicStateIndexMD(elecStateIndex);
               j++;
            }
            // time width for MD.
            if(inputTerms[j].compare(this->stringMDTimeWidth) == 0){
               double dt = atof(inputTerms[j+1].c_str()) * Parameters::GetInstance()->GetFs2AU();
               Parameters::GetInstance()->SetTimeWidthMD(dt);
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

            // MNDO
            else if(inputTerms[j].compare(this->stringTheoryMNDO) == 0){
               Parameters::GetInstance()->SetCurrentTheory(MNDO);
            }

            // AM1
            else if(inputTerms[j].compare(this->stringTheoryAM1) == 0){
               Parameters::GetInstance()->SetCurrentTheory(AM1);
            }

            // PM3
            else if(inputTerms[j].compare(this->stringTheoryPM3) == 0){
               Parameters::GetInstance()->SetCurrentTheory(PM3);
            }

            // PM3/PDG
            else if(inputTerms[j].compare(this->stringTheoryPM3PDDG) == 0){
               Parameters::GetInstance()->SetCurrentTheory(PM3PDDG);
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

   // calculate basics and check conditions
   this->CalcMolecularBasics(molecule);
   if(Parameters::GetInstance()->RequiresCIS()){
      this->CheckCisConditions(molecule);
   }
   if(Parameters::GetInstance()->RequiresMD()){
      this->CheckMdConditions();
   }

   // output conditions
   this->OutputMolecularBasics(molecule);
   this->OutputScfConditions();
   if(Parameters::GetInstance()->RequiresCIS()){
      this->OutputCisConditions();
   }
   if(Parameters::GetInstance()->RequiresMD()){
      this->OutputMdConditions();
   }
   if(0<Parameters::GetInstance()->GetIndecesMOPlot().size()){
      this->OutputMOPlotConditions();
   }

   // output inputs
   this->OutputInputTerms(inputTerms);
   cout << messageDoneParseInput;

}

void InputParser::CalcMolecularBasics(Molecule* molecule){

   molecule->CalcTotalNumberAOs();
   molecule->CalcTotalNumberValenceElectrons();
   molecule->CalcXyzCOM();
   molecule->CalcXyzCOC();
}

void InputParser::CheckCisConditions(Molecule* molecule){

   // direct CIS
   int numberOcc = molecule->GetTotalNumberValenceElectrons()/2;
   int numberVir = molecule->GetTotalNumberAOs() - numberOcc;

   // check the number of active occupied orbitals.
   if(numberOcc < Parameters::GetInstance()->GetActiveOccCIS()){
      Parameters::GetInstance()->SetActiveOccCIS(numberOcc);
   }   

   // check the number of active virtual orbitals.
   if(numberVir < Parameters::GetInstance()->GetActiveVirCIS()){
      Parameters::GetInstance()->SetActiveVirCIS(numberVir);
   }   

   // check the number of calculated excited states.
   int numberSlaterDeterminants = Parameters::GetInstance()->GetActiveOccCIS() 
                                 *Parameters::GetInstance()->GetActiveVirCIS();
   if(!Parameters::GetInstance()->IsDavidsonCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(numberSlaterDeterminants);
   }
   else{
      if(numberSlaterDeterminants < Parameters::GetInstance()->GetNumberExcitedStatesCIS()){
         Parameters::GetInstance()->SetNumberExcitedStatesCIS(numberSlaterDeterminants);
      }
      if(numberSlaterDeterminants < Parameters::GetInstance()->GetMaxDimensionsCIS()){
         Parameters::GetInstance()->SetMaxDimensionsCIS(numberSlaterDeterminants);
      }
   }
   
}

void InputParser::CheckMdConditions(){
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      int groundStateIndex = 0;
      Parameters::GetInstance()->SetElectronicStateIndexMD(groundStateIndex);
   }
}

void InputParser::OutputMolecularBasics(Molecule* molecule){

   molecule->OutputTotalNumberAtomsAOsValenceelectrons();
   molecule->OutputConfiguration();
   molecule->OutputXyzCOM();
   molecule->OutputXyzCOC();
}

void InputParser::OutputScfConditions(){

   cout << this->messageScfConditions;
   printf("%s%d\n",this->messageScfMaxIterations.c_str(),Parameters::GetInstance()->GetMaxIterationsSCF());
   printf("%s%e\n",this->messageScfRmsDensity.c_str(),Parameters::GetInstance()->GetThresholdSCF());
   printf("%s%e\n",this->messageScfDampingThresh.c_str(),Parameters::GetInstance()->GetDampingThreshSCF());
   printf("%s%e\n",this->messageScfDampingWeight.c_str(),Parameters::GetInstance()->GetDampingWeightSCF());
   printf("%s%d\n",this->messageScfDiisNumErrorVect.c_str(),Parameters::GetInstance()->GetDiisNumErrorVectSCF());
   printf("%s%e\n",this->messageScfDiisStartError.c_str(),Parameters::GetInstance()->GetDiisStartErrorSCF());
   printf("%s%e\n",this->messageScfDiisEndError.c_str(),Parameters::GetInstance()->GetDiisEndErrorSCF());
   cout << "\n";

}

void InputParser::OutputCisConditions(){
   cout << this->messageCisConditions;

   printf("%s%d\n",this->messageCisNumberActiveOcc.c_str(),Parameters::GetInstance()->GetActiveOccCIS());
   printf("%s%d\n",this->messageCisNumberActiveVir.c_str(),Parameters::GetInstance()->GetActiveVirCIS());
   printf("%s%d\n",this->messageCisNumberExcitedStates.c_str(),Parameters::GetInstance()->GetNumberExcitedStatesCIS());
   printf("%s",this->messageCisDavidson.c_str());
   if(Parameters::GetInstance()->IsDavidsonCIS()){
      printf("%s\n",this->stringYES.c_str());
      printf("%s%d\n",this->messageCisMaxIterations.c_str(),Parameters::GetInstance()->GetMaxIterationsCIS());
      printf("%s%d\n",this->messageCisMaxDimensions.c_str(),Parameters::GetInstance()->GetMaxDimensionsCIS());
      printf("%s%e\n",this->messageCisNormTolerance.c_str(),Parameters::GetInstance()->GetNormToleranceCIS());
   }
   else{
      printf("%s\n",this->stringNO.c_str());
   }

   cout << "\n";
}

void InputParser::OutputMdConditions(){
   cout << this->messageMdConditions;

   printf("%s%d\n",this->messageMdElecState.c_str(),Parameters::GetInstance()->GetElectronicStateIndexMD());
   printf("%s%d\n",this->messageMdTotalSteps.c_str(),Parameters::GetInstance()->GetTotalStepsMD());
   printf("%s%lf%s\n",this->messageMdTimeWidth.c_str(),Parameters::GetInstance()->GetTimeWidthMD()/Parameters::GetInstance()->GetFs2AU(),this->messageFs.c_str());

   cout << "\n";
}

void InputParser::OutputMOPlotConditions(){
   cout << this->messageMOPlotConditions;
   vector<int> moIndeces = Parameters::GetInstance()->GetIndecesMOPlot();
   for(int i=0; i<moIndeces.size(); i++){
      printf("%s%d\n", this->messageMOPlotIndex.c_str(), moIndeces[i]);
   }
   int* gridNum = Parameters::GetInstance()->GetGridNumberMOPlot();
   printf("%s%d %d %d\n", this->messageMOPlotGridNumber.c_str(), 
                              gridNum[XAxis], 
                              gridNum[YAxis],
                              gridNum[ZAxis]);
   double* frameLength = Parameters::GetInstance()->GetFrameLengthMOPlot();
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   printf("%s%e %e %e\n", this->messageMOPlotFrameLength.c_str(), 
                              frameLength[XAxis]/ang2AU, 
                              frameLength[YAxis]/ang2AU,
                              frameLength[ZAxis]/ang2AU);
   printf("%s%s\n", this->messageMOPlotFilePrefix.c_str(),
                      Parameters::GetInstance()->GetFileNamePrefixMOPlot().c_str());

   cout << endl;
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





