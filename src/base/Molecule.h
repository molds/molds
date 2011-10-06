#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE

using namespace std;
using namespace MolDS_base_atoms;

namespace MolDS_base{

class Molecule{
public:
   Molecule();
   ~Molecule();
   vector<Atom*>* GetAtomVect(); 
   double* GetXyzCOM();
   double* GetXyzCOC();
   void CalcXyzCOM();
   void CalcXyzCOC();
   int GetTotalNumberAOs();
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons();
   void CalcTotalNumberValenceElectrons();
   void OutputXyzCOM();
   void OutputXyzCOC();
   void OutputTotalNumberAtomsAOsValenceelectrons();
   void OutputConfiguration();
   void OutputMomenta();
   void CalcPrincipalAxes();
   void Rotate();
   void Translate();
   double GetDistanceAtoms(int atomAIndex, int atomBIndex);
   double GetDistanceAtoms(Atom* atomA, Atom* atomB);
private:
   vector<Atom*>* atomVect;
   double* xyzCOM; // x, y, z coordinates of Center of Mass;
   double* xyzCOC; // x, y, z coordinates of Center of Core;
   bool wasCalculatedXyzCOM;
   bool wasCalculatedXyzCOC;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   void CalcInertiaTensor(double** inertiaTensor, double* inertiaTensorOrigin);
   void FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments);
   void Rotate(EularAngle eularAngle, double* rotatingOrigin, RotatedObjectType rotatedObj);
   void OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments);
   void OutputInertiaTensorOrigin(double* inertiaTensorOrigin);
   void OutputRotatingConditions(RotatingType rotatingType, double* rotatingOrigin, 
                                 double* rotatingAxis, double rotatingAngle, 
                                 EularAngle rotatingEularAngles);
   void OutputTranslatingConditions(double* translatingDifference);
   string messageTotalNumberAOs;
   string messageTotalNumberAtoms;
   string messageTotalNumberValenceElectrons;
   string messageConfiguration;
   string messageConfigurationTitleAU;
   string messageConfigurationTitleAng;
   string messageMomenta;
   string messageMomentaTitleAU;
   string messageMomentaTitleGAngMolinFsin;
   string messageCOM;
   string messageCOC;
   string messageCOMTitleAU;
   string messageCOMTitleAng;
   string messageStartPrincipalAxes;
   string messageDonePrincipalAxes;
   string messagePrincipalAxes;
   string messagePrincipalAxesTitleAU;
   string messagePrincipalAxesTitleAng;
   string messageInertiaTensorOrigin;
   string messageInertiaTensorOriginTitleAU;
   string messageInertiaTensorOriginTitleAng;
   string messageStartRotate;
   string messageDoneRotate;
   string messageRotatingOrigin;
   string messageRotatingOriginTitleAU;
   string messageRotatingOriginTitleAng;
   string messageRotatingAxis;
   string messageRotatingAxisTitleAU;
   string messageRotatingAxisTitleAng;
   string messageRotatingAngle;
   string messageRotatingType;
   string messageRotatingEularAngles;
   string messageRotatingEularAnglesTitle;
   string messageStartTranslate;
   string messageDoneTranslate;
   string messageTranslatingDifference;
   string messageTranslatingDifferenceTitleAU;
   string messageTranslatingDifferenceTitleAng;
};

Molecule::Molecule(){
   this->atomVect = new vector<Atom*>;
   this->xyzCOM = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->xyzCOC = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->wasCalculatedXyzCOM = false;
   this->wasCalculatedXyzCOC = false;
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageConfiguration = "\tMolecular configration:\n";
   this->messageConfigurationTitleAU = "\t\t| i-th | atom type | x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageConfigurationTitleAng = "\t\t| i-th | atom type | x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageMomenta = "\tMomenta of each atom:\n";
   this->messageMomentaTitleAU = "\t\t| i-th | atom type | px[a.u.] | py[a.u.] | pz[a.u.] |\n ";
   this->messageMomentaTitleGAngMolinFsin = "\t\t| i-th | atom type | px | py | pz | [(g/Mol)*(angst/fs)]\n ";
   this->messageCOM = "\tCenter of Mass:\n";
   this->messageCOC = "\tCenter of Core:\n";
   this->messageCOMTitleAU = "\t\t| x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageCOMTitleAng = "\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageStartPrincipalAxes = "**********  START: Principal Axes of Inertia  **********\n";
   this->messageDonePrincipalAxes =  "**********  DONE: Principal Axes of Inertia  ***********\n\n\n";
   this->messagePrincipalAxes = "\tPrincipal Axes:\n";
   this->messagePrincipalAxesTitleAU = "\t\t| inertia moments [a.u.] | x[a.u.] | y[a.u.] | z[a.u.] | (normalized)\n";
   this->messagePrincipalAxesTitleAng = "\t\t| inertia moments [g*angust**2/mol] | x[angst.] | y[angst.] | z[angst.] | (not normalized)\n";
   this->messageInertiaTensorOrigin = "\tInertia Tensor Origin:\n";
   this->messageInertiaTensorOriginTitleAU = "\t\t| x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageInertiaTensorOriginTitleAng = "\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageStartRotate = "**********  START: Rotate molecule  **********\n";
   this->messageDoneRotate =  "**********  DONE: Rotate molecule  ***********\n\n\n";
   this->messageRotatingOrigin = "\tRotating Origin:\n";
   this->messageRotatingOriginTitleAU = "\t\t| x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageRotatingOriginTitleAng = "\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAxis = "\tRotating Axis:\n";
   this->messageRotatingAxisTitleAU = "\t\t| x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageRotatingAxisTitleAng = "\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAngle = "\tRotating Angle [degree]: ";
   this->messageRotatingType = "\tRotating Type: ";
   this->messageRotatingEularAngles = "\tRotating Eular Angles:\n";
   this->messageRotatingEularAnglesTitle = "\t\t| alpha[degree] | beta[degree] | gamma[degree] |\n";
   this->messageStartTranslate = "**********  START: Translate molecule  **********\n";
   this->messageDoneTranslate =  "**********  DONE: Translate molecule  ***********\n\n\n";
   this->messageTranslatingDifference = "\tTranslating Difference:\n";
   this->messageTranslatingDifferenceTitleAU = "\t\t| x[a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageTranslatingDifferenceTitleAng = "\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
}

Molecule::~Molecule(){
   if(this->atomVect != NULL){
      for(int i=0; i<this->atomVect->size(); i++){
         delete (*atomVect)[i];
         (*atomVect)[i] = NULL;
      }
      delete this->atomVect;
      this->atomVect = NULL;
      //cout << "atomVect deleted\n";
   }
   if(this->xyzCOM != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->xyzCOM);
      //cout << "xyzCOM deleted\n";
   }
   if(this->xyzCOC != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->xyzCOC);
      //cout << "xyzCOC deleted\n";
   }
}

vector<Atom*>* Molecule::GetAtomVect(){
   return this->atomVect;
}

double* Molecule::GetXyzCOM(){
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   return this->xyzCOM;
}

double* Molecule::GetXyzCOC(){
   if(!this->wasCalculatedXyzCOC){
      this->CalcXyzCOC();
   }
   return this->xyzCOC;
}

void Molecule::CalcXyzCOM(){
   double totalAtomicMass = 0.0;
   Atom* atom;
   double* atomicXyz;
   double atomicMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOM[j] = 0.0;
   }
      
   for(int i=0; i<this->atomVect->size(); i++){
      atom = (*this->atomVect)[i]; 
      atomicXyz = atom->GetXyz();
      atomicMass = atom->GetAtomicMass();
      totalAtomicMass += atomicMass;
      for(int j=0; j<3; j++){
         this->xyzCOM[j] += atomicXyz[j] * atomicMass;
      }
   }
   for(int i=0; i<3; i++){
      this->xyzCOM[i]/=totalAtomicMass;
   }
   this->wasCalculatedXyzCOM = true;
}

void Molecule::CalcXyzCOC(){
   double totalCoreMass = 0.0;
   Atom* atom;
   double* atomicXyz;
   double coreMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOC[j] = 0.0;
   }
      
   for(int i=0; i<this->atomVect->size(); i++){
      atom = (*this->atomVect)[i]; 
      atomicXyz = atom->GetXyz();
      coreMass = atom->GetAtomicMass() - (double)atom->GetNumberValenceElectrons();
      totalCoreMass += coreMass;
      for(int j=0; j<3; j++){
         this->xyzCOC[j] += atomicXyz[j] * coreMass;
      }
   }
   for(int i=0; i<3; i++){
      this->xyzCOC[i]/=totalCoreMass;
   }
   this->wasCalculatedXyzCOC = true;
}

int Molecule::GetTotalNumberAOs(){
   return this->totalNumberAOs;
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAOs = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAOs);
      this->totalNumberAOs += (*this->atomVect)[i]->GetValence().size();
   }
}

int Molecule::GetTotalNumberValenceElectrons(){
   return this->totalNumberValenceElectrons;
}

void Molecule::CalcTotalNumberValenceElectrons(){
   this->totalNumberValenceElectrons = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      this->totalNumberValenceElectrons += (*this->atomVect)[i]->GetNumberValenceElectrons();
   }
}

void Molecule::OutputConfiguration(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageConfiguration;
   cout << this->messageConfigurationTitleAng;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0]/ang2AU, atom->GetXyz()[1]/ang2AU, atom->GetXyz()[2]/ang2AU);
   }
   cout << "\n";

   cout << this->messageConfigurationTitleAU;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0], atom->GetXyz()[1], atom->GetXyz()[2]);
   }
   cout << "\n";

}

void Molecule::OutputMomenta(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double fs2AU = Parameters::GetInstance()->GetFs2AU();
   double gMolin2AU = Parameters::GetInstance()->GetGMolin2AU();
   double momentumUnit2AU = ang2AU*gMolin2AU/fs2AU;
   cout << this->messageMomenta;
   cout << this-> messageMomentaTitleGAngMolinFsin;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetPxyz()[0]/momentumUnit2AU, 
            atom->GetPxyz()[1]/momentumUnit2AU, 
            atom->GetPxyz()[2]/momentumUnit2AU);
   }
   cout << "\n";

   cout << this->messageMomentaTitleAU;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetPxyz()[0], atom->GetPxyz()[1], atom->GetPxyz()[2]);
   }
   cout << "\n";

}

void Molecule::OutputXyzCOM(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageCOM;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOM[0]/ang2AU,
                             this->xyzCOM[1]/ang2AU,
                             this->xyzCOM[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOM[0],
                             this->xyzCOM[1],
                             this->xyzCOM[2]);
   cout << "\n";

}

void Molecule::OutputXyzCOC(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageCOC;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOC[0]/ang2AU,
                             this->xyzCOC[1]/ang2AU,
                             this->xyzCOC[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->xyzCOC[0],
                             this->xyzCOC[1],
                             this->xyzCOC[2]);
   cout << "\n";

}

void Molecule::OutputTotalNumberAtomsAOsValenceelectrons(){
   cout << this->messageTotalNumberAtoms << this->atomVect->size() << "\n";
   cout << this->messageTotalNumberAOs << this->totalNumberAOs << "\n";
   cout << this->messageTotalNumberValenceElectrons << this->totalNumberValenceElectrons << "\n\n";
}

void Molecule::OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double gMolin2AU = Parameters::GetInstance()->GetGMolin2AU();

   cout << this->messagePrincipalAxes;
   cout << this->messagePrincipalAxesTitleAng;
   for(int i=0; i<3; i++){
      printf("\t\t%e\t%e\t%e\t%e\n",inertiaMoments[i]/gMolin2AU, 
                                    inertiaTensor[i][0]/ang2AU,
                                    inertiaTensor[i][1]/ang2AU,
                                    inertiaTensor[i][2]/ang2AU);
   }
   cout << "\n";

   cout << this->messagePrincipalAxesTitleAU;
   for(int i=0; i<3; i++){
      printf("\t\t%e\t%e\t%e\t%e\n",inertiaMoments[i], 
                                    inertiaTensor[i][0],
                                    inertiaTensor[i][1],
                                    inertiaTensor[i][2]);
   }
   cout << "\n";

}

void Molecule::OutputInertiaTensorOrigin(double* inertiaTensorOrigin){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();

   cout << this->messageInertiaTensorOrigin;
   cout << this->messageInertiaTensorOriginTitleAng;
   printf("\t\t%e\t%e\t%e\n",inertiaTensorOrigin[0]/ang2AU,
                             inertiaTensorOrigin[1]/ang2AU,
                             inertiaTensorOrigin[2]/ang2AU);
   cout << "\n";

   cout << this->messageInertiaTensorOriginTitleAU;
   printf("\t\t%e\t%e\t%e\n",inertiaTensorOrigin[0],
                             inertiaTensorOrigin[1],
                             inertiaTensorOrigin[2]);
   cout << "\n";

}

void Molecule::CalcPrincipalAxes(){

   cout << this->messageStartPrincipalAxes;

   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double inertiaTensorOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetInertiaTensorOrigin() != NULL){
      inertiaTensorOrigin[0] = Parameters::GetInstance()->GetInertiaTensorOrigin()[0];
      inertiaTensorOrigin[1] = Parameters::GetInstance()->GetInertiaTensorOrigin()[1];
      inertiaTensorOrigin[2] = Parameters::GetInstance()->GetInertiaTensorOrigin()[2];
   }

   double** inertiaTensor = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(3, 3);
   double*  inertiaMoments = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);

   try{
      this->CalcInertiaTensor(inertiaTensor, inertiaTensorOrigin);
      
      bool calcEigenVectors = true;
      MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(inertiaTensor,
                                                              inertiaMoments,
                                                              3,
                                                              calcEigenVectors);
      this->OutputPrincipalAxes(inertiaTensor, inertiaMoments);
      this->OutputInertiaTensorOrigin(inertiaTensorOrigin);
   }
   catch(MolDSException ex){
      this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);
      throw ex;
   }
   this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);

   cout << this->messageDonePrincipalAxes;
   
}

void Molecule::CalcInertiaTensor(double** inertiaTensor, double* inertiaTensorOrigin){

   Atom* atom;
   double x;
   double y;
   double z;
   double atomicMass;
   for(int a=0; a<this->atomVect->size(); a++){
      atom = (*this->atomVect)[a];
      atomicMass = atom->GetAtomicMass();
      x = atom->GetXyz()[0] - inertiaTensorOrigin[0];
      y = atom->GetXyz()[1] - inertiaTensorOrigin[1];
      z = atom->GetXyz()[2] - inertiaTensorOrigin[2];

      inertiaTensor[0][0] += atomicMass*(y*y + z*z);
      inertiaTensor[0][1] -= atomicMass*x*y;
      inertiaTensor[0][2] -= atomicMass*x*z;

      inertiaTensor[1][0] -= atomicMass*y*x;
      inertiaTensor[1][1] += atomicMass*(x*x + z*z);
      inertiaTensor[1][2] -= atomicMass*y*z;

      inertiaTensor[2][0] -= atomicMass*z*x;
      inertiaTensor[2][1] -= atomicMass*z*y;
      inertiaTensor[2][2] += atomicMass*(x*x + y*y);

   }
   
}

void Molecule::FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments){

   if(*inertiaTensor != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(inertiaTensor, 3);
      //cout << "inertiaTensor deleted\n";
   }

   if(*inertiaMoments != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(inertiaMoments);
      //cout << "inertiaMoments deleted\n";
   }

}

void Molecule::Rotate(){

   cout << this->messageStartRotate;

   // Default values are set if some conditions are not specified.
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double rotatingOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetRotatingOrigin() != NULL){
      rotatingOrigin[0] = Parameters::GetInstance()->GetRotatingOrigin()[0];
      rotatingOrigin[1] = Parameters::GetInstance()->GetRotatingOrigin()[1];
      rotatingOrigin[2] = Parameters::GetInstance()->GetRotatingOrigin()[2];
   }

   RotatingType rotatingType = Parameters::GetInstance()->GetRotatingType();
   double* rotatingAxis = Parameters::GetInstance()->GetRotatingAxis();
   EularAngle rotatingEularAngles = Parameters::GetInstance()->GetRotatingEularAngles();
   double rotatingAngle = Parameters::GetInstance()->GetRotatingAngle();

   this->OutputRotatingConditions(rotatingType, rotatingOrigin, 
                                  rotatingAxis, rotatingAngle, 
                                  rotatingEularAngles);

   // rotate
   if(rotatingType == Axis){
      EularAngle setZAxisEularAngles(rotatingAxis[0], rotatingAxis[1], rotatingAxis[2]);
      EularAngle angleAroundAxis;
      angleAroundAxis.SetAlpha(rotatingAngle);

      this->Rotate(setZAxisEularAngles, rotatingOrigin, Frame);
      this->Rotate(angleAroundAxis, rotatingOrigin, System);
      this->Rotate(setZAxisEularAngles, rotatingOrigin, System);
   }
   else if(rotatingType == Eular){
      this->Rotate(rotatingEularAngles, rotatingOrigin, System);
   }
   
   this->OutputConfiguration();
   cout << this->messageDoneRotate;
}

/***
 * rotatedObj == System: Molecule is rotated.
 * rotatedObj == Frame: De Cartesian is rotated.
 */
void Molecule::Rotate(EularAngle eularAngle, double* rotatingOrigin, RotatedObjectType rotatedObj){

   double rotatingMatrixAlpha[3][3];
   double rotatingMatrixBeta[3][3];
   double rotatingMatrixGamma[3][3];
   double inv = 1.0;
   if(rotatedObj == System){
      inv = -1.0;
   }

   CalcRotatingMatrix(rotatingMatrixAlpha, inv*eularAngle.GetAlpha(), ZAxis);
   CalcRotatingMatrix(rotatingMatrixBeta, inv*eularAngle.GetBeta(), YAxis);
   CalcRotatingMatrix(rotatingMatrixGamma, inv*eularAngle.GetGamma(), ZAxis);

   double temp1[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp1[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixGamma[k][j];
            }
            else if(rotatedObj == Frame){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixAlpha[k][j];
            }
         }
      }
   }

   double temp2[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp2[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp2[i][j] += rotatingMatrixAlpha[i][k] * temp1[k][j];
            }
            else if(rotatedObj == Frame){
               temp2[i][j] += rotatingMatrixGamma[i][k] * temp1[k][j];
            }
         }
      }
   }

   double rotatedXyz[3];
   Atom* atom;
   for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         for(int j=0; j<3; j++){
            rotatedXyz[j] = 0.0;
            for(int k=0; k<3; k++){
               rotatedXyz[j] += temp2[j][k] * (atom->GetXyz()[k] - rotatingOrigin[k]);
            }
         }
         for(int j=0; j<3; j++){
            atom->GetXyz()[j] = rotatedXyz[j] + rotatingOrigin[j];
         }
   }
}

void Molecule::OutputRotatingConditions(RotatingType rotatingType, double* rotatingOrigin, 
                                        double* rotatingAxis, double rotatingAngle, 
                                        EularAngle rotatingEularAngles){

   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double degree2Radian = Parameters::GetInstance()->GetDegree2Radian();

   // type
   printf("%s%s\n\n",this->messageRotatingType.c_str(),RotatingTypeStr(rotatingType));

   // rotating origin
   cout << this->messageRotatingOrigin;
   cout << this->messageRotatingOriginTitleAng;
   printf("\t\t%e\t%e\t%e\n\n",rotatingOrigin[0]/angst2AU,
                               rotatingOrigin[1]/angst2AU,
                               rotatingOrigin[2]/angst2AU);

   cout << this->messageRotatingOriginTitleAU;
   printf("\t\t%e\t%e\t%e\n\n",rotatingOrigin[0],
                               rotatingOrigin[1],
                               rotatingOrigin[2]);

   if(rotatingType == Axis){
      // rotating axis
      cout << this->messageRotatingAxis;
      cout << this->messageRotatingAxisTitleAng;
      printf("\t\t%e\t%e\t%e\n\n",rotatingAxis[0]/angst2AU,
                                  rotatingAxis[1]/angst2AU,
                                  rotatingAxis[2]/angst2AU);

      cout << this->messageRotatingAxisTitleAU;
      printf("\t\t%e\t%e\t%e\n\n",rotatingAxis[0],
                                  rotatingAxis[1],
                                  rotatingAxis[2]);

      // angle
      printf("%s%e\n\n",this->messageRotatingAngle.c_str(),rotatingAngle/degree2Radian);
   }
   else if (rotatingType == Eular){
      // Eular angles
      cout << this->messageRotatingEularAngles;
      cout << this->messageRotatingEularAnglesTitle;
      printf("\t\t%e\t%e\t%e\t\n\n",rotatingEularAngles.GetAlpha()/degree2Radian, 
                                    rotatingEularAngles.GetBeta()/degree2Radian,
                                    rotatingEularAngles.GetGamma()/degree2Radian);
   }

}


void Molecule::Translate(){

   cout << this->messageStartTranslate;

   double x = Parameters::GetInstance()->GetTranslatingDifference()[0];
   double y = Parameters::GetInstance()->GetTranslatingDifference()[1];
   double z = Parameters::GetInstance()->GetTranslatingDifference()[2];

   this->OutputTranslatingConditions(Parameters::GetInstance()->GetTranslatingDifference()); 

   Atom* atom;
   for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         atom->GetXyz()[0] += x;
         atom->GetXyz()[1] += y;
         atom->GetXyz()[2] += z;
   }
   
   this->wasCalculatedXyzCOM = false;
   this->CalcXyzCOM();
   this->wasCalculatedXyzCOC = false;
   this->CalcXyzCOC();

   this->OutputConfiguration();
   this->OutputXyzCOM();
   this->OutputXyzCOC();

   cout << this->messageDoneTranslate;
}

void Molecule::OutputTranslatingConditions(double* translatingDifference){

   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();

   cout << this->messageTranslatingDifference;
   cout << this->messageTranslatingDifferenceTitleAng;
   printf("\t\t%e\t%e\t%e\n\n",translatingDifference[0]/angst2AU,
                               translatingDifference[1]/angst2AU,
                               translatingDifference[2]/angst2AU);

   cout << this->messageTranslatingDifferenceTitleAU;
   printf("\t\t%e\t%e\t%e\n\n",translatingDifference[0],
                               translatingDifference[1],
                               translatingDifference[2]);

}

double Molecule::GetDistanceAtoms(int atomAIndex, int atomBIndex){
   Atom* atomA = (*this->atomVect)[atomAIndex];
   Atom* atomB = (*this->atomVect)[atomBIndex];
   return this->GetDistanceAtoms(atomA, atomB);
}

double Molecule::GetDistanceAtoms(Atom* atomA, Atom* atomB){

   double distance=0.0;
   distance = sqrt( pow(atomA->GetXyz()[0] - atomB->GetXyz()[0], 2.0)
                   +pow(atomA->GetXyz()[1] - atomB->GetXyz()[1], 2.0)
                   +pow(atomA->GetXyz()[2] - atomB->GetXyz()[2], 2.0) );
   return distance;

}





}
#endif





