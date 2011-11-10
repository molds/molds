#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE
namespace MolDS_base{

class Molecule{
public:
   Molecule();
   ~Molecule();
   std::vector<MolDS_base_atoms::Atom*>* GetAtomVect(); 
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
   double GetDistanceAtoms(MolDS_base_atoms::Atom* atomA, MolDS_base_atoms::Atom* atomB);
private:
   std::vector<MolDS_base_atoms::Atom*>* atomVect;
   double* xyzCOM; // x, y, z coordinates of Center of Mass;
   double* xyzCOC; // x, y, z coordinates of Center of Core;
   bool wasCalculatedXyzCOM;
   bool wasCalculatedXyzCOC;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   void CalcInertiaTensor(double** inertiaTensor, double* inertiaTensorOrigin);
   void FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments);
   void Rotate(MolDS_base::EularAngle eularAngle, double* rotatingOrigin, RotatedObjectType rotatedObj);
   void OutputPrincipalAxes(double** inertiaTensor, double* inertiaMoments);
   void OutputInertiaTensorOrigin(double* inertiaTensorOrigin);
   void OutputRotatingConditions(RotatingType rotatingType, double* rotatingOrigin, 
                                 double* rotatingAxis, double rotatingAngle, 
                                 MolDS_base::EularAngle rotatingEularAngles);
   void OutputTranslatingConditions(double* translatingDifference);
   std::string messageTotalNumberAOs;
   std::string messageTotalNumberAtoms;
   std::string messageTotalNumberValenceElectrons;
   std::string messageConfiguration;
   std::string messageConfigurationTitleAU;
   std::string messageConfigurationTitleAng;
   std::string messageMomenta;
   std::string messageMomentaTitleAU;
   std::string messageMomentaTitleGAngMolinFsin;
   std::string messageCOM;
   std::string messageCOC;
   std::string messageCOMTitleAU;
   std::string messageCOMTitleAng;
   std::string messageStartPrincipalAxes;
   std::string messageDonePrincipalAxes;
   std::string messagePrincipalAxes;
   std::string messagePrincipalAxesTitleAU;
   std::string messagePrincipalAxesTitleAng;
   std::string messageInertiaTensorOrigin;
   std::string messageInertiaTensorOriginTitleAU;
   std::string messageInertiaTensorOriginTitleAng;
   std::string messageStartRotate;
   std::string messageDoneRotate;
   std::string messageRotatingOrigin;
   std::string messageRotatingOriginTitleAU;
   std::string messageRotatingOriginTitleAng;
   std::string messageRotatingAxis;
   std::string messageRotatingAxisTitleAU;
   std::string messageRotatingAxisTitleAng;
   std::string messageRotatingAngle;
   std::string messageRotatingType;
   std::string messageRotatingEularAngles;
   std::string messageRotatingEularAnglesTitle;
   std::string messageStartTranslate;
   std::string messageDoneTranslate;
   std::string messageTranslatingDifference;
   std::string messageTranslatingDifferenceTitleAU;
   std::string messageTranslatingDifferenceTitleAng;
};
}
#endif





