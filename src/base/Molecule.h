#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE
namespace MolDS_base{

class Molecule{
public:
   Molecule();
   ~Molecule();
   std::vector<MolDS_base_atoms::Atom*>* GetAtomVect() const; 
   double* GetXyzCOM();
   double* GetXyzCOC();
   void CalcXyzCOM();
   void CalcXyzCOC();
   int GetTotalNumberAOs() const;
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons() const;
   void CalcTotalNumberValenceElectrons();
   void OutputXyzCOM() const;
   void OutputXyzCOC() const;
   void OutputTotalNumberAtomsAOsValenceelectrons() const;
   void OutputConfiguration() const;
   void OutputMomenta() const;
   void CalcPrincipalAxes();
   void Rotate();
   void Translate();
   double GetDistanceAtoms(int atomAIndex, int atomBIndex) const;
   double GetDistanceAtoms(const MolDS_base_atoms::Atom& atomA, 
                           const MolDS_base_atoms::Atom& atomB) const;
private:
   std::vector<MolDS_base_atoms::Atom*>* atomVect;
   double* xyzCOM; // x, y, z coordinates of Center of Mass;
   double* xyzCOC; // x, y, z coordinates of Center of Core;
   bool wasCalculatedXyzCOM;
   bool wasCalculatedXyzCOC;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   void SetMessages();
   void CalcInertiaTensor(double** inertiaTensor, 
                          double const* inertiaTensorOrigin);
   void FreeInertiaTensorMoments(double*** inertiaTensor, 
                                 double** inertiaMoments);
   void Rotate(MolDS_base::EularAngle eularAngle, 
               const double* rotatingOrigin, 
               RotatedObjectType rotatedObj);
   void OutputPrincipalAxes(double const* const* inertiaTensor, 
                            double const* inertiaMoments) const;
   void OutputInertiaTensorOrigin(double* inertiaTensorOrigin) const;
   void OutputRotatingConditions(RotatingType rotatingType, 
                                 double const* rotatingOrigin, 
                                 double const* rotatingAxis, 
                                 double rotatingAngle, 
                                 MolDS_base::EularAngle rotatingEularAngles)const;
   void OutputTranslatingConditions(double const* translatingDifference) const;
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





