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
#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE
namespace MolDS_base{

class Molecule{
public:
   Molecule();
   Molecule(const Molecule& rhs);
   Molecule& operator=(const Molecule& rhs);
   ~Molecule();
   std::vector<MolDS_base_atoms::Atom*>* GetAtomVect() const; 
   double* GetXyzCOM() const;
   double* GetXyzCOM();
   double* GetXyzCOC() const;
   double* GetXyzCOC();
   void CalcXyzCOM();
   void CalcXyzCOC();
   void CalcBasics();
   int GetTotalNumberAOs() const;
   int GetTotalNumberValenceElectrons() const;
   double GetTotalCoreMass() const;
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
   double totalCoreMass;
   void Initialize();
   void CopyInitialize(const Molecule& rhs);
   void Finalize(std::vector<MolDS_base_atoms::Atom*>** atomVect, double** xyzCOM, double**xyzCOC);
   void SetMessages();
   void CalcTotalNumberValenceElectrons();
   void CalcTotalNumberAOs();
   void CalcTotalCoreMass();
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
   std::string errorMessageGetAtomVectNull;
   std::string errorMessageGetXyzCOCNull;
   std::string errorMessageGetXyzCOMNull;
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





