#ifndef INCLUDED_ATOM
#define INCLUDED_ATOM

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include"../MallocerFreer.h"
#include"../Parameters.h"

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Atom{
private:
   void SetMessages();
   string errorMessageImuAmu;
   string errorMessageOrbitalExponent;
   string errorMessageIndoCoulombInt;
   string errorMessageIndoExchangeInt;
   string errorMessageShellType;
   string errorMessageEffectivPrincipalQuantumNumber;
protected:
   string errorMessageIndoCoreIntegral;
   string errorMessageAtomType;
   string errorMessageOrbitalType;
   double* xyz;
   AtomType atomType;
   vector<OrbitalType> valence;
   double imuAmuS;
   double imuAmuP;
   double imuAmuD;
   double bondingParameter;             // see Table 3.2 and 3.4 in J. A. Pople book
   double coreCharge;                   // = Z_A
   double effectiveNuclearChargeK;
   double effectiveNuclearChargeL;
   double effectiveNuclearChargeM;
   ShellType valenceShellType;
   int firstAOIndex;
   int GetEffectivePrincipalQuantumNumber(ShellType shellType);
   int numberValenceElectrons;
   double indoF2;                   // see (3.89) in J. A. Pople book
   double indoG1;                   // see (3.88) in J. A. Pople book
public:
   Atom();
   Atom(double x, double y, double z);
   ~Atom();
   AtomType GetAtomType();
   double* GetXyz();
   void SetXyz(double x, double y, double z);
   vector<OrbitalType> GetValence();
   double GetBondingParameter();
   double GetCoreCharge();
   int GetFirstAOIndex();
   void SetFirstAOIndex(int firstAOIndex);
   ShellType GetValenceShellType();
   int GetNumberValenceElectrons();
   double GetImuAmu(OrbitalType orbitalType);  // return 0.5*(I_mu + A_mu)
   double GetOrbitalExponent(ShellType shellType, OrbitalType orbitalType);  // (1.73) in J. A. Pople book.
   double GetIndoCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma); // (3.87) - (3.91) in J. A. Pople book.
   double GetIndoExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma); // (3.87) - (3.91) in J. A. Pople book.
   virtual double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess) = 0; // P82 - 83 in J. A. Pople book. 
};

Atom::Atom(){
   this->xyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->SetMessages();
}

Atom::Atom(double x, double y, double z){
   this->xyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->SetMessages();
   this->SetXyz(x, y, z);
}

Atom::~Atom(){
   if(xyz != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(xyz);
      xyz = NULL;
      //cout << "xyz deleted\n";
   }
   //cout << "atom deleted\n";
}

void Atom::SetMessages(){
   this->errorMessageImuAmu = "Error in base_atoms::Atom::GetImuAmu: Invalid orbitalType.\n";
   this->errorMessageOrbitalExponent = "Error in base_atoms::Atom::GetOrbitalExponent: Invalid shelltype or orbitalType.\n";
   this->errorMessageIndoCoulombInt = "Error in base_atoms::Atom::GetIndoCoulombInt: Invalid orbitalType.\n";
   this->errorMessageIndoExchangeInt = "Error in base_atoms::Atom::GetIndoExchangeInt: Invalid orbitalType.\n";
   this->errorMessageIndoCoreIntegral = "Error in base_atoms::Atom::GetIndoCoreIntegral: Invalid orbitalType.\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageShellType = "\tshell type = ";
   this->errorMessageEffectivPrincipalQuantumNumber = 
      "Error in base::Atom::GetEffectivePrincipalQuantumNumber: invalid shelltype.\n";
}

AtomType Atom::GetAtomType(){
   return this->atomType;
}

double* Atom::GetXyz(){
   return this->xyz;
}

void Atom::SetXyz(double x, double y, double z){
   xyz[0]= x * Parameters::GetInstance()->GetAngstrom2AU();
   xyz[1]= y * Parameters::GetInstance()->GetAngstrom2AU();
   xyz[2]= z * Parameters::GetInstance()->GetAngstrom2AU();
}

vector<OrbitalType> Atom::GetValence(){
   return this->valence;
}

double Atom::GetBondingParameter(){
   return this->bondingParameter;
}

double Atom::GetCoreCharge(){
   return this->coreCharge;
}

int Atom::GetFirstAOIndex(){
   return this->firstAOIndex;
}

void Atom::SetFirstAOIndex(int firstAOIndex){
   this->firstAOIndex = firstAOIndex;
}

ShellType Atom::GetValenceShellType(){
   return this->valenceShellType;
}

int Atom::GetEffectivePrincipalQuantumNumber(ShellType shellType){
   if(shellType == k){
      return 1.0;
   }
   else if(shellType == l){
      return 2.0;
   }
   else if(shellType == m){
      return 3.0;
   }
   else{
      cout << this->errorMessageEffectivPrincipalQuantumNumber;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      exit(EXIT_FAILURE);
   }
}   

int Atom::GetNumberValenceElectrons(){
   return this->numberValenceElectrons;
}


// return 0.5*(I_mu + A_mu)
double Atom::GetImuAmu(OrbitalType orbitalType){
   if(orbitalType == s){ 
      return this->imuAmuS;
   }   
   else if(orbitalType == px || orbitalType == py || orbitalType == pz ){
      return this->imuAmuP;
   }   
   else if(orbitalType == dxy || 
           orbitalType == dyz || 
           orbitalType == dzz || 
           orbitalType == dzx || 
           orbitalType == dxxyy ){
      return this->imuAmuD;
   }   
   else{
      cout << errorMessageImuAmu;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
      exit(EXIT_FAILURE);
   }   
}


// (1.73) in J. A. Pople book
double Atom::GetOrbitalExponent(ShellType shellType, OrbitalType orbitalType){
   if(shellType == k && orbitalType == s){ 
      return this->effectiveNuclearChargeK/this->GetEffectivePrincipalQuantumNumber(shellType);
   }   
   else if(shellType == l && (orbitalType == s || orbitalType == px || orbitalType == py || orbitalType == pz)){
      return this->effectiveNuclearChargeL/this->GetEffectivePrincipalQuantumNumber(shellType);
   }   
   else if(shellType == m && (orbitalType == s  || 
                              orbitalType == px || 
                              orbitalType == py || 
                              orbitalType == pz ||
                              orbitalType == dxy ||
                              orbitalType == dyz ||
                              orbitalType == dzz ||
                              orbitalType == dzx ||
                              orbitalType == dxxyy)){
      return this->effectiveNuclearChargeM/this->GetEffectivePrincipalQuantumNumber(shellType);
   }   
   else{
      cout << this->errorMessageOrbitalExponent;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
      exit(EXIT_FAILURE);
   }   
}


// (3.87) - (3.91) in J. A. Pople book.
// Indo Coulomb Interaction
double Atom::GetIndoCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma){

   double value=0.0;
   if( orbital1 == s && orbital2 == s){
      value = gamma;
   }
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){
      value = gamma;
   }
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz ) && orbital2 == s){
      value = gamma;
   }
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){
      value = gamma + 4.0*this->indoF2/25.0;
   }
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = gamma - 2.0*this->indoF2/25.0;
   }
   else{
      cout << this->errorMessageIndoCoulombInt;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      exit(EXIT_FAILURE);
   }

   return value;
}

// (3.87) - (3.91) in J. A. Pople book.
// Indo Exchange interaction
double Atom::GetIndoExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma){
   double value=0.0;

   if( orbital1 == orbital2){
      value = this->GetIndoCoulombInt(orbital1, orbital2, gamma);
   }
   else if( (orbital1 == s) && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = this->indoG1/3.0;
   }
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz) && orbital2 == s  ){
      value = this->indoG1/3.0;
   }
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = 3.0*this->indoF2/25.0;
   }
   else{
      cout << this->errorMessageIndoExchangeInt;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      exit(EXIT_FAILURE);
   }

   return value;
}

}
#endif





