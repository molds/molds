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
   string errorMessageShellType;
   string errorMessageEffectivPrincipalQuantumNumber;
   string errorMessageZindoCoreIntegral;
   double GetJss();  // Part of Eq. (13) in [BZ_1979]
   double GetJsp();  // Part of Eq. (13) in [BZ_1979]
   double GetJsd();  // Part of Eq. (13) in [BZ_1979]
   double GetJpp();  // Part of Eq. (13) in [BZ_1979]
   double GetJpd();  // Part of Eq. (13) in [BZ_1979]
   double GetJdd();  // Part of Eq. (13) in [BZ_1979]
protected:
   string errorMessageIndoCoreIntegral;
   string errorMessageZindoSCoreIntegral;
   string errorMessageAtomType;
   string errorMessageOrbitalType;
   double* xyz;
   AtomType atomType;
   vector<OrbitalType> valence;
   ShellType valenceShellType;
   int firstAOIndex;
   int numberValenceElectrons;
   double imuAmuS;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuP;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuD;                      // Table 3.4 or 3.5 in J. A. Pople book
   double bondingParameter;             // Table 3.2 and 3.4 in J. A. Pople book
   double coreCharge;                   // = Z_A in J. A. Pople book.
   double effectiveNuclearChargeK;      // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeL;      // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeMsp;    // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeMd;     // Table 1.5 in J. A. Pople book
   int GetEffectivePrincipalQuantumNumber(ShellType shellType); // Table 1.4 in J. A. Pople book
   double indoF2;                   // (3.89) in J. A. Pople book
   double indoG1;                   // (3.88) in J. A. Pople book
   double zindoF0ss;                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double zindoF0sd;                  // Table 1 in [AEZ_1986]
   double zindoF0dd;                  // Table 1 in [AEZ_1986]
   double zindoG1sp;                 // Table 3 in ref. [BZ_1979]
   double zindoF2pp;                 // Table 3 in ref. [BZ_1979]
   double zindoG2sd;                 // Table 3 in ref. [BZ_1979]
   double zindoG1pd;                 // Table 3 in ref. [BZ_1979]
   double zindoF2pd;                 // Table 3 in ref. [BZ_1979]
   double zindoG3pd;                 // Table 3 in ref. [BZ_1979]
   double zindoF2dd;                 // Table 3 in ref. [BZ_1979]
   double zindoF4dd;                 // Table 3 in ref. [BZ_1979]
   double IonPotS;   // Ionization potential, Table 4 in [BZ_1979]
   double IonPotP;   // Ionization potential, Table 4 in [BZ_1979]
   double IonPotD;   // Ionization potential, Table 4 in [BZ_1979]
   double GetZindoCoreIntegral(OrbitalType orbital, int l, int m, int n); // Eq. (13) in [BZ_1979]

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
   virtual double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory) = 0; // P82 - 83 in J. A. Pople book for INDO or Eq. (13) in [BZ_1979] for ZINDO/S
   double GetIndoF2();
   double GetIndoG1();
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
   this->errorMessageIndoCoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for INDO.\n";
   this->errorMessageZindoSCoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for ZINDO/S.\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageShellType = "\tshell type = ";
   this->errorMessageEffectivPrincipalQuantumNumber = 
      "Error in base::Atom::GetEffectivePrincipalQuantumNumber: invalid shelltype.\n";
   this->errorMessageZindoCoreIntegral = "Error in base_stoms::Atom::GetZindoCoreINtegral: Invalid orbitalType.\n";
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
                              orbitalType == pz )){
      return this->effectiveNuclearChargeMsp/this->GetEffectivePrincipalQuantumNumber(shellType);
   }   
   else if(shellType == m && (orbitalType == dxy  || 
                              orbitalType == dyz ||
                              orbitalType == dzz ||
                              orbitalType == dzx ||
                              orbitalType == dxxyy)){
      return this->effectiveNuclearChargeMd/this->GetEffectivePrincipalQuantumNumber(shellType);
   }   
   else{
      cout << this->errorMessageOrbitalExponent;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      cout << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
      exit(EXIT_FAILURE);
   }   
}


// Part of Eq. (13) in [BZ_1979]
double Atom::GetJss(){
   return this->zindoF0ss;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetJsp(){
   // F0ss = F0sp
   return this->zindoF0ss - this->zindoG1sp/6.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetJsd(){
   return this->zindoF0sd - this->zindoG2sd/10.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetJpp(){
   // F0pp = F0ss
   return this->zindoF0ss - 2.0*this->zindoF2pp/25.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetJpd(){
   // F0pd = F0sd
   return this->zindoF0sd - this->zindoG1pd/15.0 - 3.0*this->zindoG3pd/70.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetJdd(){
   return this->zindoF0dd - 2.0*(this->zindoF2dd - this->zindoF4dd)/63.0;
}

// Eq. (13) in [BZ_1979]
double Atom::GetZindoCoreIntegral(OrbitalType orbital, int l, int m, int n){
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->IonPotS 
              - this->GetJss()*(double)(l-1) 
              - this->GetJsp()*(double)m 
              - this->GetJsd()*(double)n;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->IonPotP
              - this->GetJpp()*(double)(m-1) 
              - this->GetJsp()*(double)l
              - this->GetJpd()*(double)n;
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->IonPotD
              - this->GetJdd()*(double)(n-1) 
              - this->GetJsd()*(double)l
              - this->GetJpd()*(double)m;
   }
   else{
      cout << this->errorMessageZindoCoreIntegral;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      exit(EXIT_FAILURE);
   }

   return value;
}

double Atom::GetIndoF2(){
   return this->indoF2;
}

double Atom::GetIndoG1(){
   return this->indoG1;
}



}
#endif





