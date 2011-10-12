#ifndef INCLUDED_ATOM
#define INCLUDED_ATOM

#include<stdio.h>
#include<stdlib.h>
#include<iostream>

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Atom{
public:
   Atom();
   Atom(double x, double y, double z);
   ~Atom();
   AtomType GetAtomType();
   double GetAtomicMass();
   double* GetXyz();
   void SetXyz(double x, double y, double z);
   double* GetPxyz();
   void SetPxyz(double px, double py, double pz);
   vector<OrbitalType> GetValence();
   double GetBondingParameter();
   double GetBondingParameter(TheoryType theory, OrbitalType orbital);
   double GetCoreCharge();
   int GetFirstAOIndex();
   void SetFirstAOIndex(int firstAOIndex);
   ShellType GetValenceShellType();
   int GetNumberValenceElectrons();
   double GetImuAmu(OrbitalType orbitalType);  // return 0.5*(I_mu + A_mu)
   double GetOrbitalExponent(ShellType shellType, OrbitalType orbitalType, TheoryType theory);  // See (1.73) in J. A. Pople book for CNDO, INDO, and ZINDOS. See [BT_1977] for MNDO.
   virtual double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory) = 0; // P82 - 83 in J. A. Pople book for INDO or Eq. (13) in [BZ_1979] for ZINDO/S
   double GetCoreIntegral(OrbitalType orbital, bool isGuess, TheoryType theory);
   double GetIndoF2();
   double GetIndoG1();
   double GetZindoF0ss();                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double GetZindoF0sd();                  // Table 1 in [AEZ_1986]
   double GetZindoF0dd();                  // Table 1 in [AEZ_1986]
   double GetZindoG1sp();                 // Table 3 in ref. [BZ_1979]
   double GetZindoF2pp();                 // Table 3 in ref. [BZ_1979]
   double GetZindoG2sd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoG1pd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoF2pd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoG3pd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoF2dd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoF4dd();                 // Table 3 in ref. [BZ_1979]
   double GetZindoF0ssLower();                 // Apendix in ref. [BZ_1979] 
   double GetZindoF0sdLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoF0ddLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoG1spLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoF2ppLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoG2sdLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoG1pdLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoF2pdLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoG3pdLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoF2ddLower();                 // Apendix in ref. [BZ_1979]
   double GetZindoF4ddLower();                 // Apendix in ref. [BZ_1979]
   double GetIonPot(OrbitalType orbital);
   double GetMndoAlpha(); // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoDerivedParameterD(int dIndex);    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoDerivedParameterRho(int rhoIndex);  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoElecEnergyAtom();        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoHeatsFormAtom();         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoGss();
   double GetMndoGpp();
   double GetMndoGsp();
   double GetMndoGpp2();
   double GetMndoHsp();
   double GetMndoHpp();
protected:
   string errorMessageIndoCoreIntegral;
   string errorMessageZindoSCoreIntegral;
   string errorMessageIonPot;
   string errorMessageAtomType;
   string errorMessageOrbitalType;
   double* xyz; // coordinates
   double* pxyz; // momentum. Note that this is not velocity!! 
   AtomType atomType;
   double atomicMass;  // Appendix 1 in [I_1998]
   vector<OrbitalType> valence;
   ShellType valenceShellType;
   int firstAOIndex;
   int numberValenceElectrons;
   double imuAmuS;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuP;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuD;                      // Table 3.4 or 3.5 in J. A. Pople book
   double bondingParameter;             // Table 3.2 and 3.4 in J. A. Pople book
   double bondingParameterSZindo;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double bondingParameterDZindo;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double coreCharge;                   // = Z_A in J. A. Pople book.
   double effectiveNuclearChargeK;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeL;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeMsp;    // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeMd;     // Table 1.5 in J. A. Pople book
   int GetEffectivePrincipalQuantumNumber(ShellType shellType); // Table 1.4 in J. A. Pople book
   double indoF2;                   // Table 3.6 in J. A. Pople book
   double indoG1;                   // Table 3.6 in J. A. Pople book
   double zindoF0ss;                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double zindoF0sd;                  // Table 1 in [AEZ_1986]
   double zindoF0dd;                  // Table 1 in [AEZ_1986]
   double zindoG1sp;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pp;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG2sd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG1pd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG3pd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2dd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF4dd;                 // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double ionPotS;   // Ionization potential, Table 4 in [BZ_1979]
   double ionPotP;   // Ionization potential, Table 4 in [BZ_1979]
   double ionPotD;   // Ionization potential, Table 4 in [BZ_1979]
   double mndoCoreintegralS;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoCoreintegralP;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. 
   double mndoOrbitalExponentS;      // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoOrbitalExponentP;      // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoBondingParameterS;     // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoBondingParameterP;     // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoAlpha;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoDerivedParameterD[2];    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoDerivedParameterRho[3];  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoElecEnergyAtom;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoHeatsFormAtom;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoGss;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp2;  //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoHsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double GetZindoCoreIntegral(OrbitalType orbital, int l, int m, int n); // Eq. (13) in [BZ_1979]
   double GetMndoCoreIntegral(OrbitalType orbital); // Eq. (13) in [BZ_1979]
private:
   void SetMessages();
   string errorMessageImuAmu;
   string errorMessageOrbitalExponent;
   string errorMessageShellType;
   string errorMessageEffectivPrincipalQuantumNumber;
   string errorMessageZindoCoreIntegral;
   string errorMessageMndoCoreIntegral;
   string errorMessageGetOrbitalExponentBadTheory;
   string errorMessageTheoryType;
   string errorMessageGetBondingParameterBadTheoryBadOrbital;
   string errorMessageGetDerivedParameterDBadDIndex;
   string errorMessageDIndex;
   string errorMessageGetDerivedParameterRhoBadRhoIndex;
   string errorMessageRhoIndex;
   double GetJss();  // Part of Eq. (13) in [BZ_1979]
   double GetJsp();  // Part of Eq. (13) in [BZ_1979]
   double GetJsd();  // Part of Eq. (13) in [BZ_1979]
   double GetJpp();  // Part of Eq. (13) in [BZ_1979]
   double GetJpd();  // Part of Eq. (13) in [BZ_1979]
   double GetJdd();  // Part of Eq. (13) in [BZ_1979]

};

Atom::Atom(){
   this->xyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->pxyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->SetMessages();
}

Atom::Atom(double x, double y, double z){
   this->xyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->pxyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->SetMessages();
   this->SetXyz(x, y, z);
}

Atom::~Atom(){
   if(this->xyz != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->xyz);
      //cout << "xyz deleted\n";
   }
   if(this->pxyz != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(&this->pxyz);
      //cout << "pxyz (momenta) deleted\n";
   }
   //cout << "atom deleted\n";
}

void Atom::SetMessages(){
   this->errorMessageImuAmu = "Error in base_atoms::Atom::GetImuAmu: Invalid orbitalType.\n";
   this->errorMessageOrbitalExponent = "Error in base_atoms::Atom::GetOrbitalExponent: Invalid shelltype or orbitalType.\n";
   this->errorMessageIndoCoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for INDO.\n";
   this->errorMessageZindoSCoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for ZINDO/S.\n";
   this->errorMessageMndoCoreIntegral = "Error in base_atoms::Atom::GetMndoCoreINtegral: Invalid orbitalType for MNDO.\n";
   this->errorMessageIonPot = "Error in base_atoms::Atom::GetIonPot: Invalid orbitalType.\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageShellType = "\tshell type = ";
   this->errorMessageEffectivPrincipalQuantumNumber = 
      "Error in base::Atom::GetEffectivePrincipalQuantumNumber: invalid shelltype.\n";
   this->errorMessageZindoCoreIntegral = "Error in base_atoms::Atom::GetZindoCoreINtegral: Invalid orbitalType.\n";
   this->errorMessageGetOrbitalExponentBadTheory = "Erro in base_atoms::Atom::GetOrbitalExponent: Bad theory is set.\n";
   this->errorMessageTheoryType = "Theory = ";
   this->errorMessageGetBondingParameterBadTheoryBadOrbital = "Error in base_atoms::Atom::GetBondingParameter: Bad Theory of bad orbital is set.\n";
   this->errorMessageGetDerivedParameterDBadDIndex = "Error in base_atoms::Atom::GetDerivedParameterD: Bad index for parameter D(dIndex). Only 0 and 1 are permitted.\n";
   this->errorMessageDIndex  = "dIndex = ";
   this->errorMessageGetDerivedParameterRhoBadRhoIndex = "Error in base_atoms::Atom::GetDerivedParameterRho: Bad index for parameter rho(rhoIndex). Only 0, 1, and 2 are permitted.\n";
   this->errorMessageRhoIndex = "rhoIndex = ";
}

AtomType Atom::GetAtomType(){
   return this->atomType;
}

double Atom::GetAtomicMass(){
   return this->atomicMass;
}

double* Atom::GetXyz(){
   return this->xyz;
}

double* Atom::GetPxyz(){
   return this->pxyz;
}

void Atom::SetXyz(double x, double y, double z){
   xyz[0]= x;
   xyz[1]= y;
   xyz[2]= z;
}

void Atom::SetPxyz(double px, double py, double pz){
   pxyz[0]= px;
   pxyz[1]= py;
   pxyz[2]= pz;
}

vector<OrbitalType> Atom::GetValence(){
   return this->valence;
}

double Atom::GetBondingParameter(TheoryType theory, OrbitalType orbital){

   double value = 0.0;
   if(theory == CNDO2 || theory == INDO){
      value = this->bondingParameter;
   }     
   else if(theory == ZINDOS && ( orbital == s ||
                                 orbital == px ||
                                 orbital == py ||
                                 orbital == pz ) ){
      value = this->bondingParameterSZindo;
   }
   else if(theory == ZINDOS && ( orbital == dxy ||
                                 orbital == dyz ||
                                 orbital == dzz ||
                                 orbital == dzx ||
                                 orbital == dxxyy ) ){
      value = this->bondingParameterDZindo;
   }
   else if(theory == MNDO && orbital == s){
      value = this->mndoBondingParameterS;
   }
   else if(theory == MNDO && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->mndoBondingParameterP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetBondingParameterBadTheoryBadOrbital;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << "\n";
      throw MolDSException(ss.str());
   }

   return value;

}

double Atom::GetBondingParameter(){
   return this->GetBondingParameter(CNDO2, s);
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
      stringstream ss;
      ss << this->errorMessageEffectivPrincipalQuantumNumber;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      throw MolDSException(ss.str());
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
      stringstream ss;
      ss << errorMessageImuAmu;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
      throw MolDSException(ss.str());
   }   
}

// (1.73) in J. A. Pople book
double Atom::GetOrbitalExponent(ShellType shellType, OrbitalType orbitalType, TheoryType theory){
   if(theory == CNDO2 || theory == INDO || theory == ZINDOS){
      if(shellType == k && orbitalType == s){ 
         return this->effectiveNuclearChargeK
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == l && (orbitalType == s  || 
                                 orbitalType == px || 
                                 orbitalType == py || 
                                 orbitalType == pz)){
         return this->effectiveNuclearChargeL
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == m && (orbitalType == s  || 
                                 orbitalType == px || 
                                 orbitalType == py || 
                                 orbitalType == pz )){
         return this->effectiveNuclearChargeMsp
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == m && (orbitalType == dxy  || 
                                 orbitalType == dyz ||
                                 orbitalType == dzz ||
                                 orbitalType == dzx ||
                                 orbitalType == dxxyy)){
         return this->effectiveNuclearChargeMd
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }   
   }
   else if(theory == MNDO){
      if(orbitalType == s){ 
         return this->mndoOrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->mndoOrbitalExponentP;
      }
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetOrbitalExponentBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
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
   return this->zindoF0dd - 2.0*(this->zindoF2dd + this->zindoF4dd)/63.0;
}

// Eq. (13) in [BZ_1979]
double Atom::GetZindoCoreIntegral(OrbitalType orbital, int l, int m, int n){
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->ionPotS 
              - this->GetJss()*(double)(l-1) 
              - this->GetJsp()*(double)m 
              - this->GetJsd()*(double)n;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->ionPotP
              - this->GetJpp()*(double)(m-1) 
              - this->GetJsp()*(double)l
              - this->GetJpd()*(double)n;
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->ionPotD
              - this->GetJdd()*(double)(n-1) 
              - this->GetJsd()*(double)l
              - this->GetJpd()*(double)m;
   }
   else{
      stringstream ss;
      ss << this->errorMessageZindoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetMndoCoreIntegral(OrbitalType orbital){
   double value=0.0;

   if(orbital == s){
      value = this->mndoCoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->mndoCoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageMndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetIndoF2(){
   return this->indoF2;
}

double Atom::GetIndoG1(){
   return this->indoG1;
}

double Atom::GetMndoAlpha(){
   return this->mndoAlpha;
}

double Atom::GetMndoDerivedParameterD(int dIndex){
   if(dIndex == 0 || dIndex == 1){
      return this->mndoDerivedParameterD[dIndex];
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetDerivedParameterDBadDIndex;
      ss << this->errorMessageDIndex << dIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetMndoDerivedParameterRho(int rhoIndex){
   if(rhoIndex == 0 || rhoIndex == 1 || rhoIndex == 2){
      return this->mndoDerivedParameterRho[rhoIndex];
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetDerivedParameterRhoBadRhoIndex;
      ss << this->errorMessageRhoIndex << rhoIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetMndoElecEnergyAtom(){
   return this->mndoElecEnergyAtom;
}

double Atom::GetMndoHeatsFormAtom(){
   return this->mndoHeatsFormAtom;
}

double Atom::GetMndoGss(){
   return this->mndoGss;
}

double Atom::GetMndoGpp(){
   return this->mndoGpp;
}

double Atom::GetMndoGsp(){
   return this->mndoGsp;
}

double Atom::GetMndoGpp2(){
   return this->mndoGpp2;
}

double Atom::GetMndoHsp(){
   return this->mndoHsp;
}

// see p17 in [MOPAC_1990]
double Atom::GetMndoHpp(){
   return 0.5*(this->mndoGpp - this->mndoGpp2);
}

// Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
double Atom::GetZindoF0ss(){
   return this->zindoF0ss;
}

// Table 1 in [AEZ_1986]
double Atom::GetZindoF0sd(){
   return this->zindoF0sd;
}


// Table 1 in [AEZ_1986]
double Atom::GetZindoF0dd(){
   return this->zindoF0dd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG1sp(){
   return this->zindoG1sp;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2pp(){
   return this->zindoF2pp;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG2sd(){
   return this->zindoG2sd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG1pd(){
   return this->zindoG1pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2pd(){
   return this->zindoF2pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG3pd(){
   return this->zindoG3pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2dd(){
   return this->zindoF2dd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF4dd(){
   return this->zindoF4dd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0ssLower(){
   return this->zindoF0ss;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0sdLower(){
   return this->zindoF0sd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0ddLower(){
   return this->zindoF0dd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG1spLower(){
   return this->zindoG1sp/3.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2ppLower(){
   return this->zindoF2pp/25.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG2sdLower(){
   return this->zindoG2sd/5.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG1pdLower(){
   return this->zindoG1pd/15.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2pdLower(){
   return this->zindoF2pd/35.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG3pdLower(){
   return this->zindoG3pd/245.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2ddLower(){
   return this->zindoF2dd/49.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF4ddLower(){
   return this->zindoF4dd/441.0;
}


double Atom::GetCoreIntegral(OrbitalType orbital, bool isGuess, TheoryType theory){
   return this->GetCoreIntegral(orbital, 0.0, isGuess, theory);
}

double Atom::GetIonPot(OrbitalType orbital){
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->ionPotS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->ionPotP;
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->ionPotD;
   }
   else{
      stringstream ss;
      ss << this->errorMessageIonPot;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}


}
#endif




















