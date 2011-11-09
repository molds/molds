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
   double GetOrbitalExponent(ShellType shellType, OrbitalType orbitalType, TheoryType theory);  // See (1.73) in J. A. Pople book for CNDO, INDO, and ZINDOS. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); // P82 - 83 in J. A. Pople book for INDO or Eq. (13) in [BZ_1979] for ZINDO/S. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
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
   double GetNddoAlpha(TheoryType theory); // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S for MNDO. Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoDerivedParameterD(TheoryType theory, int dIndex);    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetNddoDerivedParameterRho(TheoryType theory, int rhoIndex);  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetMndoElecEnergyAtom();        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoHeatsFormAtom();         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetNddoGss(TheoryType theory);
   double GetNddoGpp(TheoryType theory);
   double GetNddoGsp(TheoryType theory);
   double GetNddoGpp2(TheoryType theory);
   double GetNddoHsp(TheoryType theory);
   double GetNddoHpp(TheoryType theory);
   double GetNddoParameterK(TheoryType theory, int kIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterL(TheoryType theory, int lIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterM(TheoryType theory, int mIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
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
   double indoF0CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF0CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientP;       // (3.93-3.99) in J. A. Pople book
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
   int zindoL;              // see l of (13) in [BZ_1979]
   int zindoM;              // see m of (13) in [BZ_1979]
   int zindoN;              // see n (13) in [BZ_1979]
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
   double mndoDerivedParameterD[3];    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double mndoDerivedParameterRho[3];  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double mndoElecEnergyAtom;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoHeatsFormAtom;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoGss;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp2;  //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoHsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double am1CoreintegralS; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1CoreintegralP; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1OrbitalExponentS;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1OrbitalExponentP;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1BondingParameterS; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1BondingParameterP; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Alpha;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gss; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gpp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gsp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gpp2; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Hsp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1DerivedParameterD[3];    // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double am1DerivedParameterRho[3];  // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double am1ParameterK[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1ParameterL[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1ParameterM[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double pm3CoreintegralS; // Table II in ref. [S_1989].
   double pm3CoreintegralP; // Table II in ref. [S_1989].
   double pm3OrbitalExponentS;// Table II in ref. [S_1989].
   double pm3OrbitalExponentP;// Table II in ref. [S_1989].
   double pm3BondingParameterS; // Table II in ref. [S_1989].
   double pm3BondingParameterP; // Table II in ref. [S_1989].
   double pm3Alpha;// Table II in ref. [S_1989].
   double pm3Gss; // Table II in ref. [S_1989].
   double pm3Gpp; // Table II in ref. [S_1989].
   double pm3Gsp; // Table II in ref. [S_1989].
   double pm3Gpp2; // Table II in ref. [S_1989].
   double pm3Hsp; // Table II in ref. [S_1989].
   double pm3DerivedParameterD[3];    // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3DerivedParameterRho[3];  // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3ParameterK[4];// Table II in ref. [S_1989].
   double pm3ParameterL[4];// Table II in ref. [S_1989].
   double pm3ParameterM[4];// Table II in ref. [S_1989].
   double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess);
   double GetZindoCoreIntegral(OrbitalType orbital); // Eq. (13) in [BZ_1979]
   double GetMndoCoreIntegral(OrbitalType orbital); 
   double GetAm1CoreIntegral(OrbitalType orbital); 
   double GetPm3CoreIntegral(OrbitalType orbital); 
private:
   void SetMessages();
   string errorMessageImuAmu;
   string errorMessageOrbitalExponent;
   string errorMessageShellType;
   string errorMessageEffectivPrincipalQuantumNumber;
   string errorMessageZindoCoreIntegral;
   string errorMessageMndoCoreIntegral;
   string errorMessageAm1CoreIntegral;
   string errorMessagePm3CoreIntegral;
   string errorMessageGetOrbitalExponentBadTheory;
   string errorMessageTheoryType;
   string errorMessageGetBondingParameterBadTheoryBadOrbital;
   string errorMessageGetNddoAlphaBadTheory;
   string errorMessageGetNddoDerivedParameterDBadTheory;
   string errorMessageGetNddoDerivedParameterDBadDIndex;
   string errorMessageDIndex;
   string errorMessageGetNddoDerivedParameterRhoBadRhoIndex;
   string errorMessageGetNddoDerivedParameterRhoBadTheory;
   string errorMessageRhoIndex;
   string errorMessageGetNddoParameterKBadKIndex;
   string errorMessageGetNddoParameterKBadTheory;
   string errorMessageKIndex;
   string errorMessageGetNddoParameterLBadLIndex;
   string errorMessageGetNddoParameterLBadTheory;
   string errorMessageLIndex;
   string errorMessageGetNddoParameterMBadMIndex;
   string errorMessageGetNddoParameterMBadTheory;
   string errorMessageMIndex;
   string errorMessageGetNddoGssBadTheory;
   string errorMessageGetNddoGppBadTheory;
   string errorMessageGetNddoGspBadTheory;
   string errorMessageGetNddoGpp2BadTheory;
   string errorMessageGetNddoHspBadTheory;
   string errorMessageGetNddoHppBadTheory;
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
   this->errorMessageAm1CoreIntegral = "Error in base_atoms::Atom::GetAm1CoreINtegral: Invalid orbitalType for AM1.\n";
   this->errorMessagePm3CoreIntegral = "Error in base_atoms::Atom::GetPm3CoreINtegral: Invalid orbitalType for PM3.\n";
   this->errorMessageIonPot = "Error in base_atoms::Atom::GetIonPot: Invalid orbitalType.\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageShellType = "\tshell type = ";
   this->errorMessageTheoryType = "\tTheory = ";
   this->errorMessageEffectivPrincipalQuantumNumber = 
      "Error in base::Atom::GetEffectivePrincipalQuantumNumber: invalid shelltype.\n";
   this->errorMessageZindoCoreIntegral = "Error in base_atoms::Atom::GetZindoCoreINtegral: Invalid orbitalType.\n";
   this->errorMessageGetOrbitalExponentBadTheory = "Erro in base_atoms::Atom::GetOrbitalExponent: Bad theory is set.\n";
   this->errorMessageGetBondingParameterBadTheoryBadOrbital = "Error in base_atoms::Atom::GetBondingParameter: Bad Theory of bad orbital is set.\n";
   this->errorMessageGetNddoDerivedParameterDBadDIndex 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad index for parameter D(dIndex). Only 0, 1, and 2 are permitted.\n";
   this->errorMessageGetNddoDerivedParameterDBadTheory
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad theory is set.\n";
   this->errorMessageGetNddoAlphaBadTheory
      = "Error in base_atoms::Atom::GetNddoAlpha: Bad theory is set.\n";
   this->errorMessageDIndex  = "dIndex = ";
   this->errorMessageGetNddoDerivedParameterRhoBadRhoIndex 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterRho: Bad index for parameter rho(rhoIndex). Only 0, 1, and 2 are permitted.\n";
   this->errorMessageRhoIndex = "rhoIndex = ";
   this->errorMessageGetNddoDerivedParameterRhoBadTheory 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterRho: Bad thory is set.\n";
   this->errorMessageGetNddoParameterKBadKIndex 
      = "Error in base_atoms::Atom::GetNddoParameterK: Bad index for parameter K(kIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterKBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterK: Bad theory is set.\n";
   this->errorMessageKIndex  = "kIndex = ";
   this->errorMessageGetNddoParameterLBadLIndex 
      = "Error in base_atoms::Atom::GetNddoParameterL: Bad index for parameter L(lIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterLBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterL: Bad theory is set.\n";
   this->errorMessageLIndex  = "lIndex = ";
   this->errorMessageGetNddoParameterMBadMIndex 
      = "Error in base_atoms::Atom::GetNddoParameterM: Bad index for parameter M(mIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterMBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterM: Bad theory is set.\n";
   this->errorMessageKIndex  = "mIndex = ";
   this->errorMessageGetNddoGssBadTheory 
      = "Error in base_atoms::Atom::GetNddoGss: Bad theory is set.\n";
   this->errorMessageGetNddoGppBadTheory 
      = "Error in base_atoms::Atom::GetNddoGpp Bad theory is set.\n";
   this->errorMessageGetNddoGspBadTheory 
      = "Error in base_atoms::Atom::GetNddoGsp: Bad theory is set.\n";
   this->errorMessageGetNddoGpp2BadTheory 
      = "Error in base_atoms::Atom::GetNddoGpp2: Bad theory is set.\n";
   this->errorMessageGetNddoHspBadTheory 
      = "Error in base_atoms::Atom::GetNddoHsp: Bad theory is set.\n";
   this->errorMessageGetNddoHppBadTheory 
      = "Error in base_atoms::Atom::GetNddoHp: Bad theory is set.\n";
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
   else if(theory == AM1 && orbital == s){
      value = this->am1BondingParameterS;
   }
   else if(theory == AM1 && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->am1BondingParameterP;
   }
   else if(theory == PM3 && orbital == s){
      value = this->pm3BondingParameterS;
   }
   else if(theory == PM3 && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->pm3BondingParameterP;
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
   else if(theory == AM1){
      if(orbitalType == s){ 
         return this->am1OrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->am1OrbitalExponentP;
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
   else if(theory == PM3){
      if(orbitalType == s){ 
         return this->pm3OrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->pm3OrbitalExponentP;
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

// (3.93) - (3.99) in J. A. Pople book.
double Atom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess){
   double value = 0.0;
   if(orbital == s){
      value = -1.0*this->imuAmuS;
      if(!isGuess){
         value -= this->indoF0CoefficientS*gamma 
                 +this->indoG1CoefficientS*this->indoG1
                 +this->indoF2CoefficientS*this->indoF2;
      }
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->imuAmuP;
      if(!isGuess){
         value -= this->indoF0CoefficientP*gamma 
                 +this->indoG1CoefficientP*this->indoG1
                 +this->indoF2CoefficientP*this->indoF2;
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageIndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}


// Eq. (13) in [BZ_1979]
double Atom::GetZindoCoreIntegral(OrbitalType orbital){
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->ionPotS 
              - this->GetJss()*(double)(this->zindoL-1) 
              - this->GetJsp()*(double)this->zindoM
              - this->GetJsd()*(double)this->zindoN;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->ionPotP
              - this->GetJpp()*(double)(this->zindoM-1) 
              - this->GetJsp()*(double)this->zindoL
              - this->GetJpd()*(double)this->zindoN;
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->ionPotD
              - this->GetJdd()*(double)(this->zindoN-1) 
              - this->GetJsd()*(double)this->zindoL
              - this->GetJpd()*(double)this->zindoM;
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

double Atom::GetAm1CoreIntegral(OrbitalType orbital){
   double value=0.0;

   if(orbital == s){
      value = this->am1CoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->am1CoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageAm1CoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetPm3CoreIntegral(OrbitalType orbital){
   double value=0.0;

   if(orbital == s){
      value = this->pm3CoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->pm3CoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessagePm3CoreIntegral;
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

double Atom::GetNddoAlpha(TheoryType theory){
   double value = 0.0;
   if(theory == MNDO){
      value = this->mndoAlpha;
   }
   else if(theory == AM1){
      value = this->am1Alpha;
   }
   else if(theory == PM3){
      value = this->pm3Alpha;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoAlphaBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
   return value;
}

double Atom::GetNddoDerivedParameterD(TheoryType theory, int dIndex){
   if(dIndex == 0 || dIndex == 1 || dIndex == 2){
      if(theory == MNDO){
         return this->mndoDerivedParameterD[dIndex];
      }
      else if(theory == AM1){
         return this->am1DerivedParameterD[dIndex];
      }
      else if(theory == PM3){
         return this->pm3DerivedParameterD[dIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoDerivedParameterDBadDIndex;
      ss << this->errorMessageDIndex << dIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoDerivedParameterRho(TheoryType theory, int rhoIndex){
   if(rhoIndex == 0 || rhoIndex == 1 || rhoIndex == 2){
      if(theory == MNDO){
         return this->mndoDerivedParameterRho[rhoIndex];
      }
      else if(theory == AM1){
         return this->am1DerivedParameterRho[rhoIndex];
      }
      else if(theory == PM3){
         return this->pm3DerivedParameterRho[rhoIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoDerivedParameterRhoBadRhoIndex;
      ss << this->errorMessageRhoIndex << rhoIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterK(TheoryType theory, int kIndex){
   if(kIndex == 0 || kIndex == 1 || kIndex == 2 || kIndex == 3){
      if(theory == AM1){
         return this->am1ParameterK[kIndex];
      }
      else if(theory == PM3){
         return this->pm3ParameterK[kIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterKBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterKBadKIndex;
      ss << this->errorMessageKIndex << kIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterL(TheoryType theory, int lIndex){
   if(lIndex == 0 || lIndex == 1 || lIndex == 2 || lIndex == 3){
      if(theory == AM1){
         return this->am1ParameterL[lIndex];
      }
      else if(theory == PM3){
         return this->pm3ParameterL[lIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterLBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterLBadLIndex;
      ss << this->errorMessageLIndex << lIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterM(TheoryType theory, int mIndex){
   if(mIndex == 0 || mIndex == 1 || mIndex == 2 || mIndex == 3){
      if(theory == AM1){
         return this->am1ParameterM[mIndex];
      }
      else if(theory == PM3){
         return this->pm3ParameterM[mIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterMBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterMBadMIndex;
      ss << this->errorMessageMIndex << mIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetMndoElecEnergyAtom(){
   return this->mndoElecEnergyAtom;
}

double Atom::GetMndoHeatsFormAtom(){
   return this->mndoHeatsFormAtom;
}

double Atom::GetNddoGss(TheoryType theory){
   if(theory == MNDO){
      return this->mndoGss;
   }
   else if(theory == AM1){
      return this->am1Gss;
   }
   else if(theory == PM3){
      return this->pm3Gss;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGssBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGpp(TheoryType theory){
   if(theory == MNDO){
      return this->mndoGpp;
   }
   else if(theory == AM1){
      return this->am1Gpp;
   }
   else if(theory == PM3){
      return this->pm3Gpp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGppBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGsp(TheoryType theory){
   if(theory == MNDO){
      return this->mndoGsp;
   }
   else if(theory == AM1){
      return this->am1Gsp;
   }
   else if(theory == PM3){
      return this->pm3Gsp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGspBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGpp2(TheoryType theory){
   if(theory == MNDO){
      return this->mndoGpp2;
   }
   else if(theory == AM1){
      return this->am1Gpp2;
   }
   else if(theory == PM3){
      return this->pm3Gpp2;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGpp2BadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoHsp(TheoryType theory){
   if(theory == MNDO){
      return this->mndoHsp;
   }
   else if(theory == AM1){
      return this->am1Hsp;
   }
   else if(theory == PM3){
      return this->pm3Hsp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoHspBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

// see p17 in [MOPAC_1990]
double Atom::GetNddoHpp(TheoryType theory){
   if(theory == MNDO){
      return 0.5*(this->mndoGpp - this->mndoGpp2);
   }
   else if(theory == AM1){
      return 0.5*(this->am1Gpp - this->am1Gpp2);
   }
   else if(theory == PM3){
      return 0.5*(this->pm3Gpp - this->pm3Gpp2);
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoHppBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
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

double Atom::GetCoreIntegral(OrbitalType orbital, 
                             double gamma, 
                             bool isGuess, 
                             TheoryType theory){
   double value = 0.0;
   if(theory == INDO){
      value = this->GetIndoCoreIntegral(orbital, gamma, isGuess);
   }
   else if(theory == ZINDOS){
      value = this->GetZindoCoreIntegral(orbital);
   }
   else if(theory == MNDO){
      value = this->GetMndoCoreIntegral(orbital);
   }
   else if(theory == AM1){
      value = this->GetAm1CoreIntegral(orbital);
   }
   else if(theory == PM3){
      value = this->GetPm3CoreIntegral(orbital);
   }
   return value;
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




















