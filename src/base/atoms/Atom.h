#ifndef INCLUDED_ATOM
#define INCLUDED_ATOM
namespace MolDS_base_atoms{

class Atom{
public:
   Atom();
   Atom(double x, double y, double z);
   ~Atom();
   MolDS_base::AtomType GetAtomType();
   double GetAtomicMass();
   double* GetXyz();
   void SetXyz(double x, double y, double z);
   double* GetPxyz();
   void SetPxyz(double px, double py, double pz);
   std::vector<MolDS_base::OrbitalType> GetValence();
   double GetAtomicBasisValue(double x, 
                              double y, 
                              double z, 
                              int valenceIndex,
                              MolDS_base::TheoryType theory);
   double GetBondingParameter();
   double GetBondingParameter(MolDS_base::TheoryType theory, 
                              MolDS_base::OrbitalType orbital);
   double GetCoreCharge();
   int GetFirstAOIndex();
   void SetFirstAOIndex(int firstAOIndex);
   MolDS_base::ShellType GetValenceShellType();
   int GetNumberValenceElectrons();
   double GetOrbitalExponent(MolDS_base::ShellType shellType, 
                             MolDS_base::OrbitalType orbitalType, 
                             MolDS_base::TheoryType theory);  // See (1.73) in J. A. Pople book for CNDO, INDO, and ZINDOS. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
   double GetCoreIntegral(MolDS_base::OrbitalType orbital, 
                          double gamma, 
                          bool isGuess, 
                          MolDS_base::TheoryType theory); // P82 - 83 in J. A. Pople book for INDO or Eq. (13) in [BZ_1979] for ZINDO/S. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
   double GetCoreIntegral(MolDS_base::OrbitalType orbital, 
                          bool isGuess, 
                          MolDS_base::TheoryType theory);
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
   double GetZindoIonPot(MolDS_base::OrbitalType orbital);
   double GetNddoAlpha(MolDS_base::TheoryType theory); // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S for MNDO. Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoDerivedParameterD(MolDS_base::TheoryType theory, int dIndex);    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetNddoDerivedParameterRho(MolDS_base::TheoryType theory, int rhoIndex);  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetMndoElecEnergyAtom();        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoHeatsFormAtom();         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetNddoGss(MolDS_base::TheoryType theory);
   double GetNddoGpp(MolDS_base::TheoryType theory);
   double GetNddoGsp(MolDS_base::TheoryType theory);
   double GetNddoGpp2(MolDS_base::TheoryType theory);
   double GetNddoHsp(MolDS_base::TheoryType theory);
   double GetNddoHpp(MolDS_base::TheoryType theory);
   double GetNddoParameterK(MolDS_base::TheoryType theory, int kIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterL(MolDS_base::TheoryType theory, int lIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterM(MolDS_base::TheoryType theory, int mIndex);//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
protected:
   double* xyz; // coordinates
   double* pxyz; // momentum. Note that this is not velocity!! 
   MolDS_base::AtomType atomType;
   double atomicMass;  // Appendix 1 in [I_1998]
   std::vector<MolDS_base::OrbitalType> valence;
   MolDS_base::ShellType valenceShellType;
   int firstAOIndex;
   int numberValenceElectrons;
   double imuAmuS;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuP;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuD;                      // Table 3.4 or 3.5 in J. A. Pople book
   double bondingParameter;             // Table 3.2 and 3.4 in J. A. Pople book
   double coreCharge;                   // = Z_A in J. A. Pople book.
   double effectiveNuclearChargeK;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeL;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeMsp;    // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeMd;     // Table 1.5 in J. A. Pople book
   double indoF2;                   // Table 3.6 in J. A. Pople book
   double indoG1;                   // Table 3.6 in J. A. Pople book
   double indoF0CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF0CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double zindoBondingParameterS;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double zindoBondingParameterD;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double zindoF0ss;                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double zindoF0sd;        // Table 1 in [AEZ_1986]
   double zindoF0dd;        // Table 1 in [AEZ_1986]
   double zindoG1sp;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pp;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG2sd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG1pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG3pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2dd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF4dd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   int zindoL;              // see l of (13) in [BZ_1979]
   int zindoM;              // see m of (13) in [BZ_1979]
   int zindoN;              // see n (13) in [BZ_1979]
   double zindoIonPotS;   // Ionization potential, Table 4 in [BZ_1979]
   double zindoIonPotP;   // Ionization potential, Table 4 in [BZ_1979]
   double zindoIonPotD;   // Ionization potential, Table 4 in [BZ_1979]
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
   double pm3DerivedParameterD[3];    // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3DerivedParameterRho[3];  // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3ParameterK[4];// Table II in ref. [S_1989].
   double pm3ParameterL[4];// Table II in ref. [S_1989].
   double pm3ParameterM[4];// Table II in ref. [S_1989].
   double pm3Gss; // Table II in ref. [S_1989].
   double pm3Gpp; // Table II in ref. [S_1989].
   double pm3Gsp; // Table II in ref. [S_1989].
   double pm3Gpp2; // Table II in ref. [S_1989].
   double pm3Hsp; // Table II in ref. [S_1989].
private:
   std::string errorMessageIonPot;
   std::string errorMessageAtomType;
   std::string errorMessageNumberValences;
   std::string errorMessageValenceIndex;
   std::string errorMessageOrbitalType;
   std::string errorMessageOrbitalExponent;
   std::string errorMessageShellType;
   std::string errorMessageEffectivPrincipalQuantumNumber;
   std::string errorMessageCndo2CoreIntegral;
   std::string errorMessageIndoCoreIntegral;
   std::string errorMessageZindoCoreIntegral;
   std::string errorMessageMndoCoreIntegral;
   std::string errorMessageAm1CoreIntegral;
   std::string errorMessagePm3CoreIntegral;
   std::string errorMessageGetAtomicBasisValueBadValenceIndex;
   std::string errorMessageGetRealAnuglarPartAOBadValence;
   std::string errorMessageGetOrbitalExponentBadTheory;
   std::string errorMessageTheoryType;
   std::string errorMessageGetBondingParameterBadTheoryBadOrbital;
   std::string errorMessageGetNddoAlphaBadTheory;
   std::string errorMessageGetNddoDerivedParameterDBadTheory;
   std::string errorMessageGetNddoDerivedParameterDBadDIndex;
   std::string errorMessageDIndex;
   std::string errorMessageGetNddoDerivedParameterRhoBadRhoIndex;
   std::string errorMessageGetNddoDerivedParameterRhoBadTheory;
   std::string errorMessageRhoIndex;
   std::string errorMessageGetNddoParameterKBadKIndex;
   std::string errorMessageGetNddoParameterKBadTheory;
   std::string errorMessageKIndex;
   std::string errorMessageGetNddoParameterLBadLIndex;
   std::string errorMessageGetNddoParameterLBadTheory;
   std::string errorMessageLIndex;
   std::string errorMessageGetNddoParameterMBadMIndex;
   std::string errorMessageGetNddoParameterMBadTheory;
   std::string errorMessageMIndex;
   std::string errorMessageGetNddoGssBadTheory;
   std::string errorMessageGetNddoGppBadTheory;
   std::string errorMessageGetNddoGspBadTheory;
   std::string errorMessageGetNddoGpp2BadTheory;
   std::string errorMessageGetNddoHspBadTheory;
   std::string errorMessageGetNddoHppBadTheory;
   void SetMessages();
   double GetRealAnuglarPartAO(double theta, double phi, MolDS_base::OrbitalType orbital);
   double GetRadialPartAO(double dr, double orbitalExponent, MolDS_base::ShellType shell);
   int GetEffectivePrincipalQuantumNumber(MolDS_base::ShellType shellType); // Table 1.4 in J. A. Pople book
   double GetZindoJss();  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJsp();  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJsd();  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJpp();  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJpd();  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJdd();  // Part of Eq. (13) in [BZ_1979]
   double GetCndo2CoreIntegral(MolDS_base::OrbitalType orbital, double gamma, bool isGuess);
   double GetIndoCoreIntegral(MolDS_base::OrbitalType orbital, double gamma, bool isGuess);
   double GetZindoCoreIntegral(MolDS_base::OrbitalType orbital); // Eq. (13) in [BZ_1979]
   double GetMndoCoreIntegral(MolDS_base::OrbitalType orbital); 
   double GetAm1CoreIntegral(MolDS_base::OrbitalType orbital); 
   double GetPm3CoreIntegral(MolDS_base::OrbitalType orbital); 
};
}
#endif

