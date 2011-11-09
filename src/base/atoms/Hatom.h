#ifndef INCLUDED_HATOM
#define INCLUDED_HATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Hatom : public Atom {
public:
   Hatom(double x, double y, double z);
private:
};

Hatom::Hatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = H;
   this->atomicMass = 1.00794*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -12.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo = 0.0;
   this->coreCharge = 1.0;
   this->imuAmuS = 7.176*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 0.0;
   this->imuAmuD = 0.0;
   this->valenceShellType = k;
   this->effectiveNuclearChargeK = 1.2; // see P78 in J. A. Pople book
   this->effectiveNuclearChargeL = 0.0;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->numberValenceElectrons = 1;
   this->indoG1 = 0.0;
   this->indoF2 = 0.0;
   this->indoF0CoefficientS = 0.5;
   this->indoF0CoefficientP = 0.0;
   this->indoG1CoefficientS = 0.0;
   this->indoG1CoefficientP = 0.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = 0.0;
   this->zindoF0ss = 12.85 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 0.0;                 
   this->zindoF2pp = 0.0;                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 1;
   this->zindoM = 0;
   this->zindoN = 0;
   this->ionPotS = 13.06 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -11.906276 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = 0.0;         
   this->mndoOrbitalExponentS = 1.331967;      
   this->mndoOrbitalExponentP = 0.0;      
   this->mndoBondingParameterS = -6.989064 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = 0.0;     
   this->mndoAlpha = 2.544134 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -11.906276 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 52.102 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss = 12.848 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 0.0 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] = 0.0;    
   this->mndoDerivedParameterD[1] = 0.0;    
   this->mndoDerivedParameterD[2] = 0.0;    
   this->mndoDerivedParameterRho[0] = 0.5/0.4721515374;
   //this->mndoDerivedParameterRho[0] = 0.560345 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->mndoDerivedParameterRho[1] = 0.0;  
   this->mndoDerivedParameterRho[2] = 0.0;  
   this->am1CoreintegralS = -11.396427 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = 0.0;         
   this->am1OrbitalExponentS = 1.188078;      
   this->am1OrbitalExponentP = 0.0;      
   this->am1BondingParameterS = -6.173787 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = 0.0;     
   this->am1Alpha = 2.882324 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.0;    
   this->am1DerivedParameterD[2] = 0.0;    
   this->am1DerivedParameterRho[0] = 0.5/0.4721515374;
   this->am1DerivedParameterRho[1] = 0.0;  
   this->am1DerivedParameterRho[2] = 0.0;  
   this->am1ParameterK[0] = 0.122796 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.005090 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.018336 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.000000 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 2.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.20 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.80 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.10 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3CoreintegralS = -13.073321 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = 0.0;
   this->pm3OrbitalExponentS = 0.967807;      
   this->pm3OrbitalExponentP = 0.0;      
   this->pm3BondingParameterS = -5.626512 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = 0.0;
   this->pm3Alpha = 3.356386 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3Gss = 14.794208 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 0.0;
   this->pm3Gsp = 0.0;
   this->pm3Gpp2 = 0.0;   
   this->pm3Hsp = 0.0;    
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.0;    
   this->pm3DerivedParameterD[2] = 0.0;    
   this->pm3DerivedParameterRho[0] = 0.5/0.5436727936;
   this->pm3DerivedParameterRho[1] = 0.0;  
   this->pm3DerivedParameterRho[2] = 0.0;  
   this->pm3ParameterK[0] = 1.128750 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] =-1.060329 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 5.096282 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.003788 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.537465 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.570189 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
}
}
#endif
