#ifndef INCLUDED_HATOM
#define INCLUDED_HATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Hatom : public Atom {
public:
   Hatom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
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
}

double Hatom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
   double value = 0.0;

   if(theory == INDO){
      if(orbital == s){
         value = -1.0*this->imuAmuS;
         if(!isGuess){
            value -= 0.5*gamma;
         }
      }
      else{
         stringstream ss;
         ss << this->errorMessageIndoCoreIntegral;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
         throw MolDSException(ss.str());
      }
   }
   else if(theory == ZINDOS){
      if(orbital == s){
         value = this->GetZindoCoreIntegral(orbital, 1, 0, 0);
      }
      else{
         stringstream ss;
         ss << this->errorMessageZindoSCoreIntegral;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
         throw MolDSException(ss.str());
      }
   }
   else if(theory == MNDO){
      value = this->GetMndoCoreIntegral(orbital);
   }
   else if(theory == AM1){
      value = this->GetAm1CoreIntegral(orbital);
   }

   return value;
}

}
#endif
