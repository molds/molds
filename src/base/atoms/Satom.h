#ifndef INCLUDED_SATOM
#define INCLUDED_SATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Satom : public Atom {
public:
   Satom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Satom::Satom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = S;
   this->atomicMass = 32.066*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2 || 
      Parameters::GetInstance()->GetCurrentTheory() == INDO){
      this->valence.push_back(dxy);
      this->valence.push_back(dyz);
      this->valence.push_back(dzz);
      this->valence.push_back(dzx);
      this->valence.push_back(dxxyy);
   }
   this->bondingParameter = -18.150*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -14.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo =   4.0*Parameters::GetInstance()->GetEV2AU();
   this->coreCharge = 6.0;
   this->imuAmuS = 17.650*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 6.989*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.713*Parameters::GetInstance()->GetEV2AU();
   this->valenceShellType = m;
   this->effectiveNuclearChargeK = 15.70;
   this->effectiveNuclearChargeL = 11.85;
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->effectiveNuclearChargeMsp = 1.925*3.0;
      this->effectiveNuclearChargeMd = 1.731*3.0;
   }
   else{
      this->effectiveNuclearChargeMsp = 5.45;
      this->effectiveNuclearChargeMd = 5.45;
   }
   this->numberValenceElectrons = 6;
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   // the zindoF0ss for sulfer atoms are set to be equal 
   // to the one (10.09eV) in "ORCA 2.8"( http://www.thch.uni-bonn.de/tc/orca/ ).
   this->zindoF0ss = 10.09 * Parameters::GetInstance()->GetEV2AU(); 
   //this->zindoF0ss = 8.96 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 3.10 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 4.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 3.25 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG1pd = 4.31 * Parameters::GetInstance()->GetEV2AU();        
   this->zindoF2pd = 3.45 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG3pd = 2.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF2dd = 3.55 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF4dd = 2.31 * Parameters::GetInstance()->GetEV2AU();
   //this->zindoG1sp = 24807*Parameters::GetInstance()->GetKayser2AU();                 
   //this->zindoF2pp = 36600*Parameters::GetInstance()->GetKayser2AU();                 
   //this->zindoG2sd = 25972*Parameters::GetInstance()->GetKayser2AU();     
   //this->zindoG1pd = 34486*Parameters::GetInstance()->GetKayser2AU();        
   //this->zindoF2pd = 29173*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoG3pd = 20587*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoF2dd = 28411*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoF4dd = 18529*Parameters::GetInstance()->GetKayser2AU();           
   this->ionPotS = 21.11 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 12.39 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 4.11 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -72.242281 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -56.973207 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.312962;      
   this->mndoOrbitalExponentP = 2.009146;      
   this->mndoBondingParameterS = -10.761670 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -10.108433 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.478026 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -226.01239 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 66.40 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  12.88 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =   9.90 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  11.26 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  8.83 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.26 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.9189935;
   this->mndoDerivedParameterD[2] =   0.8328514;
   this->mndoDerivedParameterRho[0] = 0.5/0.4733554;
   this->mndoDerivedParameterRho[1] = 0.5/0.5544502;
   this->mndoDerivedParameterRho[2] = 0.5/0.5585244;
   this->am1CoreintegralS = -56.694056 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -48.717049 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 2.366515;      
   this->am1OrbitalExponentP = 1.667263;      
   this->am1BondingParameterS = -3.920566 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -7.905278 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.461648 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss =  11.786329 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gpp =  10.039308 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gsp =   8.663127 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gpp2 =  7.781688 * Parameters::GetInstance()->GetEV2AU();  
   this->am1Hsp =   2.532137 * Parameters::GetInstance()->GetEV2AU();   
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.6630874862;    
   this->am1DerivedParameterD[2] = 0.7345840886;    
   this->am1DerivedParameterRho[0] = 0.5/0.4331361580;
   this->am1DerivedParameterRho[1] = 0.5/0.6985577161;  
   this->am1DerivedParameterRho[2] = 0.5/0.7835093650;  
   this->am1ParameterK[0] =-0.509195 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] =-0.011863 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.012334 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 4.593691 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.865731 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 13.557336 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 0.770665 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.503313 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.009173 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
}

double Satom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
   double value = 0.0;

   if(theory == INDO){
      stringstream ss;
      ss << this->errorMessageIndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }
   else if(theory == ZINDOS){
      if(orbital == s){
         value = this->GetZindoCoreIntegral(orbital, 2, 4, 0);
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = this->GetZindoCoreIntegral(orbital, 2, 4, 0);
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
