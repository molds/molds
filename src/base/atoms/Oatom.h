#ifndef INCLUDED_OATOM
#define INCLUDED_OATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Oatom : public Atom {
public:
   Oatom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Oatom::Oatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = O;
   this->atomicMass = 15.9994*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -31.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -34.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo = 0.0;
   this->coreCharge = 6.0;
   this->imuAmuS = 25.390*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 9.111*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->valenceShellType = l;
   this->effectiveNuclearChargeK = 7.70;
   this->effectiveNuclearChargeL = 4.55;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->numberValenceElectrons = 6;
   this->indoG1 = 0.346029;
   this->indoF2 = 0.219055;
   this->zindoF0ss = 13.00 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 95298*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 55675*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->ionPotS = 32.90 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 17.28 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -99.64309 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -77.797472 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.699905;      
   this->mndoOrbitalExponentP = 2.699905;      
   this->mndoBondingParameterS = -32.688082 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -32.688082 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 3.160604 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -317.868506 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 59.559 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoDiffEnergyHeatsFormAtom = -320.451178 * Parameters::GetInstance()->GetEV2AU();
   this->mndoGss =  15.42 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  14.52 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  14.48 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 12.98 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   3.94 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.5346023927; 
   this->mndoDerivedParameterD[2] =   0.4536251725;
   this->mndoDerivedParameterRho[0] = 0.5/0.5666700426;
   this->mndoDerivedParameterRho[1] = 0.5/0.9592303457;  
   this->mndoDerivedParameterRho[2] = 0.5/0.9495760934;  
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.282894 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.240043 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.466882 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.275822 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.278628 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -97.830000 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -78.262380 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 3.108032;      
   this->am1OrbitalExponentP = 2.524039;      
   this->am1BondingParameterS = -29.272773 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -29.272773 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 4.455371 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.4988896404;    
   this->am1DerivedParameterD[2] = 0.4852321503;    
   this->am1DerivedParameterRho[0] = 0.5/0.5666700426;
   this->am1DerivedParameterRho[1] = 0.5/0.9960801167;  
   this->am1DerivedParameterRho[2] = 0.5/0.9065055775;  
   this->am1ParameterK[0] = 0.280962 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.081430 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 7.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 0.847918 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.445071 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
}

double Oatom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
   double value = 0.0;

   if(theory == INDO){
      if(orbital == s){
         value = -1.0*this->imuAmuS;
         if(!isGuess){
            value -= (this->coreCharge-0.5)*gamma - (this->coreCharge - 1.5)*this->indoG1/6.0;
         }
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = -1.0*this->imuAmuP;
         if(!isGuess){
            value -= (this->coreCharge-0.5)*gamma 
                    - this->indoG1/3.0 
                    - (this->coreCharge - 2.5)*this->indoF2*2.0/25.0;
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
