#ifndef INCLUDED_CATOM
#define INCLUDED_CATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Catom : public Atom {
public:
   Catom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Catom::Catom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = C;
   this->atomicMass = 12.0107*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -21.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -17.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo = 0.0;
   this->coreCharge = 4.0;
   this->imuAmuS = 14.051*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 5.572*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->valenceShellType = l;
   this->effectiveNuclearChargeK = 5.7;
   this->effectiveNuclearChargeL = 3.25;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->numberValenceElectrons = 4;
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   this->zindoF0ss = 11.11 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 55635*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 36375*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->ionPotS = 19.84 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 10.93 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -52.279745 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -39.205558 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 1.787537;      
   this->mndoOrbitalExponentP = 1.787537;      
   this->mndoBondingParameterS = -18.985044 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -7.934122  * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.546380 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -120.500606 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 170.89 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoDiffEnergyHeatsFormAtom = -127.910952 * Parameters::GetInstance()->GetEV2AU();
   this->mndoGss =  12.23 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  11.08 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  11.47 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  9.84 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.43 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.8074661800;
   this->mndoDerivedParameterD[2] =   0.6851577737;    
   this->mndoDerivedParameterRho[0] = 0.5/0.4494406369;  
   this->mndoDerivedParameterRho[1] = 0.5/0.6149309919;
   this->mndoDerivedParameterRho[2] = 0.5/0.6685771472;  
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.427284 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.362563 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.588660 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.430254 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.395734 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -52.028658 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -39.614239 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 1.808665;      
   this->am1OrbitalExponentP = 1.685116;      
   this->am1BondingParameterS = -15.715783 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -7.719283  * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.648274 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.8236735591;    
   this->am1DerivedParameterD[2] = 0.7268015207;    
   this->am1DerivedParameterRho[0] = 0.5/0.4494406369;
   this->am1DerivedParameterRho[1] = 0.5/0.6082783276;  
   this->am1DerivedParameterRho[2] = 0.5/0.6423370115;  
   this->am1ParameterK[0] = 0.011355 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.045924 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.020061 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] =-0.001260 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.60 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.85 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.05 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 2.65 * Parameters::GetInstance()->GetAngstrom2AU();
}

double Catom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
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
         value = this->GetZindoCoreIntegral(orbital, 2, 2, 0);
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = this->GetZindoCoreIntegral(orbital, 2, 2, 0);
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
