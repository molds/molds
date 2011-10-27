#ifndef INCLUDED_NATOM
#define INCLUDED_NATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Natom : public Atom {
public:
   Natom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Natom::Natom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = N;
   this->atomicMass = 14.00674*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -25.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -26.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo = 0.0;
   this->coreCharge = 5.0;
   this->imuAmuS = 19.316*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 7.275*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->valenceShellType = l;
   this->effectiveNuclearChargeK = 6.7;
   this->effectiveNuclearChargeL = 3.90;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->numberValenceElectrons = 5;
   this->indoG1 = 0.346029;
   this->indoF2 = 0.219055;
   this->zindoF0ss = 12.01 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 72255*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 52100*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->ionPotS = 25.69 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 14.05 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -71.932122 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -57.172319 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.255614;      
   this->mndoOrbitalExponentP = 2.255614;      
   this->mndoBondingParameterS = -20.495758 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -20.495758 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.861342 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -202.581201 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 113.00 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  13.59 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  12.98 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  12.66 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 11.59 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   3.14 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.6399036683;
   this->mndoDerivedParameterD[2] =   0.5429762678;
   this->mndoDerivedParameterRho[0] = 0.5/0.4994193177;
   this->mndoDerivedParameterRho[1] = 0.5/0.7843433156;
   this->mndoDerivedParameterRho[2] = 0.5/0.8126295047;
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.338616 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.287325 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.529751 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.337322 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.324853 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -71.860000 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -57.167581 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 2.315410;      
   this->am1OrbitalExponentP = 2.157940;      
   this->am1BondingParameterS = -20.299110 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -18.238666 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.947286 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.6433247425;    
   this->am1DerivedParameterD[2] = 0.5675527917;    
   this->am1DerivedParameterRho[0] = 0.5/0.4994193177;
   this->am1DerivedParameterRho[1] = 0.5/0.7820630445;  
   this->am1DerivedParameterRho[2] = 0.5/0.7883351388;  
   this->am1ParameterK[0] = 0.025251 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.028953 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.005806 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 2.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.50 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 2.10 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.40 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
}

double Natom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
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
         value = this->GetZindoCoreIntegral(orbital, 2, 3, 0);
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = this->GetZindoCoreIntegral(orbital, 2, 3, 0);
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
