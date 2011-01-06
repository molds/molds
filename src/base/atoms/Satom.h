#ifndef INCLUDED_SATOM
#define INCLUDED_SATOM
#include"Atom.h"

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
   this->effectiveNuclearChargeMsp = 5.45;
   this->effectiveNuclearChargeMd = 5.45;
   this->numberValenceElectrons = 6;
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   this->zindoF0ss = 8.96 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 24807*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 36600*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 25972*Parameters::GetInstance()->GetKayser2AU();     
   this->zindoG1pd = 34486*Parameters::GetInstance()->GetKayser2AU();        
   this->zindoF2pd = 29173*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoG3pd = 20587*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF2dd = 28411*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF4dd = 18529*Parameters::GetInstance()->GetKayser2AU();           
   this->IonPotS = 21.11 * Parameters::GetInstance()->GetEV2AU();
   this->IonPotP = 12.39 * Parameters::GetInstance()->GetEV2AU();
   this->IonPotD = 4.11 * Parameters::GetInstance()->GetEV2AU();
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

   return value;
}



}
#endif
