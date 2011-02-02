#ifndef INCLUDED_OATOM
#define INCLUDED_OATOM
#include"Atom.h"

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

   return value;
}



}
#endif