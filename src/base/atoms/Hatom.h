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

   return value;
}

}
#endif
