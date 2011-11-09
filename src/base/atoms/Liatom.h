#ifndef INCLUDED_LIATOM
#define INCLUDED_LIATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Liatom : public Atom {
public:
   Liatom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Liatom::Liatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = Li;
   this->atomicMass = 6.941*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = 0.0;
   this->bondingParameterDZindo = 0.0;
   this->coreCharge = 1.0;
   this->imuAmuS = 3.106*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 1.258*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->valenceShellType = l;
   this->effectiveNuclearChargeK = 2.7;
   this->effectiveNuclearChargeL = 1.3;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->numberValenceElectrons = 1;
   this->indoG1 = 0.092012;
   this->indoF2 = 0.049865;
   this->indoF0CoefficientS = 0.5;
   this->indoF0CoefficientP = 0.5;
   this->indoG1CoefficientS = 0.0;
   this->indoG1CoefficientP = -1.0/12.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = 0.0;
   this->zindoF0ss = 0.0;
   this->zindoF0sd = 0.0;        
   this->zindoF0dd = 0.0;              
   this->zindoG1sp = 20194*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF2pp = 10944*Parameters::GetInstance()->GetKayser2AU();        
   this->zindoG2sd = 0.0;                
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                
   this->zindoG3pd = 0.0;             
   this->zindoF2dd = 0.0;             
   this->zindoF4dd = 0.0;       
   this->zindoL = 1;
   this->zindoM = 0;
   this->zindoN = 0;
   this->ionPotS = 5.39 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 3.54 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
}

double Liatom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
   
   double value = 0.0;

   if(theory == INDO){
      if(orbital == s){
         value = -1.0*this->imuAmuS;
         if(!isGuess){
            value -= this->indoF0CoefficientS*gamma 
                    +this->indoG1CoefficientS*this->indoG1
                    +this->indoF2CoefficientS*this->indoF2;
         }
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = -1.0*this->imuAmuP;
         if(!isGuess){
            value -= this->indoF0CoefficientP*gamma 
                    +this->indoG1CoefficientP*this->indoG1
                    +this->indoF2CoefficientP*this->indoF2;
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
      value = this->GetZindoCoreIntegral(orbital);
   }

   return value;
}

}
#endif
