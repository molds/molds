#ifndef INCLUDED_CATOM
#define INCLUDED_CATOM
#include"Atom.h"

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Catom : public Atom {
private:
public:
   Catom(double x, double y, double z);
   double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess); // P82 - 83 in J. A. Pople book.
};

Catom::Catom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = C;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -21.0*Parameters::GetInstance()->GetEV2AU();
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
}

// P82 - 83 in J. A. Pople book.
double Catom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess){
   double value = 0.0;
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
      cout << this->errorMessageIndoCoreIntegral;
      cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      exit(EXIT_FAILURE);
   }
   return value;
}



}
#endif
