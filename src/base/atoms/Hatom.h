#ifndef INCLUDED_HATOM
#define INCLUDED_HATOM
#include"Atom.h"
#include<string>

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Hatom : public Atom {
private:
public:
   Hatom(double x, double y, double z);
   double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess); // P82 - 83 in J. A. Pople book.
};

Hatom::Hatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = H;
   this->valence.push_back(s);
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
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
}

// P82 - 83 in J. A. Pople book.
double Hatom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess){
   double value = 0.0;
   if(orbital == s){
      value = -1.0*this->imuAmuS;
      if(!isGuess){
         value -= 0.5*gamma;
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
