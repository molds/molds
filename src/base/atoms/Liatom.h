#ifndef INCLUDED_LIATOM
#define INCLUDED_LIATOM
#include<string>
#include"Atom.h"

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Liatom : public Atom {
private:
public:
   Liatom(double x, double y, double z);
   double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess); // P82 - 83 in J. A. Pople book.
};

Liatom::Liatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = Li;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
   this->coreCharge = 1.0;
   this->imuAmuS = 3.106*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 1.258*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->averageIonizationPotentialS = 5.39*Parameters::GetInstance()->GetEV2AU();
   this->averageIonizationPotentialP = 3.54*Parameters::GetInstance()->GetEV2AU();
   this->valenceShellType = l;
   this->effectiveNuclearChargeK = 2.7;
   this->effectiveNuclearChargeL = 1.3;
   this->effectiveNuclearChargeM = 0.0;
   this->numberValenceElectrons = 1;
   this->electronicConfigM = 1;
   this->electronicConfigN = 0;
   this->indoG1 = 0.092012;
   this->indoF2 = 0.049865;
}

// P82 - 83 in J. A. Pople book.
double Liatom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess){
   
   double value = 0.0;
   if(orbital == s){
      value = -1.0*this->imuAmuS;
      if(!isGuess){
         value -= 0.5*gamma;
      }
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->imuAmuP;
      if(!isGuess){
         value -= 0.5*gamma - this->indoG1/12.0;
      }
   }
   else{
      cout << this->errorMessageIndoCoreIntegral;
      cout << this->errorMessageAtomType << this->atomType << endl;
      exit(EXIT_FAILURE);
   }
   return value;
}

}
#endif
