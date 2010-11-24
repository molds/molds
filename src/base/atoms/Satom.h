#ifndef INCLUDED_SATOM
#define INCLUDED_SATOM
#include"Atom.h"

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Satom : public Atom {
private:
public:
   Satom(double x, double y, double z);
   double GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess); // P82 - 83 in J. A. Pople book.
};

Satom::Satom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = S;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->valence.push_back(dxy);
   this->valence.push_back(dyz);
   this->valence.push_back(dzz);
   this->valence.push_back(dzx);
   this->valence.push_back(dxxyy);
   this->bondingParameter = -18.150*Parameters::GetInstance()->GetEV2AU();
   this->coreCharge = 6.0;
   this->imuAmuS = 17.650*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 6.989*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.713*Parameters::GetInstance()->GetEV2AU();
   this->valenceShellType = m;
   this->effectiveNuclearChargeK = 15.70;
   this->effectiveNuclearChargeL = 11.85;
   this->effectiveNuclearChargeM = 5.45;
   this->numberValenceElectrons = 6;
   //this->indoG1 = 0.267708;
   //this->indoF2 = 0.17372;
}

// P82 - 83 in J. A. Pople book.
double Satom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess){
   double value = 0.0;
   cout << this->errorMessageIndoCoreIntegral;
   cout << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
   cout << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
   exit(EXIT_FAILURE);
   return value;
}



}
#endif
