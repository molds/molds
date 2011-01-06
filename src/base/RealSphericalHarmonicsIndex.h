#ifndef INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
#define INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
#include<string>

using namespace std;
namespace MolDS_base{


// l and m correspond to l and m in Eq.(15) in J. Comp. Chem. 20, 383(1999)
class RealSphericalHarmonicsIndex {
public:
   RealSphericalHarmonicsIndex(int l, int m);
   RealSphericalHarmonicsIndex(OrbitalType orbitalType);
   int GetL();
   int GetM();
private:
   int l;
   int m;
   RealSphericalHarmonicsIndex();
};

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(){}

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(OrbitalType orbitalType){
   string errorMessageInvalidOrbital = "Error in base::RealSphericalHarmonicIndex::RealSphericalHarmonicIndex: invalid orbitalType. Indicated orbitalType is not prepared\n";

   if(orbitalType == s){
      this->l = 0;
      this->m = 0;
   }
   else if(orbitalType == py){
      this->l = 1;
      this->m = -1;
   }
   else if(orbitalType == pz){
      this->l = 1;
      this->m = 0;
   }
   else if(orbitalType == px){
      this->l = 1;
      this->m = 1;
   }
   else if(orbitalType == dxy){
      this->l = 2;
      this->m = -2;
   }
   else if(orbitalType == dyz){
      this->l = 2;
      this->m = -1;
   }
   else if(orbitalType == dzz){
      this->l = 2;
      this->m = 0;
   }
   else if(orbitalType == dzx){
      this->l = 2;
      this->m = 1;
   }
   else if(orbitalType == dxxyy){
      this->l = 2;
      this->m = 2;
   }
   else{
      cout << errorMessageInvalidOrbital;
   }

}

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(int l, int m){
   this->l = l;
   this->m = m;
}

int RealSphericalHarmonicsIndex::GetL(){
   return this->l;
}

int RealSphericalHarmonicsIndex::GetM(){
   return this->m;
}

}
#endif
