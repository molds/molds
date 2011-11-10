#ifndef INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
#define INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
namespace MolDS_base{

// l and m correspond to l and m in Eq.(15) in J. Comp. Chem. 20, 383(1999)
class RealSphericalHarmonicsIndex {
public:
   RealSphericalHarmonicsIndex(int l, int m);
   RealSphericalHarmonicsIndex(MolDS_base::OrbitalType  orbitalType);
   int GetL();
   int GetM();
private:
   int l;
   int m;
   RealSphericalHarmonicsIndex();
};

}
#endif
