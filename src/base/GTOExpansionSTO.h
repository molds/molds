#ifndef INCLUDED_GTOEXPANSIONSTO
#define INCLUDED_GTOEXPANSIONSTO
namespace MolDS_base{

// GTOExpansionSTO is singleton
class GTOExpansionSTO{
public:
   static GTOExpansionSTO* GetInstance();
   static void DeleteInstance();
   double GetExponent(MolDS_base::STOnGType stonG, 
                      MolDS_base::ShellType shellType, 
                      OrbitalType orbitalType, 
                      int index);
   double GetCoefficient(MolDS_base::STOnGType stonG, 
                         MolDS_base::ShellType shellType, 
                         OrbitalType orbitalType, 
                         int index);

private:
   static GTOExpansionSTO* gTOExpansionSTO;
   GTOExpansionSTO();
   GTOExpansionSTO(GTOExpansionSTO&);
   void operator = (GTOExpansionSTO&);
   ~GTOExpansionSTO();

   std::string errorMessageGetCoefficientNonValidOrbital;
   std::string errorMessageGetExponentNonValidOrbital;
   std::string errorMessageOrbitalType;
   std::string errorMessageSTOnGType;
   void SetCoefficientsExponents();
   double exponents[MolDS_base::STOnGType_end]
                   [MolDS_base::ShellType_end]
                   [MolDS_base::AzimuthalType_end]
                   [6];    
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is alpha in (3) of [S_1970]. See Table I and II in [S_1970]
   double coefficients[MolDS_base::STOnGType_end]
                      [MolDS_base::ShellType_end]
                      [MolDS_base::AzimuthalType_end]
                      [6]; 
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is d in (3) of [S_1970]. See Table I and II in [S_1970]
};

}
#endif





