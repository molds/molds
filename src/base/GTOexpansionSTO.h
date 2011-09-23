#ifndef INCLUDED_GTOEXPANSIONSTO
#define INCLUDED_GTOEXPANSIONSTO

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>

using namespace std;
namespace MolDS_base{

// GTOexpansionSTO is singleton
class GTOexpansionSTO{
public:
   static GTOexpansionSTO* GetInstance();
   static void DeleteInstance();

private:
   static GTOexpansionSTO* gTOexpansionSTO;
   GTOexpansionSTO();
   GTOexpansionSTO(GTOexpansionSTO&);
   void operator = (GTOexpansionSTO&);
   ~GTOexpansionSTO();

   void SetCoefficientsExponents();
   double exponents[STOnGType_end][ShellType_end][AzimuthalType_end][6];    
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is alpha in (3) of [S_1970]. See Table I and II in [S_1970]
   double coefficients[STOnGType_end][ShellType_end][AzimuthalType_end][6]; 
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is d in (3) of [S_1970]. See Table I and II in [S_1970]

};
GTOexpansionSTO* GTOexpansionSTO::gTOexpansionSTO = NULL;

GTOexpansionSTO::GTOexpansionSTO(){
   this->SetCoefficientsExponents();
}

GTOexpansionSTO::~GTOexpansionSTO(){
}

GTOexpansionSTO* GTOexpansionSTO::GetInstance(){
   if(gTOexpansionSTO == NULL){
      gTOexpansionSTO = new GTOexpansionSTO();
   }
   return gTOexpansionSTO;
}

void GTOexpansionSTO::DeleteInstance(){
   if(gTOexpansionSTO != NULL){
      delete gTOexpansionSTO; 
   }
   gTOexpansionSTO = NULL;
}

//  see Table I and II in [S_1970]
void GTOexpansionSTO::SetCoefficientsExponents(){

   //STO-1G, k-shell
   {   
      // 1s
      exponents[STO1G][k][sAzimuthal][0] = 2.709498091*pow(10.0,-1.0);   coefficients[STO1G][k][sAzimuthal][0] = 1.0000;
   }

   //STO-1G, l-shell
   {
      // 2s
      exponents[STO1G][l][sAzimuthal][0] = 1.012151084*pow(10.0,-1.0);   coefficients[STO1G][l][sAzimuthal][0] = 1.0000;
      // 2p
      exponents[STO1G][l][pAzimuthal][0] = 1.759666885*pow(10.0,-1.0);   coefficients[STO1G][l][sAzimuthal][0] = 1.0000;
   }

   //STO-1G, m-shell
   {
      // 3s
      exponents[STO1G][m][sAzimuthal][0] = 5.296881757*pow(10.0,-1.0);   coefficients[STO1G][m][sAzimuthal][0] = 1.0000;
      // 3p
      exponents[STO1G][m][pAzimuthal][0] = 9.113614253*pow(10.0,-2.0);   coefficients[STO1G][m][sAzimuthal][0] = 1.0000;
      // 3d
      exponents[STO1G][m][dAzimuthal][0] = 1.302270363*pow(10.0,-1.0);   coefficients[STO1G][m][sAzimuthal][0] = 1.0000;
   }

   //STO-2G, k-shell
   {
      // 1s
      exponents[STO2G][k][sAzimuthal][0] = 8.518186635*pow(10.0,-1.0);   coefficients[STO2G][k][sAzimuthal][0] = 4.301284983*pow(10,-1.0); 
      exponents[STO2G][k][sAzimuthal][1] = 1.516232927*pow(10.0,-1.0);   coefficients[STO2G][k][sAzimuthal][1] = 6.789135305*pow(10,-1.0); 
   }   
/*
   //STO-2G, l-shell
   {
      // 2s
      exponents[STO2G][l][sAzimuthal][0] = 1.292278611*pow(10.0,-1.0);   coefficients[STO2G][l][sAzimuthal][0] = 7.470867124*pow(10,-1.0); 
      exponents[STO2G][l][sAzimuthal][1] = 4.908584205*pow(10.0,-2.0);   coefficients[STO2G][l][sAzimuthal][1] = 2.855980556*pow(10,-1.0); 
      // 2p
      exponents[STO2G][l][pAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO2G][l][pAzimuthal][0] = *pow(10,-.0); 
      exponents[STO2G][l][pAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO2G][l][pAzimuthal][1] = *pow(10,-.0); 
   }   

   //STO-2G, m-shell
   {
      // 3s
      exponents[STO2G][m][sAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO2G][m][sAzimuthal][0] = *pow(10,-.0); 
      exponents[STO2G][m][sAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO2G][m][sAzimuthal][1] = *pow(10,-.0); 
      // 3p
      exponents[STO2G][m][pAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO2G][m][pAzimuthal][0] = *pow(10,-.0); 
      exponents[STO2G][m][pAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO2G][m][pAzimuthal][1] = *pow(10,-.0); 
      // 3d
      exponents[STO2G][m][dAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO2G][m][dAzimuthal][0] = *pow(10,-.0); 
      exponents[STO2G][m][dAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO2G][m][dAzimuthal][1] = *pow(10,-.0); 
   }

   //STO-3G, k-shell
   {
      // 1s
      exponents[STO3G][k][sAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][k][sAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][k][sAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][k][sAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][k][sAzimuthal][2] = *pow(10.0,-.0);   coefficients[STO3G][k][sAzimuthal][2] = *pow(10,-.0); 
   }   

   //STO-3G, l-shell
   {
      // 2s
      exponents[STO3G][l][sAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][l][sAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][l][sAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][l][sAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][l][sAzimuthal][2] = *pow(10.0,-.0);   coefficients[STO3G][l][sAzimuthal][2] = *pow(10,-.0); 
      // 2p
      exponents[STO3G][l][pAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][l][pAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][l][pAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][l][pAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][l][pAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][l][pAzimuthal][1] = *pow(10,-.0); 
   }   

   //STO-3G, m-shell
   {
      // 3s
      exponents[STO3G][m][sAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][m][sAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][m][sAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][m][sAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][m][sAzimuthal][2] = *pow(10.0,-.0);   coefficients[STO3G][m][sAzimuthal][2] = *pow(10,-.0); 
      // 3p
      exponents[STO3G][m][pAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][m][pAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][m][pAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][m][pAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][m][pAzimuthal][2] = *pow(10.0,-.0);   coefficients[STO3G][m][pAzimuthal][2] = *pow(10,-.0); 
      // 3d
      exponents[STO3G][m][dAzimuthal][0] = *pow(10.0,-.0);   coefficients[STO3G][m][dAzimuthal][0] = *pow(10,-.0); 
      exponents[STO3G][m][dAzimuthal][1] = *pow(10.0,-.0);   coefficients[STO3G][m][dAzimuthal][1] = *pow(10,-.0); 
      exponents[STO3G][m][dAzimuthal][2] = *pow(10.0,-.0);   coefficients[STO3G][m][dAzimuthal][2] = *pow(10,-.0); 
   }
*/
}

}
#endif





