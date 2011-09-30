#ifndef INCLUDED_SATOM
#define INCLUDED_SATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Satom : public Atom {
public:
   Satom(double x, double y, double z);
   double GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory); 
private:
};

Satom::Satom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = S;
   this->atomicMass = 32.066*Parameters::GetInstance()->GetGMolin2AU();
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2 || 
      Parameters::GetInstance()->GetCurrentTheory() == INDO){
      this->valence.push_back(dxy);
      this->valence.push_back(dyz);
      this->valence.push_back(dzz);
      this->valence.push_back(dzx);
      this->valence.push_back(dxxyy);
   }
   this->bondingParameter = -18.150*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterSZindo = -14.0*Parameters::GetInstance()->GetEV2AU();
   this->bondingParameterDZindo =   4.0*Parameters::GetInstance()->GetEV2AU();
   this->coreCharge = 6.0;
   this->imuAmuS = 17.650*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 6.989*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.713*Parameters::GetInstance()->GetEV2AU();
   this->valenceShellType = m;
   this->effectiveNuclearChargeK = 15.70;
   this->effectiveNuclearChargeL = 11.85;
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->effectiveNuclearChargeMsp = 1.925*3.0;
      this->effectiveNuclearChargeMd = 1.731*3.0;
   }
   else{
      this->effectiveNuclearChargeMsp = 5.45;
      this->effectiveNuclearChargeMd = 5.45;
   }
   this->numberValenceElectrons = 6;
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   // the zindoF0ss for sulfer atoms are set to be equal 
   // to the one (10.09eV) in "ORCA 2.8"( http://www.thch.uni-bonn.de/tc/orca/ ).
   this->zindoF0ss = 10.09 * Parameters::GetInstance()->GetEV2AU(); 
   //this->zindoF0ss = 8.96 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 3.10 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 4.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 3.25 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG1pd = 4.31 * Parameters::GetInstance()->GetEV2AU();        
   this->zindoF2pd = 3.45 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG3pd = 2.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF2dd = 3.55 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF4dd = 2.31 * Parameters::GetInstance()->GetEV2AU();
   //this->zindoG1sp = 24807*Parameters::GetInstance()->GetKayser2AU();                 
   //this->zindoF2pp = 36600*Parameters::GetInstance()->GetKayser2AU();                 
   //this->zindoG2sd = 25972*Parameters::GetInstance()->GetKayser2AU();     
   //this->zindoG1pd = 34486*Parameters::GetInstance()->GetKayser2AU();        
   //this->zindoF2pd = 29173*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoG3pd = 20587*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoF2dd = 28411*Parameters::GetInstance()->GetKayser2AU();           
   //this->zindoF4dd = 18529*Parameters::GetInstance()->GetKayser2AU();           
   this->ionPotS = 21.11 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 12.39 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 4.11 * Parameters::GetInstance()->GetEV2AU();
}

double Satom::GetCoreIntegral(OrbitalType orbital, double gamma, bool isGuess, TheoryType theory){
   double value = 0.0;

   if(theory == INDO){
      stringstream ss;
      ss << this->errorMessageIndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }
   else if(theory == ZINDOS){
      if(orbital == s){
         value = this->GetZindoCoreIntegral(orbital, 2, 4, 0);
      }
      else if(orbital == px || orbital == py || orbital == pz){
         value = this->GetZindoCoreIntegral(orbital, 2, 4, 0);
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
