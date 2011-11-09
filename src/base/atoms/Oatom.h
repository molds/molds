#ifndef INCLUDED_OATOM
#define INCLUDED_OATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Oatom : public Atom {
public:
   Oatom(double x, double y, double z);
private:
};

Oatom::Oatom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = O;
   this->atomicMass = 15.9994*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 6.0;
   this->numberValenceElectrons = 6;
   this->valenceShellType = l;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -31.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 25.390*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 9.111*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 7.70;
   this->effectiveNuclearChargeL = 4.55;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.346029;
   this->indoF2 = 0.219055;
   this->indoF0CoefficientS = (this->coreCharge - 0.5);
   this->indoF0CoefficientP = (this->coreCharge - 0.5);
   this->indoG1CoefficientS = -1.0*(this->coreCharge - 1.5)/6.0;
   this->indoG1CoefficientP = -1.0/3.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = -2.0*(this->coreCharge - 2.5)/25.0;
   this->zindoBondingParameterS = -34.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 13.00 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 95298*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 55675*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 2;
   this->zindoM = 4;
   this->zindoN = 0;
   this->zindoIonPotS = 32.90 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 17.28 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -99.64309 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -77.797472 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.699905;      
   this->mndoOrbitalExponentP = 2.699905;      
   this->mndoBondingParameterS = -32.688082 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -32.688082 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 3.160604 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -317.868506 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 59.559 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  15.42 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  14.52 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  14.48 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 12.98 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   3.94 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.5346023927; 
   this->mndoDerivedParameterD[2] =   0.4536251725;
   this->mndoDerivedParameterRho[0] = 0.5/0.5666700426;
   this->mndoDerivedParameterRho[1] = 0.5/0.9592303457;  
   this->mndoDerivedParameterRho[2] = 0.5/0.9495760934;  
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.282894 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.240043 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.466882 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.275822 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.278628 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -97.830000 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -78.262380 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 3.108032;      
   this->am1OrbitalExponentP = 2.524039;      
   this->am1BondingParameterS = -29.272773 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -29.272773 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 4.455371 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.4988896404;    
   this->am1DerivedParameterD[2] = 0.4852321503;    
   this->am1DerivedParameterRho[0] = 0.5/0.5666700426;
   this->am1DerivedParameterRho[1] = 0.5/0.9960801167;  
   this->am1DerivedParameterRho[2] = 0.5/0.9065055775;  
   this->am1ParameterK[0] = 0.280962 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.081430 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 7.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 0.847918 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.445071 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3CoreintegralS = -86.993002 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -71.879580 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 3.796544;      
   this->pm3OrbitalExponentP = 2.389402;      
   this->pm3BondingParameterS = -45.202651 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -24.752515 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 3.217102 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3Gss = 15.755760 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 13.654016 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 10.621160 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2 = 12.40609 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp = 0.593883 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.4086173087;    
   this->pm3DerivedParameterD[2] = 0.5125738036;    
   this->pm3DerivedParameterRho[0] = 0.5/0.5790088969;
   this->pm3DerivedParameterRho[1] = 0.5/0.5299517372;  
   this->pm3DerivedParameterRho[2] = 0.5/0.8179482975;  
   this->pm3ParameterK[0] = -1.131128 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = 1.137891 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.002477 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 5.950512 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.607311 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.598395 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
}
}
#endif
