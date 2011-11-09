#ifndef INCLUDED_CATOM
#define INCLUDED_CATOM

using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

class Catom : public Atom {
public:
   Catom(double x, double y, double z);
private:
};

Catom::Catom(double x, double y, double z) : Atom(x, y, z){
   this->atomType = C;
   this->atomicMass = 12.0107*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 4.0;
   this->numberValenceElectrons = 4;
   this->valenceShellType = l;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   this->bondingParameter = -21.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 14.051*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 5.572*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 5.7;
   this->effectiveNuclearChargeL = 3.25;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   this->indoF0CoefficientS = (this->coreCharge - 0.5);
   this->indoF0CoefficientP = (this->coreCharge - 0.5);
   this->indoG1CoefficientS = -1.0*(this->coreCharge - 1.5)/6.0;
   this->indoG1CoefficientP = -1.0/3.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = -2.0*(this->coreCharge - 2.5)/25.0;
   this->zindoBondingParameterS = -17.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 11.11 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 55635*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 36375*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 2;
   this->zindoM = 2;
   this->zindoN = 0;
   this->ionPotS = 19.84 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotP = 10.93 * Parameters::GetInstance()->GetEV2AU();
   this->ionPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -52.279745 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -39.205558 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 1.787537;      
   this->mndoOrbitalExponentP = 1.787537;      
   this->mndoBondingParameterS = -18.985044 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -7.934122  * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.546380 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -120.500606 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 170.89 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  12.23 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  11.08 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  11.47 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  9.84 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.43 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.8074661800;
   this->mndoDerivedParameterD[2] =   0.6851577737;    
   this->mndoDerivedParameterRho[0] = 0.5/0.4494406369;  
   this->mndoDerivedParameterRho[1] = 0.5/0.6149309919;
   this->mndoDerivedParameterRho[2] = 0.5/0.6685771472;  
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.427284 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.362563 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.588660 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.430254 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.395734 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -52.028658 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -39.614239 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 1.808665;      
   this->am1OrbitalExponentP = 1.685116;      
   this->am1BondingParameterS = -15.715783 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -7.719283  * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.648274 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.8236735591;    
   this->am1DerivedParameterD[2] = 0.7268015207;    
   this->am1DerivedParameterRho[0] = 0.5/0.4494406369;
   this->am1DerivedParameterRho[1] = 0.5/0.6082783276;  
   this->am1DerivedParameterRho[2] = 0.5/0.6423370115;  
   this->am1ParameterK[0] = 0.011355 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.045924 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.020061 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] =-0.001260 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.60 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.85 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.05 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 2.65 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3CoreintegralS = -47.270320 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -36.266918 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 1.565085;      
   this->pm3OrbitalExponentP = 1.842345;      
   this->pm3BondingParameterS = -11.910015 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -9.802755 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 2.707807 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3Gss = 11.200708 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 10.796292 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 10.265027 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2 = 9.042566 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp = 2.290980 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.8332396384;    
   this->pm3DerivedParameterD[2] = 0.6647749859;    
   this->pm3DerivedParameterRho[0] = 0.5/0.4116151543;
   this->pm3DerivedParameterRho[1] = 0.5/0.5885706542;  
   this->pm3DerivedParameterRho[2] = 0.5/0.7647513703;  
   this->pm3ParameterK[0] = 0.050107 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = 0.050733 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.003165 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.002979 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.642214 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 0.892488 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
}
}
#endif
