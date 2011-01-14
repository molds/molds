#ifndef INCLUDED_ZINDOS
#define INCLUDED_ZINDOS

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
#include"../cndo/Cndo2.h"

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_zindo{

/***
 *  Refferences for 
 */
class ZindoS : public MolDS_cndo::Cndo2{
public:
   ZindoS();
   ~ZindoS();
   void DoesCIS();
protected:
   virtual void CalcGammaAB(double** gammaAB, Molecule* molecule);
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetFockDiagElement(Atom* atomA, int atomAIndex, 
                                     int mu, Molecule* molecule, double** gammaAB,
                                     double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                     bool isGuess);
   virtual double GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                        int mu, int nu, Molecule* molecule, double** gammaAB, double** overelap,
                                        double** orbitalElectronPopulation, bool isGuess);
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, Atom* atomA, Atom* atomB);
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              Molecule* molecule, double** fockMatrix, double** gammaAB);
private:
   double GetCoulombInt(OrbitalType orbital1, 
                        OrbitalType orbital2, 
                        Atom* atom); // Apendix in [BZ_1979]
   double GetExchangeInt(OrbitalType orbital1, 
                         OrbitalType orbital2, 
                         Atom* atom); // Apendix in [BZ_1979]
   double GetNishimotoMatagaTwoEleInt(Atom* atomA, OrbitalType orbitalA, 
                                           Atom* atomB, OrbitalType orbitalB); // ref. [MN_1957] and (5a) in [AEZ_1986]
   string errorMessageNishimotoMataga;
   string messageStartCIS;
   string messageDoneCIS;
};

ZindoS::ZindoS() : MolDS_cndo::Cndo2(){
   this->theory = ZINDOS;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //cout << "ZindoS created\n";
}

ZindoS::~ZindoS(){
   //cout << "ZindoS deleted\n";
}

void ZindoS::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in zindo::ZindoS::DoesSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in zindo::ZindoS::DoesSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in zindo::ZindoS::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in zindo::ZindoS::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_zindo::ZindoS::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_zindo::ZindoS::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageNishimotoMataga = "Error in base_zindo::ZindoS::GetNishimotoMatagaTwoEleInt: Invalid orbitalType.\n";
   this->errorMessageMolecularIntegralElement
      = "Error in zindo::ZindoS::GetMolecularIntegralElement: Non available orbital is contained.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tZINDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: ZINDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: ZINDO/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: ZINDO/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: ZINDO/S-CIS  **********\n\n\n";
}

void ZindoS::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double ZindoS::GetFockDiagElement(Atom* atomA, int atomAIndex, int mu, 
                                 Molecule* molecule, double** gammaAB,
                                 double** orbitalElectronPopulation, double* atomicElectronPopulation,
                                 bool isGuess){
   double value=0.0;
   int firstAOIndexA = atomA->GetFirstAOIndex();
   value = atomA->GetCoreIntegral(atomA->GetValence()[mu-firstAOIndexA], 
                                  isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      OrbitalType orbitalMu = atomA->GetValence()[mu-firstAOIndexA];
      for(int v=0; v<atomA->GetValence().size(); v++){
         OrbitalType orbitalLam = atomA->GetValence()[v];
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalLam, atomA);
         exchange = this->GetExchangeInt(orbitalMu, orbitalLam, atomA);
         lammda = v + firstAOIndexA;
         temp += orbitalElectronPopulation[lammda][lammda]*(coulomb - 0.5*exchange);
      }
      value += temp;
   
      temp = 0.0;
      for(int B=0; B<molecule->GetAtomVect()->size(); B++){
         if(B != atomAIndex){
            Atom* atomB = (*molecule->GetAtomVect())[B];
            for(int i=0; i<atomB->GetValence().size(); i++){
               int sigma = i + atomB->GetFirstAOIndex();
               OrbitalType orbitalSigma = atomB->GetValence()[i];
               temp += orbitalElectronPopulation[sigma][sigma]
                      *this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalSigma);
            }
            temp -= atomB->GetCoreCharge() 
                   *this->GetNishimotoMatagaTwoEleInt(atomA, s, atomB, s);
         }
      }
      value += temp;
   }

   return value;
}

double ZindoS::GetFockOffDiagElement(Atom* atomA, Atom* atomB, int atomAIndex, int atomBIndex, 
                                    int mu, int nu, Molecule* molecule, double** gammaAB, double** overlap,
                                    double** orbitalElectronPopulation, bool isGuess){
   double value = 0.0;
   OrbitalType orbitalMu = atomA->GetValence()[mu-atomA->GetFirstAOIndex()];
   OrbitalType orbitalNu = atomB->GetValence()[nu-atomB->GetFirstAOIndex()];

   double bondParameter = 0.5*(atomA->GetBondingParameter(this->theory, orbitalMu) 
                              +atomB->GetBondingParameter(this->theory, orbitalNu)); 

   if(isGuess){
      value = bondParameter*overlap[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(atomAIndex == atomBIndex){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlap[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]
                  *this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
      }
   }
   
   return value;
}

void ZindoS::CalcGammaAB(double** gammaAB, Molecule* molecule){
   // Do nothing;
}

// Apendix in [BZ_1972]
// ZINDO Coulomb Interaction
double ZindoS::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;
   
   if( orbital1 == s && orbital2 == s){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetZindoF0ssLower();
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom->GetZindoF0ssLower()
             +atom->GetZindoF2ppLower()*4.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoF0ssLower()
             -atom->GetZindoF2ppLower()*2.0;
   }   
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( orbital1 == s && ( orbital2 == dxy || 
                               orbital2 == dyz || 
                               orbital2 == dzz || 
                               orbital2 == dzx || 
                               orbital2 == dxxyy )){ 
      value = atom->GetZindoF0sdLower();
   }   
   else if( orbital2 == s && ( orbital1 == dxy || 
                               orbital1 == dyz || 
                               orbital1 == dzz || 
                               orbital1 == dzx || 
                               orbital1 == dxxyy )){ 
      value = atom->GetZindoF0sdLower();
   }
   else if( orbital1 == dzz && (orbital2 == px || orbital2==py) ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*2.0;
   }
   else if( orbital2 == dzz && (orbital1 == px || orbital1==py) ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dzz && orbital2 == pz) ||
            (orbital2 == dzz && orbital1 == pz) ){
      value = atom->GetZindoF0sdLower()
             +atom->GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == orbital2) && ( orbital1 == dxy || 
                                        orbital1 == dyz || 
                                        orbital1 == dzz || 
                                        orbital1 == dzx || 
                                        orbital1 == dxxyy )){ 
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*36.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == px) ||
            (orbital2 == dxxyy && orbital1 == px) || 
            (orbital1 == dxxyy && orbital2 == py) ||
            (orbital2 == dxxyy && orbital1 == py) || 
            (orbital1 == dxy && orbital2 == px) ||
            (orbital2 == dxy && orbital1 == px) || 
            (orbital1 == dxy && orbital2 == py) ||
            (orbital2 == dxy && orbital1 == py) ||
            (orbital1 == dzx && orbital2 == px) ||
            (orbital2 == dzx && orbital1 == px) || 
            (orbital1 == dzx && orbital2 == pz) ||
            (orbital2 == dzx && orbital1 == pz) || 
            (orbital1 == dyz && orbital2 == py) ||
            (orbital2 == dyz && orbital1 == py) || 
            (orbital1 == dyz && orbital2 == pz) ||
            (orbital2 == dyz && orbital1 == pz) ){
      value = atom->GetZindoF0sdLower()
             +atom->GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == pz) ||
            (orbital2 == dxxyy && orbital1 == pz) || 
            (orbital1 == dxy && orbital2 == pz) ||
            (orbital2 == dxy && orbital1 == pz) ||
            (orbital1 == dzx && orbital2 == py) ||
            (orbital2 == dzx && orbital1 == py) ||
            (orbital1 == dyz && orbital2 == px) ||
            (orbital2 == dyz && orbital1 == px)  ){
      value = atom->GetZindoF0sdLower()
             -atom->GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzz) ||
            (orbital2 == dxxyy && orbital1 == dzz) || 
            (orbital1 == dxy && orbital2 == dzz) ||
            (orbital2 == dxy && orbital1 == dzz)  ){
      value = atom->GetZindoF0ddLower()
             -atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*6.0;
   }
   else if( (orbital1 == dxy && orbital2 == dxxyy) ||
            (orbital2 == dxy && orbital1 == dxxyy)  ){
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*4.0
             -atom->GetZindoF4ddLower()*34.0;
   }
   else if( (orbital1 == dzx && orbital2 == dzz) ||
            (orbital2 == dzx && orbital1 == dzz) || 
            (orbital1 == dyz && orbital2 == dzz) ||
            (orbital2 == dyz && orbital1 == dzz)  ){
      value = atom->GetZindoF0ddLower()
             +atom->GetZindoF2ddLower()*2.0
             -atom->GetZindoF4ddLower()*24.0;
   }
   else if( (orbital1 == dzx && orbital2 == dxxyy) ||
            (orbital2 == dzx && orbital1 == dxxyy) || 
            (orbital1 == dzx && orbital2 == dxy) || 
            (orbital2 == dzx && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dxxyy) || 
            (orbital2 == dyz && orbital1 == dxxyy) || 
            (orbital1 == dyz && orbital2 == dxy) || 
            (orbital2 == dyz && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dzx) || 
            (orbital2 == dyz && orbital1 == dzx) ){
      value = atom->GetZindoF0ddLower()
             -atom->GetZindoF2ddLower()*2.0
             -atom->GetZindoF4ddLower()*4.0;
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageCoulombInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   
   return value;

}

// Apendix in [BZ_1972]
// ZINDO Exchange Interaction
double ZindoS::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, Atom* atom){

   double value=0.0;

   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoG1spLower();
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = atom->GetZindoG1spLower();
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom->GetZindoF2ppLower()*3.0;
   }
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( (orbital1 == s) && (orbital2 == dxy || 
                                orbital2 == dyz || 
                                orbital2 == dzz || 
                                orbital2 == dzx || 
                                orbital2 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital2 == s) && (orbital1 == dxy || 
                                orbital1 == dyz || 
                                orbital1 == dzz || 
                                orbital1 == dzx || 
                                orbital1 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital1 == px && orbital2 == dzz) ||
            (orbital2 == px && orbital1 == dzz) ||
            (orbital1 == py && orbital2 == dzz) ||
            (orbital2 == py && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()
             +atom->GetZindoG3pdLower()*18.0;
   }
   else if( (orbital1 == px && orbital2 == dxxyy) ||
            (orbital2 == px && orbital1 == dxxyy) ||
            (orbital1 == px && orbital2 == dxy) ||
            (orbital2 == px && orbital1 == dxy) ||
            (orbital1 == px && orbital2 == dzx) ||
            (orbital2 == px && orbital1 == dzx) ||
            (orbital1 == py && orbital2 == dxxyy) ||
            (orbital2 == py && orbital1 == dxxyy) ||
            (orbital1 == py && orbital2 == dxy) ||
            (orbital2 == py && orbital1 == dxy) ||
            (orbital1 == py && orbital2 == dyz) ||
            (orbital2 == py && orbital1 == dyz) ||
            (orbital1 == pz && orbital2 == dzx) ||
            (orbital2 == pz && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dyz) ||
            (orbital2 == pz && orbital1 == dyz) ){
      value = atom->GetZindoG1pdLower()*3.0
             +atom->GetZindoG3pdLower()*24.0;
   }
   else if( (orbital1 == px && orbital2 == dyz) ||
            (orbital2 == px && orbital1 == dyz) ||
            (orbital1 == py && orbital2 == dzx) ||
            (orbital2 == py && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dxxyy) ||
            (orbital2 == pz && orbital1 == dxxyy) ||
            (orbital1 == pz && orbital2 == dxy) ||
            (orbital2 == pz && orbital1 == dxy) ){
      value = atom->GetZindoG3pdLower()*15.0;
   }
   else if( (orbital1 == pz && orbital2 == dzz) ||
            (orbital2 == pz && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()*4.0
             +atom->GetZindoG3pdLower()*27.0;
   }
   else if( (orbital1 == dzz && orbital2 == dxxyy) ||
            (orbital2 == dzz && orbital1 == dxxyy) ||
            (orbital1 == dzz && orbital2 == dxy) ||
            (orbital2 == dzz && orbital1 == dxy) ){
      value = atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*15.0;
   }
   else if( (orbital1 == dzz && orbital2 == dzx) ||
            (orbital2 == dzz && orbital1 == dzx) ||
            (orbital1 == dzz && orbital2 == dyz) ||
            (orbital2 == dzz && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()
             +atom->GetZindoF4ddLower()*30.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dxy) ||
            (orbital2 == dxxyy && orbital1 == dxy) ){
      value = atom->GetZindoF4ddLower()*35.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzx) ||
            (orbital2 == dxxyy && orbital1 == dzx) ||
            (orbital1 == dxxyy && orbital2 == dyz) ||
            (orbital2 == dxxyy && orbital1 == dyz) ||
            (orbital1 == dxy && orbital2 == dzx) ||
            (orbital2 == dxy && orbital1 == dzx) ||
            (orbital1 == dxy && orbital2 == dyz) ||
            (orbital2 == dxy && orbital1 == dyz) ||
            (orbital1 == dzx && orbital2 == dyz) ||
            (orbital2 == dzx && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()*3.0
             +atom->GetZindoF4ddLower()*20.0;
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageExchangeInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   

   return value;
}

// ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt(Atom* atomA, OrbitalType orbitalA, 
                                           Atom* atomB, OrbitalType orbitalB){
   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   double gammaAA;
   if(orbitalA == s || 
      orbitalA == px ||
      orbitalA == py ||
      orbitalA == pz ){
      gammaAA = atomA->GetZindoF0ss();
   }
   /*
   else if(orbitalA == dxy ||
           orbitalA == dyz ||
           orbitalA == dzz ||
           orbitalA == dzx ||
           orbitalA == dxxyy ){
      gammaAA = atomA->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalA) << "\n";
      throw MolDSException(ss.str());
   }   

   double gammaBB;
   if(orbitalB == s || 
      orbitalB == px ||
      orbitalB == py ||
      orbitalB == pz ){
      gammaBB = atomB->GetZindoF0ss();
   }
   /*
   else if(orbitalB == dxy ||
           orbitalB == dyz ||
           orbitalB == dzz ||
           orbitalB == dzx ||
           orbitalB == dxxyy ){
      gammaBB = atomB->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomB->GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalB) << "\n";
      throw MolDSException(ss.str());
   }  

   return 1.2/( r+2.4/(gammaAA+gammaBB) );

}

void ZindoS::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, Atom* atomA, Atom* atomB){

   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);

   // see (4f) in [AEZ_1986]
   diatomicOverlap[pz][pz] *= 1.267;
   diatomicOverlap[py][py] *= 0.585;
   diatomicOverlap[px][px] *= 0.585;

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         printf("diatomicOverlap[%d][%d]=%lf\n",i,j,diatomicOverlap[i][j]);
      }
   }
   */


}

// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double ZindoS::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                          Molecule* molecule, double** fockMatrix, double** gammaAB){
   double value = 0.0;
   Atom* atomA = NULL;
   Atom* atomB = NULL;
   int firstAOIndexA;
   int firstAOIndexB;
   int numberAOsA;
   int numberAOsB;
   double gamma;
   double exchange;
   double coulomb;
   OrbitalType orbitalMu;
   OrbitalType orbitalNu;

   for(int A=0; A<molecule->GetAtomVect()->size(); A++){
      atomA = (*molecule->GetAtomVect())[A];
      firstAOIndexA = atomA->GetFirstAOIndex();
      numberAOsA = atomA->GetValence().size();

      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         orbitalMu = atomA->GetValence()[mu-firstAOIndexA];

         // CNDO term
         for(int B=0; B<molecule->GetAtomVect()->size(); B++){
            atomB = (*molecule->GetAtomVect())[B];
            firstAOIndexB = atomB->GetFirstAOIndex();
            numberAOsB = atomB->GetValence().size();

            for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
               orbitalNu = atomB->GetValence()[nu-firstAOIndexB];

               if(A==B){
                  gamma = atomA->GetZindoF0ss();
               }
               else{
                  gamma = this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
               }  

               value += gamma*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];
            }
         }

         // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
         for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
            orbitalNu = atomA->GetValence()[nu-firstAOIndexA];

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
               value += exchange*fockMatrix[moI][mu]*fockMatrix[moJ][nu]*fockMatrix[moK][nu]*fockMatrix[moL][mu];
               value += exchange*fockMatrix[moI][mu]*fockMatrix[moJ][nu]*fockMatrix[moK][mu]*fockMatrix[moL][nu];
            }

            coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

            if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                  gamma = atomA->GetZindoF0ss();
            }
            else{
               stringstream ss;
               ss << this->errorMessageMolecularIntegralElement;
               ss << this->errorMessageAtomType << AtomTypeStr(atomA->GetAtomType()) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
               throw MolDSException(ss.str());
            }   

            value += (coulomb-gamma)*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];

         }

      }
   }

   return value;
}


void ZindoS::DoesCIS(){
   cout << this->messageStartCIS;

   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberVir = this->molecule->GetTotalNumberAOs() - numberOcc;

   // check the number of active occupied orbitals.
   if(numberOcc < Parameters::GetInstance()->GetActiveOccCIS()){
      Parameters::GetInstance()->SetActiveOccCIS(numberOcc);
   }
   else{
      numberOcc = Parameters::GetInstance()->GetActiveOccCIS();
   }

   // check the number of active virtual orbitals.
   if(numberVir < Parameters::GetInstance()->GetActiveVirCIS()){
      Parameters::GetInstance()->SetActiveVirCIS(numberOcc);
   }
   else{
      numberVir = Parameters::GetInstance()->GetActiveVirCIS();
   }

   // check the number of calculated excited states.
   int numberExcitedStates = numberOcc * numberVir;
   /*
   if(numberExcitedStates < Parameters::GetInstance()->GetNumberExcitedStatesCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(numberOcc);
   }
   else{
      numberExcitedStates = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   }
   */

printf("\n\nnumber of occ orbitals: %d\n",numberOcc);
printf("number of vir orbitals: %d\n",numberVir);
printf("number of excited states: %d\n\n\n",numberExcitedStates);

   double** matrixCIS = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(numberExcitedStates, numberExcitedStates);
   double* excitedEnergies = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(numberExcitedStates);

   double value=0.0;
   int moI;
   int moA;
   int moJ;
   int moB;
   int k=0;
   int l=0;
   for(int i=0; i<numberOcc; i++){
      moI = this->molecule->GetTotalNumberValenceElectrons()/2 - i -1;      
      for(int a=0; a<numberVir; a++){
         moA = this->molecule->GetTotalNumberValenceElectrons()/2 + a;

         l=0;
         for(int j=0; j<numberOcc; j++){
            moJ = this->molecule->GetTotalNumberValenceElectrons()/2 - j -1;      
            for(int b=0; b<numberVir; b++){
               moB = this->molecule->GetTotalNumberValenceElectrons()/2 + b;

               // diagonal term
               if(k==l){
                  value = this->energiesMO[moA] - this->energiesMO[moI] 
                         +2.0*this->GetMolecularIntegralElement(moI, moA, moA, moI, 
                                                this->molecule, this->fockMatrix, NULL)
                         -    this->GetMolecularIntegralElement(moI, moI, moA, moA, 
                                                this->molecule, this->fockMatrix, NULL);

               }
               // off diagonal term (right upper)
               else if(k<l){
                  value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                                this->molecule, this->fockMatrix, NULL)
                         -    this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                this->molecule, this->fockMatrix, NULL);
               }
            
               matrixCIS[k][l] = value;
               l++;
            }
         }
         k++;
      }
   }

   bool calcEigenVectors = true;
   MolDS_mkl_wrapper::LapackWrapper::GetInstance()->Dsyevd(matrixCIS,
                                                           excitedEnergies, 
                                                           numberExcitedStates, 
                                                           calcEigenVectors);

   // output eigen energies
   for(int k=0; k<numberExcitedStates; k++){
      printf("%d-th excited energy: %e\n",k+1, excitedEnergies[k]);
   }
   cout << endl;

   if(matrixCIS != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix2d(matrixCIS, numberExcitedStates);
      matrixCIS = NULL;
      //cout << "matrixCIS deleted\n";
   }
   if(excitedEnergies != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(excitedEnergies);
      excitedEnergies = NULL;
      //cout << "exceitedEnergies deleted\n";
   }

   cout << this->messageDoneCIS;
}


}
#endif



