#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include"atoms/Atom.h"

using namespace std;
using namespace MolDS_base_atoms;

namespace MolDS_base{

class Molecule{
private:
   vector<Atom*>* atomVect;
   double* COMXyz;
   bool wasCalculatedCOMXyz;
   int totalNumberAO;
   int totalNumberValenceElectrons;
public:
   Molecule();
   ~Molecule();
   vector<Atom*>* GetAtomVect(); 
   double* GetCOMXyz();
   void CalcCOMXyz();
   int GetTotalNumberAOs();
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons();
   void CalcTotalNumberValenceElectrons();
};

Molecule::Molecule(){
   this->atomVect = new vector<Atom*>;
   this->COMXyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->wasCalculatedCOMXyz = false;
}

Molecule::~Molecule(){
   if(this->atomVect != NULL){
      for(int i=0; i<this->atomVect->size(); i++){
         delete (*atomVect)[i];
         (*atomVect)[i] = NULL;
      }
      delete this->atomVect;
      this->atomVect = NULL;
      //cout << "atomVect deleted\n";
   }
   if(this->COMXyz != NULL){
      MallocerFreer::GetInstance()->FreeDoubleMatrix1d(this->COMXyz);
      this->COMXyz = NULL;
      //cout << "COMXyz deleted\n";
   }
}

vector<Atom*>* Molecule::GetAtomVect(){
   return this->atomVect;
}

double* Molecule::GetCOMXyz(){
   if(!this->wasCalculatedCOMXyz){
      this->CalcCOMXyz();
   }
   return this->COMXyz;
}

void Molecule::CalcCOMXyz(){
   if(!this->wasCalculatedCOMXyz){
      double totalAtomicMass;
      Atom* atom;
      double* atomicXyz;
      double atomicMass;
      for(int i=0; i<this->atomVect->size(); i++){
         atom = (*this->atomVect)[i]; 
         atomicXyz = atom->GetXyz();
         atomicMass = atom->GetAtomicMass();
         totalAtomicMass += atomicMass;
         for(int j=0; j<3; j++){
            this->COMXyz[j] += atomicXyz[j] * atomicMass;
         }
      }
      for(int i=0; i<3; i++){
         this->COMXyz[i]/=totalAtomicMass;
      }
   }
   this->wasCalculatedCOMXyz = true;
}

int Molecule::GetTotalNumberAOs(){
   return this->totalNumberAO;
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAO = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAO);
      this->totalNumberAO += (*this->atomVect)[i]->GetValence().size();
   }
}

int Molecule::GetTotalNumberValenceElectrons(){
   return this->totalNumberValenceElectrons;
}

void Molecule::CalcTotalNumberValenceElectrons(){
   this->totalNumberValenceElectrons = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      this->totalNumberValenceElectrons += (*this->atomVect)[i]->GetNumberValenceElectrons();
   }
}




}
#endif





