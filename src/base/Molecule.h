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
   int totalNumberAO;
   int totalNumberValenceElectrons;
public:
   Molecule();
   ~Molecule();
   vector<Atom*>* GetAtomVect(); 
   int GetTotalNumberAOs();
   void CalcTotalNumberAOs();
   int GetTotalNumberValenceElectrons();
   void CalcTotalNumberValenceElectrons();
};

Molecule::Molecule(){
   atomVect = new vector<Atom*>;
}

Molecule::~Molecule(){
   if(atomVect != NULL){
      for(int i=0; i<atomVect->size(); i++){
         delete (*atomVect)[i];
         (*atomVect)[i] = NULL;
      }
      delete atomVect;
      atomVect = NULL;
      //cout << "mol deleted\n";
   }
}

vector<Atom*>* Molecule::GetAtomVect(){
   return this->atomVect;
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





