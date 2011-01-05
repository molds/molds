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
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   string messageTotalNumberAOs;
   string messageTotalNumberAtoms;
   string messageTotalNumberValenceElectrons;
   string messageConfiguration;
   string messageConfigurationTitleAU;
   string messageConfigurationTitleAng;
   string messageCOM;
   string messageCOMTitleAU;
   string messageCOMTitleAng;
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
   void OutputCOMXyz();
   void OutputTotalNumberAtomsAOsValenceelectrons();
   void OutputConfiguration();
};

Molecule::Molecule(){
   this->atomVect = new vector<Atom*>;
   this->COMXyz = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(3);
   this->wasCalculatedCOMXyz = false;
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageConfiguration = "\tMolecular configration:\n";
   this->messageConfigurationTitleAU = "\t\t| i-th | atom type | x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageConfigurationTitleAng = "\t\t| i-th | atom type | x [angst.] | y[angst.] | z[angst.] |\n";
   this->messageCOM = "\tCenter of Mass:\n";
   this->messageCOMTitleAU = "\t\t| x [a.u.] | y[a.u.] | z[a.u.] |\n";
   this->messageCOMTitleAng = "\t\t| x [angst.] | y[angst.] | z[angst.] |\n";
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
   return this->totalNumberAOs;
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAOs = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAOs);
      this->totalNumberAOs += (*this->atomVect)[i]->GetValence().size();
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

void Molecule::OutputConfiguration(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageConfiguration;
   cout << this->messageConfigurationTitleAng;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0]/ang2AU, atom->GetXyz()[1]/ang2AU, atom->GetXyz()[2]/ang2AU);
   }
   cout << "\n";

   cout << this->messageConfigurationTitleAU;
   for(int a=0; a<this->atomVect->size(); a++){
      Atom* atom = (*this->atomVect)[a];
      printf("\t\t%d\t%s\t%e\t%e\t%e\n",a,AtomTypeStr(atom->GetAtomType()),
            atom->GetXyz()[0], atom->GetXyz()[1], atom->GetXyz()[2]);
   }
   cout << "\n";

}

void Molecule::OutputCOMXyz(){
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   cout << this->messageCOM;
   cout << this->messageCOMTitleAng;
   printf("\t\t%e\t%e\t%e\n",this->COMXyz[0]/ang2AU,
                             this->COMXyz[1]/ang2AU,
                             this->COMXyz[2]/ang2AU);
   cout << "\n";

   cout << this->messageCOMTitleAU;
   printf("\t\t%e\t%e\t%e\n",this->COMXyz[0],
                             this->COMXyz[1],
                             this->COMXyz[2]);
   cout << "\n";

}

void Molecule::OutputTotalNumberAtomsAOsValenceelectrons(){
   cout << this->messageTotalNumberAtoms << this->atomVect->size() << "\n";
   cout << this->messageTotalNumberAOs << this->totalNumberAOs << "\n";
   cout << this->messageTotalNumberValenceElectrons << this->totalNumberValenceElectrons << "\n\n";
}

}
#endif





