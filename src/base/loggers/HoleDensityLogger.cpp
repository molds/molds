//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<string.h>
#include<math.h>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/format.hpp>
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../Uncopyable.h"
#include"../Utilities.h"
#include"../Enums.h"
#include"../EularAngle.h"
#include"../Parameters.h"
#include"../atoms/Atom.h"
#include"../Molecule.h"
#include"HoleDensityLogger.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_loggers{

HoleDensityLogger::HoleDensityLogger(const Molecule& molecule, 
                       double const* const* fockMatrix, 
                       double const* const* cisMatrix, 
                   TheoryType theory){
   this->molecule = &molecule;
   this->fockMatrix = fockMatrix;
   this->cisMatrix = cisMatrix;
   this->theory = theory;
   this->SetMessage();
}

HoleDensityLogger::HoleDensityLogger(){
}

void HoleDensityLogger::SetMessage(){
   this->errorMessageCISMatrixNULL
      = "Error in base::HolePlot::DrawDensity: CIS Matrix is NULL.\n";
   this->errorMessageFockMatrixNULL
      = "Error in base::HolePlot::DrawDensity: Fock Matrix is NULL.\n";
   this->messageCubeHeaderComment1 = "MolDS cube file (in atomic units) for Hole density.\n";
   this->messageCubeHeaderComment2 = "outer loop:x, middle loop:y, inner loop:z\n";
   this->messageSkippedElecStateIndex = "\t\tBad electronic state is skipped. The skipped electronic state: ";
   this->messageStartDensityPlot = "\t== START: Hole density plot ==\n";
   this->messageEndDensityPlot = "\t== DONE: Hole density plot ==\n\n";
   this->messageOmpElapsedTimeDensityPlot = "\tElapsed time(omp) for the Hole density plot = ";
   this->messageUnitSec = "[s].";
   this->stringCubeExtension = ".cube";
}

void HoleDensityLogger::DrawDensity(int elecStateIndex) const{
   vector<int> elecStateIndeces;
   elecStateIndeces.push_back(elecStateIndex);
   this->DrawDensity(elecStateIndeces);
}

void HoleDensityLogger::CalcOrigin(double* origin) const{
   for(int i=0; i<CartesianType_end; i++){
      origin[i] = this->molecule->GetXyzCOC()[i];
      origin[i] -= 0.5*Parameters::GetInstance()->GetFrameLengthHolePlot()[i];
   }
}

void HoleDensityLogger::DrawDensity(vector<int> elecStateIndeces) const{
   this->MatricesNullCheck();
   this->OutputLog(this->messageStartDensityPlot);
   double ompStartTime = omp_get_wtime();

   // set frame basics
   double dx=0.0, dy=0.0, dz=0.0;
   double origin[CartesianType_end] = {0.0, 0.0, 0.0};
   this->CalcGridDisplacement(&dx, &dy, &dz);
   this->CalcOrigin(origin);

   // Hole density output 
   for(int n=0; n<elecStateIndeces.size(); n++){
      // validate electronic state
      int groundState = 0;
      if(Parameters::GetInstance()->GetNumberExcitedStatesCIS() < elecStateIndeces[n] || 
         groundState == elecStateIndeces[n]){
         this->OutputLog((boost::format("%s%d\n") % this->messageSkippedElecStateIndex.c_str() 
                                                  % elecStateIndeces[n]).str()) ;
         continue;
      }

      // open the cube file
      string fileName = this->GetFileName(elecStateIndeces[n]);
      ofstream ofs(fileName.c_str());

      // output feader and molecule to the cube file
      this->OutputHeaderToFile(ofs, origin, dx, dy, dz);
      this->OutputMoleculeToFile(ofs, *this->molecule);

      // output grid data to the cube file
      int lineBreakCounter=0;
      for(int ix=0; ix<Parameters::GetInstance()->GetGridNumberHolePlot()[XAxis]; ix++){
         double x = origin[XAxis] + dx*(double)ix;
         for(int iy=0; iy<Parameters::GetInstance()->GetGridNumberHolePlot()[YAxis]; iy++){
            double y = origin[YAxis] + dy*(double)iy;
            for(int iz=0; iz<Parameters::GetInstance()->GetGridNumberHolePlot()[ZAxis]; iz++){
               double z = origin[ZAxis] + dz*(double)iz;

               double density = this->GetDensityValue(elecStateIndeces[n], x, y, z);
               ofs << (boost::format("\t%e") % density ).str();
               lineBreakCounter++;
               if(lineBreakCounter%6==0){
                  ofs << endl;
                  lineBreakCounter=0;
               }

            }
         }
      }
   }
   double ompEndTime = omp_get_wtime();
   this->OutputLog((boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeDensityPlot.c_str()
                                                 % (ompEndTime - ompStartTime)
                                                 % this->messageUnitSec.c_str()
                                                 % this->messageEndDensityPlot.c_str()).str());
}

double HoleDensityLogger::GetDensityValue(int elecStateIndex, double x, double y, double z)  const{
   double density = 0.0;
   int excitedState = elecStateIndex-1;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto) reduction(+:density)
   for(int i=0; i<numberActiveOcc; i++){
      try{
         int moI = numberActiveOcc - (i+1);
         for(int j=0; j<numberActiveOcc; j++){
            int moJ = numberActiveOcc - (j+1);
            for(int a=0; a<numberActiveVir; a++){
               int slaterDeterminatIndexIA = i*numberActiveVir + a;
               int slaterDeterminatIndexJA = j*numberActiveVir + a;
               double moIValue = this->GetMOValue(moI, *this->molecule, x, y, z);
               double moJValue = this->GetMOValue(moJ, *this->molecule, x, y, z);
               density += moIValue*this->cisMatrix[excitedState][slaterDeterminatIndexIA]
                         *moJValue*this->cisMatrix[excitedState][slaterDeterminatIndexJA];
            }
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }  
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
   return density;
}

string HoleDensityLogger::GetFileName(int elecStateIndex) const{
   int digit = 5;
   stringstream fileName;
   fileName << Parameters::GetInstance()->GetFileNamePrefixHolePlot();
   fileName << Utilities::Num2String(elecStateIndex,digit);
   fileName << this->stringCubeExtension;
   return fileName.str();
}

void HoleDensityLogger::CalcGridDisplacement(double* dx, double* dy, double* dz) const{
   *dx = Parameters::GetInstance()->GetFrameLengthHolePlot()[XAxis]
        /(double)Parameters::GetInstance()->GetGridNumberHolePlot()[XAxis];
   *dy = Parameters::GetInstance()->GetFrameLengthHolePlot()[YAxis]
        /(double)Parameters::GetInstance()->GetGridNumberHolePlot()[YAxis];
   *dz = Parameters::GetInstance()->GetFrameLengthHolePlot()[ZAxis]
        /(double)Parameters::GetInstance()->GetGridNumberHolePlot()[ZAxis];
}

double HoleDensityLogger::GetMOValue(int moIndex, const MolDS_base::Molecule& molecule, double x, double y, double z) const{
   double moValue = 0.0;
   for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
      Atom* atomA = this->molecule->GetAtom(a);
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValenceSize();
      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         double aoValue = atomA->GetAtomicBasisValue(x,
                                                     y,
                                                     z,
                                                     mu-firstAOIndexA,
                                                     this->theory);
         moValue += fockMatrix[moIndex][mu]*aoValue;
      }
   }
   return moValue;
}

void HoleDensityLogger::OutputHeaderToFile(ofstream& ofs, double const* origin, double dx, double dy, double dz) const{
   int gridNumber[CartesianType_end] = {Parameters::GetInstance()->GetGridNumberHolePlot()[XAxis], 
                                        Parameters::GetInstance()->GetGridNumberHolePlot()[YAxis],
                                        Parameters::GetInstance()->GetGridNumberHolePlot()[ZAxis]};
   char data[1000] = "";
   // output header to the cube file
   ofs << this->messageCubeHeaderComment1;
   ofs << this->messageCubeHeaderComment2;
   sprintf(data,"\t%d\t%e\t%e\t%e\n", (int)this->molecule->GetNumberAtoms(),
                                      origin[XAxis], 
                                      origin[YAxis], 
                                      origin[ZAxis]);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[XAxis], dx, 0.0, 0.0);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[YAxis], 0.0, dy, 0.0);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[ZAxis], 0.0, 0.0, dz);
   ofs << string(data);
}

void HoleDensityLogger::OutputMoleculeToFile(ofstream& ofs, const Molecule& molecule) const{
   char data[1000] = "";
   // output molecule to the cube file
   for(int a=0; a<molecule.GetNumberAtoms(); a++){
      const Atom& atomA = *molecule.GetAtom(a);
      memset(data,0,sizeof(data));
      sprintf(data,"\t%d\t%d\t%e\t%e\t%e\n", atomA.GetAtomType()+1, 
                                       atomA.GetNumberValenceElectrons(),
                                       atomA.GetXyz()[XAxis],
                                       atomA.GetXyz()[YAxis],
                                       atomA.GetXyz()[ZAxis]);
      ofs << string(data);
   }
}

void HoleDensityLogger::MatricesNullCheck() const{
   // NULL check
   if(this->cisMatrix == NULL){
      stringstream ss;
      ss << this->errorMessageCISMatrixNULL;
      throw MolDSException(ss.str());
   }
   if(this->fockMatrix == NULL){
      stringstream ss;
      ss << this->errorMessageFockMatrixNULL;
      throw MolDSException(ss.str());
   }
}

}
