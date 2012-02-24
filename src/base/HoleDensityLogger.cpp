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
#include<boost/format.hpp>
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"Utilities.h"
#include"Enums.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"Molecule.h"
#include"HoleDensityLogger.h"
using namespace std;
using namespace MolDS_base_atoms;
namespace MolDS_base{

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
   this->messageStartHoleDensityPlot = "\t== START: Hole density plot ==\n";
   this->messageEndHoleDensityPlot = "\t== DONE: Hole density plot ==\n\n";
   this->messageSkippedElecStateIndex = "\t\tBad electronic state is skipped. The skipped electronic state: ";
   this->stringCubeExtension = ".cube";
}

void HoleDensityLogger::DrawDensity(int elecStateIndex){
   vector<int> elecStateIndeces;
   elecStateIndeces.push_back(elecStateIndex);
   this->DrawDensity(elecStateIndeces);
}

void HoleDensityLogger::DrawDensity(vector<int> elecStateIndeces){
   this->MatricesNullCheck();
   this->OutputLog(this->messageStartHoleDensityPlot);
   Parameters* parameters = Parameters::GetInstance();
   int digit = 5;
   int gridNumber[CartesianType_end] = {parameters->GetGridNumberHolePlot()[XAxis], 
                                        parameters->GetGridNumberHolePlot()[YAxis],
                                        parameters->GetGridNumberHolePlot()[ZAxis]};
   double frameLength[CartesianType_end] = {parameters->GetFrameLengthHolePlot()[XAxis],
                                            parameters->GetFrameLengthHolePlot()[YAxis],
                                            parameters->GetFrameLengthHolePlot()[ZAxis]};
   double dx = frameLength[XAxis]/(double)gridNumber[XAxis];
   double dy = frameLength[YAxis]/(double)gridNumber[YAxis];
   double dz = frameLength[ZAxis]/(double)gridNumber[ZAxis];
   double origin[CartesianType_end] = {this->molecule->GetXyzCOC()[XAxis],
                                       this->molecule->GetXyzCOC()[YAxis],
                                       this->molecule->GetXyzCOC()[ZAxis]};
   for(int i=0; i<CartesianType_end; i++){
      origin[i] -= 0.5*frameLength[i];
   }

   int numberActiveOcc = parameters->GetActiveOccCIS();
   int numberActiveVir = parameters->GetActiveVirCIS();
   int numberExcitedStates = parameters->GetNumberExcitedStatesCIS();

   // file name
   string fileNamePrefix = parameters->GetFileNamePrefixHolePlot();
   char data[1000] = "";

   // Hole density output 
   int groundState = 0;
   for(int n=0; n<elecStateIndeces.size(); n++){
      if(numberExcitedStates < elecStateIndeces[n] || groundState == elecStateIndeces[n]){
         this->OutputLog((boost::format("%s%d\n") % this->messageSkippedElecStateIndex.c_str() 
                                                  % elecStateIndeces[n]).str()) ;
         continue;
      }
      int excitedState = elecStateIndeces[n]-1;

      // file open
      stringstream fileName;
      fileName << fileNamePrefix << Utilities::Num2String(elecStateIndeces[n],digit) << this->stringCubeExtension;
      ofstream ofs(fileName.str().c_str());

      this->OutputHeaderToFile(ofs, origin, gridNumber, dx, dy, dz);
      this->OutputMoleculeToFile(ofs, *this->molecule);

      // output grid data to the cube file
      int lineBreakCounter=0;
      for(int ix=0; ix<gridNumber[XAxis]; ix++){
         double x = origin[XAxis] + dx*(double)ix;
         for(int iy=0; iy<gridNumber[YAxis]; iy++){
            double y = origin[YAxis] + dy*(double)iy;
            for(int iz=0; iz<gridNumber[ZAxis]; iz++){
               double z = origin[ZAxis] + dz*(double)iz;

               double density = 0.0;
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

               memset(data,0,sizeof(data));
               sprintf(data,"\t%e",density);
               ofs << string(data);
               lineBreakCounter++;
               if(lineBreakCounter%6==0){
                  ofs << endl;
                  lineBreakCounter=0;
               }
            }
         }
      }
   }
   this->OutputLog(this->messageEndHoleDensityPlot);
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

void HoleDensityLogger::OutputHeaderToFile(ofstream& ofs, double const* origin, int const* gridNumber, double dx, double dy, double dz) const{
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
