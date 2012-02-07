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
#include"MOLogger.h"
using namespace std;
using namespace MolDS_base_atoms;
namespace MolDS_base{

MOLogger::MOLogger(const Molecule& molecule, 
                   double const* const* fockMatrix, 
                   TheoryType theory){
   this->molecule = &molecule;
   this->fockMatrix = fockMatrix;
   this->theory = theory;
   this->SetMessage();
}

MOLogger::MOLogger(){
}

void MOLogger::SetMessage(){
   this->stringCubeExtension = ".cube";
   this->messageCubeHeaderComment1 = "MolDS cube file (in atomic units).\n";
   this->messageCubeHeaderComment2 = "outer loop:x, middle loop:y, inner loop:z\n";
   this->messageStartMOPlot = "\t== START: MO Plot ==\n";
   this->messageEndMOPlot = "\t== DONE: MO Plot ==\n\n";
   this->messageSkippedMOIndex = "\t\tBad MO-index is skipped. The skipped MO-index: ";
}

void MOLogger::DrawMO(int moIndex){
   vector<int> moIndeces;
   moIndeces.push_back(moIndex);
   this->DrawMO(moIndeces);
}

void MOLogger::DrawMO(vector<int> moIndeces){
   this->OutputLog(this->messageStartMOPlot);
   Parameters* parameters = Parameters::GetInstance();
   int digit = 5;
   int gridNumber[CartesianType_end] = {parameters->GetGridNumberMOPlot()[XAxis], 
                                        parameters->GetGridNumberMOPlot()[YAxis],
                                        parameters->GetGridNumberMOPlot()[ZAxis]};
   double frameLength[CartesianType_end] = {parameters->GetFrameLengthMOPlot()[XAxis],
                                            parameters->GetFrameLengthMOPlot()[YAxis],
                                            parameters->GetFrameLengthMOPlot()[ZAxis]};
   double dx = frameLength[XAxis]/(double)gridNumber[XAxis];
   double dy = frameLength[YAxis]/(double)gridNumber[YAxis];
   double dz = frameLength[ZAxis]/(double)gridNumber[ZAxis];
   double origin[CartesianType_end] = {this->molecule->GetXyzCOC()[XAxis],
                                       this->molecule->GetXyzCOC()[YAxis],
                                       this->molecule->GetXyzCOC()[ZAxis]};
   for(int i=0; i<CartesianType_end; i++){
      origin[i] -= 0.5*frameLength[i];
   }
   // file name
   string fileNamePrefix = parameters->GetFileNamePrefixMOPlot();
   char data[1000] = "";
   // MO output 
   for(int i=0; i<moIndeces.size(); i++){
      if(this->molecule->GetTotalNumberAOs() <= moIndeces[i]){
         this->OutputLog((boost::format("%s%d\n") % this->messageSkippedMOIndex.c_str() % moIndeces[i]).str()) ;
         continue;
      }
      // file open
      stringstream fileName;
      fileName << fileNamePrefix;
      fileName << Utilities::Num2String(moIndeces[i],digit);
      fileName << this->stringCubeExtension;
      ofstream ofs(fileName.str().c_str());

      // output cube file
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

      for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
         Atom* atomA = this->molecule->GetAtom(a);
         memset(data,0,sizeof(data));
         sprintf(data,"\t%d\t%d\t%e\t%e\t%e\n", atomA->GetAtomType()+1, 
                                          atomA->GetNumberValenceElectrons(),
                                          atomA->GetXyz()[XAxis],
                                          atomA->GetXyz()[YAxis],
                                          atomA->GetXyz()[ZAxis]);
         ofs << string(data);
      }
      int lineBreakCounter=0;
      for(int ix=0; ix<gridNumber[XAxis]; ix++){
         double x = origin[XAxis] + dx*(double)ix;
         for(int iy=0; iy<gridNumber[YAxis]; iy++){
            double y = origin[YAxis] + dy*(double)iy;
            for(int iz=0; iz<gridNumber[ZAxis]; iz++){
               double z = origin[ZAxis] + dz*(double)iz;

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
                     moValue += fockMatrix[moIndeces[i]][mu]*aoValue;
                  }
               }
               memset(data,0,sizeof(data));
               sprintf(data,"\t%e",moValue);
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
   this->OutputLog(this->messageEndMOPlot);
}

}
