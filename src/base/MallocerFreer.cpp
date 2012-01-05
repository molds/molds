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
#include<math.h>
#include<string>
#include<stdexcept>
#include"MolDSException.h"
#include"Uncopyable.h"
#include"MallocerFreer.h"
using namespace std;

namespace MolDS_base{

MallocerFreer* MallocerFreer::mallocerFreer = NULL;
double MallocerFreer::currentMalloced = 0.0;
double MallocerFreer::maxMalloced = 0.0;

MallocerFreer::MallocerFreer(){
   this->errorMessageMallocFailure = "Error in base::MallocFreere: Malloc failure...\n";
   this->messageMemoryUsage = "Summary for memory usage:\n";
   this->messageMemoryCurrentHeap = "\tMax Heap: ";
   this->messageMemoryMaxHeap = "\tCurrent Heap(Leaked): ";
   this->messageMByte = " [MB].\n";
}

MallocerFreer::~MallocerFreer(){
   this->OutputMemoryUsage();
}

void MallocerFreer::OutputMemoryUsage() const{
   cout << this->messageMemoryUsage;
   cout << this->messageMemoryCurrentHeap << MallocerFreer::maxMalloced/pow(10.0,6.0) << this->messageMByte;
   cout << this->messageMemoryMaxHeap << MallocerFreer::currentMalloced/pow(10.0,6.0) << this->messageMByte;
}

MallocerFreer* MallocerFreer::GetInstance(){
   if(mallocerFreer == NULL){
      mallocerFreer = new MallocerFreer();
   }
   return mallocerFreer;
}

void MallocerFreer::DeleteInstance(){
   if(mallocerFreer != NULL){
      delete mallocerFreer;
   }
   mallocerFreer = NULL;
}

void MallocerFreer::AddCurrentMalloced(double amount){
   #pragma omp critical
   {
      MallocerFreer::currentMalloced += amount;
      if(MallocerFreer::maxMalloced < MallocerFreer::currentMalloced){
         MallocerFreer::maxMalloced = MallocerFreer::currentMalloced;
      }
   }
}

void MallocerFreer::SubtCurrentMalloced(double amount){
   #pragma omp critical
   MallocerFreer::currentMalloced -= amount;
}

}
