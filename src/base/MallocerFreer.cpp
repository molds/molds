#include<stdio.h>
#include<stdlib.h>
#include<iostream>
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
   this->messageMemoryUsage = "Memory summary related to temporary arraies (Heap reagion in a node).\n";
   this->messageMemoryUsageCurrent = "\tMax malloced: ";
   this->messageMemoryUsageMax = "\tCurrent malloced: ";
   this->messageKByte = " [kb].\n";
}

MallocerFreer::~MallocerFreer(){
   this->OutputMemoryUsage();
}

void MallocerFreer::OutputMemoryUsage() const{
   cout << this->messageMemoryUsage;
   cout << this->messageMemoryUsageCurrent << MallocerFreer::maxMalloced/1000.0 << this->messageKByte;
   cout << this->messageMemoryUsageMax << MallocerFreer::currentMalloced/1000.0 << this->messageKByte;
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
