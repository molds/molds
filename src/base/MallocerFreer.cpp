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
   this->messageMemoryUsage = "Memory summary related to temporary arraies in a node...\n";
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

int* MallocerFreer::MallocIntMatrix1d(int size1) const{
   int* matrix;  
   matrix = new int[size1];
   MallocerFreer::AddCurrentMalloced((double)(size1*sizeof(int)));
   this->InitializeIntMatrix1d(matrix, size1);
   return(matrix); 
}

void MallocerFreer::InitializeIntMatrix1d(int* matrix, int size1) const{
   for(int i=0;i<size1;i++){
      matrix[i] = 0;
   }
}

void MallocerFreer::FreeIntMatrix1d(int** matrix, int size1) const{
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*sizeof(int)));
   *matrix = NULL;
}

double* MallocerFreer::MallocDoubleMatrix1d(int size1) const{
   double* matrix;  
   matrix = new double[size1];
   MallocerFreer::AddCurrentMalloced((double)(size1*sizeof(double)));
   this->InitializeDoubleMatrix1d(matrix, size1);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix1d(double* matrix, int size1) const{
   for(int i=0;i<size1;i++){
      matrix[i] = 0.0;
   }
}

void MallocerFreer::FreeDoubleMatrix1d(double** matrix, int size1) const{
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*sizeof(double)));
   *matrix = NULL;
}

double** MallocerFreer::MallocDoubleMatrix2d(int size1, int size2) const{
   double** matrix;  
   matrix = new double*[size1];
   if(matrix==NULL){
      throw MolDSException(this->errorMessageMallocFailure);
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double[size2];
      if (matrix[i]==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
   }
   MallocerFreer::AddCurrentMalloced((double)(size1*size2*sizeof(double)));
   this->InitializeDoubleMatrix2d(matrix, size1, size2);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix2d(double** matrix, int size1, int size2) const{
   for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
         matrix[i][j] = 0.0;
      }
   }
}

void MallocerFreer::FreeDoubleMatrix2d(double*** matrix, int size1, int size2) const{
   int i=0;
   for(i=0;i<size1;i++){
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*size2*sizeof(double)));
   *matrix = NULL;
}

double*** MallocerFreer::MallocDoubleMatrix3d(int size1, int size2, int size3) const{
   double*** matrix;  
   matrix = new double**[size1];
   if(matrix==NULL){
      throw MolDSException(this->errorMessageMallocFailure);
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double*[size2];
      if(matrix[i]==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int j=0;j<size2;j++){
         matrix[i][j] = new double[size3];
         if(matrix[i][j]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
      }
   }
   MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*sizeof(double)));
   this->InitializeDoubleMatrix3d(matrix, size1, size2, size3);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix3d(double*** matrix, int size1, int size2, int size3) const{
   for(int i=0;i<size1;i++) {
      for(int j=0;j<size2;j++){
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = 0.0;
         }
      }
   }
}

void MallocerFreer::FreeDoubleMatrix3d(double**** matrix, int size1, int size2, int size3) const{
   int i=0, j=0;
   for (i=0;i<size1;i++) {
      for (j=0;j<size2;j++) {
         delete [] (*matrix)[i][j];
      }
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*sizeof(double)));
   *matrix = NULL;
}

double**** MallocerFreer::MallocDoubleMatrix4d(int size1, int size2, int size3, int size4) const{
   double**** matrix;  
   matrix = new double***[size1];
   if(matrix==NULL){
      throw MolDSException(this->errorMessageMallocFailure);
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double**[size2];
      if(matrix[i]==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int j=0;j<size2;j++){
         matrix[i][j] = new double*[size3];
         if(matrix[i][j]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = new double[size4];
            if(matrix[i][j][k]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
         }
      }
   }
   MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*sizeof(double)));
   this->InitializeDoubleMatrix4d(matrix, size1, size2, size3, size4);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3, int size4) const{
   for(int i=0;i<size1;i++) {
      for(int j=0;j<size2;j++){
         for(int k=0;k<size3;k++){
            for(int l=0;l<size4;l++){
               matrix[i][j][k][l] = 0.0;
            }
         }
      }
   }
}

void MallocerFreer::FreeDoubleMatrix4d(double***** matrix, int size1, int size2, int size3, int size4) const{
   int i=0, j=0, k=0;
   for (i=0;i<size1;i++) {
      for (j=0;j<size2;j++) {
         for (k=0;k<size3;k++) {
            delete [] (*matrix)[i][j][k];
         }
         delete [] (*matrix)[i][j];
      }
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*sizeof(double)));
   *matrix = NULL;
}

double***** MallocerFreer::MallocDoubleMatrix5d(int size1, int size2, int size3, 
                                                int size4, int size5) const{
   double***** matrix;  
   matrix = new double****[size1];
   if(matrix==NULL){
      throw MolDSException(this->errorMessageMallocFailure);
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double***[size2];
      if(matrix[i]==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int j=0;j<size2;j++){
         matrix[i][j] = new double**[size3];
         if(matrix[i][j]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = new double*[size4];
            if(matrix[i][j][k]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int l=0;l<size4;l++){
               matrix[i][j][k][l] = new double[size5];
               if(matrix[i][j][k][l]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
            }
         }
      }
   }
   MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*size5*sizeof(double)));
   this->InitializeDoubleMatrix5d(matrix, size1, size2, size3, size4, size5);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix5d(double***** matrix, int size1, int size2, int size3, 
                                                                 int size4, int size5) const{
   for(int i=0;i<size1;i++) {
      for(int j=0;j<size2;j++){
         for(int k=0;k<size3;k++){
            for(int l=0;l<size4;l++){
               for(int m=0;m<size5;m++){
                  matrix[i][j][k][l][m] = 0.0;
               }
            }
         }
      }
   }
}

void MallocerFreer::FreeDoubleMatrix5d(double****** matrix, int size1, int size2, int size3,
                                                            int size4, int size5) const{
   for (int i=0;i<size1;i++) {
      for (int j=0;j<size2;j++) {
         for (int k=0;k<size3;k++) {
            for (int l=0;l<size4;l++) {
               delete [] (*matrix)[i][j][k][l];
            }
            delete [] (*matrix)[i][j][k];
         }
         delete [] (*matrix)[i][j];
      }
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*size5*sizeof(double)));
   *matrix = NULL;
}

double****** MallocerFreer::MallocDoubleMatrix6d(int size1, int size2, int size3, 
                                                 int size4, int size5, int size6) const{
   double****** matrix;  
   matrix = new double*****[size1];
   if(matrix==NULL){
      throw MolDSException(this->errorMessageMallocFailure);
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double****[size2];
      if(matrix[i]==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int j=0;j<size2;j++){
         matrix[i][j] = new double***[size3];
         if(matrix[i][j]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = new double**[size4];
            if(matrix[i][j][k]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int l=0;l<size4;l++){
               matrix[i][j][k][l] = new double*[size5];
               if(matrix[i][j][k][l]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
               for(int m=0;m<size5;m++){
                  matrix[i][j][k][l][m] = new double[size6];
                  if(matrix[i][j][k][l][m]==NULL){
                     throw MolDSException(this->errorMessageMallocFailure);
                  }
               }
            }
         }
      }
   }
   MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*size5*size6*sizeof(double)));
   this->InitializeDoubleMatrix6d(matrix, size1, size2, size3, size4, size5, size6);
   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix6d(double****** matrix, int size1, int size2, int size3, 
                                                                  int size4, int size5, int size6) const{
   for(int i=0;i<size1;i++) {
      for(int j=0;j<size2;j++){
         for(int k=0;k<size3;k++){
            for(int l=0;l<size4;l++){
               for(int m=0;m<size5;m++){
                  for(int n=0;n<size6;n++){
                     matrix[i][j][k][l][m][n] = 0.0;
                  }
               }
            }
         }
      }
   }
}

void MallocerFreer::FreeDoubleMatrix6d(double******* matrix, int size1, int size2, int size3,
                                                             int size4, int size5, int size6) const{
   for (int i=0;i<size1;i++) {
      for (int j=0;j<size2;j++) {
         for (int k=0;k<size3;k++) {
            for (int l=0;l<size4;l++) {
               for (int m=0;m<size5;m++) {
                  delete [] (*matrix)[i][j][k][l][m];
               }
               delete [] (*matrix)[i][j][k][l];
            }
            delete [] (*matrix)[i][j][k];
         }
         delete [] (*matrix)[i][j];
      }
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*size5*size6*sizeof(double)));
   *matrix = NULL;
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
