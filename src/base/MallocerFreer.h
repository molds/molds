#ifndef INCLUDED_MALLOCERFREER
#define INCLUDED_MALLOCERFREER

#include<stdio.h>
#include<iostream>
#include<math.h>
#include <stdlib.h>

using namespace std;

namespace MolDS_base{

// MallocerFreer is singleton
class MallocerFreer{
private:
   MallocerFreer();
   MallocerFreer(MallocerFreer&);
   void operator = (MallocerFreer&);
   ~MallocerFreer();
   static MallocerFreer* mallocerFreer;
public:
   static MallocerFreer* GetInstance();
   static void DeleteInstance();

   // int matrix
   // 1d
   int* MallocIntMatrix1d(int size1);
   void InitializeIntMatrix1d(int* matirx, int size1);
   void FreeIntMatrix1d(int* matrix);

   // real number matrix
   // 1d
   double* MallocDoubleMatrix1d(int size1);
   void InitializeDoubleMatrix1d(double* matrix, int size1);
   void FreeDoubleMatrix1d(double* matrix);
   // 2d
   double** MallocDoubleMatrix2d(int size1, int size2);
   void InitializeDoubleMatrix2d(double** matrix, int size1, int size2);
   void FreeDoubleMatrix2d(double** matrix, int size1);
   // 4d
   double**** MallocDoubleMatrix4d(int size1, int size2, int size3, int size4);
   void InitializeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3, int size4);
   void FreeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3);

};
MallocerFreer* MallocerFreer::mallocerFreer = NULL;

MallocerFreer::MallocerFreer(){
}

MallocerFreer::~MallocerFreer(){
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

int* MallocerFreer::MallocIntMatrix1d(int size1){
   
   int* matrix;  
   matrix = new int[size1];

   this->InitializeIntMatrix1d(matrix, size1);

   return(matrix); 
}

void MallocerFreer::InitializeIntMatrix1d(int* matrix, int size1){
   for(int i=0;i<size1;i++){
      matrix[i] = 0;
   }
}

void MallocerFreer::FreeIntMatrix1d(int* matrix){
   delete [] matrix;
}

double* MallocerFreer::MallocDoubleMatrix1d(int size1){
   
   double* matrix;  
   matrix = new double[size1];

   this->InitializeDoubleMatrix1d(matrix, size1);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix1d(double* matrix, int size1){
   for(int i=0;i<size1;i++){
      matrix[i] = 0.0;
   }
}

void MallocerFreer::FreeDoubleMatrix1d(double* matrix){
   delete [] matrix;
}

double** MallocerFreer::MallocDoubleMatrix2d(int size1, int size2){
   
   double** matrix;  

   matrix = new double*[size1];
   if(matrix==NULL){
      exit(EXIT_FAILURE); // malloc failure
   }

   for(int i=0;i<size1;i++) {
      matrix[i] = new double[size2];
      if (matrix[i]==NULL){
         exit(EXIT_FAILURE); //malloc failure
      }
   }

   this->InitializeDoubleMatrix2d(matrix, size1, size2);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix2d(double** matrix, int size1, int size2){
   for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
         matrix[i][j] = 0.0;
      }
   }
}

void MallocerFreer::FreeDoubleMatrix2d(double** matrix, int size1){

   int i=0;
   for(i=0;i<size1;i++){
      delete [] matrix[i];
   }
   delete [] matrix;
}

double**** MallocerFreer::MallocDoubleMatrix4d(int size1, int size2, int size3, int size4){
   
   double**** matrix;  

   matrix = new double***[size1];
   if(matrix==NULL){
      exit(EXIT_FAILURE); //malloc failure
   }
   for(int i=0;i<size1;i++) {
      matrix[i] = new double**[size2];
      if(matrix[i]==NULL){
         exit(EXIT_FAILURE); //malloc failure
      }
      for(int j=0;j<size2;j++){
         matrix[i][j] = new double*[size3];
         if(matrix[i][j]==NULL){
           exit(EXIT_FAILURE); //malloc failure
         }
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = new double[size4];
            if(matrix[i][j][k]==NULL){
               exit(EXIT_FAILURE); //malloc failure
            }
         }
      }
   }

   this->InitializeDoubleMatrix4d(matrix, size1, size2, size3, size4);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3, int size4){
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

void MallocerFreer::FreeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3){

   int i=0, j=0, k=0;

   for (i=0;i<size1;i++) {
      for (j=0;j<size2;j++) {
         for (k=0;k<size3;k++) {
            delete [] matrix[i][j][k];
         }
         delete [] matrix[i][j];
      }
      delete [] matrix[i];
   }
   delete [] matrix;

}

}
#endif
