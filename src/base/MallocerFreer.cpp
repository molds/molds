#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<stdexcept>
#include"MolDSException.h"
#include"MallocerFreer.h"
using namespace std;

namespace MolDS_base{

MallocerFreer* MallocerFreer::mallocerFreer = NULL;

MallocerFreer::MallocerFreer(){
   this->errorMessageMallocFailure = "Error in base::MallocFreere: Malloc failure...\n";
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

void MallocerFreer::FreeIntMatrix1d(int** matrix){
   delete [] *matrix;
   *matrix = NULL;
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

void MallocerFreer::FreeDoubleMatrix1d(double** matrix){
   delete [] *matrix;
   *matrix = NULL;
}

double** MallocerFreer::MallocDoubleMatrix2d(int size1, int size2){
   
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

void MallocerFreer::FreeDoubleMatrix2d(double*** matrix, int size1){

   int i=0;
   for(i=0;i<size1;i++){
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   *matrix = NULL;
}

double*** MallocerFreer::MallocDoubleMatrix3d(int size1, int size2, int size3){
   
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

   this->InitializeDoubleMatrix3d(matrix, size1, size2, size3);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix3d(double*** matrix, int size1, int size2, int size3){
   for(int i=0;i<size1;i++) {
      for(int j=0;j<size2;j++){
         for(int k=0;k<size3;k++){
            matrix[i][j][k] = 0.0;
         }
      }
   }
}

void MallocerFreer::FreeDoubleMatrix3d(double**** matrix, int size1, int size2){

   int i=0, j=0;

   for (i=0;i<size1;i++) {
      for (j=0;j<size2;j++) {
         delete [] (*matrix)[i][j];
      }
      delete [] (*matrix)[i];
   }
   delete [] *matrix;
   *matrix = NULL;

}

double**** MallocerFreer::MallocDoubleMatrix4d(int size1, int size2, int size3, int size4){
   
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

void MallocerFreer::FreeDoubleMatrix4d(double***** matrix, int size1, int size2, int size3){

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
   *matrix = NULL;

}

double***** MallocerFreer::MallocDoubleMatrix5d(int size1, int size2, int size3, 
                                                int size4, int size5){
   
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

   this->InitializeDoubleMatrix5d(matrix, size1, size2, size3, size4, size5);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix5d(double***** matrix, int size1, int size2, int size3, 
                                                                 int size4, int size5){
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
                                                             int size4){

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
   *matrix = NULL;

}

double****** MallocerFreer::MallocDoubleMatrix6d(int size1, int size2, int size3, 
                                                 int size4, int size5, int size6){
   
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

   this->InitializeDoubleMatrix6d(matrix, size1, size2, size3, size4, size5, size6);

   return(matrix); 
}

void MallocerFreer::InitializeDoubleMatrix6d(double****** matrix, int size1, int size2, int size3, 
                                                                  int size4, int size5, int size6){
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
                                                             int size4, int size5){

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
   *matrix = NULL;

}

}
