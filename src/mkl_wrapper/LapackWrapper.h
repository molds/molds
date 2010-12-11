#ifndef INCLUDED_LAPACKWRAPPER
#define INCLUDED_LAPACKWRAPPER
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<string>
#include<stdlib.h>
#include"mkl.h"

using namespace std;
using namespace MolDS_base;

namespace MolDS_mkl_wrapper{

// LapackWrapper is singleton
class LapackWrapper{
private:
   LapackWrapper();
   LapackWrapper(LapackWrapper&);
   void operator = (LapackWrapper&);
   ~LapackWrapper();
   static LapackWrapper* lapackWrapper;

   string errorMessageDsyevdInfo;
   string errorMessageDsyevdSize;
public:
   static LapackWrapper* GetInstance();
   static void DeleteInstance();
   int Dsyevd(double** matrix, double* eigenValues, int size, bool calcEigenVectors);
};
LapackWrapper* LapackWrapper::lapackWrapper = NULL;

LapackWrapper::LapackWrapper(){
   this->errorMessageDsyevdInfo = "Error in mkl_wrapper::LapackWrapper::Dsyevd: info != 0\n";
   this->errorMessageDsyevdSize = "Error in mkl_wrapper::LapackWrapper::Dsyevd: size of matirx < 1\n";
}

LapackWrapper::~LapackWrapper(){
}

LapackWrapper* LapackWrapper::GetInstance(){
   if(lapackWrapper == NULL){
      lapackWrapper = new LapackWrapper();
   }
   return lapackWrapper;
}

void LapackWrapper::DeleteInstance(){
   if(lapackWrapper != NULL){
      delete lapackWrapper;
   }
   lapackWrapper = NULL;
}


/***
 *
 * return notice
 *    i-th eigen value is eigenValues[i].
 *    i-th eigen vector is (matirx[i][0], matirx[i][1], matirx[i][2], ....).
 *
 * ***/
int LapackWrapper::Dsyevd(double** matrix, double* eigenValues, int size, bool calcEigenVectors){
   int info = 0;
   int k = 0;
   int lwork;
   int liwork;
   char job;
   char uplo = 'U';
   int lda = size;
   double* convertedMatrix;
   double* work;
   int* iwork;

   // set job type
   if(calcEigenVectors){
      job = 'V';
   }
   else{
      job = 'N';
   }

   // calc. lwork and liwork
   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDsyevdSize;
      throw MolDSException(ss.str());
   }
   else if(size == 1){
      lwork = 1;
      liwork = 1;
   }
   else if(1 < size && job == 'N'){
      lwork = 2*size + 1;
      liwork = 2;
   }
   else{
      // calc. k
      double temp = log((double)size)/log(2.0);
      if( (double)((int)temp) < temp ){
         k = (int)temp + 1;
      }
      else{
         k = (int)temp;
      }
      lwork = 3*size*size + (5+2*k)*size + 1;
      liwork = 5*size + 3;
   }


   // malloc
   work = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(lwork);
   iwork = MallocerFreer::GetInstance()->MallocIntMatrix1d(liwork);
   convertedMatrix = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(size*size);

   for(int i = 0; i < size; i++){
      for(int j = i; j < size; j++){
         convertedMatrix[i+j*size] = matrix[i][j];
      }
   }

   // call Lapack
   dsyevd(&job, &uplo, &size, convertedMatrix, &lda, eigenValues, work, &lwork, iwork, &liwork, &info);

   for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
         matrix[i][j] = convertedMatrix[j+i*size];  //i-th row is i-th eigen vector
         //matrix[j][i] = convertedMatrix[j+i*size];  //i-th column is i-th eigen vector
      }
   }

   for(int i=0;i<size;i++){
      if(matrix[i][0]<0){
         for(int j=0;j<size;j++){
            matrix[i][j]*=-1.0;
         }
      }
   }   
   //printf("size=%d lwork=%d liwork=%d k=%d info=%d\n",size,lwork,liwork,k,info);

   // free
   MallocerFreer::GetInstance()->FreeDoubleMatrix1d(convertedMatrix);
   convertedMatrix = NULL;
   MallocerFreer::GetInstance()->FreeDoubleMatrix1d(work);
   work = NULL;
   MallocerFreer::GetInstance()->FreeIntMatrix1d(iwork);
   iwork = NULL;
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDsyevdInfo;
      throw MolDSException(ss.str());
   }
   return info;
}



}
#endif
