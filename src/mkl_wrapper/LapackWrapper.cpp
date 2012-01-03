#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<stdexcept>
#include"mkl.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"LapackWrapper.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_mkl_wrapper{
LapackWrapper* LapackWrapper::lapackWrapper = NULL;

LapackWrapper::LapackWrapper(){
   this->calculatedDsysvBlockSize = false;
   this->dsysvBlockSize = 64;
   this->errorMessageDsyevdInfo = "Error in mkl_wrapper::LapackWrapper::Dsyevd: info != 0: info = ";
   this->errorMessageDsyevdSize = "Error in mkl_wrapper::LapackWrapper::Dsyevd: size of matirx < 1\n";
   this->errorMessageDsysvInfo = "Error in mkl_wrapper::LapackWrapper::Dsysv: info != 0\n";
   this->errorMessageDsysvSize = "Error in mkl_wrapper::LapackWrapper::Dsysv: size of matirx < 1\n";
}

LapackWrapper::~LapackWrapper(){
}

LapackWrapper* LapackWrapper::GetInstance(){
   if(lapackWrapper == NULL){
      lapackWrapper = new LapackWrapper();
      //cout << "LapackWrapper created.\n\n" << endl;
   }
   return lapackWrapper;
}

void LapackWrapper::DeleteInstance(){
   if(lapackWrapper != NULL){
      delete lapackWrapper;
      //cout << "LapackWrapper deleted\n\n" << endl;
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
   double* tempEigenValues;
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
   work = (double*)mkl_malloc( sizeof(double)*lwork, 16 );
   iwork = (int*)mkl_malloc( sizeof(int)*liwork, 16 );
   convertedMatrix = (double*)mkl_malloc( sizeof(double)*size*size, 16 );
   tempEigenValues = (double*)mkl_malloc( sizeof(double)*size, 16 );

   for(int i = 0; i < size; i++){
      for(int j = i; j < size; j++){
         convertedMatrix[i+j*size] = matrix[i][j];
      }
   }

   // call Lapack
   dsyevd(&job, &uplo, &size, convertedMatrix, &lda, tempEigenValues, work, &lwork, iwork, &liwork, &info);

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

   for(int i = 0; i < size; i++){
      eigenValues[i] = tempEigenValues[i];
   }
   //printf("size=%d lwork=%d liwork=%d k=%d info=%d\n",size,lwork,liwork,k,info);

   // free
   mkl_free(work);
   mkl_free(iwork);
   mkl_free(convertedMatrix);
   mkl_free(tempEigenValues);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDsyevdInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}

/***
 *
 * Slove matrix*X=b, then we get X by this method.
 * The X is stored in b.
 *
 */
int LapackWrapper::Dsysv(double const* const* matrix, double* b, int size){
   int info = 0;
   int lwork;
   char uplo = 'U';
   int lda = size;
   int ldb = size;
   int nrhs = 1;
   double* convertedMatrix;
   double* work;
   double* tempB;
   int* ipiv;

   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDsysvSize;
      throw MolDSException(ss.str());
   }

   // malloc
   ipiv = (int*)mkl_malloc( sizeof(int)*2*size, 16 );
   convertedMatrix = (double*)mkl_malloc( sizeof(double)*size*size, 16 );
   tempB = (double*)mkl_malloc( sizeof(double)*size, 16 );

   for(int i = 0; i < size; i++){
      for(int j = i; j < size; j++){
         convertedMatrix[i+j*size] = matrix[i][j];
      }
   }
   for(int i = 0; i < size; i++){
      tempB[i] = b[i];
   }

   // calc. lwork
   #pragma omp critical
   {
      if(!this->calculatedDsysvBlockSize){
         lwork = -1;
         double tempWork[3]={0.0, 0.0, 0.0};
         dsysv(&uplo, &size, &nrhs, convertedMatrix, &lda, ipiv, tempB, &ldb, tempWork, &lwork, &info);
         this->calculatedDsysvBlockSize = true;
         this->dsysvBlockSize = tempWork[0]/size;
      }
   }
   info = 0;
   lwork = this->dsysvBlockSize*size;
   work = (double*)mkl_malloc( sizeof(double)*lwork, 16 );

   // call Lapack
   dsysv(&uplo, &size, &nrhs, convertedMatrix, &lda, ipiv, tempB, &ldb, work, &lwork, &info);
   for(int i = 0; i < size; i++){
      b[i] = tempB[i];
   }

   // free
   mkl_free(convertedMatrix);
   mkl_free(ipiv);
   mkl_free(work);
   mkl_free(tempB);
  
   if(info != 0){
      cout << "info=" << info << endl;
      stringstream ss;
      ss << errorMessageDsysvInfo;
      throw MolDSException(ss.str());
   }
   return info;
}
}
