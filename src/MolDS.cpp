#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<string>
#include<time.h>
#include<list>
#include<vector>
#include"base/Enums.h"
#undef INCLUDED_ENUMS
#define RENUMSTR_BODY 1  
#include"base/Enums.h"
#include"base/Molecule.h"
#include"base/atoms/Atom.h"
#include"base/atoms/Hatom.h"
#include"base/atoms/Catom.h"
#include"base/atoms/Liatom.h"
#include"base/MallocerFreer.h"
#include"base/InputParser.h"
#include"base/Utilities.h"
#include"base/EularAngle.h"
#include"base/Parameters.h"
#include"cndo/Cndo2.h"
#include"indo/Indo.h"
#include"mkl_wrapper/LapackWrapper.h"



using namespace std;
using namespace MolDS_base;
using namespace MolDS_mkl_wrapper;


int main(){
   // Welcome Message
   time_t startTime;
   struct tm *ltm;
   char s[50];
   clock_t startTick = clock();
   time(&startTime);
   ltm = localtime(&startTime);
   fmttm(s, ltm);
   cout << "\n\n     >>>>>  Welcome to the MolDS world at " << s << "  <<<<<\n\n\n";

   // declare
   InputParser::GetInstance();
   Molecule* molecule = new Molecule();
   MallocerFreer::GetInstance();
   Parameters::GetInstance();
   LapackWrapper::GetInstance();


   // Parse input
   InputParser::GetInstance()->Parse(molecule);

   // cndo2
   /*
   MolDS_cndo::Cndo2* cndo2 = new MolDS_cndo::Cndo2();
   cndo2->SetMolecule(molecule);
   cndo2->DoesSCF();
   delete cndo2;
   */

   // indo
   MolDS_indo::Indo* indo = new MolDS_indo::Indo();
   indo->SetMolecule(molecule);
   indo->DoesSCF();
   delete indo;



   /*** test lapack ***
   {
   int size = 3;
   double** matrix;
   double* eigenValues;
   eigenValues = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(size);
   matrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(size,size);

   matrix[0][0] = 3.0;  matrix[0][1] = 2.0;  matrix[0][2] = 1.0;
   matrix[1][0] = 2.0;  matrix[1][1] = 4.0;  matrix[1][2] = 1.0;
   matrix[2][0] = 1.0;  matrix[2][1] = 1.0;  matrix[2][2] = 6.0;

   double oldMatrix[3][3];
   for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
         oldMatrix[i][j] = matrix[i][j];
      }
   }
      

   LapackWrapper::GetInstance()->Dsyevd(matrix, eigenValues, size, true);

   double diag[3][3];
   for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
         diag[i][j] = 0.0;
         for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
               diag[i][j] += matrix[i][l] * oldMatrix[l][k] * matrix[j][k];
            }
         }
         printf("diag[%d][%d] = %lf\t",i,j,diag[i][j]);
      }
      cout << "\n";
   }

   cout << "eigen:1 " << eigenValues[0] << "\n";
   cout << matrix[0][0] << "\t" << matrix[0][1] << "\t" << matrix[0][2] << "\n";

   cout << "eigen:2 " << eigenValues[1] << "\n";
   cout << matrix[1][0] << "\t" << matrix[1][1] << "\t" << matrix[1][2] << "\n";

   cout << "eigen:3 " << eigenValues[2] << "\n";
   cout << matrix[2][0] << "\t" << matrix[2][1] << "\t" << matrix[2][2] << "\n";

   MallocerFreer::GetInstance()->FreeDoubleMatrix2d(matrix,size);
   MallocerFreer::GetInstance()->FreeDoubleMatrix1d(eigenValues);
   }
   /*** end( test lapack ) ***/

   /*** test lapack ***
   {
   int size = 8;
   double** matrix;
   double* eigenValues;
   eigenValues = MallocerFreer::GetInstance()->MallocDoubleMatrix1d(size);
   matrix = MallocerFreer::GetInstance()->MallocDoubleMatrix2d(size,size);

   matrix[0][0] = -0.114142;  
   matrix[0][1] = 0.2;
   matrix[0][4] = -0.209441;  
   matrix[0][6] = 0.168121;

   matrix[1][1] = -0.046230;  
   matrix[1][5] = -0.153291;  
   matrix[2][2] = -0.046230;
   matrix[2][4] = 0.168121;  
   matrix[2][6] = -0.054602;  

   matrix[3][3] = -0.046230;
   matrix[3][7] = -0.153291;

   matrix[4][4] = -0.114142;
   matrix[5][5] = -0.046230;
   matrix[6][6] = -0.046230;
   matrix[7][7] = -0.046230;


   double oldMatrix[size][size];
   for(int i=0;i<size;i++){
      for(int j=i;j<size;j++){
         oldMatrix[i][j] = matrix[i][j];
      }
      for(int j=0;j<i;j++){
         oldMatrix[i][j] = matrix[j][i];
      }
   }
      

   LapackWrapper::GetInstance()->Dsyevd(matrix, eigenValues, size, true);

   double diag[size][size];
   for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){
         diag[i][j] = 0.0;
         for(int k=0;k<size;k++){
            for(int l=0;l<size;l++){
               diag[i][j] += matrix[i][l] * oldMatrix[l][k] * matrix[j][k];
            }
         }
         printf("diag[%d][%d] = %lf\t",i,j,diag[i][j]);
      }
      cout << "\n";
   }


   MallocerFreer::GetInstance()->FreeDoubleMatrix2d(matrix,size);
   MallocerFreer::GetInstance()->FreeDoubleMatrix1d(eigenValues);
   }
   /*** end( test lapack ) ***/



   //Free 
   LapackWrapper::DeleteInstance(); 
   Parameters::DeleteInstance();
   MallocerFreer::DeleteInstance();
   delete molecule;
   InputParser::DeleteInstance();


   // Farewell Message
   clock_t endTick = clock();
   double consumedTime = (double)(endTick - startTick)/(double)CLOCKS_PER_SEC;
   cout << "\n\n     >>>>>  The MolDS finished normally!  <<<<<\n";
   cout <<     "     >>>>>  Consumed time (CPU time): " << consumedTime << "[s].  <<<<<\n";
   cout <<     "     >>>>>  See you.  <<<<<\n\n\n";

   return 0;
}













