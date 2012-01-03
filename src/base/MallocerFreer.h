#ifndef INCLUDED_MALLOCERFREER
#define INCLUDED_MALLOCERFREER
namespace MolDS_base{

// MallocerFreer is singleton
class MallocerFreer: private Uncopyable{
public:
   static MallocerFreer* GetInstance();
   static void DeleteInstance();

   // int matrix
   // 1d
   int* MallocIntMatrix1d(int size1);
   void InitializeIntMatrix1d(int* matirx, int size1);
   void FreeIntMatrix1d(int** matrix, int size1);

   // real number matrix
   // 1d
   double* MallocDoubleMatrix1d(int size1);
   void InitializeDoubleMatrix1d(double* matrix, int size1);
   void FreeDoubleMatrix1d(double** matrix, int size1);
   // 2d
   double** MallocDoubleMatrix2d(int size1, int size2);
   void InitializeDoubleMatrix2d(double** matrix, int size1, int size2);
   void FreeDoubleMatrix2d(double*** matrix, int size1, int size2);
   // 3d
   double*** MallocDoubleMatrix3d(int size1, int size2, int size3);
   void InitializeDoubleMatrix3d(double*** matrix, int size1, int size2, int size3);
   void FreeDoubleMatrix3d(double**** matrix, int size1, int size2, int size3);
   // 4d
   double**** MallocDoubleMatrix4d(int size1, int size2, int size3, int size4);
   void InitializeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3, int size4);
   void FreeDoubleMatrix4d(double***** matrix, int size1, int size2, int size3, int size4);
   // 5d
   double***** MallocDoubleMatrix5d(int size1, int size2, int size3, 
                                    int size4, int size5);
   void InitializeDoubleMatrix5d(double***** matrix, int size1, int size2, int size3, 
                                                     int size4, int size5);
   void FreeDoubleMatrix5d(double****** matrix, int size1, int size2, int size3,
                                                 int size4, int size5);
   // 6d
   double****** MallocDoubleMatrix6d(int size1, int size2, int size3, 
                                     int size4, int size5, int size6);
   void InitializeDoubleMatrix6d(double****** matrix, int size1, int size2, int size3, 
                                                      int size4, int size5, int size6);
   void FreeDoubleMatrix6d(double******* matrix, int size1, int size2, int size3,
                                                 int size4, int size5, int size6);
private:
   MallocerFreer();
   ~MallocerFreer();
   static MallocerFreer* mallocerFreer;
   static double currentMalloced;
   static double maxMalloced;
   static void AddCurrentMalloced(double amount);
   static void SubtCurrentMalloced(double amount);
   std::string errorMessageMallocFailure;
   std::string messageMemoryUsage;
   std::string messageMemoryUsageCurrent;
   std::string messageMemoryUsageMax;
   std::string messageKByte;
   void OutputMemoryUsage() const;


};
}
#endif
