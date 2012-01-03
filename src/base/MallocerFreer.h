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
   int* MallocIntMatrix1d(int size1) const;
   void InitializeIntMatrix1d(int* matirx, int size1) const;
   void FreeIntMatrix1d(int** matrix, int size1) const;

   // real number matrix
   // 1d
   double* MallocDoubleMatrix1d(int size1) const;
   void InitializeDoubleMatrix1d(double* matrix, int size1) const;
   void FreeDoubleMatrix1d(double** matrix, int size1) const;
   // 2d
   double** MallocDoubleMatrix2d(int size1, int size2) const;
   void InitializeDoubleMatrix2d(double** matrix, int size1, int size2) const;
   void FreeDoubleMatrix2d(double*** matrix, int size1, int size2) const;
   // 3d
   double*** MallocDoubleMatrix3d(int size1, int size2, int size3) const;
   void InitializeDoubleMatrix3d(double*** matrix, int size1, int size2, int size3) const;
   void FreeDoubleMatrix3d(double**** matrix, int size1, int size2, int size3) const;
   // 4d
   double**** MallocDoubleMatrix4d(int size1, int size2, int size3, int size4) const;
   void InitializeDoubleMatrix4d(double**** matrix, int size1, int size2, int size3, int size4) const;
   void FreeDoubleMatrix4d(double***** matrix, int size1, int size2, int size3, int size4) const;
   // 5d
   double***** MallocDoubleMatrix5d(int size1, int size2, int size3, 
                                    int size4, int size5) const;
   void InitializeDoubleMatrix5d(double***** matrix, int size1, int size2, int size3, 
                                                     int size4, int size5) const;
   void FreeDoubleMatrix5d(double****** matrix, int size1, int size2, int size3,
                                                 int size4, int size5) const;
   // 6d
   double****** MallocDoubleMatrix6d(int size1, int size2, int size3, 
                                     int size4, int size5, int size6) const;
   void InitializeDoubleMatrix6d(double****** matrix, int size1, int size2, int size3, 
                                                      int size4, int size5, int size6) const;
   void FreeDoubleMatrix6d(double******* matrix, int size1, int size2, int size3,
                                                 int size4, int size5, int size6) const;
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
