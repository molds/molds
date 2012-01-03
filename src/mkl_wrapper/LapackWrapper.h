#ifndef INCLUDED_LAPACKWRAPPER
#define INCLUDED_LAPACKWRAPPER
namespace MolDS_mkl_wrapper{
// LapackWrapper is singleton
class LapackWrapper: private MolDS_base::Uncopyable{
public:
   static LapackWrapper* GetInstance();
   static void DeleteInstance();
   int Dsyevd(double** matrix, double* eigenValues, int size, bool calcEigenVectors);
   int Dsysv(double const* const* matrix, double* b, int size);
private:
   LapackWrapper();
   ~LapackWrapper();
   static LapackWrapper* lapackWrapper;
   bool calculatedDsysvBlockSize;
   int dsysvBlockSize;
   std::string errorMessageDsyevdInfo;
   std::string errorMessageDsyevdSize;
   std::string errorMessageDsysvInfo;
   std::string errorMessageDsysvSize;
};
}
#endif
