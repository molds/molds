#ifndef INCLUDED_MALLOCERFREER
#define INCLUDED_MALLOCERFREER
namespace MolDS_base{

// MallocerFreer is singleton
class MallocerFreer: private Uncopyable{
public:
   static MallocerFreer* GetInstance();
   static void DeleteInstance();

   //1d
   template<typename T> T* Malloc(int size1) const{
      T* matrix;  
      matrix = new T[size1];
      MallocerFreer::AddCurrentMalloced((double)(size1*sizeof(T)));
      this->Initialize<T>(matrix, size1);
      return(matrix); 
   }

   template<typename T> void Initialize(T* matrix, int size1) const{
      for(int i=0;i<size1;i++){
         matrix[i] = (T)0;
      }
   }

   template<typename T> void Free(T** matrix, int size1) const{
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced((double)(size1*sizeof(T)));
      *matrix = NULL;
   }

   //2d
   template<typename T> T** Malloc(int size1, int size2) const{
      T** matrix;  
      matrix = new T*[size1];
      if(matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         matrix[i] = new T[size2];
         if (matrix[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
      }
      MallocerFreer::AddCurrentMalloced((double)(size1*size2*sizeof(T)));
      this->Initialize<T>(matrix, size1, size2);
      return(matrix); 
   }

   template<typename T> void Initialize(T** matrix, int size1, int size2) const{
      for(int i=0;i<size1;i++){
         for(int j=0;j<size2;j++){
            matrix[i][j] = (T)0.0;
         }
      }
   }

   template<typename T> void Free(T*** matrix, int size1, int size2) const{
      int i=0;
      for(i=0;i<size1;i++){
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced((double)(size1*size2*sizeof(T)));
      *matrix = NULL;
   }

   // 3d
   template<typename T> T*** Malloc(int size1, int size2, int size3) const{
      T*** matrix;  
      matrix = new T**[size1];
      if(matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         matrix[i] = new T*[size2];
         if(matrix[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            matrix[i][j] = new T[size3];
            if(matrix[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
         }
      }
      MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*sizeof(T)));
      this->Initialize<T>(matrix, size1, size2, size3);
      return(matrix); 
   }

   template<typename T> void Initialize(T*** matrix, int size1, int size2, int size3) const{
      for(int i=0;i<size1;i++) {
         for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
               matrix[i][j][k] = (T)0.0;
            }
         }
      }  
   }

   template<typename T> void Free(T**** matrix, int size1, int size2, int size3) const{
      int i=0, j=0;
      for (i=0;i<size1;i++) {
         for (j=0;j<size2;j++) {
            delete [] (*matrix)[i][j];
         }
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*sizeof(T)));
      *matrix = NULL;
   }

   //4d
   template<typename T> T**** Malloc(int size1, int size2, int size3, int size4) const{
      T**** matrix;  
      matrix = new T***[size1];
      if(matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         matrix[i] = new T**[size2];
         if(matrix[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            matrix[i][j] = new T*[size3];
            if(matrix[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               matrix[i][j][k] = new T[size4];
               if(matrix[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*sizeof(T)));
      this->Initialize<T>(matrix, size1, size2, size3, size4);
      return(matrix); 
   }

   template<typename T> void Initialize(T**** matrix, int size1, int size2, int size3, int size4) const{
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

   template<typename T> void Free(T***** matrix, int size1, int size2, int size3, int size4) const{
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
      MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*sizeof(T)));
      *matrix = NULL;
   }

   //5d
   template<typename T> T***** Malloc(int size1, int size2, int size3, int size4, int size5) const{
      T***** matrix;  
      matrix = new T****[size1];
      if(matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         matrix[i] = new T***[size2];
         if(matrix[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            matrix[i][j] = new T**[size3];
            if(matrix[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               matrix[i][j][k] = new T*[size4];
               if(matrix[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
               for(int l=0;l<size4;l++){
                  matrix[i][j][k][l] = new T[size5];
                  if(matrix[i][j][k][l]==NULL){
                     throw MolDSException(this->errorMessageMallocFailure);
                  }
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*size5*sizeof(T)));
      this->Initialize<T>(matrix, size1, size2, size3, size4, size5);
      return(matrix); 
   }

   template<typename T> void Initialize(T***** matrix, int size1, int size2, int size3, int size4, int size5) const{
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

   template<typename T> void Free(T****** matrix, int size1, int size2, int size3, int size4, int size5) const{
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
      MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*size5*sizeof(T)));
      *matrix = NULL;
   }

   //6d
   template<typename T> T****** Malloc(int size1, int size2, int size3, int size4, int size5, int size6) const{
      T****** matrix;  
      matrix = new T*****[size1];
      if(matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         matrix[i] = new T****[size2];
         if(matrix[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            matrix[i][j] = new T***[size3];
            if(matrix[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               matrix[i][j][k] = new T**[size4];
               if(matrix[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
               for(int l=0;l<size4;l++){
                  matrix[i][j][k][l] = new T*[size5];
                  if(matrix[i][j][k][l]==NULL){
                     throw MolDSException(this->errorMessageMallocFailure);
                  }
                  for(int m=0;m<size5;m++){
                     matrix[i][j][k][l][m] = new T[size6];
                     if(matrix[i][j][k][l][m]==NULL){
                        throw MolDSException(this->errorMessageMallocFailure);
                     }
                  }
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced((double)(size1*size2*size3*size4*size5*size6*sizeof(T)));
      this->Initialize<T>(matrix, size1, size2, size3, size4, size5, size6);
      return(matrix); 
   }

   template<typename T> void Initialize(T****** matrix, int size1, int size2, int size3, int size4, int size5, int size6) const{
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

   template<typename T> void Free(T******* matrix, int size1, int size2, int size3, int size4, int size5, int size6) const{
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
      MallocerFreer::SubtCurrentMalloced((double)(size1*size2*size3*size4*size5*size6*sizeof(T)));
      *matrix = NULL;
   }
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
