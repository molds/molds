#ifndef INCLUDED_EULARANGLE
#define INCLUDED_EULARANGLE
#include<math.h>
#include<string>


namespace MolDS_base{

class EularAngle{
public:
   EularAngle();
   EularAngle(double x, double y, double z);
   EularAngle(double* angles);
   double GetAlpha();
   double GetBeta();
   double GetGamma();
   void SetAlpha(double alpha);
   void SetBeta(double beta);
   void SetGamma(double gamma);
private:
   string errorMessageInvalidXYZ;
   double alpha;
   double beta;
   double gamma;
   void SetMessage();
};

EularAngle::EularAngle(){
   this->SetMessage();
   this->alpha = 0.0;
   this->beta = 0.0;
   this->gamma = 0.0;
}

EularAngle::EularAngle(double x, double y, double z){
   this->SetMessage();
   double r = 0.0;

   // calc. beta
   if(x==0.0 && y==0.0 && z==0.0){
      stringstream ss;
      ss << this->errorMessageInvalidXYZ;
      throw MolDSException(ss.str());
   }

   r = sqrt( pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) );
   this->beta = acos(z/r);

   // calc. alpha
   if(x==0.0 && y==0.0){
      this->alpha = 0.0;
   }
   else{ 
      r = sqrt( pow(x, 2.0) + pow(y, 2.0) );
      this->alpha = atan2(y/r, x/r);
   }

   // set gamma
   this->gamma = 0.0;
   
}

EularAngle::EularAngle(double* angles){
   this->SetMessage();
   this->alpha = angles[0];
   this->beta = angles[1];
   this->gamma = angles[2];
}

double EularAngle::GetAlpha(){
   return this->alpha;
}

double EularAngle::GetBeta(){
   return this->beta;
}

double EularAngle::GetGamma(){
   return this->gamma;
}

void EularAngle::SetAlpha(double alpha){
   this->alpha = alpha;
}

void EularAngle::SetBeta(double beta){
   this->beta = beta;
}

void EularAngle::SetGamma(double gamma){
   this->gamma = gamma;
}

void EularAngle::SetMessage(){
   this->errorMessageInvalidXYZ="Error in base::EularAngle: Invalid coordinates. x=y=z=0.\n";
}



}
#endif
