#ifndef INCLUDED_EULARANGLE
#define INCLUDED_EULARANGLE
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
   std::string errorMessageInvalidXYZ;
   double alpha;
   double beta;
   double gamma;
   void SetMessage();
};

}
#endif
