#ifndef INCLUDED_MATHUTILITIES
#define INCLUDED_MATHUTILITIES
namespace MolDS_base{
// n!
int Factorial(int n);
// nCk
int Conbination(int n, int k);
// max
template <typename T> T Max(T a, T b);
// min
template <typename T> T min(T a, T b);
// rotating matrix
void CalcRotatingMatrix(double matrix[][3], double sita, CartesianType cartesianType);
}
#endif
 
