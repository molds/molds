#include<string>
#include<stdexcept>
#include"MolDSException.h"
using namespace std;
namespace MolDS_base{
MolDSException::MolDSException(string cause) : domain_error(cause){
}
}





