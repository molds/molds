#ifndef INCLUDED_MOLDSEXCEPTION
#define INCLUDED_MOLDSEXCEPTION

#include <stdexcept>

using namespace std;

namespace MolDS_base{

class MolDSException : public domain_error {
public:
   MolDSException(string cause);
private:
};

MolDSException::MolDSException(string cause) : domain_error(cause){


}

}
#endif





