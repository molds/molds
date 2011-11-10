#ifndef INCLUDED_MOLDSEXCEPTION
#define INCLUDED_MOLDSEXCEPTION
namespace MolDS_base{
class MolDSException : public std::domain_error {
public:
   MolDSException(std::string cause);
private:
};
}
#endif





