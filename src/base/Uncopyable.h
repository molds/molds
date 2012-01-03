#ifndef INCLUDED_UNCOPYABLE
#define INCLUDED_UNCOPYABLE
namespace MolDS_base{

class Uncopyable{
public:
protected:
   Uncopyable(){};
   ~Uncopyable(){};
private:
   Uncopyable(const Uncopyable&);
   Uncopyable& operator = (const Uncopyable&);
};
}
#endif
