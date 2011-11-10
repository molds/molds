#ifndef INCLUDED_PM3
#define INCLUDED_PM3
namespace MolDS_pm3{

/***
 *  Main Refferences for PM3 are [S_1989, S_1989-2, S_1991, S_2004, S_2007]
 */
class Pm3 : public MolDS_am1::Am1{
public:
   Pm3();
   virtual ~Pm3();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
private:
};

}
#endif



