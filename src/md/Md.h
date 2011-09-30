#ifndef INCLUDED_MD
#define INCLUDED_MD

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_md{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class Md : public MolDS_cndo::Cndo2{
public:
   Md();
   ~Md();
private:
   MolDS_zindo::ZindoS* zindoS;
   vector<TheoryType> enableTheoryTypes;
   void SetMessages();
   void SetEnableTheoryTypes();
};

Md::Md(){
   //cout << "Md created\n";
}

Md::~Md(){
   //cout << "Md\n";
}

void Md::SetMessages(){
}

void Md::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
}


}
#endif



