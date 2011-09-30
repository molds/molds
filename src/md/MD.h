#ifndef INCLUDED_MD
#define INCLUDED_MD

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_md{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class MD : public MolDS_cndo::Cndo2{
public:
   MD();
   MD(MolDS_zindo::ZindoS* zindoS);
   ~MD();
private:
   MolDS_zindo::ZindoS* zindoS;
   vector<TheoryType> enableTheoryTypes;
   void SetMessages();
   void SetEnableTheoryTypes();
};

MD::MD(){
   cout << "MD created 1\n";
}

MD::MD(MolDS_zindo::ZindoS* zindoS){
   this->zindoS = zindoS;
   cout << "MD created 2\n";
}

MD::~MD(){
   //cout << "MD\n";
}

void MD::SetMessages(){
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
}


}
#endif



