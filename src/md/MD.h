#ifndef INCLUDED_MD
#define INCLUDED_MD

using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_md{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class MD{
public:
   MD();
   ~MD();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
private:
   MolDS_cndo::Cndo2* cndo;
   vector<TheoryType> enableTheoryTypes;
   string errorMessageNotEnebleTheoryType;
   string errorMessageTheoryType;
   void CheckEnableTheoryType(TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
};

MD::MD(){
   this->cndo = NULL;
   this->SetEnableTheoryTypes();
   this->SetMessages();
   //cout << "MD created \n";
}

MD::~MD(){
   //cout << "MD deleted\n";
}

void MD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in md::MD::CheckEnableTheoryType: Non available theory is set.\n";
}

void MD::SetTheory(MolDS_cndo::Cndo2* cndo){
   // check enable electonic theory
   this->CheckEnableTheoryType(cndo->GetTheoryType());
   this->cndo = cndo;
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
}

void MD::CheckEnableTheoryType(TheoryType theoryType){
   bool isEnable = false;
   for(int i=0; i<this->enableTheoryTypes.size();i++){
      if(theoryType == this->enableTheoryTypes[i]){
         isEnable = true;
         break;
      }
   }
   if(!isEnable){
      stringstream ss;
      ss << this->errorMessageNotEnebleTheoryType;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      throw MolDSException(ss.str());
   }
}

}
#endif



