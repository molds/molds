#ifndef INCLUDED_MD
#define INCLUDED_MD
namespace MolDS_md{

/***
 *  Velocty Verlet is used here.
 */
class MD{
public:
   MD();
   ~MD();
   void SetTheory(MolDS_cndo::Cndo2* cndo);
   void DoesMD();
private:
   std::string messageinitialConditionMD;
   std::string messageStartMD;
   std::string messageEndMD;
   std::string messageStartStepMD;
   std::string messageEndStepMD;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreKineticEnergy;
   std::string messageCoreRepulsionEnergy;
   std::string messageElectronicEnergy;
   std::string messageTotalEnergy;
   std::string messageErrorEnergy;
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   MolDS_cndo::Cndo2* cndo;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
   void OutputEnergies(double initialEnergy);
   double OutputEnergies();
};

}
#endif



