#ifndef INCLUDED_PM3PDDG
#define INCLUDED_PM3PDDG
namespace MolDS_pm3{

/***
 *  Main Refferences for PM3/PDDG are [RCJ_2002, BGRJ_2003, and BGJ_2003]
 */
class Pm3Pddg : public MolDS_pm3::Pm3{
public:
   Pm3Pddg();
   virtual ~Pm3Pddg();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double CalcDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB);
   //virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
   //                                                     int atomBIndex, 
   //                                                     MolDS_base::CartesianType axisA);
private:
};

}
#endif



