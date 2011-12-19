#ifndef INCLUDED_AM1
#define INCLUDED_AM1
namespace MolDS_am1{

/***
 *  Main Refferences for AM1 are [DZHS_1985, DY_1990]
 */
class Am1 : public MolDS_mndo::Mndo{
public:
   Am1();
   virtual ~Am1();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA) const;
   virtual void CalcHFProperties();
   virtual void OutputHFResults(double const* const* fockMatrix, 
                                double const* energiesMO, 
                                double const* atomicElectronPopulation, 
                                const MolDS_base::Molecule& molecule) const;
private:
};

}
#endif



