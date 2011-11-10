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
   virtual double CalcDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB);
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA);
   virtual void CalcHFProperties();
   virtual void OutputHFResults(double** fockMatrix, 
                                double* energiesMO, 
                                double* atomicElectronPopulation, 
                                MolDS_base::Molecule* molecule);
private:
};

}
#endif



