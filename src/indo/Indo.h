#ifndef INCLUDED_INDO
#define INCLUDED_INDO
namespace MolDS_indo{

/***
 *  Refferences for Indo are [PB_1970] and [PS_1966].
 */
class Indo : public MolDS_cndo::Cndo2{
public:
   Indo();
   virtual ~Indo();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetFockDiagElement(MolDS_base_atoms::Atom* atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     MolDS_base::Molecule* molecule, 
                                     double** gammaAB,
                                     double** orbitalElectronPopulation, 
                                     double* atomicElectronPopulation,
                                     double****** twoElecTwoCore,
                                     bool isGuess);
   virtual double GetFockOffDiagElement(MolDS_base_atoms::Atom* atomA, 
                                        MolDS_base_atoms::Atom* atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
                                        int mu, 
                                        int nu, 
                                        MolDS_base::Molecule* molecule, 
                                        double** gammaAB, 
                                        double** overelap,
                                        double** orbitalElectronPopulation,
                                        double****** twoElecTwoCore,
                                        bool isGuess);
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              MolDS_base::Molecule* molecule, 
                                              double** fockMatrix, 
                                              double** gammaAB);
private:
   double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                        MolDS_base::OrbitalType orbital2, 
                        double gamma, 
                        MolDS_base_atoms::Atom* atom); // Indo Coulomb Interaction, (3.87) - (3.91) in J. A. Pople book.
   double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                         MolDS_base::OrbitalType orbital2, 
                         double gamma, 
                         MolDS_base_atoms::Atom* atom); // Indo Exchange Interaction, (3.87) - (3.91) in J. A. Pople book.
};

}
#endif



