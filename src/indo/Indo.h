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
   virtual double GetFockDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     const MolDS_base::Molecule& molecule, 
                                     double const* const* gammaAB,
                                     double const* const* orbitalElectronPopulation, 
                                     double const* atomicElectronPopulation,
                                     double const* const* const* const* const* const* twoElecTwoCore,
                                     bool isGuess) const;
   virtual double GetFockOffDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                        const MolDS_base_atoms::Atom& atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
                                        int mu, 
                                        int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overelap,
                                        double const* const* orbitalElectronPopulation,
                                        double const* const* const* const* const* const* twoElecTwoCore,
                                        bool isGuess) const;
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
private:
   double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                        MolDS_base::OrbitalType orbital2, 
                        double gamma, 
                        const MolDS_base_atoms::Atom& atom) const; // Indo Coulomb Interaction, (3.87) - (3.91) in J. A. Pople book.
   double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                         MolDS_base::OrbitalType orbital2, 
                         double gamma, 
                         const MolDS_base_atoms::Atom& atom) const; // Indo Exchange Interaction, (3.87) - (3.91) in J. A. Pople book.
};

}
#endif



