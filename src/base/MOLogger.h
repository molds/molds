#ifndef INCLUDED_MOLOGGER
#define INCLUDED_MOLOGGER
namespace MolDS_base{

class MOLogger{
public:
   MOLogger(const MolDS_base::Molecule& molecule, double const* const* fockMatrix, MolDS_base::TheoryType theory);
   void DrawMO(int moIndex);
   void DrawMO(std::vector<int> moIndeces);
private:
   std::string stringCubeExtension;
   std::string messageCubeHeaderComment1;
   std::string messageCubeHeaderComment2;
   std::string messageStartMOPlot;
   std::string messageEndMOPlot;
   std::string messageSkippedMOIndex;
   MOLogger();
   MolDS_base::Molecule const* molecule;
   double const* const* fockMatrix;
   MolDS_base::TheoryType theory;
   void SetMessage();
};

}
#endif
