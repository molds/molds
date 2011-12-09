#ifndef INCLUDED_INPUT_PARSER
#define INCLUDED_INPUT_PARSER
namespace MolDS_base{

// InputParser is singleton
class InputParser{
public:
   static InputParser* GetInstance();
   static void DeleteInstance();
   void Parse(Molecule* molecule);
private:
   static InputParser* inputParser;
   InputParser();
   InputParser(InputParser&);
   void operator = (InputParser&);
   ~InputParser();
   std::string messageStartParseInput;
   std::string messageDoneParseInput;
   std::string messageTotalNumberAOs;
   std::string messageTotalNumberAtoms;
   std::string messageTotalNumberValenceElectrons;
   std::string messageInputTerms;
   std::string messageScfConditions;
   std::string messageScfMaxIterations;
   std::string messageScfRmsDensity;
   std::string messageScfDampingThresh;
   std::string messageScfDampingWeight;
   std::string messageScfDiisNumErrorVect;
   std::string messageScfDiisStartError;
   std::string messageScfDiisEndError;
   std::string messageCisConditions;
   std::string messageCisNumberActiveOcc;
   std::string messageCisNumberActiveVir;
   std::string messageCisNumberExcitedStates;
   std::string messageCisDavidson;
   std::string messageCisNormTolerance;
   std::string messageCisMaxIterations;
   std::string messageCisMaxDimensions;
   std::string messageMdConditions;
   std::string messageMdTotalSteps;
   std::string messageMdElecState;
   std::string messageMdTimeWidth;
   std::string messageMOPlotConditions;
   std::string messageMOPlotIndex;
   std::string messageMOPlotGridNumber;
   std::string messageMOPlotFrameLength;
   std::string messageMOPlotFilePrefix;
   std::string messageFs;
   std::string stringYES;
   std::string stringNO;
   std::string stringSpace;
   std::string stringCommentOut;
   std::string stringTheory;
   std::string stringTheoryEnd;
   std::string stringTheoryCNDO2;
   std::string stringTheoryINDO;
   std::string stringTheoryZINDOS;
   std::string stringTheoryMNDO;
   std::string stringTheoryAM1;
   std::string stringTheoryPM3;
   std::string stringTheoryPM3PDDG;
   std::string stringTheoryPrincipalAxes;
   std::string stringTheoryTranslate;
   std::string stringTheoryRotate;
   std::string stringTheoryNONE;
   std::string stringGeometry;
   std::string stringGeometryEnd;
   std::string stringScf;
   std::string stringScfEnd;
   std::string stringScfMaxIter;
   std::string stringScfRmsDensity;
   std::string stringScfDampingThresh;
   std::string stringScfDampingWeight;
   std::string stringScfDiisNumErrorVect;
   std::string stringScfDiisStartError;
   std::string stringScfDiisEndError;
   std::string stringMO;
   std::string stringMOPlot;
   std::string stringMOPlotEnd;
   std::string stringMOPlotGridNumber;
   std::string stringMOPlotFrameLength;
   std::string stringMOPlotFilePrefix;
   std::string stringInertiaTensor;
   std::string stringInertiaTensorEnd;
   std::string stringInertiaTensorOrigin;
   std::string stringRotate;
   std::string stringRotateEnd;
   std::string stringRotatingOrigin;
   std::string stringRotatingAxis;
   std::string stringRotatingAngle;
   std::string stringRotatingAngles;
   std::string stringRotatingType;
   std::string stringRotatingTypeAxis;
   std::string stringRotatingTypeEularAngle;
   std::string stringTranslate;
   std::string stringTranslateEnd;
   std::string stringTranslatingDifference;
   std::string stringCIS;
   std::string stringCISEnd;
   std::string stringCISActiveOcc;
   std::string stringCISActiveVir;
   std::string stringCISNStates;
   std::string stringCISDavidson;
   std::string stringCISMaxIter;
   std::string stringCISMaxDimensions;
   std::string stringCISNormTolerance;
   std::string stringMD;
   std::string stringMDEnd;
   std::string stringMDTotalSteps;
   std::string stringMDElecState;
   std::string stringMDTimeWidth;
   void CalcMolecularBasics(Molecule* molecule);
   void CheckCisConditions(Molecule* molecule);
   void CheckMdConditions();
   void OutputMolecularBasics(Molecule* molecule);
   void OutputScfConditions();
   void OutputCisConditions();
   void OutputMdConditions();
   void OutputMOPlotConditions();
   void OutputInputTerms(std::vector<std::string> inputTerms);
   bool IsCommentOut(std::string str);
   std::vector<std::string> GetInputTerms();
};

}
#endif





