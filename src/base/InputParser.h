//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#ifndef INCLUDED_INPUT_PARSER
#define INCLUDED_INPUT_PARSER
namespace MolDS_base{

// InputParser is singleton
class InputParser: private Uncopyable{
public:
   static InputParser* GetInstance();
   static void DeleteInstance();
   void Parse(Molecule* molecule) const;
private:
   static InputParser* inputParser;
   InputParser();
   ~InputParser();
   void SetMessages();
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
   std::string messageMemoryConditions;
   std::string messageMemoryLimitHeap;
   std::string messageMemoryMB;
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
   std::string stringMemory;
   std::string stringMemoryEnd;
   std::string stringMemoryLimitHeap;
   std::string stringMD;
   std::string stringMDEnd;
   std::string stringMDTotalSteps;
   std::string stringMDElecState;
   std::string stringMDTimeWidth;
   void CalcMolecularBasics(Molecule* molecule) const;
   void CheckCisConditions(const Molecule& molecule) const;
   void CheckMdConditions() const;
   void OutputMolecularBasics(Molecule* molecule) const;
   void OutputScfConditions() const;
   void OutputMemoryConditions() const;
   void OutputCisConditions() const;
   void OutputMdConditions() const;
   void OutputMOPlotConditions() const;
   void OutputInputTerms(std::vector<std::string> inputTerms) const;
   bool IsCommentOut(std::string str) const;
   std::vector<std::string> GetInputTerms() const;
};

}
#endif





