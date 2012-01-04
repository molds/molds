#//************************************************************************//
#// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
#//                                                                        // 
#// This file is part of MolDS.                                            // 
#//                                                                        // 
#// MolDS is free software: you can redistribute it and/or modify          // 
#// it under the terms of the GNU General Public License as published by   // 
#// the Free Software Foundation, either version 3 of the License, or      // 
#// (at your option) any later version.                                    // 
#//                                                                        // 
#// MolDS is distributed in the hope that it will be useful,               // 
#// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
#// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
#// GNU General Public License for more details.                           // 
#//                                                                        // 
#// You should have received a copy of the GNU General Public License      // 
#// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
#//************************************************************************//
CC = icpc
LIBS32 = -lmkl_intel      -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LIBS64 = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
CFLAGS = -O0 -openmp -openmp-report2
EXENAME = MolDS.out
DEPFILE = obj/objfile.dep

ALL_CPP_FILES = MolDS.cpp base/MolDSException.cpp mkl_wrapper/LapackWrapper.cpp base/Utilities.cpp base/Enums.cpp base/MathUtilities.cpp base/MallocerFreer.cpp base/EularAngle.cpp base/Parameters.cpp base/atoms/Atom.cpp base/atoms/Hatom.cpp base/atoms/Liatom.cpp base/atoms/Catom.cpp base/atoms/Natom.cpp base/atoms/Oatom.cpp base/atoms/Satom.cpp base/Molecule.cpp base/InputParser.cpp base/GTOExpansionSTO.cpp base/RealSphericalHarmonicsIndex.cpp base/MOLogger.cpp cndo/Cndo2.cpp indo/Indo.cpp zindo/ZindoS.cpp mndo/Mndo.cpp am1/Am1.cpp pm3/Pm3.cpp pm3/Pm3Pddg.cpp md/MD.cpp 
ALL_HEAD_FILES = base/MolDSException.h base/Uncopyable.h mkl_wrapper/LapackWrapper.h base/Utilities.h base/Enums.h base/MathUtilities.h base/MallocerFreer.h base/EularAngle.h base/Parameters.h base/atoms/Atom.h base/atoms/Hatom.h base/atoms/Liatom.h base/atoms/Catom.h base/atoms/Natom.h base/atoms/Oatom.h base/atoms/Satom.h base/Molecule.h base/InputParser.h base/GTOExpansionSTO.h base/RealSphericalHarmonicsIndex.h base/MOLogger.h cndo/Cndo2.h indo/Indo.h zindo/ZindoS.h mndo/Mndo.h am1/Am1.h pm3/Pm3.h pm3/Pm3Pddg.h md/MD.h 
ALL_OBJ_FILES = obj/MolDS.o obj/MolDSException.o obj/LapackWrapper.o obj/Utilities.o obj/Enums.o obj/MathUtilities.o obj/MallocerFreer.o obj/EularAngle.o obj/Parameters.o obj/Atom.o obj/Hatom.o obj/Liatom.o obj/Catom.o obj/Natom.o obj/Oatom.o obj/Satom.o obj/Molecule.o obj/InputParser.o obj/GTOExpansionSTO.o obj/RealSphericalHarmonicsIndex.o obj/MOLogger.o obj/Cndo2.o obj/Indo.o obj/ZindoS.o obj/Mndo.o obj/Am1.o obj/Pm3.o obj/Pm3Pddg.o obj/MD.o 

$(EXENAME): depend $(ALL_OBJ_FILES)
ifeq ($(INTEL), 64)
	$(CC) -o $@ $(ALL_OBJ_FILES) $(LIBS64) $(CFLAGS) 
else
	$(CC) -o $@ $(ALL_OBJ_FILES) $(LIBS32) $(CFLAGS) 
endif

ifeq ($(wildcard $(DEPFILE)), $(DEPFILE))
include $(DEPFILE)
endif

depend: $(ALL_HEAD_FILES) $(ALL_CPP_FILES)
ifeq ($(INTEL), 64)
	$(CC) -MM $(ALL_CPP_FILES) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -o $$@ $$< $$(LIBS64) $$(CFLAGS) -c /g' > $(DEPFILE)
else
	$(CC) -MM $(ALL_CPP_FILES) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -o $$@ $$< $$(LIBS32) $$(CFLAGS) -c /g' > $(DEPFILE)
endif

.PHONY: clean
clean:
	rm -f $(ALL_OBJ_FILES) $(EXENAME) $(DEPFILE)

all: clean $(EXENAME)
