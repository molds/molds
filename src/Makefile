CC = icpc
LIBS32 = -lmkl_intel      -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LIBS64 = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
CFLAGS = -O3 -openmp -openmp-report2
DEPFILE = obj/objfile.dep
EXENAME = MolDS.out

ALL_CPP_FILES = MolDS.cpp base/MolDSException.cpp mkl_wrapper/LapackWrapper.cpp base/Utilities.cpp base/Enums.cpp base/MathUtilities.cpp base/MallocerFreer.cpp base/EularAngle.cpp base/Parameters.cpp base/atoms/Atom.cpp base/atoms/Hatom.cpp base/atoms/Liatom.cpp base/atoms/Catom.cpp base/atoms/Natom.cpp base/atoms/Oatom.cpp base/atoms/Satom.cpp base/Molecule.cpp base/InputParser.cpp base/GTOExpansionSTO.cpp base/RealSphericalHarmonicsIndex.cpp cndo/Cndo2.cpp indo/Indo.cpp zindo/ZindoS.cpp mndo/Mndo.cpp am1/Am1.cpp pm3/Pm3.cpp md/MD.cpp 
ALL_OBJ_FILES = obj/MolDS.o obj/MolDSException.o obj/LapackWrapper.o obj/Utilities.o obj/Enums.o obj/MathUtilities.o obj/MallocerFreer.o obj/EularAngle.o obj/Parameters.o obj/Atom.o obj/Hatom.o obj/Liatom.o obj/Catom.o obj/Natom.o obj/Oatom.o obj/Satom.o obj/Molecule.o obj/InputParser.o obj/GTOExpansionSTO.o obj/RealSphericalHarmonicsIndex.o obj/Cndo2.o obj/Indo.o obj/ZindoS.o obj/Mndo.o obj/Am1.o obj/Pm3.o obj/MD.o 

$(EXENAME) : $(ALL_OBJ_FILES)
ifeq ($(INTEL), 64)
	$(CC) -o $@ $(ALL_OBJ_FILES) $(LIBS64) $(CFLAGS) 
else
	$(CC) -o $@ $(ALL_OBJ_FILES) $(LIBS32) $(CFLAGS) 
endif

depend : $(DEPFILE)

$(DEPFILE) : $(ALL_CPP_FILES)
ifeq ($(INTEL), 64)
	$(CC) -MM $(ALL_CPP_FILES) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -o $$@ $$< $$(LIBS64) $$(CFLAGS) -c /g' > $(DEPFILE)
else
	$(CC) -MM $(ALL_CPP_FILES) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -o $$@ $$< $$(LIBS32) $$(CFLAGS) -c /g' > $(DEPFILE)
endif

ifeq ($(wildcard $(DEPFILE)), $(DEPFILE))
include $(DEPFILE)
endif

.PHONY : clean
clean :
	rm -f $(ALL_OBJ_FILES)
	rm -f $(EXENAME)
	rm -f $(DEPFILE)


