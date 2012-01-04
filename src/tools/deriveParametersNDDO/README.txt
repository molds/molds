compile for 32bit OS:

icpc ../../base/MolDSException.cpp ../../base/Utilities.cpp ../../base/Enums.cpp ../../base/MathUtilities.cpp ../../base/MallocerFreer.cpp ../../base/EularAngle.cpp ../../base/Parameters.cpp deriveParametersNDDO.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O0 -openmp -openmp-report2
