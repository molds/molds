compile for 32bit OS:

icpc ../../base/MolDSException.cpp ../../base/Enums.cpp ../../base/MathUtilities.cpp paramYZ.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O0/../base/MathUtilities.cpp paramYZ.cpp -lmkl_intel -lmkl_intel_thread -lmkl_c
