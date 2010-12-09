for 32 bit
   $icc MolDS.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

for 64 bit
   $icc MolDS.cpp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
