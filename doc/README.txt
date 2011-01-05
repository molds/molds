
Compile: 
   for 32 bit
   $icc MolDS.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

   for 64 bit
   $icc MolDS.cpp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread


Carring Out:
   $./a.out < input.in

Principal Axes (Diagonalizing the inertia tensor):
   Write "principal-axes" in theory-directive.

   E.g. 
         THEORY
            principal-axes
         THEORY_END
   
   -options
    option is "origin" only for setting the origin of the inertia tensor.
    options are written in inertia-directive in angstrom unit.
    Center of mass is used as origin when the "origin" is not set.

    E.g.
         INERTIA
            origin 1.2 2.3 3.4
         INERTIA_END

