
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

Rotate Molecule:
   Write "rotate" in theory-directive.

   E.g. 
         THEORY
            rotate
         THEORY_END

   -options
    "type", "origin", "axis", "angle" and "angles" are prepared as options.
    These options are written in rotate-directive. Examples are shown below.

    "type" indicates whether the rotating is carring out around a axis or acording to Euler angles.
    Default value of "type" is axis.

    "origin" indicates the origin of the rotation in angstrom unit.
    Default value of "origin" is center of mass.

    "axis" indicates a axis around which the rotation is carried out in angstrom unit.
    Default value of "axis" is z-axis.
    This option is valid only for "type" set as axis.

    "angle" indicates angle for the rotation around the "axis" in degree unit.
    Default value of "angle" is 0.
    This option is valid only for "type" set as axis.

    "angles" indicates Euler angles for the rotation in degree unit.
    Default values of "angles" are 0, 0, and 0.
    This option is valid only for "type" set as Euler angles.
   

   E.g. for "type" set as axis
         ROTATE
            type axis
            origin 1.0 2.0 3.0 
            axis 3.0 4.0 5.0 
            angle 30
         ROTATE_END

   E.g. for "type" set as Euler angles
         ROTATE
            type eular_angle
            angles 15 25 35
         ROTATE_END



