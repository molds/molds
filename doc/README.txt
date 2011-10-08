
Compile: 
   for 32 bit
   $icc MolDS.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3
   $icc MolDS.cpp -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3 -openmp -openmp-report2

   for 64 bit
   $icc MolDS.cpp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
   $icc MolDS.cpp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3 -openmp -openmp-report2


Carring Out:
   $./a.out < input.in

Capabilities:

            | HF  | CIS | MD(gs) | MD(es) |
   ---------|-----|-----|--------|--------|
   CNDO2    | OK  | --  | --     | --     |
   ---------|-----|-----|--------|--------|
   INDO     | OK  | --  | --     | --     |
   ---------|-----|-----|--------|--------|
   ZINDO/S  | OK  | OK  | OK     | --     |
   ---------|-----|-----|--------|--------|
   MNDO     | Sch | Sch | Sch    | Sch    |

      "OK", "Sch", and "--" mean available, shceduled, and non-scheduled methods, respectively.
      MD(gs) and MD(es) mean Molecular Dynamics on ground and excited states, respectively. 

How to Write Input-files:

   Comment Out:
      Lines starting with "//" or "#" in input-files are treated as comments.


   SCF:
      Write "cndo/2", "indo", or "zindo/s" in theory-directive.

      E.g. 
         THEORY
            indo 
         THEORY_END
   
      -options
       "max_iter", "rms_density", "damping_thresh", "damping_weight", 
       "diis_num_error_vect", "diis_start_error", and "diis_end_error" are prepared as options.

       Default value of "max_iter" is 100.
       Default value of "rms_density" is 10**(-8.0).
       Default value of "damping_thresh" is 1.
       Default value of "damping_weight" is 0.8.
       Default value of "diis_num_error_vect" is 5.
       Default value of "diis_start_error" is 0.01.
       Default value of "diis_end_error" is 10**(-8.0).

       E.g.
         SCF
            max_iter 200
            rms_density 0.00000001
            damping_thresh 0.1
            damping_weight 0.7
            diis_num_error_vect 6
            diis_start_error 0.01
            diis_end_error 0.00000001
         SCF_END


   CIS:
      Write CIS-directive.

      E.g.
         CIS
            (options)
         CIS_END
   
      -options
       "davidson", "active_occ", "active_vir", "max_iter", "max_dim", "norm_tol", 
       and "nstates" are prepared as options.

       "davidson" should be set as "yes" or "no". 
       Default value of "davidson" is "yes".

       "active_occ" ("active_vir") is set to the number of occupied (virtual) orbitals
       if user set "active_occ" ("active_vir") to be greater than 
       the number of occupied (virtual) orbitals. 
       Default value of "active_occ" is 10. Default value of "active_vir" is 10.

       "nstates" is valid for the Davidson algorithm only, 
       hence "nstates" is set to "active_occ*active_vir" 
       in direct CIS algorithm (without the Davidson algorithem). 
       Default value of "nstates" is 5 for the Davidson algorithem.

       "max_iter" is valid for the Davidson algorithm only. 
       This option means the number of times of Davidson roop. 
       Default value of "max_iter" is 100.

       "max_dim" is valid for the Davidson algorithm only. 
       This option means the number of slater determinans used by expansion of the excited states. 
       Note that Hartree-Fock state (groudn state) is not included in the "max_dim".
       Default value of "max_dim" is 100.

       "norm_tol" is valid for the Davidson algorithm only. 
       This option means the max tolerance for the norm of the residual vectors.
       Default value of "norm_tol" is 10**(-6.0).

       E.g.
         CIS
            davidson no
            active_occ 2
            active_vir 2
            nstates 1000
            max_iter 100
            max_dim 100
            norm_tol 0.000001
         CIS_END


   MD (Molecular dynamics):
      Write MD-directive.

      E.g.
         MD 
            (options)
         MD_END
  
      -options
       "total_steps", "electronic_state", 
       and "dt" are prepared as options.

       "electronic_state" means the electronic eigen state 
       on which the system runs.
       Default value of "electronic_state" is 0. That is, 
       electronic ground state.

       Default value of "total_steps" is 10. 

       "dt" should be set in femto-second.
       Default value of "dt" is 0.1[fs].


   Principal Axes (Diagonalizing the inertia tensor):
      Write "principal_axes" in theory-directive.

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


   Translate Molecule:
      Write "translate" in theory-directive.

      E.g. 
         THEORY
            translate 
         THEORY_END

      -options
       "difference" indicates difference for the translation in angstrom unit.
       This option is written in translate-directive.
       Default values are 0, 0, and 0.

       E.g. 
         TRANSLATE
            difference 12 30 45
         TRANSLATE_END



