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
==============================================================================
COMPILE(using GNUmake): 
   In the "src" directory in the MolDS package.

   for 32 bit
   $ make depend INTEL=32
   $ make INTEL=32

   for 64 bit
   $ make depend INTEL=64
   $ make INTEL=64

   The compile succeeded if you could fine "MolDS.out" in the "src" directory. 
   Use "$ make clean" when you wanna clean the compilation.

COMPILE(primitive method): 
   In the "src" directory in the MolDS package.
   for 32 bit
   $ icc <MolDS.cpp and all cpp-files> -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3 -openmp -openmp-report2

   for 64 bit
   $ icc <MolDS.cpp and all cpp-files> -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O3 -openmp -openmp-report2

==============================================================================
CARRY OUT MolDS:
   In the "src" directory in the MolDS package.
   $ ./MolDS.out < input.in

==============================================================================
SAMPLE and TEST
   See files in "test" directories for sample files.
   In the "test" directory, *.in files are input files, then *.dat files are
   associated output files. For test calculations, carry out below ruby-script 
   in the "test" directory. This scripb will finished in a few minutes with big 
   output(a few thausands lines).

   $ ruby Test_Of_MolDS.rb

==============================================================================
CAPABILITIES:

   Electronic state and molecular dynamics:
            | HF  | CIS | MD(gs) | MD(es) | MC(gs) | MC(es) | 
   ---------|-----|-----|--------|--------|--------|--------|
   CNDO2    | OK  | --  | --     | --     | Sch    | --     |
   ---------|-----|-----|--------|--------|--------|--------|
   INDO     | OK  | --  | --     | --     | Sch    | --     |
   ---------|-----|-----|--------|--------|--------|--------|
   ZINDO/S  | OK  | OK  | OK     | --     | Sch    | Sch    |
   ---------|-----|-----|--------|--------|--------|--------|
   MNDO     | OK  | OK  | OK     | OK     | Sch    | Sch    |
   ---------|-----|-----|--------|--------|--------|--------|
   AM1      | OK  | OK  | OK     | OK     | Sch    | Sch    |
   ---------|-----|-----|--------|--------|--------|--------|
   PM3      | OK  | OK  | OK     | OK     | Sch    | Sch    |
   ---------|-----|-----|--------|--------|--------|--------|
   PM3/PDDG | OK  | OK  | OK     | OK     | Sch    | Sch    |

      "OK", "Sch", and "--" mean available, shceduled, and non-scheduled methods, respectively.
      MD(gs) and MD(es) mean Born-Oppenheimer Molecular Dynamics on ground and excited states, respectively. 

   Elements:
   CNDO2    | H, Li, C, N, O, and S
   INDO     | H, Li, C, N, and O
   ZINDO/S  | H, C, N, O, and S
   MNDO     | H, C, N, O, and S 
   AM1      | H, C, N, O, and S 
   PM3      | H, C, N, O, and S 
   PM3/PDDG | H, C, N, O, and S 

==============================================================================
HOW TO WRITE INPUT:

   <Comment Out>
      Lines starting with "//" or "#" in input-files are treated as comments.

   <SCF>
      Write "cndo/2", "indo", "zindo/s", "mndo", "am1", 
      "pm3/pddg", or "pm3/pddg" in theory-directive.
      MNDO only supports (can calculate) Heats of formation.

      E.g. 
         THEORY
            indo 
         THEORY_END
   
      -options
       "max_iter", "rms_density", "damping_thresh", "damping_weight", 
       "diis_num_error_vect", "diis_start_error", and "diis_end_error" are prepared as options.

       The default value of "max_iter" is 100.
       The default value of "rms_density" is 10**(-8.0).
       The default value of "damping_thresh" is 1.
       The default value of "damping_weight" is 0.8.
       The default value of "diis_num_error_vect" is 5.
       The default value of "diis_start_error" is 0.01.
       The default value of "diis_end_error" is 10**(-8.0).

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
   
   <MEMORY>
      For settings of memory usage, write options in memory-directive.

      E.g.
         MEMORY
            (options)
         MEMORY_END

      -options
       "limit_heap" is only prepared. Note that this limitation is not 
       the exact limitation of heap usage. Please consider this option as a rough limitation.
       The value of this option should be written with the MByte unit.
       The default value is 256 [MB].

      E.g.
         MEMORY
            limit_heap 512
         MEMORY_END

   <MO Plot>
      writ MO plot directive

      E.g.
         MOPlot 
            (options)
         MOPlot_END

      -options
       "mo", "grid_number", "frame_length", and "file_prefix" are prepared.

       "mo" is index of the molcular orbital. mo=0 means the lowest energy MO.
       The default value of "mo" is not set.

       "grid_number" is the grid number of the frame in xyz-coordinates.
       The default values are 25, 25, and 25 for x, y, and z coordinates, respectively.

       "frame_length" is the length of the frame of each coordinate.
       The default values are 10, 10, and 10[angst.] for x, y, and z coordinates.

       "file_prefix" is a prefix of the file name to which the MO is written.
       The default values is "MO_".

      E.g.
         MOPLOT
            mo 5
            mo 8
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix MOPlot_
         MOPLOT_END

   <CIS>
      Write CIS-directive.

      E.g.
         CIS
            (options)
         CIS_END
   
      -options
       "davidson", "active_occ", "active_vir", "max_iter", "max_dim", "norm_tol", 
       and "nstates" are prepared as options.

       "davidson" should be set as "yes" or "no". 
       The default value of "davidson" is "yes".

       "active_occ" ("active_vir") is set to the number of occupied (virtual) orbitals
       if user set "active_occ" ("active_vir") to be greater than 
       the number of occupied (virtual) orbitals. 
       The default value of "active_occ" is 10. The default value of "active_vir" is 10.

       "nstates" is valid for the Davidson algorithm only, 
       hence "nstates" is set to "active_occ*active_vir" 
       in direct CIS algorithm (without the Davidson algorithem). 
       The default value of "nstates" is 5 for the Davidson algorithem.

       "max_iter" is valid for the Davidson algorithm only. 
       This option means the number of times of Davidson roop. 
       The default value of "max_iter" is 100.

       "max_dim" is valid for the Davidson algorithm only. 
       This option means the number of slater determinans used by expansion of the excited states. 
       Note that Hartree-Fock state (groudn state) is not included in the "max_dim".
       The default value of "max_dim" is 100.

       "norm_tol" is valid for the Davidson algorithm only. 
       This option means the max tolerance for the norm of the residual vectors.
       The default value of "norm_tol" is 10**(-6.0).

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


   <MD (Molecular dynamics)>
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
       The default value of "electronic_state" is 0. That is, 
       electronic ground state is default.

       The default value of "total_steps" is 10. 

       "dt" means the time width of molecular dynamics.
       "dt" should be set in femto-second.
       The default value of "dt" is 0.1[fs].

      E.g.
         MD
            total_steps 50
            electronic_state 0
            dt 0.05
         MD_END

   <MC (Monte Carlo)>
      Write MC-directive. The canonical sampling is only implemented.

      E.g.
         MC 
            (options)
         MC_END
  
      -options
       "total_steps", "electronic_state", "temperature"
       and "dr" are prepared as options.

       The default value of "total_steps" is 10. 

       "electronic_state" means the electronic eigen state 
       on which the system walks.
       The default value of "electronic_state" is 0. That is, 
       electronic ground state is default.

       "temperature" means the temperature in the MC sampling.
       The default value of "temeprture" is 300[K].

       "dr" means the max absolute displacement (step) width of each Cartesian coordinate.
       Namely, the actual displacement in the MC is in the range [-dr, dr).
       "dr" should be set in angstrom unit.
       The default value of "dr" is 0.1[angstrom].

      E.g.
         MC
            total_steps 50
            electronic_state 0
            dr 0.05
         MC_END

   <Principal Axes (Diagonalizing the inertia tensor)>
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

   <Rotate Molecule>
      Write "rotate" in theory-directive.

      E.g. 
         THEORY
            rotate
         THEORY_END

      -options
       "type", "origin", "axis", "angle" and "angles" are prepared as options.
       These options are written in rotate-directive. Examples are shown below.

       "type" indicates whether the rotating is carring out around a axis or acording to Euler angles.
       The default value of "type" is axis.

       "origin" indicates the origin of the rotation in angstrom unit.
       The default value of "origin" is center of mass.

       "axis" indicates a axis around which the rotation is carried out in angstrom unit.
       The default value of "axis" is z-axis.
       This option is valid only for "type" set as axis.

       "angle" indicates angle for the rotation around the "axis" in degree unit.
       The default value of "angle" is 0.
       This option is valid only for "type" set as axis.

       "angles" indicates Euler angles for the rotation in degree unit.
       The default values of "angles" are 0, 0, and 0.
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


   <Translate Molecule>
      Write "translate" in theory-directive.

      E.g. 
         THEORY
            translate 
         THEORY_END

      -options
       "difference" indicates difference for the translation in angstrom unit.
       This option is written in translate-directive.
       The default values are 0, 0, and 0.

       E.g. 
         TRANSLATE
            difference 12 30 45
         TRANSLATE_END



