#!/usr/bin/ruby -w

class TesterOmp
   @@surfixDat = ".dat"
   @@surfixInp = ".in"
   @@tempFile = "temp.dat"
   @@moldsBin = "../src/MolDS.out"
   @@command = "command: "
   def doesTestOmp(prefix, mklNumThreads, ompNumThreads)
      setPrefix(prefix)
      ENV["MKL_NUM_THREADS"] = mklNumThreads
      ENV["OMP_NUM_THREADS"] = ompNumThreads
      system("echo MPI:no")
      system("echo MKL_NUM_THREADS:$MKL_NUM_THREADS")
      system("echo OMP_NUM_THREADS:$OMP_NUM_THREADS")
      print @@command + @moldsCommand + "\n"
      system(@moldsCommand)
      print @@command + @diffCommand + "\n"
      system(@diffCommand)
      system("echo '\n\n'")
   end
   def setPrefix(prefix)
      @inputFile = prefix + @@surfixInp
      @outputFile = prefix + @@surfixDat
      @moldsCommand = @@moldsBin + " < " + @inputFile + " > " + @@tempFile 
      @diffCommand = "diff " + @outputFile + " " + @@tempFile
   end
   private :setPrefix
end

system("echo ")
system("echo '*****************************************'")
system("echo '***                                   ***'")
system("echo '***                                   ***'")
system("echo '***       Start Test for MolDS        ***'")
system("echo '***                                   ***'")
system("echo '***                    Power by Ruby  ***'")
system("echo '*****************************************\n\n'")

testerOmp = TesterOmp.new

system("echo '-------------------------------------------'")
system("echo '----------   Test of CNDO2/HF     ---------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of INDO/HF    -----------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_indo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_indo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet     ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2O <<<\n'")
prefix = "h2o_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet  ---------'")
system("echo '----------  With Davidson for the CIS  ---------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of ZINDO/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of MNDO/HF     ----------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet      ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet      ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of MNDO/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet-force --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet-force --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of AM1/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet       ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet       ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of AM1/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet-force  --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet-force  --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of PM3/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet       ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet       ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of PM3/HF-Force  --------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet-force  --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet-force  --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '---------- Test of PM3/PDDG/HF ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/CIS-singlet  ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/CIS-singlet  ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/HF-Force  ---------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/PDDG/CIS-singlet-force  ----'")
system("echo '---------  Without Davidson for the CIS    --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/PDDG/CIS-singlet-force  ----'")
system("echo '---------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("rm -rf temp.dat")
