#!/usr/bin/ruby -w

class TesterOmp
   @@surfixDat = ".dat"
   @@surfixInp = ".in"
   @@tempFile = "temp.dat"
   @@moldsBin = "../src/MolDS.out"
   def doesTestOmp(prefix, mklNumThreads, ompNumThreads)
      setPrefix(prefix)
      ENV["MKL_NUM_THREADS"] = mklNumThreads
      ENV["OMP_NUM_THREADS"] = ompNumThreads
      system("echo MPI:no")
      system("echo MKL_NUM_THREADS:")
      system("echo $MKL_NUM_THREADS")
      system("echo OMP_NUM_THREADS:")
      system("echo $OMP_NUM_THREADS")
      print @moldsCommand + "\n"
      system(@moldsCommand)
      print @diffCommand + "\n"
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
system("echo '***                power by ruby      ***'")
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

system("rm -rf temp.dat")
