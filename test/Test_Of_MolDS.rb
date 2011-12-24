#!/usr/bin/ruby -w

system("echo ")
system("echo '*****************************************'")
system("echo '***                                   ***'")
system("echo '***                                   ***'")
system("echo '***       Start Test for MolDS        ***'")
system("echo '***                                   ***'")
system("echo '***                power by ruby      ***'")
system("echo '*****************************************\n\n'")

system("echo '-------------------------------------------'")
system("echo '----------   Test of CNDO2/HF     ---------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_cndo2.in > temp.dat")
system("diff ch4_cndo2.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_cndo2.in > temp.dat")
system("diff ch4_cndo2.dat temp.dat")
system("echo '\n\n'")

system("echo '---------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet     ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_zindos_directCIS_singlet.in > temp.dat")
system("diff ch4_zindos_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_zindos_directCIS_singlet.in > temp.dat")
system("diff ch4_zindos_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '\t\t\t>>> C2H6 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < c2h6_zindos_directCIS_singlet.in > temp.dat")
system("diff c2h6_zindos_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < c2h6_zindos_directCIS_singlet.in > temp.dat")
system("diff c2h6_zindos_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet  ---------'")
system("echo '----------  With Davidson for the CIS  ---------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_zindos_davidsonCIS_singlet.in > temp.dat")
system("diff ch4_zindos_davidsonCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_zindos_davidsonCIS_singlet.in > temp.dat")
system("diff ch4_zindos_davidsonCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '\t\t\t>>> C2H6 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < c2h6_zindos_davidsonCIS_singlet.in > temp.dat")
system("diff c2h6_zindos_davidsonCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < c2h6_zindos_davidsonCIS_singlet.in > temp.dat")
system("diff c2h6_zindos_davidsonCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '-------------------------------------------'")
system("echo '----------   Test of MNDO/HF     ----------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_mndo.in > temp.dat")
system("diff ch4_mndo.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_mndo.in > temp.dat")
system("diff ch4_mndo.dat temp.dat")
system("echo '\n\n'")

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet      ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_mndo_directCIS_singlet.in > temp.dat")
system("diff ch4_mndo_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_mndo_directCIS_singlet.in > temp.dat")
system("diff ch4_mndo_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '-------------------------------------------'")
system("echo '----------   Test of AM1/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_am1.in > temp.dat")
system("diff ch4_am1.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_am1.in > temp.dat")
system("diff ch4_am1.dat temp.dat")
system("echo '\n\n'")

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet       ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_am1_directCIS_singlet.in > temp.dat")
system("diff ch4_am1_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_am1_directCIS_singlet.in > temp.dat")
system("diff ch4_am1_directCIS_singlet.dat temp.dat")
system("echo '\n\n'")

system("echo '-------------------------------------------'")
system("echo '----------   Test of PM3/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_pm3.in > temp.dat")
system("diff ch4_pm3.dat temp.dat")
system("echo '\n\n'")

system("echo MPI:no")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_pm3.in > temp.dat")
system("diff ch4_pm3.dat temp.dat")
system("echo '\n\n'")

system("rm -rf temp.dat")
