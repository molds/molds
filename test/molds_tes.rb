#!/usr/bin/ruby -w

system("echo ")
system("echo '*****************************************'")
system("echo '***                                   ***'")
system("echo '***                                   ***'")
system("echo '***       Start Test for MolDS        ***'")
system("echo '***                                   ***'")
system("echo '***                power by ruby      ***'")
system("echo '*****************************************'")
system("echo ")
system("echo ")


system("echo ------------------------------------------")
system("echo ---------- START: CH4 with CNDO2 ---------")
system("echo ------------------------------------------")
system("echo MPI_NODES:")
system("echo 1")
ENV["MKL_NUM_THREADS"] = "1"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "1"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_cndo2.in > temp.dat")
system("diff temp.dat ch4_cndo2.dat")
system("echo ")
system("echo ")

system("echo MPI_NODES:")
system("echo 1")
ENV["MKL_NUM_THREADS"] = "2"
system("echo MKL_NUM_THREADS:")
system("echo $MKL_NUM_THREADS")
ENV["OMP_NUM_THREADS"] = "2"
system("echo OMP_NUM_THREADS:")
system("echo $OMP_NUM_THREADS")
system("../src/MolDS.out < ch4_cndo2.in > temp.dat")
system("diff temp.dat ch4_cndo2.dat")
system("echo ------------------------------------------")
system("echo ---------- END: CH4 with CNDO2 -----------")
system("echo ------------------------------------------")


system("rm -rf temp.dat")
