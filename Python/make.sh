#This script makes all of the lammps input files and run files for HPC

num=5

for i in {1..50}
do
	seed=$RANDOM
	sed -e "s/<RANDOM>/${seed}/g" ./Template-FeNiCr.in > FeNiCr-${seed}.in
	sed -e "s/<RANDOM>/${seed}/g" ./script > script-${seed}
	qsub script-${seed}
done


############
date
