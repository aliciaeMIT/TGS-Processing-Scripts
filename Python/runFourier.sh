#This script makes all of the lammps input files and run files for HPC


printf '#Peak frequency presented: run (Random Seed) - Peak Frequency \n\n' > ./PeakFreq.txt
printf '#SAW speed presented: run (Random Seed) - SAW speed \n\n' > ./SAWv.txt
printf '#Acoustic Damping constant presented: run (Random Seed) - Tau \n\n' > ./Tau.txt
printf '#Thermal Diffusivity presented: run (Random Seed) - k (sigma_k)\n\n' > ./K.txt

mkdir Transforms
mkdir Signals
mkdir AcousDampFits
mkdir TempProfs

count=1
for i in {0..32767}
do
	if [ -a "COM-${i}.txt" ]
	then
		echo "${count} - ${i}"
		sed -e "s/COM.txt/COM-${i}.txt/g" ./Fourier_MD.py > Fourier_MD-${i}.py
		sed -i -e "s/TEMP.txt/TEMP-${i}.txt/g" ./Fourier_MD-${i}.py
		sed -i -e "s/DATA.txt/DATA-${i}.txt/g" ./Fourier_MD-${i}.py
		sed -i -e "s/Signal.png/Signal-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/AcousDampFit.png/AcousDampFit-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/Transform.png/Transform-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/TempProf.png/TempProf-${i}.png/g" ./Fourier_MD-${i}.py
		python Fourier_MD-${i}.py

		#Peak Frequency file
		printf "${count} (${i}):  " >> ./PeakFreq.txt
		sed '5q;d' DATA-${i}.txt >> ./PeakFreq.txt

		#Saw Speed file
		printf "${count} (${i}):  " >> ./SAWv.txt
		sed '9q;d' DATA-${i}.txt >> ./SAWv.txt

		#Acoustic Damping Constant file
		printf "${count} (${i}):  " >> ./Tau.txt
		sed '11q;d' DATA-${i}.txt >> ./Tau.txt

		#Thermal Diffusivity file
		printf "${count} (${i}):  " >> ./K.txt
		sed '13q;d' DATA-${i}.txt >> ./K.txt

		rm Fourier_MD-${i}.py-e
		mv Signal-${i}.png ./Signals/
		mv Transform-${i}.png ./Transforms/
		mv AcousDampFit-${i}.png ./AcousDampFits/
		mv TempProf-${i}.png ./TempProfs/

		let "count=count+1"
	fi
done


############
date
