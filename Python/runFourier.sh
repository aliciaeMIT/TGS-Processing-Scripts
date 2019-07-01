#This script makes all of the lammps input files and run files for HPC


printf '#Peak frequency presented: run (Random Seed) - Peak Frequency \n\n' > ./PeakFreq.txt
printf '#SAW speed presented: run (Random Seed) - SAW speed \n\n' > ./SAWv.txt
printf '#Acoustic Damping constant presented: run (Random Seed) - Tau \n\n' > ./Tau.txt
printf '#Thermal Diffusivity presented: run (Random Seed) - k (sigma_k)\n\n' > ./K.txt
printf '#SRO parameter presented: run (Random Seed) - alpha \n\n' > ./Alpha.txt

mkdir Transforms
mkdir Signals
mkdir AcousDampFits
mkdir TempProfs
mkdir SROparams
mkdir DataFiles

count=1
for i in {0..32767}
do
	if [ -a "COM-${i}.txt" ]
	then
		echo "${count} - ${i}"

        # TGS Analysis script
		sed -e "s/COM.txt/COM-${i}.txt/g" ./Fourier_MD.py > Fourier_MD-${i}.py
		sed -i -e "s/TEMP.txt/TEMP-${i}.txt/g" ./Fourier_MD-${i}.py
		sed -i -e "s/DATA.txt/DATA-${i}.txt/g" ./Fourier_MD-${i}.py
		sed -i -e "s/Signal.png/Signal-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/AcousDampFit.png/AcousDampFit-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/Transform.png/Transform-${i}.png/g" ./Fourier_MD-${i}.py
		sed -i -e "s/TempProf.png/TempProf-${i}.png/g" ./Fourier_MD-${i}.py
		python Fourier_MD-${i}.py

        #Short Range Order
        sed -e "s/rdf.txt/rdf-${i}.txt/g" ./SRO.py > SRO-${i}.py
        sed -i -e "s/SRO_params.txt/SRO_params-${i}.txt/g" ./SRO-${i}.py
        python SRO-${i}.py

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

        #Total SRO File
        printf "${count} (${i}):  " >> ./Alpha.txt
        sed '28q;d' SRO_params-${i}.txt >> ./Alpha.txt

		let "count=count+1"
	fi
done

rm Fourier_MD-*.py-e
rm SRO-*.py-e
mv Signal-*.png ./Signals/
mv Transform-*.png ./Transforms/
mv AcousDampFit-*.png ./AcousDampFits/
mv TempProf-*.png ./TempProfs/
mv SRO_params-*.txt ./SROparams
mv DATA-*.txt ./DataFiles

############
date
