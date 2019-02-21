#!/bin/sh

## Grab a random pair of POS/NEG data files from the TGS sync directory

cd /home/hereiam/projects/TGS-Random-Generator

filename_pos="$(find /home/hereiam/Documents/Sync/Machines/TGS/LSAW_Data/ -type f -name "*POS-1.txt" | sort -R | tail -1)"
filename_neg="$(echo $filename_pos | rev | cut -c 10- | rev)""NEG-1.txt"
# echo $filename_pos
# echo $filename_neg

cp $filename_pos POS.txt
cp $filename_neg NEG.txt

## Extract the rough grating spacing from the filename

temp=${filename_pos%%um-*}
grating="$(echo $temp | rev | cut -c -4 | rev)"
# echo $temp
echo Grating = $grating

## Extract interesting bits of info to overlay on the final image

StudyName="$(grep -nrH 'Study Name' POS.txt)"
StudyName="$(echo $StudyName | cut -c 22- | tr -d '\r')"
echo $StudyName
SampleName="$(grep -nrH 'Sample Name' POS.txt)"
SampleName="$(echo $SampleName | cut -c 23- | tr -d '\r')"
echo $SampleName
Date="$(grep -nrH 'Date' POS.txt)"
Date="$(echo $Date | cut -c 16- | tr -d '\r')"
echo $Date
Time="$(grep -nrH 'Time' POS.txt | head -n1)"
Time="$(echo $Time | cut -c 16- | tr -d '\r')"
echo $Time

## Create a MATLAB .m file to process the two data files

sed 's;<GRATING>;'$grating';' GoGo-Template.m > GoGo.m
sed -i "s;<TIMESTAMP>;\"Date: "$Date", Time: "$Time", \\\lambda="$grating"\\\mu m\";" GoGo.m
sed -i "s;<SAMPLESTAMP>;\"Study: "$StudyName", Sample: "$SampleName"\";" GoGo.m
sed -i 's/_/\\_/g' GoGo.m

/usr/local/MATLAB/R2018b/bin/./matlab -nodisplay -nodesktop -r "GoGo"

## Copy all the image files to a web directory for hosting

cp *.png /var/www/html/Data/TGS/

## Delete any MATLAB crash files, if they exist

rm /home/hereiam/matlab_crash_dump*
