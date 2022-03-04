#!/bin/sh

#  Script.sh
#  
#
#  Created by Laura Mulvey on 15.07.21.
#

####
# step 1: move all rec files into one folder
mkdir rec_files
files=$(ls $1)
counter=1
for i in $files

do
    cd $1/$i
   cp parasites_READY.trees.ale.uml_rec $counter.ale.uml_rec
   
cd ../../

mv $1/$i/$counter.ale.uml_rec rec_files

let counter=counter+1

done

####
# step 2: extract all lines


counter=1
line=1017
mkdir brn
while [ $line -le 1049 ];
do
rec=$(ls rec_files)
for i in $rec
do
sed -n $line'p' rec_files/$i >> $counter.txt
done
mv $counter.txt brn
let line=line+1
let counter=counter+1

done


#######
# step 3 add them together and get the average

files=$(ls brn)
for i in $files
do
echo $i >> names.txt
done

mkdir aver
Rscript average.R

rm -r brn


############
# step 4 take out last line of each file and put it into 1
aver=$(ls -1 aver | sort -V)

for i in $aver
do
#echo $i
sed -n '1002p' aver/$i >> final.txt
done

rm -r rec_files
rm names.txt
rm -r aver
